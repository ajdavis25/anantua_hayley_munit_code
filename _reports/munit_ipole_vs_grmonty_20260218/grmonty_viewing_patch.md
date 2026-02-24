# Optional GRMONTY Core Patch: Runtime Cone Filter During Photon Accumulation

Status: **design only, not applied**.

Reason: postprocessing cone tool is complete and validated. A core accumulation filter is feasible but touches hot-path photon accounting in `record_super_photon()` and should be compiled/tested in a controlled branch before enabling in production runs.

## Goal
Add runtime options to keep current behavior by default (`all_sky`) and optionally restrict accumulated photons to a viewing cone (`cone`) around `thetacam_deg` with half-angle `cone_half_angle_deg`.

## Scope (minimal)
- Parse three new parfile keys in `igrmonty/src/par.c`:
  - `view_mode` (`0=all_sky`, `1=cone`)
  - `thetacam_deg`
  - `cone_half_angle_deg`
- Add these fields to `Params` in `igrmonty/src/par.h`.
- In `igrmonty/model/iharm/model.c`, inside `record_super_photon()`:
  - compute photon BL `th` (already available),
  - convert to folded angle used by existing theta bins,
  - reject photon if outside selected cone in cone mode.
- Write view settings to HDF5 output metadata under `/params/view/*`.

## Proposed diff (not applied)
```diff
diff --git a/igrmonty/src/par.h b/igrmonty/src/par.h
--- a/igrmonty/src/par.h
+++ b/igrmonty/src/par.h
@@
   double jet_thetae;
   double jet_ne_mult;
+
+  // viewing restriction (optional)
+  int view_mode;                // 0=all_sky (default), 1=cone
+  double thetacam_deg;          // camera theta in degrees
+  double cone_half_angle_deg;   // cone half-angle in degrees
 
   char loaded;
 } Params;

diff --git a/igrmonty/src/par.c b/igrmonty/src/par.c
--- a/igrmonty/src/par.c
+++ b/igrmonty/src/par.c
@@
   params->jet_thetae = 0.0;
   params->jet_ne_mult = 1.0;
+
+  params->view_mode = 0;
+  params->thetacam_deg = 17.0;
+  params->cone_half_angle_deg = 180.0;
@@
     read_param(line, "jet_thetae", &(params->jet_thetae), TYPE_DBL);
     read_param(line, "jet_ne_mult", &(params->jet_ne_mult), TYPE_DBL);
+
+    read_param(line, "view_mode", &(params->view_mode), TYPE_INT);
+    read_param(line, "thetacam_deg", &(params->thetacam_deg), TYPE_DBL);
+    read_param(line, "cone_half_angle_deg", &(params->cone_half_angle_deg), TYPE_DBL);
 
     // set model parameters
     try_set_radiation_parameter(line);

diff --git a/igrmonty/model/iharm/model.c b/igrmonty/model/iharm/model.c
--- a/igrmonty/model/iharm/model.c
+++ b/igrmonty/model/iharm/model.c
@@
 static double jet_thetae = 0.0;
 static double jet_ne_mult = 1.0;
+static int view_mode = 0;
+static double view_thetacam_deg = 17.0;
+static double view_cone_half_angle_deg = 180.0;
+
+static inline double fold_theta_deg(double theta_deg)
+{
+  double th = fmod(theta_deg, 360.0);
+  if (th < 0.0) th += 360.0;
+  if (th > 180.0) th = 360.0 - th;
+  if (th > 90.0) th = 180.0 - th;
+  return th;
+}
+
+static inline int photon_in_view_cone(double th_bl_rad)
+{
+  if (view_mode == 0) return 1;
+  double th_deg = th_bl_rad * 180.0 / M_PI;
+  double th_fold = fold_theta_deg(th_deg);
+  double tc = fold_theta_deg(view_thetacam_deg);
+  double a = view_cone_half_angle_deg;
+  if (a < 0.0) a = 0.0;
+  double lo = tc - a;
+  double hi = tc + a;
+  if (lo < 0.0) lo = 0.0;
+  if (hi > 90.0) hi = 90.0;
+  return (th_fold >= lo && th_fold <= hi);
+}
@@
   bl_coord(ph->X, &r, &th);
+  if (!photon_in_view_cone(th))
+    return;
@@
     jet_thetae = params->jet_thetae;
     jet_ne_mult = params->jet_ne_mult;
+    view_mode = params->view_mode;
+    view_thetacam_deg = params->thetacam_deg;
+    view_cone_half_angle_deg = params->cone_half_angle_deg;
@@
   h5io_add_data_int(fid, "/params/N_THBINS", N_THBINS);
+  h5io_add_group(fid, "/params/view");
+  h5io_add_data_int(fid, "/params/view/mode", view_mode);
+  h5io_add_data_dbl(fid, "/params/view/thetacam_deg", view_thetacam_deg);
+  h5io_add_data_dbl(fid, "/params/view/cone_half_angle_deg", view_cone_half_angle_deg);
```

## Expected behavior
- Backward compatible defaults (`view_mode=0`) keep current outputs unchanged.
- Cone mode (`view_mode=1`) zeroes out contributions from photons outside cone before theta-bin accumulation.
- Outputs remain schema-compatible (`/output/nuLnu`, `/output/dOmega` still present).

## Example parfile additions
```text
view_mode 1
thetacam_deg 17.0
cone_half_angle_deg 10.0
```

## Risks
- `record_super_photon()` is performance-critical; extra logic may impact runtime.
- Existing theta folding means this is still a folded-theta cone, not full camera/FOV image-plane equivalence.
- If users set inconsistent `thetacam_deg` conventions (17 vs 163), folding hides sign/hemisphere information.

## Rollback plan
- Remove/ignore `view_*` keys by resetting defaults and deleting the cone gate in `record_super_photon()`.
- Existing outputs from all-sky mode remain directly comparable to current pipeline.

## Recommendation
Keep this core patch optional. For current analysis/comparisons, use the no-recompile postprocess utility:
`igrmonty/tools/viewing_cone_postprocess.py`.
