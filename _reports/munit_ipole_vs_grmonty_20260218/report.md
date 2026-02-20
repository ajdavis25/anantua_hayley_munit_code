# Munit Scaling and IPOLE-vs-GRMONTY Consistency Report

Generated: 20260218

## Executive Summary

This report audits how IPOLE and GRMONTY outputs in `/work/vmo703/` depend on `M_unit`, whether matching the 230 GHz bin is a valid cross-code strategy, and whether an IPOLE-like viewing restriction can be imposed on GRMONTY without recompiling.

Key conclusions:

1. IPOLE in this pipeline is a camera-based image calculation at a single observer setup: one `thetacam`, fixed image resolution (`nx`,`ny`), fixed image-plane FOV (`fovx_dsource`,`fovy_dsource`), and fixed observing frequency (`freqcgs=2.28e11 Hz`).
2. GRMONTY outputs are angle-binned spectra at infinity (theta bins, azimuth-averaged), not finite-FOV images. The tuned pipeline computes 4pi-averaged flux from `/output/nuLnu` and `/output/dOmega` and then applies `D=16.8 Mpc`.
3. In this implementation, `RHO_unit` and `Ne_unit` scale linearly with `M_unit`, `B_unit` scales as `M_unit^(1/2)`, and `Thetae_unit` is independent of `M_unit` (for the selected electron-model branches).
4. Empirically, reliable 230 GHz slopes `d ln F_nu / d ln M_unit` are generally >2 in available multi-Munit GRMONTY sets (median near 2.68). This is strongly nonlinear and frequency-dependent.
5. Matching only the 230 GHz bin is not sufficient to ensure broader SED equivalence: mm normalization can match while high-frequency Compton-dominated behavior and bolometric output differ significantly.
6. IPOLE-like FOV restriction does not exist in current GRMONTY outputs. The closest no-recompile option is post-processing by theta-bin selection (or coarse cone over theta bins). True camera/FOV matching requires code changes.

Concrete Munit offset example from paired tables: SANE RBETA a=-0.50 t=4000 pos0: IPOLE M_unit=1.064e+29, GRMONTY best=1.373e+29, delta=29.0%.

## Data Inventory and Reproducibility

- IPOLE HDF5 runs scanned: 580 (`ipole_outputs/M87/**/*.h5`).
- GRMONTY HDF5 runs scanned: 39 (`igrmonty_outputs/m87/**/*.h5`, parsable spectrum names only).
- IPOLE 230 GHz flux stats (Ftot): min=0.00293 Jy, median=0.498 Jy, max=9.76 Jy.
- GRMONTY 230 GHz postprocessed flux stats (4pi): min=0.0391 Jy, median=0.502 Jy, max=7.73 Jy.
- All generated artifacts are in this report folder (`tables/`, `plots/`, `notes.json`).
- Bracket/bisection/tuning scripts identified in-tree: `igrmonty/auto_munit_bracket.py`, `igrmonty/tune_munit_once.py`, `scratch/munit_finder.py`, `sgrA/auto_munit_finder_sgrA.py`, `munit_finder/optimizer.py`.

## 1) What IPOLE Restricts (Viewing Angle / FOV)

Evidence from parfiles and HDF5 headers:

- `ipole_pars/run_Ma+0.94_5000_RBETA_pos0.par` sets `thetacam 17`, `nx 160`, `ny 160`, `fovx_dsource 160`, `fovy_dsource 160`, `dsource 16.9e6`, `freqcgs 230.e9`.
- `runIPOLE.par` shows the same pipeline concept with `thetacam 163.0`, `nx=320`, `ny=320`, `fov=160.0`, `freqcgs=228000000000.0`.
- Output headers (e.g. `ipole_outputs/M87/betacrit/Ma+0.94_6000_CRITBETA_f0.5_bc0.1_1.00.h5`) store camera values directly in `header/camera/*` and scalar flux in `Ftot`.

Interpretation for this pipeline:

- Primary restriction is **single observer inclination** (`thetacam`) plus **finite image-plane FOV** (`fovx_dsource`,`fovy_dsource`) and fixed pixel grid (`nx`,`ny`).
- No evidence was found for a separate post-hoc acceptance cone in IPOLE outputs; restriction is applied by camera geometry at ray tracing time.
- Notable convention detail: parfiles often show 17 deg, while many output headers show 163 deg. This is consistent with a north/south camera-angle convention flip in wrapper usage.

## 2) What GRMONTY Outputs Represent

Code-level evidence:

- `igrmonty/src/main.c` header comment: spectrum is computed vs **latitudinal viewing angle**, averaged over azimuth.
- `igrmonty/model/iharm/model.h` defines `N_THBINS 18`.
- `igrmonty/model/iharm/model.c` bins escaped photons by BL theta in `record_super_photon()` and writes `/output/nuLnu` with shape `[N_TYPEBINS, N_EBINS, N_THBINS]` and `/output/dOmega`.
- `igrmonty/model/iharm/model.c` computes total luminosity by integrating over theta bins with `dOmega` weighting (`L += nuLnu * dOmega * dlE / (4pi)`).
- `igrmonty/auto_munit_bracket.py` then converts to observer flux at distance by integrating angular bins and dividing by `4pi D^2` (with `D=16.8 Mpc`).

Therefore, in this pipeline GRMONTY output used for Munit tuning is **not a camera flux for one line of sight**; it is a postprocessed angular integration over theta bins (effectively 4pi-averaged flux).

Important implementation mismatch to note:

- `igrmonty/notebooks/scrap.ipynb` example cells sum only `nuLnu[0,:,:]` over theta with `np.sum(..., axis=-1)` and then apply `1/(4pi D^2)` directly.
- `igrmonty/auto_munit_bracket.py` instead performs `dOmega`-weighted angular integration and includes all spectral components before converting to flux.

So notebook quicklooks are useful diagnostics but are not guaranteed to match the auto-tuning flux definition.

## 3) Expected Munit Scaling in This Implementation

From `igrmonty/model/iharm/model.c`:

- `RHO_unit = M_unit / L_unit^3`  =>  `RHO_unit ∝ M_unit` (MBH fixed).
- `Ne_unit = RHO_unit / (m_p + m_e)`  =>  `Ne_unit ∝ M_unit`.
- `B_unit = c * sqrt(4pi * RHO_unit)`  =>  `B_unit ∝ M_unit^(1/2)`.
- `Thetae_unit` is set by electron-model branch and constants (not by `M_unit`) in the common RBETA/CRITBETA branches used here.

Implications (qualitative):

- Electron scattering optical depth roughly scales with density (`tau_es ~ n_e * sigma_T * l`), so increases with `M_unit` if geometry is fixed.
- Synchrotron characteristic frequency shifts with magnetic field (`nu_c ~ B * theta_e^2`), so `nu_c` shifts upward as `M_unit^(1/2)` when `theta_e` is nearly unchanged.
- Compton importance increases with optical depth and therefore can change rapidly with `M_unit`; high-frequency tail sensitivity is expected.

Log-printed unit sanity check (same model, different `M_unit`):

- `igrmonty/logs/SANE_RBETA_a-0.5_t4000_pos0_trial01.log`: `M_unit=1.0124e29`, `B_unit=1136.69`, `Ne_unit=6.83597e7`, `RHO_unit=1.14402e-16`, `Thetae_unit=244.82`.
- `igrmonty/logs/SANE_RBETA_a-0.5_t4000_pos0_trial02.log`: `M_unit=1.56322e29`, `B_unit=1412.46`, `Ne_unit=1.05552e8`, `RHO_unit=1.76645e-16`, `Thetae_unit=244.82`.

This measured change matches the expected trends (`Ne_unit`,`RHO_unit` linear in `M_unit`; `B_unit` ~ sqrt scaling; `Thetae_unit` unchanged).

## 4) Empirical Munit Scaling Measurements

Measured from existing GRMONTY HDF5 spectra (no reruns), using the same 4pi conversion approach as `auto_munit_bracket.py`.

- Detailed slopes are in `tables/scaling_summary.csv`.
- Reliable subsets (>=3 Munit points and Munit dynamic range >=1.3) show 230 GHz slopes typically above 2.
- Slopes vary strongly with frequency and model state, confirming that one-frequency matching does not preserve full-spectrum behavior.

## 5) Is Matching 230 GHz Alone Sufficient?

Short answer: **No** for full physical comparability.

Frequency-anchor alignment detail:

- IPOLE outputs here are fixed at `freqcgs = 2.28e11 Hz` (228 GHz).
- GRMONTY tuning scripts target `230e9 Hz`, but with current log-spaced bins the nearest bin used in existing spectra is typically `2.381685551976e11 Hz` (~238.17 GHz), as recorded in `tables/scaling_summary.csv` (`freq_bin_center_hz` for the GHz anchors).

This means the practical anchor is not exactly identical unless the GRMONTY binning is adjusted or interpolation is used.

Valid use case:

- If the only requirement is to normalize each model to the same 230 GHz flux point, matching that bin is acceptable as a scalar calibration target.

What remains different even after 230 GHz matching:

- Frequency-dependent slope means nearby bands (e.g., 86/345 GHz) do not necessarily co-match.
- High-frequency emission (Compton-sensitive) changes disproportionately with Munit in several groups.
- Bolometric luminosity and radiative efficiency can differ even if one mm anchor matches.
- Cross-code mismatch remains if one code uses camera-L.o.S. flux (IPOLE) and the other uses 4pi-averaged postprocessed flux (current GRMONTY pipeline).

## 6) IPOLE vs GRMONTY Pairing at 0.5 Jy

- Pair table: `tables/flux_match_pairs.csv`.
- IPOLE side uses curated near-0.5 Jy entries from `data/summary_table.csv`.
- GRMONTY side uses best-near-target entries from `data/munits_tuning_history.csv` by `(state,model,spin,dump,pos)`.

- Median Munit ratio (GRMONTY / IPOLE) across available pairs: 0.701.
- Several paired cases show O(10-50%) Munit offsets despite both being tuned near 0.5 Jy at mm frequency.

## 7) Can GRMONTY Enforce IPOLE-like Viewing Restriction?

Current capability without code changes:

- GRMONTY already stores theta-binned spectra (`N_THBINS=18`), so postprocessing can select one theta bin (nearest to IPOLE inclination) or a coarse theta-cone average.
- This is still azimuth-averaged and does not impose finite image-plane/FOV selection.

Not available from current outputs:

- Per-photon direction catalogs are not saved in the HDF5 outputs, so arbitrary acceptance-cone cuts cannot be retroactively applied at higher angular resolution than stored bins.
- Camera-plane/FOV cuts analogous to IPOLE are not represented in GRMONTY outputs.

Optional minimal code-change path (not applied here):

1. Add optional output of escaped-photon direction bins with finer theta/phi (or photon list) before compression into current spectra.
2. Add postprocessing that applies user-specified angular cone and then computes Fnu at distance.
3. For true IPOLE-like comparison, add a camera/ray-plane module (substantially larger change than cone cut).

## 8) Recommendations

1. Keep 230 GHz matching only as an initial normalization, then compare at least 86/230/345 GHz plus one high-frequency diagnostic metric.
2. For IPOLE-vs-GRMONTY comparisons, align **both** distance and angular definition: use a GRMONTY theta-bin proxy near IPOLE inclination rather than 4pi-integrated flux alone.
3. If precise viewing-angle comparability is required, add a lightweight escaped-photon angular output (or finer angular histogram) to GRMONTY; this is the minimal feasible extension short of full camera imaging.
4. Track resolution explicitly (IPOLE data in this tree include both 160 and 320 px camera products) when assembling comparison subsets.

## Assumptions / Limitations

- GRMONTY 230 GHz conversion used `D=16.8 Mpc` following `igrmonty/auto_munit_bracket.py` and `igrmonty/notebooks/scrap.ipynb`.
- IPOLE files are single-frequency products (`header/freqcgs=2.28e11 Hz`), so full IPOLE SED-vs-Munit comparison is not available from these outputs alone.
- Some historical tuning rows reference spectrum paths that no longer exist; pairing uses recorded CSV values for best-near-target Munit where needed.

## Artifact Index

- `tables/run_manifest.csv`
- `tables/scaling_summary.csv`
- `tables/flux_match_pairs.csv`
- `plots/` (SED overlays, ratio plots, IPOLE image comparison)
- `notes.json`

## Optional (Not Applied): Minimal GRMONTY Patch Plan For Viewing-Cone Output

Goal: preserve current outputs, and add an extra optional output that enables post-hoc cone cuts around a target inclination.

Proposed files and edits:

1. `igrmonty/src/par.h`
- Add optional parameters:
  - `int save_escape_dirs` (0/1)
  - `int N_PHI_BINS_OUT` (e.g., default 36)
  - `int N_TH_BINS_OUT` (e.g., default 180)

2. `igrmonty/src/par.c`
- Parse these keys from `.par` files:
  - `save_escape_dirs`
  - `N_PHI_BINS_OUT`
  - `N_TH_BINS_OUT`

3. `igrmonty/model/iharm/model.c`
- In `record_super_photon()`:
  - Compute `(theta, phi)` at escape (already has BL theta; add phi from `ph->X[3]` normalized to `[0,2pi)`).
  - If `save_escape_dirs=1`, accumulate energy histogram into `escape_nuLnu[phi_bin][theta_bin][E_bin]`.
- In output writer section (`report_spectrum`):
  - Write optional datasets:
    - `/output_escape/nuLnu_phi_theta` (N_TYPEBINS x N_EBINS x N_PHI x N_TH)
    - `/output_escape/theta_edges`
    - `/output_escape/phi_edges`

4. Postprocessing script (new, optional)
- Add `igrmonty/postprocess_cone_flux.py`:
  - Inputs: `--spec path.h5 --theta0_deg 17 --half_angle_deg 2 --distance_mpc 16.8`
  - Integrates only bins within cone and returns `Fnu(ν)` and `F230`.

Example command sequence (after patch + rebuild only GRMONTY, not IPOLE):

```bash
cd /work/vmo703/igrmonty
make -j
./grmonty -par /work/vmo703/igrmonty/logs/SANE_RBETA_a-0.5_t4000_pos0_trial04.par
python3 /work/vmo703/igrmonty/postprocess_cone_flux.py \
  --spec /work/vmo703/igrmonty_outputs/m87/spectrum_Sa-0.5_4000_RBETA_pos0_trial04.h5 \
  --theta0_deg 17 --half_angle_deg 2 --distance_mpc 16.8
```

Notes:
- This remains an angular-cone spectrum restriction, not a full finite image-plane camera/FOV ray trace.
- For exact IPOLE equivalence, a camera/image-plane module is required (larger development effort).
