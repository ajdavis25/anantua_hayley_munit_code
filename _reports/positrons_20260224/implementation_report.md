# GRMONTY Positron Implementation Report (2026-02-24)

## Scope
This report summarizes the positron-capability changes integrated into `igrmonty` (GRMONTY), including physics conventions, code updates, tests, and a controlled pair-scaling sweep.

## Locked Physics/Units Decisions

### 1) `Ne_unit` / `Ne` meaning (code-side convention)
- `Ne` passed into emissivity/opacity routines is treated as baseline ion-associated density `n_i`.
- Pair content is parameterized by `positron_ratio = f` with:
- `n_+ = f n_i`
- `n_- = (1 + f) n_i`
- `n_lep = n_- + n_+ = (1 + 2f) n_i`
- Anchors:
- `src/jnu_mixed.c:52-66`
- `src/radiation.c:168-175`
- `src/main.c:108-114`

### 2) `kappa_es` units and scattering opacity form
- `kappa_es(ν, Θe)` is used as mass scattering opacity in `cm^2 g^-1`.
- Invariant scattering opacity is locked to:
- `alpha_inv_scatt = ν * kappa_es(ν, Θe) * (n_lep * m_p)`
- Unit check documented in code:
- `cm^2 g^-1 * g cm^-3 = cm^-1`, then multiplied by `ν` for invariant form.
- Anchor:
- `src/radiation.c:164-175`

### 3) Bremsstrahlung under pairs
- Minimum required behavior implemented:
- Electron-ion term scales as `n_i * n_lep`.
- Electron-electron same-sign term scales as `n_-^2 + n_+^2`.
- Anchors:
- `src/jnu_mixed.c:223-227`
- `src/jnu_mixed.c:233-240` (alternate branch also uses `n_i * n_lep`)

### 4) Downstream convention to avoid pair double counting
- GRMONTY now applies pair lepton scaling internally from `positron_ratio`.
- Downstream scripts should keep `M_unit` / density baryonic and must **not** pre-scale by `(1 + 2f)`.
- Anchors:
- `README:50-63`
- `src/main.c:110-114`
- `auto_munit_bracket.py:13-15,950-953`

## Code Changes (What Was Modified)

### Runtime parameters + CLI plumbing
- Added `positron_ratio` to runtime parameters.
- Added CLI/param switch for test execution without source edits:
- `run_tests` parameter
- `--run-tests`, `--run_tests`, and `--run-tests=<0|1>` support
- Files:
- `src/par.h:19,51`
- `src/par.c:12,33,48-63,89,123,125`
- `template.par:6,12`

### Main composition state and validation
- Added global `positron_ratio` state.
- Validates nonnegative finite value at startup.
- Prints composition convention and warning against external pair pre-scaling.
- Optional `run_all_tests()` execution path.
- Files:
- `src/decs.h:151,178`
- `src/main.c:74,95-115`

### Radiative microphysics updates
- Scattering opacity uses total leptons (`n_lep`) via `positron_ratio`.
- Synchrotron absorptivity paths scale with `Ne_lep`.
- Files:
- `src/radiation.c:164-175,188-190,238,265,296+`

### Emissivity and brems updates
- Added pair-density helpers (`n_lep`, `n_-`, `n_+`).
- Updated synchrotron and integrated emissivities to use total radiating leptons.
- Implemented pair-aware brems scaling described above.
- File:
- `src/jnu_mixed.c:47-66,171-240,248+,363+`

### Output metadata additions
- Writes positron ratio into output HDF metadata:
- iharm: `/params/electrons/positron_ratio` and `/params/electrons/positronRatio`
- riaf/sphere: `/params/positron_ratio`
- Files:
- `model/iharm/model.c:1430-1432`
- `model/riaf/model.c:646`
- `model/sphere/model.c:546`

### Documentation and tooling notes
- Added positron convention and usage notes to docs/templates.
- Added script-side warning for nonzero pair ratio.
- Files:
- `README:6-7,24,50-63`
- `template.par:6,12`
- `auto_munit_bracket.py:13-15,950-953`

### Build stability fix
- Prevented OpenMP/TLS + LTO linker breakage by disabling LTO when OpenMP is enabled.
- File:
- `makefile:56-60`

## Tests and Validation Hooks Added

### Internal regression switch
- `./grmonty --run-tests` now runs `run_all_tests()` and exits cleanly.
- Added pair-scaling test hook in test suite.
- File:
- `src/tests.c:6,16,19,109,133-135`

### Runtime check observed
- Command run: `./grmonty --run-tests`
- Observed exit code: `42` (expected test-suite exit path)
- Log includes: `testing pair scaling hooks (positron_ratio)`

## Controlled Pair Sweep Results
Sweep requested for one dump with `positron_ratio = 0, 0.5, 1.0`.

- Input dump: `/work/vmo703/grmhd_dump_samples/Ma-0.5_4000.h5`
- Spectra output dir: `/work/vmo703/igrmonty_outputs/m87/pair_sweep_20260224/`
- Collated table: `/work/vmo703/_reports/positrons_20260224/tables/pair_sweep_20260224.csv`
- Sweep logs/parfiles: `/work/vmo703/igrmonty/logs/pair_sweep_20260224/`

| positron_ratio | flux_jy | flux_scale_vs_pos0 | tau_scatt_mean_scale_vs_pos0 | tau_scatt_weighted_scale_vs_pos0 |
|---:|---:|---:|---:|---:|
| 0.0 | 0.2851732247 | 1.0000000000 | 1.0000000000 | 1.0000000000 |
| 0.5 | 0.4689221430 | 1.6443414120 | 2.0570899488 | 2.1256508463 |
| 1.0 | 0.8394235969 | 2.9435568423 | 2.9743584649 | 3.1878686993 |

Interpretation:
- Scattering depth trends increase with total lepton content as expected.
- Flux increases monotonically with pair fraction in this setup.
- Exact factors are not strictly linear due to transfer/absorption and Monte Carlo effects, but directionality is correct.

## Backward Compatibility
- Default `positron_ratio = 0.0` preserves baseline electron-ion behavior.
- New features are opt-in.
- Existing workflows remain valid if they do not externally pre-scale density for pairs.

## Changed Files Snapshot
- `README`
- `auto_munit_bracket.py`
- `makefile`
- `model/iharm/model.c`
- `model/riaf/model.c`
- `model/sphere/model.c`
- `src/decs.h`
- `src/jnu_mixed.c`
- `src/main.c`
- `src/par.c`
- `src/par.h`
- `src/radiation.c`
- `src/tests.c`
- `template.par`

## QC Figures (Generated 2026-02-24)
- Plot script:
- `/work/vmo703/_reports/positrons_20260224/scripts/make_pair_sweep_plots.py`
- Figures:
- `/work/vmo703/_reports/positrons_20260224/figures/pair_sed_overlay_nuLnu_Lsun.png`
- `/work/vmo703/_reports/positrons_20260224/figures/pair_scaling_relative.png`
- `/work/vmo703/_reports/positrons_20260224/figures/pair_scaling_absolute.png`

Checks covered:
- SED overlay compares `f_+ = 0, 0.5, 1.0` for the same dump/seed using angle-averaged `nuLnu`.
- Relative-scaling plot checks flux and scattering-depth scale factors vs `f_+=0`.
- Absolute-metrics plot checks raw flux and scattering optical-depth trends vs pair fraction.
