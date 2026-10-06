# MWL SED check — existing m87 corpus vs observed core anchors

Generated 2026-08-06 by `_qa/mwl_sed_check.py` + `_qa/mwl_sed_figures.py`.

Motivation: committee remark that Doppler boosting of the M87 jet is verified
in radio/optical/X-ray. The same observing campaigns give quantitative anchors
that constrain the wJET supplement. This check applies them to the spectra
already on disk (it does NOT cover the pending Phase-4 fixed-supplement reruns).

## Anchors

- **230 GHz**: 0.5 Jy — the pipeline tuning target (`auto_munit_bracket.py`
  `F_TARGET_JY`), consistent with the EHT 2017 compact-flux scale (~0.5–1 Jy).
- **2–10 keV core**: L_X = (4.4 ± 0.1)e40 erg/s — quasi-simultaneous
  Chandra+NuSTAR, 2017 Apr 12/14, core dominant over HST-1 in the low state
  (EHT MWL WG 2021, ApJL 911 L11; arXiv:2104.06855).

## Method

Conventions reused verbatim from `igrmonty/tools/viewing_cone_postprocess.py`
(nuLnu in L_sun; θ folded to [0°,90°], 18 × 5° bins; dΩ-weighted):
`baseline_4pi` = tuner-style all-sky average; `cone_renorm_4pi_equiv` =
isotropic equivalent seen from 17° ± 10° (bins 1–4, the P4.3 convention;
163° folds to 17°). L_X = ∫ nuLnu dlnν over 2–10 keV (9 grid bins, all
populated in every run). 48 canonical spectra (2 stale files were also
checked on the first pass, then quarantined — see below).

## Results

1. **X-ray gate: PASS for all 48 canonical runs** (`fig1_lx_gate.png`).
   L_X(17°)/L_obs by family (min–max over 2 spins × 3 dumps):
   MAD Crit-β 5e-4–1.8e-3 · MAD R-β 0.004–0.075 ·
   SANE Crit-β wJET 0.0015–0.217 · SANE R-β wJET 7e-4–0.198.
   Worst headroom 0.22× (Sa−0.5_5000 CRITBETAwJET pos1) — 4.6× under the
   bound. pos1 (pair-loaded) runs are ×3.3 (median) more X-ray-luminous than
   matched pos0; scattered (Compton) fraction of total L reaches 0.50 in the
   hottest wJET pos1 runs. The only OVERSHOOT on the first pass was the
   stale `..._5000_RBETAwJET_rh160_pos1_trial01.h5` (untuned iteration, 6 Jy
   at 230 GHz, 7.6× the X-ray bound); it and `..._CRITBETA_pos0_TEST.h5`
   were moved to `test_scrap/` on 2026-08-06, and the CSVs/figures were
   regenerated over the clean 48-file top level.

2. **Tuning-frame gap** (`fig2_f230_frame_gap.png`): every run is tuned on the
   4π-averaged flux (0.49–0.57 Jy achieved), but the 17°-cone flux the
   observer actually measures is lower in every family: MAD Crit-β 1.3–2.2×,
   MAD R-β 1.4–6.6×, SANE −0.5 wJET 1.3–3.1×, **SANE +0.94 wJET 4–25×**
   (0.02–0.13 Jy). Consequences: (a) comparisons to the EHT compact flux are
   family-inconsistently normalized; (b) an IPOLE image at i=163° tuned to the
   same target will land on a different M_unit than these grmonty values,
   most severely for Sa+0.94 wJET; (c) retuning in the cone frame would raise
   M_unit by ≈ gap^(1/P) (tuner p_factor P≈2) and raise L_X superlinearly
   (Compton roughly ∝ M_unit²⁺) — the worst family's headroom would shrink
   from ~4.6× toward ~2×, and Sa+0.94 would lose 1–2 orders of its margin.

3. **M_unit provenance split**: the canonical corpus on disk is the **July
   2026 retune** under the wjet-era code (tuning history 2026-05-01 →
   2026-07-28; e.g. Sa−0.5_5000 RBETAwJET pos1: warm-start at the frozen
   7.147e28 gave 5.56 Jy on 05-01, converged at 2.805e28 → 0.506 Jy on
   07-11). `data/final_grmonty_paper.csv` is the **April 2026 pre-supplement
   freeze**; its M_units differ from the files on disk by 0.15–1.95×.
   Any paper table should cite the July values recorded inside the h5 files
   (`params/M_unit`), not the April CSV.

## Caveats

- θ bins are folded about the equator: each cone value averages the
  approaching and receding funnel. True approaching-side (boosted) fluxes are
  somewhat higher than reported here; only unfolded binning or IPOLE at 163°
  can quantify that.
- This corpus predates the supplement B-scaling fix — Phase-4 reruns will
  move these numbers. The check is push-button to rerun on new outputs.
- pos1 spectra predate the PP-1 e⁺e⁻ brems fix (commit 6bd8db9): NR-branch
  channel only, subdominant for these models.
- Ns = 1e6 per run; THETAE_MAX = 1000 active in all runs.

## Era separation for future (fixed-code) runs

Added 2026-08-06 so post-fix production never mixes with this corpus:

- `auto_munit_bracket.py --run-subdir NAME` writes spectra to
  `igrmonty_outputs/m87/NAME/` and par/log files to `igrmonty/logs/NAME/`
  (both auto-created). Since 2026-08-06 the default (no flag) is an auto
  date-stamped `run_YYYY-MM-DD/`; under `run_auto_munit.slurm` the date is
  the array job's SUBMISSION date, so all tasks of one campaign share a
  folder even when they start on different days. `--run-subdir .` restores
  the old top-level behavior. To `--resume` a campaign started on an
  earlier day, pass its folder name explicitly.
- `run_auto_munit.slurm` passes it through:
  `sbatch --export=ALL,AUTO_MUNIT_RUN_SUBDIR=postfix_2026-08 run_auto_munit.slurm`
- The tuning-history CSV stays shared (each row records its spec/par/log
  paths, so eras remain traceable there).
- `wjet_single_par.slurm` runs take their spectrum path from the .par file
  itself — put the era subdir in the par's output path for those.
- The pre-fix corpus (the 48 canonical files) was archived to
  `igrmonty_outputs/m87/pre_bugfix_2026-07/` on 2026-08-06 (user decision:
  not paper-ready; Phase-4 fixed-code reruns supersede it). The matching
  56 .par + 56 .log files were archived the same day to
  `igrmonty/logs/nonpaper_ready/`. It remains the
  warm-start source for Phase-4 M_units (`params/M_unit` inside each h5).
  Historical CSVs (munits_tuning_history, wjet_provenance_report) still
  record the old top-level paths; `tag_wjet_provenance.py` globs
  recursively so it can be rerun anytime. Only legacy `tune_munit_once.py`
  or an explicit `--run-subdir .` still writes to the (now-empty) top
  level. A root-level `--resume` no longer finds these spectra.

## Files

- `mwl_sed_check.csv` (per file × mode), `mwl_sed_verdicts.csv` (gates),
  `plots/mwl_sed/fig{1,2,3}*.png`.
- Rerun: `/work/vmo703/ipole_venv/bin/python _qa/mwl_sed_check.py [--dir CORPUS_DIR]`
  then `... mwl_sed_figures.py [--dir CORPUS_DIR]`. Default `--dir` is the
  archived `pre_bugfix_2026-07/` corpus; point it at a new era subdir (e.g.
  `postfix_2026-08/`) after Phase-4 runs land. Output CSV/figure names are
  fixed, so a rerun overwrites the previous era's copies.
