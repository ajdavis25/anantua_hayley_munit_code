# Phase-4 production staging — 2026-09-09

Staged while the legacy A/B arm (job 808101, anantuabhg) finishes. Everything
here is submit-ready EXCEPT the two model decisions listed at the bottom.

## What changed

1. **`/work/vmo703/data/final_grmonty_paper.csv`** (tuner warm-start table):
   - `MunitUsed_pos0/pos1` for all 24 rows replaced with the July-2026 tuned
     values read from each archived spectrum's `/params/M_unit`
     (`igrmonty_outputs/m87/pre_bugfix_2026-07/`, 48 files, all matched 1:1).
     April values were stale by 0.15–2.0x. Backup: `final_grmonty_paper.csv.bak_apr30`
     here (and the April version is committed in git).
   - New informational columns `bias_jul_pos0/1`: each branch's fitted bias
     from July (MADs 0.05–0.17, SANE wJET 1.6–10.2).
   - **Convergence: 48/48.** The Jul-21 freeze said 42/48, but post-freeze
     retries (history rows through Jul 28) converged the rest — every warm
     start is a converged value. (The history CSV mixes schemas: 19-col 2025
     header, 27-col current rows; parse positionally by width.)
   - **6 new rows (row_id 24–29): MAD RBETAwJET rh80**, spins +0.94/−0.5 x
     dumps 4000/5000/6000 — the go-forward paper configs. Warm start
     7.443e24 = P31 M_unit (1.3488e25) rescaled by the measured paper-literal
     F230(4pi) = 1.642 Jy (tuner power-law, p=2); pos1 = pos0 x 0.587
     (median July MAD pos1/pos0).

2. **`igrmonty/auto_munit_bracket.py`** (uncommitted, with the rest):
   - `constant_beta_paper_literal` is now a first-class wJET parameter:
     `WJET_DEFAULTS` sets it to **1** (paper-literal P_B = B²/8π per
     Anantua+2020 / Emami+2021 eq. 25; see
     `igrmonty/docs/2026-09-08_12pi_provenance_verdict.md`), overridable
     per-row via CSV column or `--constant-beta-paper-literal`.
   - It flows into every generated wJET par (verified:
     `sample_par_row24_pos0.par`) and into the tuning-history CSV (appended
     as the LAST column so older-era positional parsers stay aligned).
   - **Consequence:** SANE wJET reruns also run paper-literal now. Their
     supplement is nearly extinct (38 zones in the rh160 census), so the
     effect is ~nil, but it is a deliberate change vs. the July era —
     one equation everywhere.

## Launch (when the user gives the go)

    cd /work/vmo703/igrmonty
    sbatch --array=0-59%5 run_auto_munit.slurm

(The file's `#SBATCH --array=0-47%5` is overridden on the CLI; 60 tasks =
30 rows x pos0/pos1, verified against the launcher's own expansion logic.
Partitions already prefer `anantuabhg`; fit_bias 1 from bias 0.05 with
backoff — the tested July machinery.)

**Binary (added 09-15, launch blocker fixed):** the repo-root `./grmonty`
(Jun 11) predates `constant_beta_paper_literal` and SILENTLY IGNORES the key
— wJET branches would run legacy physics while their pars claim
paper-literal. The launcher now defaults `AUTO_MUNIT_GRMONTY_BIN` to
`scratch/p4_staging/build/grmonty_p4` (rebuilt 09-15 from the freeze commit
— githash 65835d9, clean, MD5 5acf470f0a4b7f5be58d18703186f49e; smoke-tested:
init prints "constant-beta P_B form: paper-literal B^2/8pi" on the sample
par) and refuses to launch if that binary is missing. The Jun-11 binary
remains untouched. Freeze commits: igrmonty 65835d9, outer repo 33cd9c9.
LAUNCHED 2026-09-15.

## Rescue (2026-09-25) — original array 812893 cancelled at 9/60 usable

Three failure classes, all diagnosed (docs/2026-09-25_critbeta_floor_decision.md
+ agenda decision 6): Crit-β floor 1e-3 killed photon generation (24 tasks);
bias-guard-5 backoff spiral budget-killed MAD RBETA (11+4 tasks); MAD wJET
needed multicore + hit a fitter-path anomaly (fixed-bias workaround; fitter
bug still open — see below). Rescue binary **grmonty_p4b** (githash 388a18d
clean, MD5 26b96ee331e008993b2650ce8a067041; crit-β validated 71,831 photons
where p4 made 0). Freeze commits: igrmonty 46c1118+388a18d, outer f6efd0a+7515c2d.

Rescue arrays (all -p anantuabhg, subdir run_2026-09-15, binary p4b):
- **818690 p4r_crit** — tasks 24-47%5, fresh (floor fixed), 5d wall.
- **818691 p4r_resume** — tasks 2,3,11-23%5, --resume reuses on-disk trials,
  abort_ratio 100, 5d wall.
- **818692 p4r_madwjet** — tasks 48-59%3, 16 cpus/task, fit_bias 0 bias 0.05
  (A/B-proven; bypasses the fitter anomaly), 7d wall.
- **818706 p4_qa_watch** — QA watchdog (14d): runs the QA harness on every
  landing final spectrum + logs bad task terminals.
  Live feed: `tail -f igrmonty_outputs/m87/_qa/p4_qa_live.log`

OPEN (non-blocking): grmonty's internal bias fitter measures ratio=0 on some
MAD configs while the identical par at fit_bias 0 generates and scatters
normally. Confirmed scope (09-27): ALL MAD CRITBETA (both spins — killed the
12 p4r_crit tasks 36-47, resubmitted fixed-bias as **820422 p4r_madcrit**,
bias 0.05 = July's fitted value) and MAD wJET a+0.94; MAD RBETA and all SANE
fit fine. Root cause unisolated — group/debt list. Also: IPOLE needs BOTH
one-line gates (paper-literal P_B and crit_floor 3e-2) before P4.3 imaging.

RESULT (09-27, first QA sweep of 26 finals): all F230 anchors and provenance
PASS. The 12 X-ray-gate FAILs are physics, not mechanics: **every MAD wJET
branch tuned to 0.5 Jy overproduces the 2017 core X-ray at 4pi** (a-0.5:
1.05-2.3e41 = 2.4-5.2x over; a+0.94: 5.7e40-2.5e41 = 1.3-5.7x). In the 17-deg
observer cone a+0.94 drops to 1.3-5.1e40 (5/6 under the gate) but a-0.5 stays
1.6-3.3x OVER in-frame -- retrograde MAD wJET at beta_e0=0.1/rh80 is
X-ray-excluded in both frames; prograde survives only in the observer frame
(tuning-frame decision now has teeth). Under the legacy 12pi form all of
these would be ~9x higher = instantly excluded. beta_e0 scan (agenda
decision 5; Emami best-bet is 1e-2, ours 0.1) is the natural next knob.

## Still waiting on ashton's calls (NOT staged)

- **β-arm scope** (`jet_beta_cut`, currently 0.1 = status quo): the β ≤ 0.1
  override arm carries ~94% of SANE wJET emission (mid-latitude disk/wind
  painted as jet). Tighten to (σ AND β) / lower the cut / keep. A change
  lands as a CSV column or WJET_DEFAULTS edit — one line, then relaunch.
- **σ trust boundary on MAD** (Ryan+2018 §3.3 precedent): status quo = emit
  everywhere. No emission-cut knob exists in grmonty wJET modes today
  (σ ≥ 10 carries only ~1.5% of MAD proxy); implementing one is a small code
  change if sheath-only wins.
- **ipole**: the same one-line paper-literal gate must go into the ipole
  working copy before the P4.3 imaging comparison (spectra unaffected).

## Verification artifacts

- `sample_par_row24_pos0.par` — full production par for the first new row.
- Staging log: rerun `stage_phase4.py` ONLY after restoring the backup
  (it appends notes/rows; not idempotent).
