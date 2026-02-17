# GRMONTY QA Report

Generated: 2026-02-17 00:17:22 UTC
Base directory: `/work/vmo703/igrmonty_outputs/m87`
Selection mode: `crossref-logs` + `canonical-only`
Logs directory: `/work/vmo703/igrmonty/logs`
Log configs discovered: **28**
Filtered out non-canonical filenames: **7**
Filtered out not matched to log configs: **0**
Log configs with no canonical `.h5` found: **3**

Missing canonical outputs for log configs:
- `('M', -0.5, 4000, 'RBETA', 'OFF', 0)`
- `('M', -0.5, 4000, 'RBETA', 'ON', 0)`
- `('M', 0.94, 4000, 'RBETA', 'ON', 0)`

## Overall Counts

- Total `.h5` outputs found: **25**
- Open failures: **0**
- Hard sanity failures: **0**
- Suspicious (soft) failures: **0**
- Pass sanity checks (`looks good`): **25**
- Total non-pass files: **0**

Artifacts:
- Inventory CSV: `/work/vmo703/igrmonty_outputs/m87/_qa/grmonty_output_inventory.csv`
- Sanity CSV: `/work/vmo703/igrmonty_outputs/m87/_qa/grmonty_h5_sanity_report.csv`
- Pairwise CSV: `/work/vmo703/igrmonty_outputs/m87/_qa/grmonty_pairwise_report.csv`
- Plot directory: `/work/vmo703/igrmonty_outputs/m87/_qa/plots`

## Top 10 Worst Offenders

- No offenders. All files passed current checks.

## Summary by Heating Model

| Group | Total | Pass | Fail |
|---|---:|---:|---:|
| CRITBETA | 11 | 11 | 0 |
| RBETA | 14 | 14 | 0 |

## Summary by WJET Flag

| Group | Total | Pass | Fail |
|---|---:|---:|---:|
| OFF | 17 | 17 | 0 |
| ON | 8 | 8 | 0 |

## Summary by State (MAD/SANE)

| Group | Total | Pass | Fail |
|---|---:|---:|---:|
| MAD | 8 | 8 | 0 |
| SANE | 17 | 17 | 0 |

## Summary by Spin

| Group | Total | Pass | Fail |
|---|---:|---:|---:|
| -0.5 | 14 | 14 | 0 |
| 0.94 | 11 | 11 | 0 |

## Pairwise Comparative Checks

- Heating pairs (`RBETA` vs `CRITBETA`) evaluated: **10**
- Suspicious heating pairs: **0**
- WJET OFF/ON pairs evaluated: **8**
- Suspicious WJET pairs: **0**
- WJET behavior summary (based on median ON/OFF ratio): ON>OFF in **8** pairs, ON<OFF in **0** pairs.

## Common Problem Patterns

- No failure patterns detected under current sanity thresholds.

