# Confirm-run overrides to the Phase-4 pass/fail table

The pass/fail table and factorial report raw campaign values. Where a gate
verdict sat near the line on lean Compton statistics, dedicated confirm runs
(same tuned M_unit, higher fixed bias, independent seed, outside the
production folder) supersede the campaign value. Record:

## Sa-0.5_5000_CRITBETAwJET_rh20_bc1_f0.5_pos1  —  FAIL → **PASS (confirmed 2026-10-03)**

| run | bias | seed | scatter ratio | L_X(2–10 keV, 4π) | ×gate |
|---|---|---|---|---|---|
| production (campaign) | 0.05 | -1 | 0.0016 | 6.71e40 | 1.52 FAIL |
| confirm (job 825417, 33 min @16c) | 0.5 | 43 | 0.0100 | 1.93e40 | **0.44 PASS** |

Verdict: the campaign value was low-bias Monte-Carlo noise (the SANE Crit-β
X-ray band is scattering-starved at bias 0.05; a single high-weight packet
dominated the band). The confirm at 10× sampling passes with margin.

**Consequences:**
- The SANE+wJET column is fully clean → every gate failure in the campaign
  lives in {MAD} × {+wJET}. The factorial's original design claim is restored
  exactly.
- The "positrons tip marginal cells over" pattern now rests solely on the
  MAD CRITBETAwJET pos1 branches (1.1–1.5× over, also bias-0.05 statistics,
  though with f_compton ≈ 0.33–0.48 their sampling is far healthier than this
  SANE case). RECOMMENDED before the paper quotes that cell: same 30-minute
  confirm treatment for the five MAD CRITBETAwJET pos1 branches (and the
  near-boundary pos0 at 0.89×). The 12 MAD RBETAwJET failures are robust —
  margins 2.4–5.7× with strong Compton sampling.
- Confirm spectra live in scratch/p4_staging/confirm/ (never in the
  production run folder, so the campaign record stays raw).
