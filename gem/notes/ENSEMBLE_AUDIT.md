# Audit of the uncertainty ensemble (FIG_ensemble) — retracted as delivered

## Verified (foundation is sound)
- Activity function identical to `build_etcgem_tpc.py`: peak = exp(-(T-Topt)²/2σ²),
  death = 1/(1+exp((T-Tm)/w)), floor 1e-6, cost = base/activity.
- Calibration constants read from `etcgem_calib.json` (sig 10.2314, w 8.7885, P 0.35697,
  SCALE 2.07412) — same source as the original pipeline.
- Protein-pool bound held at P0 throughout; not mutated inside the growth call.
- Zero-jitter reproduces `etcgem_tpc_pred.csv` exactly for all four species at 34/40/44 °C.

## Defect 1 — wrong sigma for Tm
Delivered run used σ = 12 °C for BOTH Topt and Tm. The predictor's own README gives
**Seq2Topt RMSE = 12.26 °C, Seq2Tm RMSE = 7.57 °C**. Tm noise was inflated ~60 %.
Re-run with correct σ changes the numbers only slightly (auris median µ44 0.0334 → 0.0364;
joint 0.05 → 0.04; rank 0.39 → 0.36) but the delivered figure is built on wrong inputs.

## Defect 2 — Panel A interpretation is invalid (Jensen effect)
Median predicted µ(44 °C) vs noise level (independent per-enzyme, Tm scaled 7.57/12.26):

| σ(Topt) | auris | haemulonii | duobushaemulonii | parapsilosis |
|---|---|---|---|---|
| 0.00 | 0.4474 | 0.3833 | 0.2843 | 0.4496 |
| 2.00 | 0.4374 | 0.3615 | 0.2729 | 0.4407 |
| 4.00 | 0.3700 | 0.3588 | 0.2574 | 0.3996 |
| 6.00 | 0.3042 | 0.2661 | 0.2165 | 0.3586 |
| 8.00 | 0.2460 | 0.1517 | 0.1473 | 0.2544 |
| 12.26 | 0.0529 | 0.0165 | 0.0131 | 0.0313 |

Growth declines monotonically for every species. Activity enters the pool constraint as
1/a, so by Jensen's inequality zero-mean noise raises E[1/a] and depresses predicted
growth **by construction**. This is a property of the noise, not a finding about the model.
**The claim that "Fig 4 Panel A should be qualified" is withdrawn** — Panel A's original
statement about the point-estimate model was correct.

## Defect 3 — the headline 0/100 is confounded
At σ = 12 all species collapse below the permissive floor, so "relatives fail at 40 °C AND
all species retain growth at 34 °C" cannot occur regardless of separability. The statistic
measures noise magnitude, not whether the model can separate species.

## Never run — the actual separability test
Independent per-enzyme noise averages toward no separation. The structure that *could*
create separation is a **systematic per-species offset** (orthologues are near-identical,
so predictor errors are correlated within a species). I asserted that correlated errors
"would reduce separation further" — untested, and probably backwards.

## What survives
*C. auris* ranks above every relative at 44 °C in ~30–45 % of draws at every σ from 2 to
12.26, vs 25 % by chance — a weak, stable tendency consistent with its marginally warmer
predicted Topt. This does not depend on the defects above.

## Corrected design (not yet run)
1. σ(Topt) = 12.26, σ(Tm) = 7.57.
2. Add the **systematic per-species offset** ensemble as the separability test.
3. Replace the confounded joint statistic: **condition on viable draws** (all species retain
   permissive growth), then ask whether separation occurs within that subset.
4. Drop or reframe Panel A; the Jensen effect is a diagnostic, not a result.
5. Report Monte-Carlo error; N = 100 gives only ~±3 %.
