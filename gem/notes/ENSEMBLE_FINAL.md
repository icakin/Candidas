# Uncertainty ensemble — corrected and re-run (supersedes the retracted first attempt)

## What was fixed relative to the retracted version
1. Correct sigmas: **Seq2Topt RMSE 12.26 °C, Seq2Tm RMSE 7.57 °C** (predictor README).
   The first attempt used 12 °C for both.
2. Added the **systematic per-species offset** model — the error structure that CAN create
   species separation, because orthologues are near-identical so predictor errors are
   correlated within a species. The first attempt only used independent per-enzyme noise,
   which averages toward no separation.
3. The separation question is now asked **only among viable draws** (all four species retain
   growth at 34 °C). This removes the confound that sank the first attempt, where large
   independent noise collapsed everything and made the joint condition unreachable for
   reasons unrelated to separability (Jensen effect: activity enters the pool as 1/a, so
   zero-mean noise depresses growth by construction).
4. Two seeds, 500 draws pooled, Wilson confidence intervals reported.

## Sanity check
Zero-offset reproduces `etcgem_tpc_pred.csv` exactly for all four species at 34/40/44 °C
in both runs.

## Result (500 draws, seeds 2024 and 99, pooled)
- Viable draws (all species µ(34 °C) ≥ 0.20): **148 / 500**
- Observed phenotype (auris grows at 44 °C AND all three relatives fail at 40 °C)
  among viable draws: **0 / 148**, 95 % CI **0 – 2.5 %**
- Viable draws reaching the ≈16 °C Topt separation the counterfactual says is required:
  **3 / 148** — and none of those three reproduces the phenotype.
- *C. auris* ranks above every relative at 44 °C: **141 / 500 = 0.282**, 95 % CI
  **0.244 – 0.323** — the interval contains chance (0.25).
- Across all 500 draws only one produced the phenotype at all, and it required a **+40.4 °C**
  systematic Topt separation (>3σ) — and that draw was not viable at 34 °C.

Sensitivity to the permissive floor (seed 2024, N = 300): 0/95 at µ34 ≥ 0.20,
0/63 at ≥ 0.30, 0/33 at ≥ 0.40.

## Interpretation
Separation and viability are mutually exclusive under realistic predictor error. A
systematic offset large enough to separate the species also pushes cells below permissive
growth, so the model cannot produce *C. auris*'s high-temperature phenotype while remaining
a living cell anywhere in the credible parameter space. This is consistent with, and
independent of, the counterfactual result (≈16 °C Topt / ≈32 °C Tm separation required).

## Correction to an earlier claim
I previously reported that *C. auris* "ranks top at 44 °C in 39–42 % of draws vs 25 %
chance — a weak but real tendency." That came from the independent-noise ensemble. Under
the appropriate systematic-error model the estimate is 0.282 with CI 0.244–0.323, which
**includes chance**. The weak-tendency claim is withdrawn.

## Limitations (stated, not glossed)
- The systematic model assumes perfectly correlated error within a species; reality lies
  between that and the independent case. Both extremes were tested; neither reproduces the
  phenotype in a viable model.
- kcat uncertainty is not included in this run.
- 0/148 means <2.5 % at 95 % confidence, not zero.
- Nothing here bears on Fig 4 Panel A; the earlier "Panel A should be qualified" claim was
  itself an error (Jensen effect) and remains withdrawn.
