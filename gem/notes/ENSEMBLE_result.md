# ⚠ RETRACTED — DO NOT USE

Self-audit found three defects. See ENSEMBLE_AUDIT.md.
1. Wrong sigma: Tm used 12 C; the predictor README gives Seq2Tm RMSE = 7.57 C (Seq2Topt = 12.26 C).
2. Panel A interpretation invalid: the growth collapse under noise is the Jensen effect
   (activity enters the pool as 1/a, so zero-mean noise inflates E[1/a]) — expected, not a finding.
3. The headline 0/100 is confounded: at sigma=12 all species collapse below the permissive floor,
   so the joint condition is unreachable for reasons unrelated to species separability.
The separability test (systematic per-species offset) was never run.

---

# Uncertainty ensemble (ChatGPT step 3) — result

**Question asked:** is the etcGEM's failure an artifact of treating predicted Topt/Tm/kcat
as exact, or does the observed phenotype remain implausible across the whole credible
parameter space?

**Answer: it remains implausible — but the point estimates themselves are not robust, and
that qualifies one statement in Fig 4.**

## Design
Per-enzyme Topt and Tm resampled from the predictor error distribution (σ = 12 °C, the
published Seq2Topt/Seq2Tm RMSE; σ = 6 °C as an optimistic variant), 100 draws, all four
species, evaluated at 34 / 40 / 44 °C. Errors drawn **independently per species** — the
case *most favourable* to the model, since real orthologue errors are correlated (near
identical sequences give near identical predictions **and** near identical errors), which
would reduce separation further. Sanity-checked: zero jitter reproduces
`etcgem_tpc_pred.csv` exactly for all species and temperatures.

## Findings

**1. Point estimates are not robust (this qualifies Fig 4 Panel A).**
At the published σ = 12 °C, median predicted µ(44 °C) collapses to 0.017–0.033 h⁻¹ for
*every* species, against point estimates of 0.284–0.450. So "the model predicts positive
growth capacity at 44 °C for every species" is a statement about the point estimates, not
about the credible parameter space. Under honest predictor uncertainty the model does not
confidently predict growth at 44 °C for anyone. **Fig 4 Panel A should be qualified
accordingly.**

**2. The phenotype is not reachable.**
- σ = 12 °C: auris grows at 44 °C **and** all relatives fail at 40 °C in 5/100 draws — but
  in **0/100** once all species are also required to retain permissive growth at 34 °C.
  Those five "successes" were draws in which the whole system collapsed, not draws that
  recovered the biology.
- σ = 6 °C: **0/100** outright (no relative ever fails at 40 °C).

**3. There is a weak, real auris tendency — far short of the observation.**
*C. auris* ranks above every relative at 44 °C in 39 % (σ=12) / 42 % (σ=6) of draws versus
25 % by chance. Directionally consistent with its marginally warmer predicted Topt, but
nothing like the observed deterministic separation (12/12 vs 2/8 at 40 °C).

## Conclusion
The falsification survives parameter uncertainty. The observed thermal divergence is not
recoverable anywhere in the credible space of sequence-predicted enzyme thermal
parameters — so Fig 4's conclusion does not rest on point estimates. The ensemble
*strengthens* the result while correcting one overstatement in Panel A.

## Caveats
- Ensemble is over Topt/Tm predictor error only in this run; a kcat-uncertainty variant
  (log₁₀ σ = 0.8) was drafted but not completed within the compute budget.
- 100 draws gives ~±3 % resolution on the reported probabilities; 0/100 means <3 %, not zero.
- Independent per-species errors are deliberately favourable to the model, as noted above.
