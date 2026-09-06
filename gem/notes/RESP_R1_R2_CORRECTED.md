# Respiration modelling (R1 + R2) — CORRECTED, with retraction

**This supersedes RESP_R1_result.md and RESP_R2_result.md. Both earlier conclusions were
wrong, for one reason: the respiration predictor contained the prediction target.**

## What went wrong
The predictor I used, `respiration_C_per_C_h` (and equivalently `R_O2_mg_cell_min`), is
**exactly** total O₂ consumed divided by the biomass integral:

    R_O2_mg_cell_min == C_tot_O2_mg_per_L / biomass_integral_cells_min_per_L   (identity to 1e-16)

The biomass integral is computed from the *fitted growth curve* (r, K, N₀). So the
"respiration rate" has growth in its denominator. Faster-growing wells automatically get a
smaller respiration value. Two consequences:

1. **The "uncoupling" was a division artifact, not biology.** Within species×temperature
   strata: total O₂ vs growth rate r = **−0.07** (essentially unrelated); biomass integral
   (the divisor) vs r = **+0.63**; O₂÷biomass vs r = **−0.68** (39/39 strata negative). The
   universal negative relationship is manufactured by the division.
2. **The R1 gate was target leakage.** Re-run with growth-independent oxygen measures
   (leave-one-species-out, per-well growth MAE, standardisation from training folds only):

   | predictor | LOSO MAE |
   |---|---|
   | temperature only | 0.1418 |
   | + total O₂ consumed | 0.1417 (no gain) |
   | + O₂ per inoculated cell·min | 0.2442 (worse) |
   | + O₂ ÷ biomass (circular) | 0.1214 ← the only "improvement", and it is invalid |

   **Oxygen adds no information about growth beyond temperature.**

A second, smaller error is also corrected: the originally reported temp+O₂ MAE of 0.109
came from standardising the O₂ feature separately within train and test folds (a scaling
leak). Done correctly that number is 0.1214 — and it is invalid anyway, per above.

## Retractions
- **R1 retracted.** "Raw O₂ adds out-of-fold growth-prediction information beyond
  temperature (0.142→0.109, −23%)" is **false**. Corrected: 0.1418 → 0.1417, no gain.
- **R1's "respiration rises while growth falls = uncoupling" retracted** as a biological
  claim. Per-biomass respiration rises at high temperature partly *because* growth (hence
  the divisor) falls.
- **R2 is void, not informative.** Maintenance-energetics models did fail badly (LOSO MAE
  0.19–0.34 vs 0.1214 across linear, threshold and exponential maintenance, multi-start),
  and the unconstrained fit ran to unphysical *negative* maintenance. But that failure
  tells us nothing about biology: the model was being asked to explain an artifact. The
  earlier conclusion that "respiration and growth uncouple in a direction energetics
  cannot represent" is withdrawn.

## What actually stands
- Total O₂ consumed is essentially uninformative about growth rate at fixed temperature and
  species, in these processed data.
- Therefore the whole respiration-constrained programme (R1–R6: oxygen constraint,
  maintenance, proteome reallocation, inverse accounting) **has no independent respiratory
  signal to work with as the data are currently derived.** It should not proceed on this
  basis.
- Nothing here touches Figure 3 (thermal phenotype/discordance) or Figure 4 (the
  counterfactual falsification), which use growth data only and are unaffected.

## What would be needed to revisit respiration properly
A specific respiration rate that is **not** normalised by a growth-derived quantity — e.g.
the O₂ depletion slope over a defined window divided by an independently measured cell
density at that time (direct counts/OD at the window), or a per-well O₂ rate normalised by
inoculum. That requires going back to the raw O₂ traces and the plate records, not the
derived table.

## Process note
This was found only by redoing the analysis from first principles and checking how each
variable is constructed. The lesson is the one ChatGPT flagged in advance — "do not use a
quantity derived from growth to predict that same growth" — and I violated it by assuming
`respiration_C_per_C_h` was an independent biomass-normalised measurement.
