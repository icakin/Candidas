# PHASE 3 VERDICT — predictive test, FROZEN
# Predictors: ORIGINAL HAND-TRIMMED fits restricted to 22-38 °C. Taxon identity never used.
# Outcomes  : frozen v4 states at 40-44 °C and the interval-censored upper thermal limit.
# Validation: leave-one-ISOLATE-out, complete isolate withheld, 20 isolates.
# Ambiguity : both bounds reported; ambiguous outcomes never reallocated.

## PREDICTORS (per isolate, from 22-38 °C only)
  growth      : r at 38 °C ; growth optimum from a quadratic on log r over 22-38 °C
  metabolism  : log(K/r) at 38 °C ; slope of log(K/r) over 34-38 °C
  M1 growth-only | M2 metabolism-only | M3 combined, with the best single predictor of each
  type selected INSIDE each training fold (never on the held-out isolate).

## RESULT — metabolic predictors DO NOT improve held-out prediction

Δ held-out log predictive score against growth-only (positive = metabolism helps):

  outcome                                   M2 - M1     M3 - M1
  confirmed growth at 40 °C, lower bound    -0.1726     -0.0564
  confirmed growth at 40 °C, upper bound    -0.0126     -0.0637
  confirmed growth at 42 °C, lower bound    +0.0488     +0.0754
  confirmed growth at 42 °C, upper bound    -0.0533     +0.0042
  confirmed growth at 44 °C, upper bound    -0.0478     -0.0128

Metabolism improves prediction in ONE of five outcomes (42 °C, lower bound) and is worse or
indistinguishable in the other four. Growth-only reaches AUC 1.000 for confirmed growth at
40 °C under the lower bound; no metabolic model beats it anywhere except at 42 °C lower.

Upper thermal limit, mean absolute error on the 18 uncensored isolates:
  growth-only 1.10 °C | metabolism-only 1.15 °C | combined 1.17 °C
Metabolism is 0.06 °C WORSE, combined 0.07 °C worse. No model predicted either censored
isolate (Clade1_2070, Hae_1724) at >= 44 °C: both were under-predicted at ~43.2-43.8 °C.

## NAMED ISOLATES (as required by the plan)
  Hae_1724     observed >44   growth 43.4   metabolism 43.5   combined 43.8
               — every model under-predicts the outlier relative; none identifies it as
                 exceptional from 22-38 °C data alone.
  para_2052    observed 40    growth 41.0   metabolism 41.0   combined 41.1
  Clade2_2073  observed 38    growth 40.5   metabolism 38.0   combined 41.0
               — the one C. auris isolate that fails early. ONLY the metabolism model
                 predicts it correctly; the growth model over-predicts by 2.5 °C.
  Clade4_2079  observed 42    growth 45.3   metabolism 43.8   combined 45.1
               — the largest single error in the set, over-predicted by 3.3 °C by growth.

## A REDUNDANCY FOUND, NOT A RESULT
"Entry into respiration without confirmed growth at 40 °C" is EXACTLY the complement of
"confirmed growth at 40 °C, upper bound" — every isolate has either >= 3 growth-or-ambiguous
wells or >= 3 respiration-only wells, with none in between. It is not an independent outcome
and is not reported as a separate test.

## CONCLUSION
Metabolic predictors measured below fever do NOT improve held-out prediction of failure at
40-44 °C beyond ordinary growth measurements. Per the prespecified stopping rule, this line
is closed: no other temperatures, transformations or predictor combinations were searched.

The one suggestive exception is recorded without being pursued: at 42 °C under the lower
bound, metabolism adds +0.049 and the combined model +0.075 log predictive score, and the
metabolic model alone correctly predicts the early failure of Clade2_2073, the single
C. auris isolate that behaves like a relative. With 20 isolates this is not evidence of an
early-warning signal; it is a hypothesis for a larger panel.
