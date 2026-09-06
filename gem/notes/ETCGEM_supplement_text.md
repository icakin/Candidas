# etcGEM — supplementary text (boundary/falsification analysis)

## MAIN-TEXT PARAGRAPH (one paragraph; points to Supplementary Figure SX)
To test whether the striking mismatch between phylogenomic proximity and thermal
phenotype is explained by metabolic-enzyme sequence, we built an enzyme- and
temperature-constrained genome-scale model (etcGEM) for each species, parameterized
from sequence-predicted turnover numbers, enzyme molecular weights, and per-enzyme
optimal and melting temperatures, with a shared protein pool and shared thermal
hyperparameters calibrated on C. auris and frozen (Supplementary Fig. SX; Methods).
Under leave-one-species-out cross-validation the model did not transfer: for every
held-out species a temperature-specific cross-species baseline predicted growth as
well or better than the genome-parameterized model (Fig. SX-A), genome-specific
thermal parameters added no out-of-fold advantage over an identical shared thermal
profile, both models were worse than a trivial median-optimum predictor for the
thermal optimum (Fig. SX-C), and oxygen flux was non-identifiable so respiratory
partitioning could not be predicted (Fig. SX-D). The model reproduced the
calibration species and a generic mesophilic response but systematically
over-predicted growth beyond each species' observed loss-of-growth boundary
(Fig. SX-B). We therefore did not pursue the etcGEM as a predictive model.

## MAXIMUM DEFENSIBLE SENTENCE (use verbatim if you want the single strongest line)
An auris-calibrated, genome-parameterized etcGEM failed to outperform a
temperature-specific cross-species baseline for any held-out species and left
respiratory flux non-identifiable; the observed thermal divergence was not
recoverable from metabolic-network structure and current sequence-based
enzyme-property predictions, implicating factors outside — or below the resolution
of — the present model, including condition-specific regulation and proteome
allocation, membrane and cell-wall physiology, stress responses, or localized
sequence effects.

## SUPPLEMENTARY FIGURE LEGEND
**Supplementary Figure SX. An auris-calibrated genome-parameterized etcGEM does not
transfer to predict interspecific thermal performance.**
(A) Out-of-fold growth mean absolute error (pre-boundary temperatures, defined
independently of growth as temperatures below the first point where <80% of wells
remained growth-positive) for the etcGEM versus a temperature-specific baseline
(held-out species predicted from the mean of the training species at each
temperature); C. auris is calibration (in-sample) and shown for reference only. For
every transfer species the baseline is as good or better. (B) Out-of-fold predicted
(line) and measured (points) thermal-performance curves; grey marks the observed
loss-of-growth region, where the capacity model over-predicts growth. (C) Thermal-
optimum and whole-curve error relative to a trivial baseline (dashed line = 1);
genome-specific and identical-profile-null models are indistinguishable, and both
exceed 1 (worse than trivial) for the thermal optimum. (D) Oxygen-flux variability
analysis at 90% of optimal growth (C. auris): the feasible qO2/mu range spans 3.5-
19.5 (5.6x), i.e. respiration is non-identifiable owing to respiro-fermentative
freedom. Measured zeros are zero-inclusive realized-growth summaries, not proof of
cell death; steady-state flux-balance models represent feasible growth capacity and
contain no damage, repair or loss-of-growth process.

## RESCUE ATTEMPTS (preregistered, both null — document exhaustion)
Two principled rescues were tried before abandoning genome-only prediction:
(1) Per-species pool calibrated at a single non-stressful reference (30 C) and
prediction of the NORMALIZED thermal response — this removes baseline growth-scale
differences. Out-of-fold it did not identify C. auris as the high-temperature
survivor: predicted normalized growth at 42 C was 0.80 (auris), 0.78 (para), 0.77
(hae), 0.73 (duo) — indistinguishable — while only auris grew (measured 0.94; the
others 0). (2) A localized essential-enzyme reserve test: the melting temperatures
of the ACTIVE, flux-carrying enzymes at 42 C were not higher in C. auris (median
53.8 C) than in the clade (para 54.3 C, duo 53.7 C, hae 53.4 C). Both rescues being
null indicates the exhaustion of defensible genome-only approaches with current
sequence-based predictors, not a limitation of a single model implementation.

## TRANSPORT/RESPIRATION CORRECTNESS NOTE
An audit showed that in these reconstructions ~89% of active internal flux at high
temperature runs through reactions with no gene-protein rule — including glucose and
amino-acid uptake, oxygen/phosphate/CO2 transport, and cytochrome c oxidase — so the
membrane-transport and respiratory layer is thermally unconstrained. Adding a
sequence-derived thermal cost to these reactions lowered predicted high-temperature
growth (e.g. C. auris mu at 44 C 0.45 -> 0.20) but did NOT create species
discrimination, because these transporter/complex sequences are conserved across the
clade like the metabolic enzymes. Thus the discriminating signal is absent at every
layer the model can represent (metabolic enzymes, proteome composition, active-enzyme
Tm, and membrane transport/respiration).

## FRAMING (the diagnostic conclusion — use this)
This analysis does not show that thermal persistence is fundamentally unpredictable;
it shows precisely which information is insufficient. Thermal persistence is
estimable, but not from conserved metabolic-enzyme sequences alone: it would require
either a genome-wide predictor trained across many fungi with known thermal limits,
or direct measurement of temperature-resolved regulation, proteome allocation and
stress physiology. The etcGEM's value here is diagnostic — it localizes the missing
information to processes outside the metabolic-enzyme layer.

## METHODS SENTENCE (add to methods)
etcGEM performance was evaluated by leave-one-species-out cross-validation against
(i) a temperature-specific cross-species baseline, (ii) an identical shared thermal
profile (global-median Topt/Tm), and (iii) trivial median-optimum and mean-curve
baselines; shared hyperparameters were fitted only within training folds and never
per species. Oxygen exchange was left unconstrained and respiration was assessed by
flux-variability analysis at 90% of optimal growth without fitting to measured
respiration.

## WHAT NOT TO CLAIM (internal note)
- Do NOT say the divergence "lies in" regulation/membranes — we excluded only what
  this model represents and what current predictors resolve. Use "implicating
  factors outside or below the resolution of the model."
- Do NOT report the r=0.67 / auris-inclusive statistic as performance.
- Frame as a preregistered boundary/falsification analysis, supplement + one
  paragraph, not a capstone.
