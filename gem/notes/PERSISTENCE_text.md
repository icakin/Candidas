# Persistence decomposition — positive main-figure result (revised, rigorous)

## FIGURE TITLE
C. auris differs by retaining detectable growth across the febrile boundary — not in
the conditional growth rate or respiratory cost of a competent cell.

## PRIMARY RESULT (the defensible statistical claim)
At 40 C the fraction of wells with detectable growth was 0.98 (59/60) for C. auris
versus 0.33 (5/15) C. haemulonii, 0.00 (0/10) C. duobushaemulonii and 0.27 (4/15)
C. parapsilosis. Predefined contrast, C. auris vs the pooled sampled relatives:
risk difference +0.38 (Fisher exact p = 6x10^-7); at 42 C, +0.53 (p = 3x10^-8).
Because there is a single auris species, this is a contrast conditional on the
sampled taxa, and wells are not independent species replicates.
A pooled-slope binomial logistic places the persistence midpoint (T50) at >44 C for
C. auris (right-censored, above the tested range) versus ~40 C for the others; T50 is
reported only secondarily because a shared slope cannot capture the sharp non-auris
drop-off — the model-free 40/42 C fractions above are the primary evidence.

## CONDITIONAL PERFORMANCE (what does NOT differ)
Among growth-positive wells in the permissive range (persistence >= 0.8), growth rate
broadly overlaps across species (~0.5-0.8 h^-1, 28-34 C). Respiration-to-growth rises
from 36 to 40 C in every surviving isolate (shared direction): median fold-change
1.4x across 12 C. auris isolates, and ~3.5x in the single surviving C. haemulonii
(1724) and C. parapsilosis (2052) isolates (n=1 each — descriptive only, no species
mean implied). So the febrile respiratory cost is general, not unique to one species;
the surviving non-auris isolates may pay more, but this cannot be generalized.

## FIGURE LEGEND
**Figure X. Thermal divergence among these species is a growth-persistence boundary.**
(A) Fraction of inoculated wells with detectable growth vs temperature: per-species
means with 95% Wilson intervals, faint per-isolate proportions behind, and a
pooled-slope logistic guide; shaded band = 40 C, where C. auris (0.98) vastly exceeds
the sampled relatives (0-0.33; risk difference +0.38, Fisher p<10^-6). C. auris T50 is
right-censored (>44 C). (B) Growth rate among growth-positive wells only; solid where
persistence >= 0.8 (permissive range, where the comparison is unbiased), open faded
points = survivor-only subsets at higher temperature (not connected). Conditional
growth broadly overlaps in the permissive range. (C) Isolate-level respiration/growth
fold-change from 36 to 40 C: every surviving isolate increases (dashed line = no
change); C. auris (n=12) modestly (median 1.4x), the lone surviving non-auris isolates
more (n=1 each). n = isolates per species; non-auris species have 2-3 isolates.

## HOW IT PAIRS WITH THE etcGEM (supplement)
This decomposition is the positive result; the etcGEM is the supplementary
sufficiency test showing that genome-parameterized metabolic/transport/respiratory
capacity is not sufficient to reproduce this persistence boundary.

## HONEST CAVEATS
- Non-auris species rest on 2-3 isolates; the robust claim is C. auris vs the sampled
  relatives, not a precise ranking among the three non-auris species.
- "Detectable growth" is a finite-assay-window observation, not a viability assay.
- Confirm that QC/hardware-failed wells are coded missing, not as no-growth.
- Conditional and R/G values at the highest temperatures reflect few survivors.
