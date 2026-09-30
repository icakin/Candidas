# etcGEM counterfactual falsification figure — title, caption, verdict (4-panel, strengthened)

## TITLE
Within the etcGEM, predicted metabolic-enzyme thermal properties are insufficient to
reproduce the observed thermal divergence

## THE QUESTION (two-sided sufficiency test)
The baseline model already lets every species grow at 44 °C, so "how much must C. auris
rise" is meaningless. Instead: what minimum interspecies difference in a model parameter
reproduces BOTH sides of the phenotype — C. auris still growing at 44 °C AND the median
relative failing by 40 °C — while permissive (≤34 °C) growth stays above a floor
(0.7 × baseline µ at 34 °C)? Tested for UNIFORM proteome-wide shifts and for TARGETED
(sparse) shifts found by a dynamic search.

## CAPTION
**Figure X | Neither uniform nor targeted shifts in predicted metabolic-enzyme
thermostability reproduce the observed thermal divergence.** (A) The calibrated etcGEM
predicts detectable growth at 44 °C for all four species (bars = model-predicted µ);
observed isolate-level detection at 44 °C is annotated (n/N: C. auris 10/12,
C. haemulonii 1/3, C. duobushaemulonii 0/2, C. parapsilosis 0/3). Bars and n/N are
different quantities (model prediction vs observed isolate counts). (B) Imposing a
uniform Tm shift on all enzymes of an illustrative relative (C. duobushaemulonii; the
other relatives behave the same): 40 °C growth falls below detection only near a −32 °C
shift, but permissive 34 °C growth drops below its floor far earlier — no shift satisfies
both. (C) Interspecies separation predicted from sequence versus the minimum required by
the model (dumbbells; fold-gap annotated); the requirement is stable across detection
thresholds 0.03–0.10 h⁻¹. (D) Dynamic sparse test: predicted 40 °C growth as bottleneck
enzymes are shifted −8 °C in Tm, added best-first with marginal effects recomputed after
every addition (so bottlenecks exposed by flux rerouting are caught); growth plateaus far
above detection for every relative. A beam search over all pairs and triples of the
strongest candidates likewise never approaches detection (lowest pair/triple µ at 40 °C
≈ 0.30–0.55 vs 0.05 detection). Values are model-implied minima, not biological Tm/Topt
estimates; a preregistered sufficiency/falsification test, not causal inference; 44 °C is
assay-censored (10/12 C. auris isolates still growing there).

## MECHANISM FEASIBILITY (for caption/supplement)
| Mechanism (metabolic-enzyme layer) | Required interspecies change | 40 °C failure achieved? | 34 °C growth preserved? | Both satisfied? |
|---|---|---|---|---|
| Uniform Tm shift | ≈32 °C (≈62× predicted) | only at ≈32 °C | No | No |
| Uniform Topt shift | ≈16 °C (≈16× predicted) | only at ≈16 °C | No | No |
| Targeted/sparse (dynamic, ≤8 °C, ≤30 enzymes + beam pairs/triples) | not reached | No | yes | No |
| Flat protein-pool capacity | no feasible value | only by also killing 34 °C | No | No |

## WHY (rerouting, documented)
Shifting the single strongest bottleneck enzyme's Tm barely lowers 40 °C growth
(e.g. C. duobushaemulonii: 0.435 → 0.40 h⁻¹); flux redistributes (mitochondrial
ATP-synthase and alternative dehydrogenase reactions adjust), so no small enzyme set
forms a controlling bottleneck. Because the growth optimum itself is unique, alternative
optimal flux distributions do not change these µ values.

## GATE VERDICT (passes on uniform and dynamic-targeted)
Required ≫ predicted (≈16–62×), consistent across all relatives, both parameters, and
all detection thresholds; the uniform change that could work is non-physical (destroys
permissive growth); and no dynamically-searched targeted set (single, greedy-reranked, or
beam pairs/triples, ≤8 °C) reproduces the boundary while preserving permissive growth.
This is a strong search, not a formal impossibility proof (it does not exhaustively
enumerate all enzyme combinations or arbitrarily large mixed Tm+Topt changes), but the
convergent negative across uniform and targeted regimes supports the title.

## PLACEMENT
Main-figure-quality (Fig 4) if the Results explicitly introduce bulk metabolic-enzyme
thermostability as a tested biological hypothesis; otherwise a substantially stronger
etcGEM supplement than the prior transfer/LOOCV figure (can replace or accompany
FIG_SUPP_etcGEM).

## NEXT (gated secondary, only if desired)
Respiration-constrained inverse model as a strict go/no-go: raw O₂ only (never O₂/growth
or CUE), compared against a temperature-plus-O₂ regression; promote only if it beats that
baseline. Otherwise stop modelling — this falsification figure is the computational result.
