# Figure 3 — consistent febrile growth in C. auris carries an increased respiratory cost

Reading path: phylogenetic relationship -> species-level 40 °C count (12/12, 1/3,
0/2, 1/3) -> isolate-level well counts -> conditional R/G fold change.

## Panels
- **a** Phylogenomic context. Tree from 600 single-copy orthologs / 289,787 aa,
  ML (LG+Γ), rooted on the Debaryomycetaceae. The four *C. auris* clade reference
  genomes are collapsed to one compact tip labelled "4 clade references collapsed"
  (their internal spread is only ≤0.013 subs/site, so the tip is deliberately NOT
  drawn spanning the 12 isolate rows). Grey tips = reference genomes with no
  phenotype in this study. Full ten-genome tree with supports: supplement.
- **b** Persistence. Wells with detectable growth, of 5 per isolate, at 36–44 °C.
  Raw counts — no model. Neutral grey ramp for 1–5, pale red for 0 (blue is
  reserved as C. auris' identity colour).
- **c** (fused into panel b as a numeric column) R/G fold change, 37→40 °C.
  **Point-estimate fold changes — ratios of marginal posterior medians, NOT
  posterior medians of the ratio (no joint draws). Caption wording:** "All 14
  isolates with detectable growth at 40 °C had higher model-estimated R/G at 40
  than at 37 °C (point-estimate fold changes, 1.21–1.94×)." Hierarchical
  estimates with credible intervals: Fig. 2c.  Fold change
  shown ONLY for isolates with detectable growth at 40 °C. Isolates without it are
  marked as such rather than assigned an extrapolated ratio.
  NOT a second gate — it has no pass/fail boundary and cannot compare species
  with one surviving comparator isolate each.

## Numbers
- 7 taxa, **20 isolates** (23 in the dataset; 3 *C. glabrata* held out via
  EXCLUDE_GROUPS pending a re-run — assay noise, not a fitting problem).
- Persistence at 40 °C: **14/20 isolates**, 12 of them *C. auris*.
  *C. auris* 12/12 · *C. haemulonii* 1/3 (1724) · *C. parapsilosis* 1/3 (2052) ·
  *C. duobushaemulonii* 0/2.  At 42 °C: 13/20 (12 *C. auris* + 1724). At 44 °C: 11/20.
- Conditional cost among those 14: R/G higher at 40 than 37 in **14/14**, with
  non-overlapping marginal 95 % CrIs in 14/14 (strict, correlation-ignoring test).
  *C. auris* 1.21–1.76× · *C. haemulonii* 1724 1.51× · *C. parapsilosis* 2052 1.94×.
  **Cost does not separate the species** — 1724 sits inside the *C. auris* range.
- The 38→40 °C step is where the separation appears; at 36–38 °C almost every
  isolate grows. 40 °C is therefore not a cherry-picked threshold.

## Wording constraints
- "fever threshold: 40 °C" marks the 37→40 °C transition used throughout. Do NOT
  call 40–44 °C "the febrile range" — 42–44 °C is supra-febrile stress.
- "detectable growth", never survival/death — no viability assay was done.
- Cost is **conditional** on growth being estimable at 40 °C.
- Do NOT call the cost stage a "gate", and do NOT claim cost separates species.
- Persistence is a difference in **prevalence/consistency across sampled isolates**,
  not exclusive possession of thermotolerance — 1724 and 2052 are counterexamples.
- Target sentence: "Among the 20 analysed isolates, *C. auris* alone maintained
  detectable growth at 40 °C across every isolate, whereas all isolates that grew
  at 40 °C had higher estimated R/G than at 37 °C, separating consistent
  febrile-temperature persistence from its conditional respiratory cost."
- Discussion-level: "passing the thermal filter does not make fever free."
- No family-level boundary: *C. parapsilosis* is in a different family yet
  phenocopies the haemulonii complex. The pattern maps to the *C. auris* lineage.
- Do not generalise to "the haemulonii complex" — *C. pseudohaemulonii* was not
  phenotyped.

## Consequence for the species-level numbers in Fig 2
Species-level fever cost pools isolates that never grow at 40 °C, so those taxa's
values rest partly on model extrapolation:
- *C. haemulonii* 1.95× pools 1724 (1.51, measured) with 1768 (2.12) and 1769
  (2.26), neither of which grows at 40 °C.
- *C. duobushaemulonii* 6.13× is **entirely** extrapolated — 0/10 wells at 40 °C.
- *C. parapsilosis* 1.92× rests on 2052 alone among the three isolates.
This must be stated in the Fig 2 legend / Results.
