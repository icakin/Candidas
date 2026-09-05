# What sets the thermal limit in emerging *Candida* pathogens?

*Draft case for support. Every figure quoted below is sourced in the notes at the end.*

## Background

*Candida auris* emerged through the nearly simultaneous appearance of genetically distinct
clades on three continents and is classified as a WHO critical-priority fungal pathogen.
Among its most striking traits is its ability to grow at temperatures at which its closest
relatives lose detectable growth. Our preliminary work quantifies that divergence, measures
quantifies its inferred respiratory cost relative to growth, and establishes by three
complementary routes that the
divergence is not explained by coarse gene content, by sequence-predicted metabolic-enzyme
properties, or by current metabolic models. The result is a
well-posed question about which cellular process actually sets the thermal limit in these
organisms, and answering it is the aim of this proposal.

## Preliminary result 1: growth and respiration decouple, and carbon efficiency peaks below host temperature

We assayed 20 isolates across seven lineages by dissolved-oxygen respirometry at twelve
temperatures from 22 to 44 °C, five replicate wells per isolate per temperature: the four
*C. auris* clades (three isolates each), *C. haemulonii* (three), *C. duobushaemulonii*
(two) and *C. parapsilosis* (three). Fitting the *C. auris* clades and *C. parapsilosis*
hierarchically in a Bayesian framework, with isolates nested within clades, we find that
growth rate rises to a clade-level optimum near 33 to 36 °C and then declines, while per-cell respiration shows no detectable
turnover between 22 and 44 °C: it increases monotonically across the measured range, with leave-one-out cross-validation
favouring a monotonic Arrhenius description over a unimodal one. The two processes are
governed by different thermal sensitivities. The apparent activation energy of growth is
0.83 to 1.14 eV against 0.30 to 0.52 eV for respiration, with the difference credibly
greater than zero in all five fitted groups.

Carbon-use efficiency is therefore itself unimodal, peaking early and then collapsing. The
position of that economic optimum is identifiable without any assumption that sets absolute
scale: writing G = r·q and R = c·K, the ratio c/q drops out of the optimum entirely, so it
follows from the fitted growth rate and volumetric consumption rate alone. This requires no
value for the per-cell carbon quota q, but it does require q not to vary with temperature,
an assumption Aim 2 tests directly; if cells follow the usual temperature-size rule and are
smaller when warm, the optimum moves cooler still. On that scale-free basis the optimum sits
at 26.1 to 31.5 °C, between 5.5 and 10.9 °C below body temperature, with posterior
probability greater than 0.9999 of lying below 37 °C in every fitted group.
Host-relevant temperatures therefore lie on the declining side of the in vitro
carbon-efficiency curve for every species we measured.

## Preliminary result 2: a fever imposes a quantifiable, taxon-specific carbon cost

At host-relevant temperatures in vitro, every taxon lies beyond its carbon-efficiency
optimum, and the size of that penalty differs. At 40 °C the respiratory cost per unit
growth reaches 1.34-fold for *C. auris* Clade IV but 3.13-fold for *C. parapsilosis*. These
costs are each normalised to the taxon's own cheapest temperature, so the absolute value of
the adopted per-cell carbon quota, including the four-fold larger quota assigned to
*C. parapsilosis*, cancels before any cross-taxon comparison is made; the same is true of
the 37 to 40 °C changes and the pairwise contrasts below. They depend only on the quota
being constant within a taxon across temperature, the assumption Aim 2 tests. The quota's
absolute value enters the level of per-cell respiration and carbon-use efficiency, not any
comparison reported in this section. Across
the clinically direct step from 37 to 40 °C, *C. parapsilosis* retains only 58% of its
growth rate while its cost per unit growth rises 1.95-fold; Clade IV retains 89% at a cost
increase of 1.25-fold. All four *C. auris* clades credibly pay less at fever than
*C. parapsilosis*, with cost ratios of 0.43 to 0.72, and that contrast survives every
treatment of the starting cell number. We do not report a finer ranking among the *C. auris*
clades because the data do not support one.

Across the five fitted groups, the cost of a fever falls as the environmental growth optimum rises.
With five fitted groups this is suggestive rather than established, and the sign rather than the
magnitude is the interpretable quantity. The pattern nevertheless suggests that
environmental thermal selection contributes to tolerance of the mammalian host.

## Preliminary result 3: the thermotolerance divergence is robust at isolate level

The *haemulonii*-complex relatives, assayed on the same plates and at the same twelve
temperatures, provide the comparison that matters for thermotolerance itself.
In the prespecified isolate-level contrast at 40 °C, all twelve *C. auris* isolates retain
detectable growth while only three of eight isolates of the three closest relatives do (risk
difference 0.63, 95% CI 0.22 to 0.86; Fisher *P* = 0.004). At 44 °C, the ceiling of the
assay, the separation persists at 10 of 12 against 1 of 8. Isolate-level thermal performance
curves separate *C. auris* from the relatives with a common-language effect size of 0.95:
that is, a randomly chosen *C. auris* isolate outperforms a randomly chosen relative isolate
95% of the time. The divergence is therefore not an artefact of pooling or of one tolerant
strain.

## Preliminary result 4: the genome does not predict the phenotype

Three complementary analyses converge here.

First, within *C. auris* the four clades differ up to two-fold in peak growth rate and
measurably in fever cost, yet median whole-proteome amino-acid identity across clades I to
IV is 99.4 to 100%. Broad coding-sequence divergence is unlikely to explain the physiology,
although a small number of causal coding or regulatory variants remains possible.

Second, across ten annotated proteomes we built a phylogeny from 2,931 single-copy
orthologs, a 600-gene supermatrix with full support at every node, and asked whether gene
content tracks thermal phenotype. Copy number of the canonical thermal machinery, HSP90,
small heat-shock proteins, the trehalose pathway, calcineurin, fatty-acid desaturase and
ergosterol biosynthesis, is invariant across all ten genomes. The single exception,
alternative oxidase, separates families rather than thermal phenotypes. Coarse
proteome-composition indices appear to correlate with thermal optimum until the
pseudoreplication is corrected: collapsing the four near-identical *C. auris* clades to one
phylogenetic point removes the association, and the index that appears predictive separates
the two families completely, independent of thermal optimum. Most tellingly, those four
clades differ by 2.7 °C in carbon-use-efficiency optimum while being 25 times less divergent
than the nearest species split. Thermal performance is evolutionarily labile and is not
predicted by coarse gene content.

Third, the available transcriptional evidence does not resolve what is different. In the
public clade-resolved RNA-seq dataset we identified, generated under a single controlled
condition, central-metabolic enzyme genes
are not coordinately lower in the slower-growing clade (median log2 fold-change −0.005
against +0.031 for the genome background, Wilcoxon *P* = 0.93). Being compositional, from a
single isolate at one temperature, that dataset does not support a coordinated relative
down-allocation of metabolic transcripts, but cannot resolve differences in absolute
output.

## Preliminary result 5: metabolic modelling cannot reproduce the divergence either

We constructed enzyme- and temperature-constrained genome-scale models for four species from
their respective genomes, with per-enzyme catalytic and thermal parameters predicted from
protein sequence, calibrating three shared parameters on *C. auris* alone and then
freezing them. We tested ten model formulations and constraints spanning enzyme melting
temperature, catalytic optimum, proteome capacity, maintenance energetics and oxygen
availability. Only two temperature-sloped formulations, maintenance demand and proteome
reallocation, could in principle reproduce both sides of the pattern by preserving
*C. auris* growth while driving its relatives below detection, and both required values far
outside empirical bounds: a maintenance slope of β = 8.1 per °C, compared with the measured
slope of 0.241 per °C, and
loss of 87% of the metabolic proteome by 40 °C against reported heat-shock reallocations of
10 to 30%. Every temperature-independent formulation failed structurally, removing as much
capacity at permissive temperatures as at restrictive ones.

In the baseline thermal-enzyme formulation, reproducing the observed pattern required
approximately 33 °C of interspecies separation in enzyme melting temperature, against
approximately 0.4 °C in deduplicated paired sequence predictions and 1.6 °C in a measured
whole-proteome comparison of two thermally divergent congeneric yeasts; after deduplication
no relative showed a resolvable difference from *C. auris* in predicted catalytic optimum at
all. Holding each species out and refitting every parameter on the other three leaves the
failure intact: the model still predicts 79 to 85% of its own peak growth at temperatures
where the observed thermal-performance curves of the relatives fall to or below detection.

## The puzzle this exposes

The models predict cellular failure, defined as the temperature at which maximum feasible
growth falls below the experimental detection threshold, at 52.7 to 54.5 °C across the four
species. The relatives begin losing detectable growth around 40 °C and *C. auris* around
45 °C. A missing cellular failure process therefore becomes limiting 9 to 14 °C below the
model-predicted limit, before the bulk metabolic-enzyme unfolding regime the model encodes.

This gap is not an artefact of the model's known defects. The three we identified all
inflate capacity, so a sceptic could argue the predicted limit is spuriously high; we
therefore corrected all three together, constraining the respiratory quotient by tightening
the amino-acid pool, blocking the energy-generating cycle, and setting the enzyme-unfolding
width to its physically plausible value with the growth scale refitted to hold *C. auris*
peak growth. The predicted limit moved by less than one degree, to 53.4 to 54.1 °C, leaving
the gap at +9 °C for *C. auris* and +13 to 14 °C for the relatives. The reason is
structural: temperature enters the model only as a capacity throttle keyed to enzyme
melting temperature, so the predicted limit tracks the median enzyme melting point (53 to
54 °C) regardless of any metabolic correction. Metabolism was never what set the limit,
which is why the limit must be measured directly. Identifying the process that does set it
is the central aim of this proposal.

## Aims

**Aim 1. Establish the uncensored thermal limit and the temporal order of failure.**
*C. auris* was still growing at the ceiling of our assay, so its true limit is unknown, and
that number provides the reference point for all subsequent mechanistic comparisons. We will
extend the thermal range across the isolate panel, then determine the order in which
membrane integrity, proteostatic capacity and mitochondrial energy conversion fail in the
same isolates at the same temperatures. Establishing which system fails first will identify
the most plausible upstream event and prioritise the causal perturbations in Aim 2.

**Aim 2. Test causality, and quantify regulatory and bioenergetic allocation.** Correlated
failure is not causal failure, so we will perturb each system independently and ask which
perturbation moves the thermal limit, giving priority to whichever system Aim 1 identifies
as failing earliest and most selectively between *C. auris* and its relatives. In the same experiments we will resolve the
regulatory-versus-bioenergetic question that our own analysis identifies as unresolved:
absolute, spike-in-calibrated proteomics under these thermal conditions, together with
direct measurement of the carbon quota across temperature, which currently enters the
absolute carbon scale as an assumption. The four *C. auris* clades are the natural
experiment for this, differing measurably in thermal phenotype while being genomically
near-identical, so any allocation difference tracking that phenotype is unlikely to be
confounded by gene content. The mitochondrial arm is further supported by preliminary
evidence of *C. auris*-specific gene gains and losses enriched in mitochondrial assembly
functions, together with the uncoupling of respiration from growth documented above.

**Aim 3. Discover adaptive routes without presupposing them.** Replicated experimental
evolution of the relatives under progressively increasing temperature provides a
hypothesis-independent arm. The evolved lines will test whether increased thermotolerance
can arise without coordinated shifts in enzyme thermal properties, and convergence among
independent lines will nominate a mechanism. Nomination is not demonstration: evolution need
not replay the route *C. auris* itself took, so candidate changes will be validated by
direct reconstruction or perturbation before any causal claim is made.

## Significance

Thermotolerance may help environmental fungi cross the mammalian thermal barrier and remain
viable within febrile hosts. Together with the exceptional environmental persistence of
*C. auris*, that capacity magnifies its threat in healthcare settings. Our preliminary work
rules out two obvious explanations within the resolution of our data: broad differences in
thermal-response gene content, and sequence-predicted thermostability of metabolic enzymes.
What remains is a specific 13 °C window and a defined set of cellular failure processes. By combining
comparative physiology, causal perturbation, absolute proteomics and experimental evolution,
this project will identify the cellular process that sets that limit, and provide both a
framework and a candidate cellular vulnerability for understanding the emergence of
thermotolerant fungal pathogens more generally.

---

## Sources for every number above

*Internal audit table. In the submitted proposal, convert the external benchmarks
(Lockhart et al. on emergence, the WHO priority list, Walunjkar et al. on the measured
congeneric benchmark) into formal citations and drop the file paths.*

| Claim | Source |
|---|---|
| 20 isolates, 7 lineages, 12 temperatures, 5 wells each | `results/tables/Oxygen_Data_Filtered.csv` |
| q is a fixed per-taxon constant (1320 fg C/cell; 5280 for *C. parapsilosis*), so G/q = r exactly and panel 1a's h⁻¹ label is exact | `results/tables/derived_N0_R_results_with_carbon.csv`; `scripts/14_main_figures.R` lines 388–396 |
| growth optimum 33–36 °C; respiration monotonic; LOO favours Arrhenius | Fig. 1a,b |
| E_G 0.83–1.14 eV, E_R 0.30–0.52 eV, difference credible in all five | Fig. 1c |
| scale-free CUE optimum 26.1–31.5 °C, 5.5–10.9 °C below 37, P > 0.9999; assumes T-independent q | Fig. 1d,e and the c/q identity |
| fever cost 1.34× (Clade IV) vs 3.13× (*C. parapsilosis*) at 40 °C | Fig. 2a |
| 37→40 °C: growth retained 0.58 vs 0.89; cost ×1.95 vs ×1.25 | Fig. 2c |
| all four clades pay less than *C. parapsilosis*, ratios 0.43–0.72 | Fig. 2d |
| fever cost falls as growth optimum rises (suggestive, n = 5) | Fig. 2b |
| 12/12 vs 3/8 isolates at **40 °C** (prespecified); RD 0.63 [0.22, 0.86]; Fisher P = 0.0036; CLES 0.95 | Fig. 3; `gem/FIG3_caption.md` |
| 10/12 vs 1/8 isolates at 44 °C (assay ceiling) | `gem/counts_at_44.json` |
| whole-proteome identity 99.4–100% across clades I–IV | manuscript, coding-sequence section |
| GSE165762 median log2FC −0.005 vs +0.031, Wilcoxon P = 0.93 | manuscript, same section |
| 2,931 orthologs; 600-gene supermatrix; all nodes 1.00; 25× divergence | `phylo/phylogenomics_results.md` |
| thermal-machinery copy number invariant; AOX the family-tracking exception | same |
| composition correlation collapses on phylogenetic collapse | same |
| ten formulations; β 8.1 vs 0.241; 87% vs 10–30% | `gem/FIG4_etcgem_caption.md` |
| 33 °C required; 0.4 °C deduplicated; 1.6 °C measured benchmark | Fig. 4B,C; `gem/paired_dedup_audit.py` |
| no resolvable Topt difference after deduplication | `gem/paired_dedup_audit.py` |
| held-out 79–85% of peak at observed failure temperature | `gem/loocv_nested.py` |
| model-predicted limit 52.7–54.5 °C; 53.4–54.1 °C after correcting all three defects | `gem/thermal_sensitivity.py`; `gem/gap_robustness.py` |
| spike-in proteomics + carbon quota as the required next experiment | manuscript, Conclusion |

## Two things to settle before submission

**Nomenclature.** Recent revisions place the *auris*/*haemulonii* group in *Candidozyma*
(2024) and *C. parapsilosis* in *Lodderomyces* (2026). The manuscript retains *Candida*
throughout for continuity with the clinical literature. State this once, or a referee will
raise it.

**Personnel and facilities.** Aims 1 to 3 are all wet-lab. Name the people, the platforms
and the collaborators, particularly for the proteomics and the evolution experiment.
