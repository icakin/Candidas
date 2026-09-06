# etcGEM growth-TPC pipeline (Candidas) — v3

> **v3.** PI directive: build an etcGEM using the project's genomic data to
> estimate growth. Objective is now **species growth TPCs**, not CUE. Task 0 is
> RESOLVED (see below) and clears the enzyme-pool worry. Design updated to a
> leave-one-species-out prediction test so "genome estimates growth" is a real,
> non-circular claim.

## Directive & honest framing
etcGEM (the etcYeast method) uses the genome for network + enzyme identities and
**calibrates** thermal parameters against measured growth data — it does not
predict TPCs blind from sequence (seq->Topt is ~12 C RMSE, far coarser than the
2-4 C interspecies contrast). So: genome IS the input (reactions, enzymes, kcat,
Topt priors); the measured TPCs pin the thermal layer. Frame the result as
"a genome-scale mechanistic model parameterised from genomic data and calibrated
to our TPCs," never "the genome alone predicts the curves."

**The strong, non-circular test — leave-one-species-out:** calibrate the thermal
framework on species with best support (auris = iRV973, parapsilosis = iDC1003),
then PREDICT a held-out species' TPC (haemulonii or duobushaemulonii, built by
orthology transfer) from its genome and test vs the measured curve. The held-out
curve never enters calibration -> genuine genome-driven growth estimate.

## v4 — second review: conditional GO, non-negotiables locked

**Task 0b — RESOLVED (by necessity): `r` is the operational growth rate.**
No independent OD/cell-count time-courses exist (only inoculum counts in otu_inoc.csv),
so `r` is the sole growth-rate source. Calibrate against r*60; state the etcGEM
predicts "growth as defined by the O2 method." The growth/activation conflation and
the fast-doubling point are inherited, manuscript-wide limitations to STATE, not new
risks — no cleaner ground truth is available or required. Original note follows:

**(context) is `r` the right KIND of rate?** `r` is the
O2-acceleration rate that (per the project's methods paper) conflates population
growth with per-cell activation — the identifiability wall. An etcGEM predicts
pure biomass-specific growth. Options: (1) accept `r` as the operational growth
rate and state the etcGEM predicts "growth as defined by the O2 method" (consistent
with the manuscript); (2) if OD / cell-count time-courses exist, derive true
population mu from them and calibrate to THAT. Also: mu~0.5-0.9/h implies 46-83 min
doubling — fast; sanity-check that raw curves actually double on that timescale
(iRV973 measured auris ~0.11/h in defined medium at 30 C).

**Locked design changes (non-negotiable):**
- THERMAL IDENTIFIABILITY: 2 curves cannot fit per-enzyme thermal params (etcYeast
  had 2,292, still underdetermined). Calibrate only a LOW-DIM shared hyperparameter
  set: global bias+variance on sequence-predicted Topt; shared/class dCp‡; low-dim
  NGAM(T); one protein-pool scale. NEVER hundreds of per-enzyme values. Keep a
  parameter ENSEMBLE; propagate to held-out predictive intervals (not one best fit).
- HELD-OUT mu_max RULE: if the held-out species' mu_max scales its uptake, the
  claim is NORMALIZED shape mu(T)/mu_max, NOT absolute TPC. Primary endpoint =
  normalized shape + Topt + upper decline, conditional on ONE benign-T anchor.
  Absolute mu(T) = stretch, using ONLY training-derived uptake scaling. Held-out
  curve / mu_max / respiration untouched for the primary test.
- HOLDOUT DESIGN: train iRV973 + iDC1003 (frozen), predict BOTH haemulonii &
  duobushaemulonii. Not "train-3-test-1" (only 1 prediction, uses a draft species
  in calibration). n small -> mechanistic case study, not powered CV; state so.
- POOL BINDING: verify, don't infer from mu. Report pool utilisation + shadow price
  per T; relax pool 10-20% and show growth moves; compare GEM/ecGEM/etcGEM TPCs;
  audit growth-carrying flux for enzyme-free bypasses.
- BENCHMARKS: must beat (a) taxon-blind mean-of-training TPC and (b) a
  SEQUENCE-ABLATED etcGEM (homogenised thermal priors). If ablation barely changes
  held-out predictions, the genome isn't providing the signal. Evaluate full-curve
  error, Topt, T50-decline, max growth T, interval coverage - not a 4-point correlation.
- MEDIUM (see medium_YMS.csv v2): POOLED AA C/N budget (not 20 caps); single carbon
  path sucrose->glc+fru via invertase with transport cost (no free glc/fru/maltose);
  malt carbon as separate budgeted scenario; 3 preregistered scenarios (min / proxy
  / rich). Robust only if Topt & collapse persist within one 2 C assay step.
- BOTTLENECK ROBUSTNESS: for any enzyme claimed to set the ceiling, check if its
  product is importable from the proxy, remove that import, vary the AA budget,
  test alternate optima; require recurrence across medium + reconstruction ensembles
  before any enzyme-specific mechanistic claim.

## Calibration target (Task 0 output)
`results/tables/gem_growth_tpc_target.csv` — median mu (per HOUR) vs T per species.
This is what the etcGEM fits. auris mu_max ~0.80/h @36C; hae/duo/para 0.5-0.85/h.

# Among-species metabolic-model pipeline (Candidas) — v2

> **v2 (post-review).** An external methods review corrected several points in v1;
> the changes are folded in below and flagged **[REVIEW]**. Net effect: the
> reconstruction method changed, the "sequence predicts CUE" loop was restructured
> to remove circularity, two potential hard blockers were identified, and a
> lower-risk alternative (a focused mitochondrial energy model) is now the
> recommended fallback. **Do not begin six reconstructions until Task 0 and the
> medium question are resolved.**

## Headline framing [REVIEW — restructured]
An ecGEM can compute a CUE-like quantity, but if carbon uptake AND O2 AND CO2 are
all constrained, CUE is imposed, not predicted — that is not validation. The loop
"sequence → ecGEM → predicted CUE → measured CUE" is only sound if the validated
quantity is left UNCONSTRAINED. Use the constraint hierarchy in Step 4. Honest
claims available: the model *reproduces* the ranking under shared independent
constraints; *identifies which assumptions are required*; *predicts withheld
exchange fluxes*. Not: "sequence predicts observed CUE."

Goal: build one model per **species**, estimate growth and the growth/respiration
partition, and test against the respirometry data. Interspecies only (auris clades
collapse to one — identical genomes). Isothermal (37 °C); temperature is a stretch
goal only.

Species (phenotype available → validation targets): C. auris (12 iso),
C. haemulonii (3), C. duobushaemulonii (2), C. parapsilosis (3).
Non-phenotyped, for anchoring: C. albicans, C. pseudohaemulonii.

## Curated models available [REVIEW — more than v1 listed]
- C. auris: **iRV973** (FEMS Yeast Res 2023).
- C. parapsilosis: **iDC1003** (curated, 4-compartment).
- C. albicans: **iRV781**.
So 2 of 4 phenotyped species have curated models → reconstruction is
harmonise-and-transfer, NOT de novo carving.

---

## Task 0 — unit reconciliation — RESOLVED [v3]
**`r` is the specific growth rate per MINUTE, not per hour** (07_oxygen_fits.R:957,
`growth_fgC_h = r * cell_carbon_fg * MIN_TO_H`, MIN_TO_H=60; whole O2 model runs in
minutes). So **mu[/h] = r * 60**. The "58-fold discrepancy" was exactly this 60x
min->h factor: 0.0128/min * 1320 fg * 60 = 1013 fg C/cell/h (matches the data).
Consequences:
- real growth rates are mu ~0.5-0.9 /h, NOT 0.013/h -> **the enzyme pool WILL bind;
  the etcGEM does not collapse to a plain GEM. ChatGPT's slack-pool worry is void
  (it rested on the same units error).**
- growth_fgC_h is genuinely fg C/cell/h; gDW bridge unchanged and correct.
- RQ is ASSUMED = 1 (config.R:435), respiration C = O2 * RQ * 12/32 — CO2 is NOT
  independently measured. Fine for growth-TPC calibration (needs mu only); a caveat
  only if carbon balance / CO2 is ever used.
- NOTE: `growth_C_per_C_h` = mu[/h], NOT CUE. (CUE = growth_fgC/(growth+resp).)

### (former hard-stop notes, now closed)
Reviewer independently confirmed the discrepancy: under balanced growth
G_C = mu * m_C,cell = 0.013 * 1320 ≈ 17.2 fg C/cell/h, NOT ~1000 — a ~58-fold gap.
A 0.013/h rate is a ~53 h doubling time, implausibly slow for these yeasts at
37 °C, which is strong evidence that `r` is the O2-model acceleration rate, not
the specific growth rate mu. Resolve before ANY flux bound is built.
gDW bridge (confirmed correct): gDW/cell = cell_carbon_fg / f_C * 1e-15.
Also pin: whether r = mu; cell-count basis (buds/clumps/dead); dry mass & carbon
at the same phase/temperature; O2 units + background correction; RQ = qCO2/qO2;
organic C secretion (ethanol/glycerol/acetate); balanced vs batch-transition growth.
CUE definitions to report all three: apparent G_C/(G_C+CO2_C) [matches experiment];
uptake yield G_C/C_in; carbon recovery (G_C+C_out)/C_in.

## Step 1 — reconstruction [REVIEW — method corrected]
**CarveMe / gapseq are bacterially-oriented and strip eukaryotic compartments —
NOT appropriate defaults for compartmentalized Candida metabolism.** Instead:
1. Harmonise curated **iRV973** (auris) + **iDC1003** (parapsilosis) into one
   common Candida template (shared namespace, compartments, directionality,
   biomass, GAM/NGAM/P/O), then transfer reactions to haemulonii &
   duobushaemulonii by orthology, preserving mito/peroxisome/cytosol localization.
2. Or a fungal-specific pipeline (CarveFungi) + Candida curation.
Comparability needs identical namespace, compartments, directionality, biomass
conventions, harmonised GAM/NGAM/P/O, identical gap-fill rules, and explicit
gene-supported vs gap-filled tracking. "Same automated method" alone just imposes
the same errors everywhere. Benchmark the automated pipeline against iRV973 /
iDC1003 (what does it miss / falsely add?). **Do not claim species-unique
capabilities from gap-filled reactions** — require gene/orthology/localization
(+ phenotype) evidence; iRV973's authors found many "unique" ECs were artifacts.

## Step 2 — medium — SPECIFIED (complex; proxy required) [v3]
**YMS = 4 g/L yeast extract + 4 g/L malt extract + 4 g/L sucrose.** Complex/
undefined: the extracts are amino-acid/peptide/vitamin/mineral/sugar mixtures, not
encodable quantitatively from the recipe. Encoded in `gem/medium_YMS.csv`:
- sucrose = defined C source (+ glc/frc from hydrolysis); malt adds maltose (extra
  unquantified C, flagged).
- 20 amino acids + standard vitamins as a documented YE proxy; **amino-acid uptake
  CAPPED low so extracts supplement, not replace, biosynthesis** (free AA otherwise
  let FBA bypass synthesis and inflate growth).
- minerals / NH4 / O2 open (not limiting).
- absolute uptake SCALE fixed during calibration to a species' measured mu_max;
  temperature then modulates.
**Why not fatal for the TPC objective:** the temperature SHAPE (Topt, collapse) is
enzyme-driven; the medium sets only the baseline scale, which calibration absorbs.
Complex medium degrades absolute-growth prediction (not claimed) far more than
temperature-response prediction (the actual objective).
**Required robustness check:** sensitivity analysis over the amino-acid uptake cap
— it sets how much biosynthesis bypass is allowed, which shifts the rate-limiting
enzymes and therefore the predicted TPC. Report TPC ranking robustness across caps.
Exchange IDs in the CSV are BiGG hints; map to each model's namespace during
harmonisation (iRV973 / iDC1003 may differ).

### (former blocker note)
Exchange bounds are RATES, not concentrations. "20 g/L glucose" does not fix
uptake in mmol/gDW/h. If YMS contains yeast extract / peptone / malt (undefined),
it cannot be encoded quantitatively — the optimiser will pick a favourable
nutrient mix and inflate growth & CUE. Need: defined C/N sources, measured/
defensible uptake bounds, AA/vitamin availability, O2 transfer, pH, lipid/sterol
assumptions, by-product measurements. Rich YMS with no uptake data → qualitative
feasibility/sensitivity only, NOT quantitative CUE validation.

## Step 3 — enzyme constraint [REVIEW — pool + kcat cautions]
Fixed pool = THREE assumptions (equal total protein/gDW, equal metabolic-protein
fraction, equal saturation), not one. Fitting the pool to auris makes auris a
CALIBRATION species that can no longer validate. Preferred, descending:
(1) measure protein/gDW (+ proteome allocation) per species;
(2) external protein priors with a propagated range;
(3) shared-pool + species-specific sensitivity;
(4) fit to an EXTERNAL control (C. albicans), leaving all 4 targets free.
DLKcat: fills orders of magnitude, but poor generalization below ~60% identity to
training and no assay-condition encoding — do not treat as resolving subtle
ortholog differences. Run shared ortholog-kcat baseline + species DLKcat + kcat
uncertainty ensembles, and **test whether the protein pool is even binding**:
if true growth ~0.013/h, capacity may be slack and the ecGEM collapses to the GEM
→ no interspecies separation. Check this early.

## Step 4 — constraints [REVIEW — hierarchy, leave validated quantity free]
| mode | constrain | test (left free) |
|---|---|---|
| predictive | medium + measured C uptake | growth, O2, CO2, CUE |
| respiratory validation | growth + C uptake | O2/CO2, by-products |
| yield validation | C uptake + O2 or CO2 | growth, organic secretion |
| diagnostic/calibration | growth + uptake + O2 + CO2 | infer NGAM, P/O, allocation |
Last row is reconciliation, NOT validation.
**Respiratory-chain structure is decisive:** iRV973 imports GAM (yeast), NGAM
(C. tropicalis), P/O (S. cerevisiae) — these inherited values largely SET CUE.
Models need explicit coupled reactions for cytochrome respiration, AOX,
alternative NADH dehydrogenases, proton translocation, ATP synthase, and leak.
**SHAM/azide partitioning directly constrains phosphorylating vs
non-phosphorylating O2 flux and is likely more informative here than DLKcat.**

## Step 5 — validation [REVIEW — ranking, not correlation; nested levels]
n=4: pre-specify the observed ranking (auris > para > duo ≈ hae), report absolute
errors vs experimental uncertainty, and test whether predicted intervals are
NARROWER than the observed 0.21 CUE contrast and whether the ranking is robust
across plausible NGAM/P/O/uptake/pool values. Compare nested levels:
(1) common stoichiometric model; (2) species reaction content; (3) enzyme
constraints; (4) explicit alternative-respiration architecture; (5) species
physiological inputs — to see which layer carries explanatory power. Auxotrophy /
substrate-use claims valid only if NOT used for gap-fill. Reconstruction controls:
iRV973, iDC1003 (more relevant than albicans as an implementation check).

## Step 6 — temperature (STRETCH, hypothesis-only) [REVIEW — agreed]
etcYeast used Bayesian inference against organism-specific phenotype data, not
clean sequence-only prediction. Seq-based enzyme-Topt is ~12 °C RMSE — far coarser
than the 2–4 °C interspecies contrast, and Topt alone is insufficient (need
T-dependent activity, stability, abundance, maintenance). Hypothesis generation only.

## The three worries — reviewer verdicts [REVIEW]
1. Conserved core metabolism: real problem. On rich medium, harmonised models
   predict near-identical yields. Enzyme layer separates species ONLY if capacity
   binds, kcat differences are real, and allocation differs / is specified —
   otherwise apparent separation is DLKcat/annotation/gap-fill artifact. Given the
   comparative-genomics result, regulation/maintenance/respiratory coupling are
   more plausible differentiators than reaction presence/absence.
2. CUE as an ecModel output: POTENTIAL CUE yes; observed CUE from sequence alone
   no. Report a feasible CUE envelope (FVA/sampling), not one optimum — growth
   maximization returns max yield, not realized physiology.
3. Pilot: auris+haemulonii is a POOR first test (only auris curated; networks near-
   identical). Better pilot: harmonise curated auris + parapsilosis, verify units &
   carbon balance, run CUE sensitivity WITHOUT using measured CUE, check ranking
   robustness vs the 0.21 contrast, THEN transfer the template to haemulonii.

## GO / NO-GO [REVIEW]
Do NOT start six reconstructions until: unit discrepancy resolved; carbon uptake
available as a RATE; empirical vs modeled CUE definitions matched; plausible
NGAM/P/O do not reverse the ranking; the enzyme pool actually binds; predicted CUE
uncertainty < the biological contrast; no validation variable used for fit/gap-fill.

**Fallback if uptake/RQ/secretion data are unavailable (recommended):** drop
"sequence predicts CUE" and build a focused **mitochondrial energy model** —
cytochrome respiration + AOX + alternative NADH dehydrogenases + P/O + maintenance
— constrained by the SHAM/azide partitioning data. Closer to the biological
question (growth–respiration decoupling, the para/AOX story), weeks not months, and
not hostage to inherited genome-scale energetic assumptions.

---
## What runs where
Reconstruction / GECKO / solvers: Mac or cluster (DIAMOND, MILP solver, MATLAB or
py-ecModel). Constraint prep / unit bridging / validation plots: plain Python/R.
