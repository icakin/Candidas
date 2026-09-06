# etcGEM for interspecies Candida thermal prediction — full method & result (for review)

Goal: build an enzyme- and temperature-constrained genome-scale metabolic model
(etcGEM) that uses each species' genome to predict its growth thermal-performance
curve (TPC), and test whether it reproduces the measured interspecies differences
across four Candida species: *C. auris*, *C. haemulonii*, *C. duobushaemulonii*
(the Candidozyma clade) and the outgroup *C. parapsilosis*.

## 1. Base metabolic models
- *C. auris*: curated iRV973. Its GPRs used B8441 GenBank locus tags (CJI97*) that
  exist in no protein sequence DB, so we re-keyed all genes to the B9J08/RefSeq
  proteome (XP_ accessions) via KofamScan KO annotation. Result: fully
  sequence-linked model.
- *C. parapsilosis*: curated iDC1003 (BioModels), native KEGG namespace.
- *C. haemulonii*, *C. duobushaemulonii*: no curated model exists. Built drafts by
  projecting the auris scaffold and assigning each metabolic reaction's GPR from
  the TARGET species' own genes by evidence (KO->KEGG-reaction primary, KO->EC
  secondary); reactions with no species evidence were kept but GPR-cleared and
  flagged; ALL foreign auris gene IDs removed. (Complex-AND structure not
  reconstructed — curated models have only 4-6 transporter AND-complexes, no
  metabolic complexes, so all GPRs treated as isozyme OR.)
- Medium: YMS (4 g/L yeast extract + 4 g/L malt extract + 4 g/L sucrose) encoded
  as a glucose carbon proxy + pooled amino-acid carbon budget + biotin auxotrophy.
  All four models grow on YMS.

## 2. Enzyme (GECKO/sMOMENT) layer
- Molecular weight per gene from its protein sequence (100% coverage).
- kcat per reaction predicted from sequence + substrate SMILES by DLKcat (ran in
  sandbox). 1062-1313 reactions per species.
- Single-pool sMOMENT constraint: sum_r (MW_r / kcat_r) * v_r <= P.
- Pool-binding test (auris): plain FBA gives 2.05/h; the ecGEM at a literature
  protein budget P~0.25-0.36 g/gDW gives ~0.76/h, matching measured mu_max ~0.80/h.
  => the pool constraint is load-bearing (not slack); the ecGEM is in the right
  regime. This is the key precondition and it passes.

## 3. Temperature (etc) layer
- For every enzyme, predicted optimal temperature Topt and melting temperature Tm
  from sequence, using Seq2Topt / Seq2Tm (ESM2-8M embeddings + attention head;
  reported held-out RMSE ~12 C for Topt). Predicted for all 3,042 enzymes
  (auris 704, hae 679, duo 662, para 997).
- Temperature-dependent activity per enzyme:
    act_g(T) = exp(-(T-Topt_g)^2 / (2 sig^2)) * 1/(1+exp((T-Tm_g)/w))
  (Gaussian peak at the predicted Topt, sigmoid denaturation cutoff at the
  predicted Tm.) Effective kcat_r(T) = kcat_r * act_g(T); pool cost = MW_r/kcat_r(T).
- Shape hyperparameters {sig, w, P} are SHARED across species and calibrated ONCE
  on the measured *C. auris* TPC (Nelder-Mead, multi-start), then FROZEN to predict
  the other three species with no further fitting. So any interspecies TPC
  difference comes ONLY from the sequence-predicted per-enzyme {Topt, Tm}.
  Calibrated: sig~10.2 C, w~8.8, P~0.36; one global growth-scale from the auris fit.
- TPC = sweep T, solve FBA with the pool constraint at each T, record max growth.

## 4. Measured TPCs and a data-quality correction
- Measured growth rates are from O2 respirometry (r per minute, x60 -> /h), per
  isolate x 5 replicate wells, 22-44 C.
- IMPORTANT: the original per-T summary took the median over SURVIVING wells only.
  At high T most wells are dead, so that median is over a heat-tolerant handful and
  is inflated — e.g. haemulonii declined normally 30->36 C (n=15 wells) then
  "peaked" at 42 C on only 5 surviving wells (survivorship artifact).
- We recomputed the TPCs counting dead/invalid wells as zero growth (correct for a
  thermal-TOLERANCE curve). haemulonii then peaks at 30 C, not 42 C. Corrected
  measured optima span 30-36 C (6 C), not 32-42 C (10 C). We recalibrated on the
  corrected auris curve.

## 5. Result
Predicted optima (fine 0.5 C grid): auris 34.5, hae 34.5, duo 34.0, para 33.5 C
— span ~1 C (auris warmest, para coolest: directionally correct). Measured optima:
auris 36, hae 30, duo 32, para 32 C — span 6 C. (On the 2 C experimental grid all
predictions round to 34 C.)
- The model reproduces absolute growth rate (pool binds at measured mu_max) and the
  general mesophilic TPC shape, and nails the calibration species.
- It does NOT reproduce the interspecies differences: predicted optima compress to
  ~1 C, and it does not capture auris's distinctive high-T survival (auris grows to
  44 C; the other three die at 38 C).

## 6. Why — mechanism, verified
- The predicted Topts are NOT uniform: within a species they span ~62 C (21-83 C),
  so the predictor is strongly sequence-sensitive; all 3,042 input sequences are
  distinct (no duplicated input).
- What coincides is the species MEDIAN, because the species share the same enzymes
  as orthologs at 86% median identity. The SAME enzyme in auris vs hae differs by a
  median of only 1.8 C in predicted Topt, vs a 62 C range across different enzymes.
  The difference tracks divergence (orthologs >92% identity: 1.2 C; <79%: 2.1 C),
  confirming the predictor responds to sequence — there is just too little
  between-species sequence difference to yield more than ~1-2 C per enzyme.
- Tm is likewise conserved: medians within 0.4 C, and even the 1st-percentile Tm
  (death-relevant) within 0.5 C across all four species.

## 7. Attempts to legitimately widen the predicted spread (all reported honestly)
1. Flux-control weighting (weight thermal scaling by rate-controlling enzymes at
   high T): does NOT help — the rate-controlling core enzymes are the MOST
   conserved, so the flux-weighted Topts compress further (33.9-34.7 C).
2. Alternative whole-proteome OGT proxy (IVYWREL composition, Zeldovich 2007):
   ranks auris highest and parapsilosis lowest (directionally consistent with
   auris's heat tolerance) but the spread is 0.59 pct-points (~3-4 C, prokaryote-
   calibrated, noisy) — weak corroboration, not a calibrated fix.
3. We did NOT fit per-species parameters to the measured curves (that would be
   circular / non-predictive). The compression is robust to every legitimate fix.

## 7b. Validation audits (added after external review)
- FREE-REACTION AUDIT: at the optimum, ~88% of internal flux-sum runs through
  uncosted reactions, but these are legitimately enzyme-free (53% multi-compartment
  transport, 47% spontaneous/no-gene); ZERO gene-associated reactions were left
  uncosted (no mis-costing). The pool still binds hard (auris 0.39 with pool vs
  2.05 without, 5.3x), so growth is gated by the costed enzymes that carry the
  thermal layer. Drafts have the same free-flux profile as curated auris (no
  draft-specific bypass). Reversible reactions charged by |v| (no negative cost).
- SIGNAL-INJECTION / RECOVERY: shifting every auris enzyme Topt+Tm by +2/+4/+6 C
  moves the predicted CELLULAR optimum by +2.0/+4.0/+6.0 C (1:1). => the shared
  response architecture faithfully propagates enzyme thermal differences; the flat
  predictions come from flat INPUTS, not from architectural compression. This
  directly refutes the "shared smooth response predetermines the null" objection.
- OUT-OF-FOLD DECOMPOSITION (the key benchmark; supersedes the earlier in-sample
  null test). Leave-one-species-out, genome-specific vs identical-profile null vs
  trivial baselines:
    * optimum MAE:  trivial median-optimum 1.5 C  <  genome 2.5 C  =  null 2.5 C
      -> the etcGEM does NOT beat a trivial "predict the median optimum" baseline,
         and genome-specific = null (no thermal-sequence advantage).
    * whole-curve RMSE:  null 0.170  ~=  genome 0.183  <  trivial mean-curve 0.226
      -> the ecGEM DOES capture broad TPC structure (beats trivial curve), but the
         genome-specific thermal layer adds NOTHING over a shared profile out-of-fold.
    * NOTE: the earlier "genome beats null 0.251 vs 0.324" was IN-SAMPLE (auris-
      calibrated); under proper LOOCV the advantage vanishes. This is the honest
      number and it is a cleaner (near-null) result, not a partial success.
- LOOCV (all 4 folds): mean |optimum error| 2.5 C, and every fold predicts
  substantial growth (mu 0.44-0.68) PAST each species' measured no-growth temperature
  (death-cliff false positives). Panel D quantifies this (hae, para at 40 C).
- REALISTIC SIGNAL-INJECTION (answers "the uniform 1:1 test is tautological"):
  shifting Topt/Tm by +6 C for ONLY the top 10% highest-enzyme-demand reactions
  moves the cellular optimum +3 C (top 25% -> +4 C). So the model IS sensitive to
  realistic, localized thermal adaptation in the rate-controlling enzymes; the flat
  prediction reflects absent sequence signal, not model insensitivity.
- CORRELATED PREDICTOR-ERROR SENSITIVITY: with ortholog-error correlation rho=0 the
  auris optimum 90% interval is 6 C wide; at rho=0.5 it is 20 C wide (24-44 C). Under
  any plausible correlated error the optimum is essentially unconstrained - so the
  species optima are unresolvable, whichever error model is assumed.
- UNCERTAINTY PROPAGATION: resampling each enzyme's Topt/Tm within predictor error
  (sigma_Topt=12, sigma_Tm=7 C; independent per-enzyme jitter, the correct model for
  the interspecies DIFFERENCE since systematic bias is common-mode) gives predicted-
  optimum 90% CIs that are WIDE and FULLY OVERLAPPING: auris 34 [32-40], hae 34
  [32-40], duo 34 [30-38], para 36 [32-38]. => within predictor resolution the
  species optima are statistically indistinguishable; the ~1 C point-spread is noise.
- COMPLEX-GPR SENSITIVITY: recomputing with enzyme-complex AND rules (pool cost =
  SUM of subunit MW) instead of isozyme-OR (mean MW) leaves the optima compressed
  (OR spread 2 C, AND spread 0 C). Robust to GPR complex structure.
- COMMON-NETWORK CONTROL: running each species' enzyme profile on the IDENTICAL
  auris scaffold (orthology-mapped) still compresses the optima to 2 C - so the
  compression is a SEQUENCE effect, not a network-reconstruction artifact of the
  projected drafts.
- MEDIUM SENSITIVITY: optima spread stays 1-2 C under restricted / base / rich
  amino-acid availability. kcat MEDIAN-vs-MAX over substrates: spread stays 1-2 C.
  Robust to medium and kcat aggregation.

## 8. Honest conclusion (final - out-of-fold / model-boundary framing)
Across leave-one-species-out tests, a genome-parameterized metabolic model
reproduced the growth-rate SCALE (after shared calibration) and broad mesophilic
temperature dependence, and beat a trivial mean-curve baseline on whole-curve RMSE
(0.18 vs 0.23) - so metabolic capacity transfers partially across species. But it
FAILED OUT-OF-FOLD to recover species-specific thermal optima (optimum MAE 2.5 C,
WORSE than a trivial median-optimum predictor at 1.5 C) and upper-temperature
growth (death-cliff false positives), and the genome-specific enzyme-sequence
thermal layer added NO out-of-fold advantage over an identical-profile null.
Thus, within the resolution of current enzyme-property predictors, metabolic-enzyme
sequence contains limited-but-insufficient information to explain the observed
thermal divergence; the missing predictive information is CONSISTENT WITH
CONTRIBUTIONS FROM unresolved enzyme-level effects and/or regulation, proteome
allocation, membrane physiology and cellular stress responses (implicated, not
identified). We predicted (did not measure) enzyme thermostability, and a ~12 C-RMSE
predictor cannot establish biological equality at 1-6 C, so this is a cross-
validated model-boundary / falsification result, NOT a causal localization.

Recommended one-line claim (per reviewer): "Current genome-parameterized metabolic
capacity transfers partially across species but fails at the species-specific
thermal boundary."
NOTE: we predicted enzyme thermostability, we did not measure it; and a ~12 C-RMSE
predictor cannot establish biological equality at the 1-6 C scale. So this is a
model-boundary / falsification result, NOT a demonstrated causal localization.

## 9. Known limitations / points for critique
- Only one outgroup (parapsilosis) — the "distance" axis is effectively 2-point.
- Two of four models are drafts (hae, duo) without complex-AND GPRs.
- The thermal model shape is a shared 3-parameter phenomenology, not the full
  Li-2021 thermodynamic form (no ΔG_folding/ΔCp per enzyme).
- Seq2Topt/Seq2Tm resolution (~12 C RMSE) may itself mask sub-degree real
  differences — so "not encoded in sequence" is partly "not resolvable by this
  predictor"; both point to no usable per-species signal.
- Death cliff (organismal upper limit) is not modelled from the least-stable
  essential enzyme; median/percentile Tm did not separate species anyway.

## 10. How performance could be genuinely improved (not by fitting the answer)
Within-clade prediction has a low ceiling because the signal is not in the
metabolic-enzyme sequences; the real gains are in adding biology or scope:
- Death-cliff via ESSENTIAL enzymes: identify growth-essential reactions by FBA
  single-deletion, then set the organismal upper limit from the least-stable
  essential enzyme's Tm rather than the median. (We checked bulk/percentile Tm —
  conserved — but not the essential subset specifically; worth one clean test.)
- Expression/abundance weighting: weight each enzyme's pool cost by measured
  transcript/protein abundance (RNA-seq or proteomics at 2-3 temperatures). If
  auris up-regulates heat-stable isoenzymes or chaperones, this would capture a
  real, non-sequence signal the current model omits.
- Add the missing biology the result points to: membrane-lipid remodelling,
  heat-shock/chaperone capacity — these are prime candidates for where the thermal
  divergence lives; even a coarse chaperone-capacity term could help.
- Broaden the phylogenetic range: add 1-2 species at intermediate distance
  (C. albicans, C. lusitaniae — genomes already in hand) to turn the 2-point
  distance axis into a real distance-vs-accuracy curve; the framework works better
  between clades than within one.
- Rigorous performance metric: leave-one-species-out cross-validation (calibrate
  on 3, predict the 4th, rotate) instead of calibrating only on auris — reports an
  honest generalisation number rather than an in-sample fit.
- Better thermal inputs: a fungal-fine-tuned Topt/Tm predictor, or experimental
  Tm for a handful of key enzymes, would lower the ~12 C predictor noise floor —
  though it cannot manufacture a difference that isn't in the sequence.
The one thing NOT to do: per-species tuning of pool size / thermal offset to match
the measured curves — that converts prediction into curve-fitting and voids the
result.
