# etcGEM — publication-ready text (Figure, Results, Methods)

## FIGURE LEGEND

**Figure X. A genome-parameterized enzyme- and temperature-constrained metabolic
model transfers partially across species but fails at the species-specific thermal
boundary.**
(**A**) Measured growth thermal-performance curves (TPCs) for *C. auris*,
*C. haemulonii*, *C. duobushaemulonii* and *C. parapsilosis* from O2 respirometry
(median specific growth rate; wells with no detectable growth counted as zero).
(**B**) Leave-one-species-out (out-of-fold) etcGEM predictions (lines) versus
measured values (points); each predicted curve is the fold in which that species
was held out. (**C**) Out-of-fold benchmark decomposition. For thermal optimum
(mean absolute error, °C) the genome-specific model (2.5 °C) equals an
identical-thermal-profile null (2.5 °C) and is worse than a trivial
median-optimum predictor (1.5 °C); for whole-curve RMSE the ecGEM (genome 0.18,
null 0.17) beats a trivial mean-curve baseline (0.23), i.e. it captures broad TPC
structure but the sequence-derived thermal layer adds nothing over a shared
profile. (**D**) Death-cliff false positives: at temperatures where growth is
measured to be zero, the out-of-fold model still predicts substantial growth
(0.5–0.7 h⁻¹). Shared thermal-shape hyperparameters were calibrated once and frozen;
only sequence-predicted per-enzyme optimal (Topt) and melting (Tm) temperatures
differ between species.

## RESULTS (paragraph)

To test whether interspecies thermal divergence is encoded in metabolic-enzyme
sequence, we built enzyme- and temperature-constrained genome-scale models
(etcGEMs) for the four species, parameterizing every reaction with a sequence-
predicted turnover number (DLKcat) and its enzyme's sequence-predicted optimal and
melting temperatures (Seq2Topt/Seq2Tm), and constraining flux by a shared protein
pool. After shared calibration the model reproduced the growth-rate scale and the
broad mesophilic temperature dependence, and in leave-one-species-out
cross-validation it beat a trivial mean-curve baseline on whole-curve error
(RMSE 0.18 vs 0.23), indicating that genome-parameterized metabolic capacity
transfers partially across species. However, the model failed to recover the
species-specific features that distinguish the TPCs: out-of-fold thermal-optimum
error (2.5 °C) was no better than an identical-thermal-profile null (2.5 °C) and
worse than a trivial median-optimum predictor (1.5 °C), and the model produced
large false-positive growth predictions (0.5–0.7 h⁻¹) at temperatures where no
growth was observed. Predicted per-enzyme thermostability was clade-conserved
(median paired ortholog ΔTopt = 1.9 °C; 1st-percentile Tm within 0.5 °C across
species), so predicted optima collapsed onto a common ~34 °C (≈1 °C spread) against
a measured spread of 6 °C; propagating the predictor's error made the predicted
optima statistically unresolvable (90% intervals ≥6 °C wide, fully overlapping).
The failure was not an artifact of the model: a signal-injection test confirmed the
architecture propagates real enzyme-level shifts to the cellular optimum
(shifting the top 10% highest-demand enzymes by +6 °C moved the optimum +3 °C), and
the compression was robust to network reconstruction (identical-scaffold control),
GPR complex structure, growth medium, and kcat aggregation. Thus, within the
resolution of current enzyme-property predictors, metabolic-enzyme sequence carries
limited-but-insufficient information to explain the observed thermal divergence — a
result consistent with contributions from unresolved enzyme-level effects and/or
regulation, proteome allocation, membrane physiology and cellular stress responses
(implicated, not identified), and consistent with the phylogenomic near-identity of
the Candidozyma clade.

## METHODS (paragraph)

**Enzyme- and temperature-constrained models.** Base genome-scale models were
iRV973 (*C. auris*, re-keyed from B8441 to the B9J08/RefSeq proteome by KofamScan KO
annotation) and iDC1003 (*C. parapsilosis*); *C. haemulonii* and
*C. duobushaemulonii* were reconstructed by projecting the auris scaffold and
assigning each metabolic reaction's gene–protein–reaction rule from the target
species' own genes by KO→reaction (primary) and KO→EC (secondary) evidence, with
all foreign genes removed and unsupported reactions retained but flagged. Each
gene-associated reaction was assigned a molecular weight (per protein sequence) and
a turnover number kcat predicted from sequence and substrate SMILES by DLKcat, and
flux was constrained by a single sMOMENT protein pool Σ_r (MW_r/kcat_r)·v_r ≤ P.
The temperature layer scaled each enzyme's effective kcat by act(T) =
exp(−(T−Topt)²/2σ²)·[1+exp((T−Tm)/w)]⁻¹, where Topt and Tm were predicted per
enzyme from sequence (Seq2Topt/Seq2Tm, ESM2-8M embeddings). The three shared
hyperparameters (σ, w, P) plus one global growth scale were calibrated once and
frozen; the only species-specific inputs were the sequence-predicted per-enzyme
Topt and Tm. Growth thermal-performance curves were computed by solving the
constrained model across temperature. **Evaluation.** Measured TPCs were the median
specific growth rate per temperature with non-growing wells scored as zero. Model
performance was assessed by leave-one-species-out cross-validation against an
identical-thermal-profile null (global-median Topt/Tm) and trivial baselines
(median measured optimum; mean of the other species' measured curves), reporting
raw and maximum-normalized RMSE, thermal-optimum error, and predicted growth at the
first observed no-growth temperature. Robustness was assessed by predictor-error
propagation (Topt σ=12 °C, Tm σ=7 °C; independent and ρ=0.5 correlated),
signal-injection recovery, an identical-scaffold common-network control, isozyme-OR
versus complex-AND GPR treatment, three amino-acid medium scenarios, and median
versus maximum kcat aggregation. Analyses used cobrapy with a GLPK solver.

## ONE-LINE TAKE
Genome-parameterized metabolic capacity transfers partially across species but
fails at the species-specific thermal boundary; the differentiating information is
not in metabolic-enzyme sequence.
