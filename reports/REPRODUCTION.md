# REPRODUCTION — Candidas, C1

Everything below was produced by re-running this repository end to end from the
committed raw data, on the machine recorded in `env/versions.json`
(macOS 26.2, Apple silicon, 16 cores; R 4.5.2; Python 3.9.6; rstan 2.32.7;
etc-GEM engine pinned to `7d32383`). Nothing in `results/`, `data/` or the
committed etc-GEM outputs was modified — everything went to `runs/C1_reproduction/` and
`runs/C1/supp_data/`.

This is a **reproduction report**, not a science report. It changes no analysis
decision. Where a number does not come back, it says so, with the number.

---

## VERDICT

**No — not today, not as it stands, and the failure is specific rather than
general.** After the execution fixes in this pass (which change no analysis
decision), a third party can now clone the repository, run one command, and
regenerate the whole pipeline unattended in **49 minutes**. The **respirometry
half reproduces**: every table from scripts 02–07 comes back IDENTICAL or
NUMERICALLY EQUIVALENT (worst disagreement 2.4 × 10⁻¹⁰), and every load-bearing
Bayesian number — E_G, E_R, ΔE, T_opt(CUE), the 37 °C and 40 °C carbon taxes,
all 10 pairwise fever contrasts, and which 8 of them resolve — comes back within
**0.4 % relative**, well inside its own credible interval. The **etc-GEM half
does not reproduce at all**: the committed effective-capacity values
(kcat_scale 0.614 / 0.410 / 0.622 / 0.710) regenerate as **0.833 / 0.561 /
0.851 / 0.945**, a systematic rescaling by 1.33–1.37×. The cause is identified,
not merely observed: the committed `calibration_r2.csv` and the `calibrated`
column of `calibration_curves.csv` were built with the **previous single-point
anchor** (`ANCHOR_MU = 1.1 @ 313.15 K`) while the code in the repository now uses
`0.95 @ 307.15 K`, and 1.089/0.7955 = **1.369** is exactly the observed rescale.
Re-evaluating the committed parameters under the old anchor reproduces the
committed R² to **|Δ| < 5 × 10⁻¹⁶ on all four clades**; under the current anchor
it does not (R² falls to 0.09–0.46). Beyond that, `outputs/supp_data/` is not a
single build at all: four of its thirteen files are byte-identical to
`supp_data_original_backup/` while the rest are not, and
`figure3_data/model_curves.csv` uses a temperature grid (32 points, 1 °C,
18–49 °C) that the current code cannot emit (53 points, 0.5 °C, 18–44 °C). What
is missing, precisely: (i) a record of which anchor produced the committed
etc-GEM outputs and a re-run of every etc-GEM file under one anchor; (ii)
`capacity_isolates.csv`, `capacity_bootstrap_draws.csv`, `capacity_corr.csv` and
`figure3_data/` regenerated from the same build as everything else they are
compared against; (iii) seeds on the four unseeded stochastic stages listed in
§3. Until (i) and (ii) are done, every absolute capacity number in the paper —
0.61/0.41/0.62/0.71, the 12 isolate capacities, ICC1 = 90 %, and the two
cross-phenotype correlations — rests on files whose provenance the repository
cannot establish. Their *ordering* and *relative* structure do reproduce.

---

## 1. File-by-file numeric comparison

Classification, as specified: **IDENTICAL** (max |Δ| < 1e-10) / **NUMERICALLY
EQUIVALENT** (< 1e-6) / **STOCHASTIC-BUT-CONSISTENT** (differs, but every 95 %
interval in the file overlaps its partner) / **DIVERGENT** (an interval does not
overlap) / **NOT REPRODUCED**.

Raw JSON with per-column detail: `reports/tools/compare_tables.json`,
`compare_etcgem.json`, `compare_fig3.json`.

### 1a. `results/` vs `runs/C1_reproduction/` (R pipeline)

| file | rows | max abs Δ | max rel Δ | worst column | verdict |
|---|---:|---:|---:|---|---|
| `Oxygen_Curve_Code_Key.csv` | 1080 | 0 | 0 | – | IDENTICAL |
| `Oxygen_Trimmed_Series_Metadata.csv` | 1080 | 6.37e-12 | 9.27e-15 | main_run_score | IDENTICAL |
| `activation_energy_summary.csv` | 18 | 8.88e-16 | 3.64e-15 | respiration_E_eV | IDENTICAL |
| `arrhenius_growth_fgC_h_coefs.csv` | 18 | 2.67e-15 | 3.64e-15 | alpha | IDENTICAL |
| `arrhenius_respiration_fgC_h_coefs.csv` | 18 | 3.55e-15 | 2.61e-15 | alpha | IDENTICAL |
| `cue_quadratic_fit_coefs.csv` | 18 | 4.62e-14 | 5.15e-14 | T_opt_C | IDENTICAL |
| `diagnostic_hot_temps_growth_check.csv` | 111 | 0 | 0 | – | IDENTICAL |
| `failed_fits_at_bound.csv` | 17 | 0 | 0 | – | IDENTICAL |
| `group_lookup_with_delta.csv` | 1080 | 0 | 0 | – | IDENTICAL |
| `removed_points.csv` | 72 | 0 | 0 | – | IDENTICAL |
| `replication_summary_by_temperature.csv` | 204 | 0 | 0 | – | IDENTICAL |
| `sharpe_schoolfield_K_O2_rate_coefs.csv` | 72 | 0 | 0 | – | IDENTICAL |
| `sharpe_schoolfield_growth_fgC_h_coefs.csv` | 72 | 0 | 0 | – | IDENTICAL |
| `sharpe_schoolfield_respiration_fgC_h_coefs.csv` | 72 | 0 | 0 | – | IDENTICAL |
| `summary_hot_temps_growth_check.csv` | 26 | 0 | 0 | – | IDENTICAL |
| **`bayes_resp_arr_summary.csv`** | 46 | 1.73e-11 | 2.73e-15 | ess_tail | **IDENTICAL** — see note |
| `derived_N0_R_results_with_carbon.csv` | 991 | 2.36e-10 | 2.61e-12 | delta_aicc_curv | NUMERICALLY EQUIVALENT |
| `fit_coefficients_long.csv` | 3240 | 2.36e-10 | 8.57e-11 | delta_aicc_curv | NUMERICALLY EQUIVALENT |
| `fit_coefficients_wide.csv` | 1063 | 2.36e-10 | 1.40e-11 | delta_aicc_curv | NUMERICALLY EQUIVALENT |
| `fit_metrics.csv` | 1080 | 2.36e-10 | 8.57e-11 | delta_aicc_curv | NUMERICALLY EQUIVALENT |
| `bayes_E_resp_minus_growth.csv` | 5 | 5.99e-03 | 1.44e-02 | lo | STOCHASTIC-BUT-CONSISTENT |
| `bayes_clade_contrasts.csv` | 40 | 1.66e-01 | 4.69e-01 | lo | STOCHASTIC-BUT-CONSISTENT |
| `bayes_clade_params.csv` | 5 | 1.87e-02 | 6.98e-03 | growth_Th_K | STOCHASTIC-BUT-CONSISTENT |
| `bayes_cue_by_clade.csv` | 600 | 3.59e-03 | 1.75e-02 | hi | STOCHASTIC-BUT-CONSISTENT |
| `bayes_cue_by_isolate.csv` | 1800 | 1.42e-03 | 1.45e-02 | hi | STOCHASTIC-BUT-CONSISTENT |
| `bayes_growth_ss_summary.csv` | 72 | 3.40e+03 | 6.02e-01 | ess_bulk | STOCHASTIC-BUT-CONSISTENT — see note |
| `carbon_tax_contrasts_fever.csv` | 4 | 1.04e-02 | 1.95e-02 | lo | STOCHASTIC-BUT-CONSISTENT |
| `carbon_tax_curves.csv` | 445 | 4.31e-01 | 4.32e-02 | hi | STOCHASTIC-BUT-CONSISTENT |
| `carbon_tax_group.csv` | 20 | 4.32e-02 | 1.04e-02 | tax | STOCHASTIC-BUT-CONSISTENT |
| `carbon_tax_isolate.csv` | 30 | 1.94e-02 | 6.23e-03 | hi | STOCHASTIC-BUT-CONSISTENT |
| `fig_contrasts.csv` | 10 | 1.04e-02 | 1.95e-02 | lo | STOCHASTIC-BUT-CONSISTENT |
| `fig_values.csv` | 5 | 2.22e-01 | 1.44e-02 | g40_hi | STOCHASTIC-BUT-CONSISTENT |
| `bayes_resp_model_comparison_loo.csv` | 2 | 8.90e-05 | 1.23e-05 | looic | (no interval) — see note |
| `manual_fit_windows.csv`, `plot_exclude_points.csv`, `otu_cell_sizes.csv`, `otu_inoc.csv`, `otu_names.csv` | – | – | – | – | **NOT REGENERATED — BY DESIGN** |

Notes.

* **`bayes_resp_arr_summary.csv` is exactly reproduced.** Every posterior
  quantity — `mean`, `sd`, `q2.5`, `q50`, `q97.5`, `rhat` — differs by **exactly
  0.0**. The only non-zero column is `ess_bulk`/`ess_tail` at ~1e-11, which is
  floating-point summation order in the ESS estimator, not sampling. This is the
  proof that **brms with `seed = 1234` is bit-deterministic on this machine**:
  given bit-identical input data, the sampler returns bit-identical draws.
* **`bayes_growth_ss_summary.csv`'s 3.4e+03 is in `ess_bulk`**, an MCMC
  diagnostic, not an estimate. Every parameter's posterior mean, sd and
  quantiles agrees to ≤ 1.9e-2, and every 95 % interval overlaps.
  Why this model and not the respiration one? See §3.
* **`bayes_resp_model_comparison_loo.csv`** carries no lo/hi pair, so the rule
  cannot classify it. Its worst disagreement is 8.9e-5 on `looic` (1.2e-5
  relative) against an SE of ~30 — numerically equivalent in substance.
* The five **NOT REGENERATED** files are the committed outputs of the four
  click-driven apps. They are **inputs**, treated as data; regenerating them
  would require a human to re-make hand judgements. See §5.
* `runs/C1_reproduction/tables/` additionally contains three large intermediates
  (`Oxygen_All_Long.csv`, `Oxygen_Data_Filtered.csv`,
  `Oxygen_Data_Smoothed_Trimmed.csv`) that `.gitignore` keeps out of the repo.
* **A figure-inventory difference worth recording.** `results/figures/` has 83
  files, `runs/C1_reproduction/figures/` 78. Every difference is accounted for. Only in
  `results/`: **7 stale top-level copies** of the manuscript figures
  (`FIG1_decoupling.png`, `FIG2_the_bill.png`, `FIG_MODEL.png`,
  `FIG_MODEL_SUPP.png`, `FIG_MODEL_SUPP_validation.png`,
  `FIG_MODEL_SUPP_consistency.png`, `FIG_model_schematic.png`) sitting beside
  the real ones in `figures/manuscript/`. **No current script writes them** —
  12, 13 and 14 all write to `figures/manuscript/`. They are leftovers from an
  older layout, and a reader cannot tell which copy is current. Not deleted
  here (that would be a change to `results/`). Only in `runs/C1_reproduction/`: the two
  large diagnostic PDFs `per_series_fits.pdf` and
  `oxygen_trimming_diagnostics.pdf`.

### 1b. `outputs/supp_data/` vs `runs/C1/supp_data/` (etc-GEM)

| file | rows | max abs Δ | max rel Δ | worst column | verdict |
|---|---:|---:|---:|---|---|
| `calibration_measured.csv` | 48 | 0 | 0 | – | IDENTICAL (pass-through of the measured data) |
| `capacity_corr.csv` | 1 | 3.74e-02 | 3.93e-02 | r_lo | STOCHASTIC-BUT-CONSISTENT |
| `calibration_r2.csv` | 4 | **1.03e+01** | 1.44e+00 | r2_apriori | **DIVERGENT** |
| `calibration_curves.csv` | 228 | 1.36e-01 | 1.46e-01 | calibrated | **DIVERGENT** |
| `predicted_vs_observed.csv` | 48 | 1.36e-01 | 1.63e+00 | predicted | **DIVERGENT** |
| `capacity_isolates.csv` | 12 | 1.82e-01 | 2.13e-01 | capacity | **DIVERGENT** |
| `capacity_bootstrap.csv` | 4 | 2.64e-01 | 2.84e-01 | cap_hi | **DIVERGENT** |
| `capacity_bootstrap_draws.csv` | 6000 | 2.73e-01 | 4.27e-01 | capacity | **DIVERGENT** |
| `identifiability_optima.csv` | 12 | 2.35e-01 | 4.52e-01 | opt_value | **DIVERGENT** |
| `identifiability_profiles.csv` | 300 | 2.52e-01 | 8.46e-01 | rmse | **DIVERGENT** |
| `bayes_summary.csv` | 12 | 2.93e-01 | 6.86e-01 | hi | **DIVERGENT** |
| `bayes_predictive.csv` | 116 | 1.72e-01 | 1.78e-01 | hi | **DIVERGENT** |
| `bayes_posterior_samples.csv` | 10000 | 2.36e+00 | 2.00e+00 | dTopt | **DIVERGENT** |
| `bayes_chain_raw.npz` | – | – | – | – | new in C1 (raw emcee chains; no committed counterpart) |

### 1c. `outputs/figure3_data/` vs `runs/C1/figure3_data/`

| file | rows (committed → C1) | max abs Δ | verdict |
|---|---|---:|---|
| `measured.csv` | 48 → 48 | 0 | IDENTICAL |
| `capacity.csv` | 4 → 4 | 1.72e-01 | **DIVERGENT** |
| `model_curves.csv` | **128 → 212** | – | **SHAPE MISMATCH** |
| `model_uniform.csv` | **32 → 53** | – | **SHAPE MISMATCH** |

The shape mismatch is itself a provenance finding. The committed
`model_curves.csv` holds **32 temperatures per clade at 1.0 °C spacing, spanning
18–49 °C**. The current `stage_curves()` writes `np.arange(18, 50, 0.5)` filtered
to `18 ≤ T ≤ 44`, i.e. **53 temperatures at 0.5 °C spacing, 18–44 °C**. The
committed file cannot have been produced by the code in this repository.

---

## 2. The manuscript's load-bearing numbers, committed vs regenerated

Extracted with `reports/tools/headline_numbers.py`; JSON in
`reports/tools/headline_committed.json` and `headline_C1.json`.

### 2a. Bayesian thermal economics — **reproduces**

| taxon | quantity | committed | regenerated | rel Δ |
|---|---|---:|---:|---:|
| Clade1 | E_G (eV) | 0.913145 | 0.913614 | 5.1e-04 |
| Clade1 | E_R (eV) | 0.519262 | 0.519262 | **0.0** |
| Clade1 | ΔE = E_G − E_R | 0.394630 | 0.394884 | 6.4e-04 |
| Clade1 | T_opt(CUE) °C | 31.2659 | 31.2690 | 9.9e-05 |
| Clade1 | tax37 | 1.297758 | 1.298296 | 4.1e-04 |
| Clade1 | tax40 | 1.899335 | 1.899009 | 1.7e-04 |
| Clade2 | E_G | 1.144072 | 1.142917 | 1.0e-03 |
| Clade2 | E_R | 0.414078 | 0.414078 | **0.0** |
| Clade2 | ΔE | 0.730652 | 0.729453 | 1.6e-03 |
| Clade2 | T_opt(CUE) | 32.9723 | 32.9817 | 2.8e-04 |
| Clade2 | tax37 | 1.301866 | 1.300853 | 7.8e-04 |
| Clade2 | tax40 | 2.248168 | 2.246613 | 6.9e-04 |
| Clade3 | E_G | 0.946514 | 0.946846 | 3.5e-04 |
| Clade3 | E_R | 0.401281 | 0.401281 | **0.0** |
| Clade3 | ΔE | 0.544282 | 0.546534 | 4.1e-03 |
| Clade3 | T_opt(CUE) | 32.7806 | 32.7700 | 3.2e-04 |
| Clade3 | tax37 | 1.169894 | 1.170666 | 6.6e-04 |
| Clade3 | tax40 | 1.606811 | 1.608426 | 1.0e-03 |
| Clade4 | E_G | 1.019261 | 1.019220 | 4.0e-05 |
| Clade4 | E_R | 0.300647 | 0.300647 | **0.0** |
| Clade4 | ΔE | 0.718466 | 0.717902 | 7.9e-04 |
| Clade4 | T_opt(CUE) | 34.1041 | 34.1104 | 1.9e-04 |
| Clade4 | tax37 | 1.073349 | 1.073057 | 2.7e-04 |
| Clade4 | tax40 | 1.341536 | 1.341000 | 4.0e-04 |
| para | E_G | 0.829153 | 0.829667 | 6.2e-04 |
| para | E_R | 0.363160 | 0.363160 | **0.0** |
| para | ΔE | 0.467647 | 0.467663 | 3.5e-05 |
| para | T_opt(CUE) | 31.3114 | 31.2975 | 4.4e-04 |
| para | tax37 | 1.596521 | 1.596734 | 1.3e-04 |
| para | tax40 | 3.134770 | 3.127264 | 2.4e-03 |

Every ΔE remains credibly non-zero, and every T_opt(CUE) remains below 37 °C, in
both runs. **Classification: STOCHASTIC-BUT-CONSISTENT** (worst relative
difference 4.1e-3, every 95 % interval overlapping).

E_R is bit-identical for all five taxa because the respiration Arrhenius model
receives bit-identical input data and is seeded (§3).

### 2b. The 10 pairwise fever contrasts — **8 of 10 resolve, in both**

Committed (ratio [95 % CrI], `credible`):

| a | b | ratio | 95 % CrI | resolves? |
|---|---|---:|---|---|
| Clade4 | para | 0.4289 | 0.3285 – 0.5524 | ✔ |
| Clade3 | para | 0.5140 | 0.3902 – 0.6660 | ✔ |
| Clade1 | para | 0.6069 | 0.4534 – 0.7937 | ✔ |
| Clade2 | para | 0.7177 | 0.5338 – 0.9567 | ✔ |
| Clade1 | Clade2 | 0.8456 | 0.6788 – 1.0423 | ✘ |
| Clade1 | Clade3 | 1.1820 | 0.9842 – 1.4200 | ✘ |
| Clade3 | Clade4 | 1.1975 | 1.0320 – 1.3886 | ✔ |
| Clade2 | Clade3 | 1.3979 | 1.1429 – 1.7111 | ✔ |
| Clade1 | Clade4 | 1.4146 | 1.1954 – 1.6683 | ✔ |
| Clade2 | Clade4 | 1.6740 | 1.3856 – 2.0263 | ✔ |

Regenerated: **the same 8 of 10 resolve, the same two do not** (Clade1–Clade2,
Clade1–Clade3). Whole-file agreement 1.04e-2 absolute / 1.95e-2 relative
(worst column `lo`). **STOCHASTIC-BUT-CONSISTENT.**

### 2c. Spearman ρ (environmental optimum vs fever cost)

| statistic | committed | regenerated |
|---|---|---|
| ρ on the 5 posterior medians | −1.0000 | −1.0000 |
| posterior median ρ printed by 12, run against the committed fits | **−1.00 [−1.00, −0.70]** | −1.00 [−1.00, −0.70] |
| **value printed in the manuscript** | **−0.90 [−1.00, −0.30]** | — |

⚠ **The manuscript's ρ does not match the committed results.**
`_body_main.qmd` line 143 states "Spearman ρ = −0.90, 95 % CrI −1.00 to −0.30".
Running `12_main_figures.R` against the **committed** `results/rds/*.rds` prints
`rho = -1.00 [-1.00, -0.70]`, in two independent runs. The manuscript value is
stale relative to the results it cites. **NOT REPRODUCED** as written.

### 2d. etc-GEM — **does not reproduce**

| clade | R² committed | R² C1 | kcat_scale committed | kcat_scale C1 | C1 / committed |
|---|---:|---:|---:|---:|---:|
| I | 0.955597 | 0.955121 | 0.613644 | **0.833030** | 1.3575 |
| II | 0.942169 | 0.942169 | 0.409658 | **0.560790** | 1.3689 |
| III | 0.965218 | 0.964427 | 0.621805 | **0.851270** | 1.3690 |
| IV | 0.963019 | **0.877050** | 0.710416 | **0.945418** | 1.3308 |

The R² of the fit reproduces for I–III (Δ ≤ 8e-4). Clade IV's drops from 0.963 to
0.877 — the differential evolution converged to a worse optimum, and this is the
clade whose emcee posterior also presses its upper prior bound (§3).

**The capacity axis is rescaled by a near-constant factor.** The anchor ratio is
1.0890 / 0.7955 = **1.3690**; Clades II and III rescale by 1.3689 and 1.3690 — the
anchor ratio to four decimals. Clades I (1.3575) and IV (1.3308) differ from it
only by DE convergence noise. This is not drift — it is a units change. **DIVERGENT.**

### 2e. The 12 isolate capacities

| clade | isolate | committed | regenerated |
|---|---|---:|---:|
| I | Clade1_2068 | 0.670282 | 0.842693 |
| I | Clade1_2069 | 0.641290 | 0.813702 |
| I | Clade1_2070 | 0.699273 | 0.862021 |
| II | Clade2_2071 | 0.510471 | 0.638101 |
| II | Clade2_2072 | 0.471816 | 0.589782 |
| II | Clade2_2073 | 0.365513 | 0.464152 |
| III | Clade3_2074 | 0.678886 | 0.860934 |
| III | Clade3_2075 | 0.688550 | 0.860934 |
| III | Clade3_2076 | 0.669222 | 0.841606 |
| IV | Clade4_2077 | 0.784251 | 0.945418 |
| IV | Clade4_2078 | 0.822907 | 0.993737 |
| IV | Clade4_2079 | 0.735932 | 0.906763 |

**DIVERGENT** (max |Δ| 0.182). Note that in the regenerated file Clade3_2074 and
Clade3_2075 land on the *same* grid point — the 120-point capacity grid
(step 0.0096639) is too coarse to separate them.

Worse, the committed values are not on the same footing as the rest of
`supp_data/`: they are integer multiples of the grid step away from the **backup**
clade kcat_scale, not the current one (§4.3).

### 2f. Downstream of capacity

| statistic | committed | regenerated |
|---|---:|---:|
| ICC1 (between-clade variance) | 0.89894 (90 %) | 0.90862 (91 %) |
| F(3,8) | 27.684 | 30.831 |
| p | 1.414e-04 | 1.00e-04 |
| r(capacity, fever tax at 40 °C) | **−0.8609** (p = 3.24e-04) | **−0.8474** (p = 5.03e-04) |
| r(capacity, respiration 32–36 °C) | **−0.6432** (p = 2.41e-02) | **−0.6201** (p = 3.15e-02) |
| r(capacity, peak growth) [consistency check] | +0.9620 | +0.9626 |

The manuscript's **r = −0.86** and **r = −0.64** are reproduced exactly from the
committed files (−0.8609, −0.6432). Under the regenerated capacities they become
−0.8474 and −0.6201 — the *conclusions* survive the rescaling, because a Pearson
correlation is invariant to a linear rescale of one variable and the rescale is
only approximately linear. **STOCHASTIC-BUT-CONSISTENT in substance, but resting
on DIVERGENT inputs.**

⚠ Their **confidence intervals are not reproducible at all**: `12_main_figures.R`
computes them with `replicate(5000, sample.int(n, n, TRUE))` and **no seed** (§3).

### 2g. 11_capacity_expression.R (GEO GSE165762)

Runs from the **committed** count matrix; confirmed it does **not** download
(`Using committed GEO count matrix (no download)` in the log). Regenerated
result: 921 metabolic enzyme genes tested; median log2FC metabolic **−0.005**
vs genome background **+0.031**; Wilcoxon W = 2 064 624, **p = 0.928**; 49.6 %
of metabolic genes higher in Clade I — i.e. the capacity-as-expression
hypothesis is *not* supported, which is the committed conclusion.

Against `data/expression/capacity_expression_CladeI_vs_II.csv`: same 5 396
genes, all matched, **max |ΔlogFC| = 8.58e-04**. **NUMERICALLY EQUIVALENT** at
the level the conclusion depends on (the residual is edgeR/limma version drift,
not sampling — edgeR is deterministic). The regenerated copy went to
`runs/C1_reproduction/expression/`; `data/` was not written to.

---

## 3. Seeds and determinism

| stage | RNG | seed set? | where | bit-reproducible on re-run? |
|---|---|---|---|---|
| **08 brms — respiration Arrhenius** | Stan | ✔ `BAYES_SEED <- 1234` | `08_bayesian_models.R:132` | **YES — proven.** Every posterior quantity 0.0 different across runs. |
| **08 brms — growth Sharpe-Schoolfield** | Stan | ✔ same seed | same | **YES given identical data.** It differs here only because its input `se_y` differs in the last bits (below). |
| **08 brms — respiration SS (for LOO)** | Stan | ✔ same seed | same | as above |
| **etc-GEM `calibrate` — differential evolution** | numpy via scipy | ✔ `seed = 3` | `generate_model_data.py` | **PARTLY.** Clades II, III, IV came back bit-identical across two runs; Clade I did not — `dTopt = −0.3785` (746 objective evaluations) vs `−0.3578` (778). |
| **etc-GEM `bayes` — emcee sampler** | emcee internal | ✘ **NOT SEEDED** | `emcee.EnsembleSampler(...)` takes no `random_state` | **NO.** Only the *initial walker positions* are seeded (`default_rng(3)`); the chain itself is not. |
| **etc-GEM `bayes` — posterior thinning** | numpy | ✔ `default_rng(7)`, `default_rng(11)` | `stage_bayes` | yes, but it thins a non-reproducible chain |
| **etc-GEM `boot` — isolate bootstrap** | numpy | ✔ `default_rng(3)` | `stage_boot` | yes, given the same DE fit — which is not reproducible |
| **etc-GEM `boot` — correlation bootstrap** | numpy | ✔ same generator | `stage_boot` | as above |
| **10_carbon_tax.R — posterior subsample** | R | ✘ **NOT SEEDED** | `sample.int(nrow(p), 1500)`, `10_carbon_tax.R:217` | **NO — measured, see below** |
| **12_main_figures.R — posterior subsample ×3** | R | ✘ **NOT SEEDED** | `sample.int(nrow(p), 1200)` ×3 | figure curves only; the written tables are unaffected |
| **12_main_figures.R — correlation CI bootstrap** | R | ✘ **NOT SEEDED** | `replicate(B, sample.int(n, n, TRUE))`, `12_main_figures.R:1149` | **NO.** The r = −0.86 / −0.64 CIs printed on Fig 2 are not reproducible. |
| **11_capacity_expression.R** | – | n/a | – | yes (edgeR is deterministic) |

`set.seed()` appears **nowhere** in the analysis scripts.

**Cheap stages run twice, as specified.**

* **10_carbon_tax.R**, run twice against the *same* committed brms fits:
  `carbon_tax_group.csv`, `carbon_tax_isolate.csv` and
  `carbon_tax_contrasts_fever.csv` came back **byte-identical**, but
  **`carbon_tax_curves.csv` differed by up to 0.823 absolute / 4.92 % relative**
  — purely from the unseeded 1500-draw subsample.
* **12_main_figures.R**, run twice the same way: `fig_values.csv` and
  `fig_contrasts.csv` came back **byte-identical**, and — importantly — **also
  byte-identical to the committed versions**. So every table 12 writes is a
  deterministic function of the brms fits. Its unseeded sampling touches only
  the plotted curves.
* **etc-GEM calibrate**, run twice: reproducible for three clades of four.
  `differential_evolution` is seeded, but its objective is an LP solve, and
  degenerate LPs do not return bit-identical vertices. Once one objective value
  differs in the last bits the search path diverges, which is exactly what the
  differing evaluation counts for Clade I (746 vs 778) show. A seed is not
  enough when the objective itself is not a pure function.
* **etc-GEM apriori**, run twice: `apriori_Topt` came back **30/30/30/31** and
  **30/30/32/32** — a 2 °C swing on a quantity that is clade-independent by
  construction (§4.5).

**Whole-pipeline determinism, measured on two complete `run_all.sh` runs:**
39 of the 40 CSVs in `runs/C1_reproduction/tables/` are **byte-identical** between runs —
including every brms output — and the one exception is `carbon_tax_curves.csv`.
The etc-GEM side is not reproducible at all. Full numbers in **§8**.

**Why the growth model moved and the respiration model did not.** `resp_dat`
depends only on `respiration_fgC_h`, which is bit-identical between committed and
C1. `growth_dat` additionally carries `se_y = SE/Estimate` from
`fit_coefficients_long.csv`, whose `SE` column differs by up to **1.04e-16**
(3.3e-13 relative). A seeded Stan chain given data that differs in the last bit
takes a different trajectory. The 1e-2 disagreement in the growth posterior is
the amplification of a 1e-16 input perturbation, not sampler noise.

---

## 4. Known inconsistencies — verified or refuted by direct test

Nothing here is fixed. Test scripts: `reports/tools/known_inconsistencies.py`
and `reports/tools/anchor_and_ceiling_test.py`; JSON alongside them.

### 4.1 Clade IV's `calibrated` curve exceeds its own `apriori` plateau — **CONFIRMED, and explained**

`outputs/supp_data/calibration_curves.csv`:

| clade | apriori peak | calibrated peak | kcat_scale | exceeds? |
|---|---:|---:|---:|---|
| I | 0.795520 @ 31.0 °C | 0.804287 @ 34.5 °C | 0.6136 | yes, by 0.008767 |
| II | 0.795520 @ 31.0 °C | 0.531135 @ 34.0 °C | 0.4097 | no (−0.264385) |
| III | 0.795520 @ 31.0 °C | 0.814734 @ 35.0 °C | 0.6218 | yes, by 0.019214 |
| IV | 0.795520 @ 31.0 °C | **0.931441 @ 36.0 °C** | 0.7104 | **yes, by 0.135921** |

A single build cannot do this — and we did not have to assume it. Rebuilding the
model and evaluating both curves from one build:

| anchor | a-priori peak | calibrated peak I / II / III / IV |
|---|---:|---|
| current, `0.95 @ 307.15 K` | 0.795520 | 0.5875 / 0.3880 / 0.5952 / 0.6804 — **all below** |
| older, `1.1 @ 313.15 K` | **1.089000** | **0.8043 / 0.5311 / 0.8147 / 0.9314** — all below |

The old-anchor calibrated peaks are **exactly** the committed ones. So:

> The `apriori` column of `calibration_curves.csv` was built with the **current**
> anchor (ceiling 0.7955); the `calibrated` column in the **same file** was built
> with the **old** anchor (ceiling 1.0890). Clade IV's 0.9314 never violated any
> ceiling — it is below 1.0890. The file splices two builds.

### 4.2 The anchor — **CONFIRMED**

R² of the committed parameters against the measured clade means, recomputed:

| clade | R² committed | R² under **current** anchor (0.95 @ 307.15 K) | R² under **older** anchor (1.1 @ 313.15 K) | \|Δ\| vs committed, old anchor |
|---|---:|---:|---:|---:|
| I | 0.955597 | **0.0897** | 0.955597 | < 5e-16 |
| II | 0.942169 | **0.4636** | 0.942169 | < 5e-16 |
| III | 0.965218 | **0.1545** | 0.965218 | < 5e-16 |
| IV | 0.963019 | **0.2349** | 0.963019 | < 5e-16 |

RMSE matches to the same precision (0.037627 / 0.034870 / 0.034643 / 0.043100 in
both). The committed `calibration_r2.csv` was produced with `ANCHOR_MU = 1.1 @
313.15 K` — the value the surviving comment in `_supp_common.py` still names
("Was an arbitrary 1.1 at 40 C") — while the code now uses `0.95 @ 307.15 K`.
This single fact explains the whole capacity rescale in §2d.

### 4.3 `capacity_isolates.csv` is from the backup build — **CONFIRMED**

* `sha256` of `supp_data/capacity_isolates.csv` and
  `supp_data_original_backup/capacity_isolates.csv`:
  `907dd418ea2f1c301791328d0f476a5d15ddd72bb7b6cc81ec44555b3fa10f47` —
  **byte-identical**.
* Capacity grid `np.linspace(0.20, 1.35, 120)` → step **0.00966386554**.
* Offsets from the clade kcat_scale, in grid steps:

| clade | isolate | vs **BACKUP** kcat | vs **CURRENT** kcat |
|---|---|---:|---:|
| I | Clade1_2068 | **0** | 5.861 |
| I | Clade1_2069 | **−3** | 2.861 |
| I | Clade1_2070 | **+3** | 8.861 |
| II | Clade2_2071 | **+7** | 10.432 |
| II | Clade2_2072 | **+3** | 6.432 |
| II | Clade2_2073 | **−8** | −4.568 |
| III | Clade3_2074 | **0** | 5.907 |
| III | Clade3_2075 | **+1** | 6.907 |
| III | Clade3_2076 | **−1** | 4.907 |
| IV | Clade4_2077 | **+1** | 7.640 |
| IV | Clade4_2078 | **+5** | 11.640 |
| IV | Clade4_2079 | **−4** | 2.640 |

**12 / 12 integer against the backup; 0 / 12 against the current.**

Does the C1 re-run change the downstream numbers? Yes, but modestly:

| statistic | committed | regenerated |
|---|---:|---:|
| ICC1 | 0.89894 | 0.90862 |
| F(3,8) | 27.684 | 30.831 |
| p | 1.414e-04 | 1.00e-04 |
| r(capacity, fever tax 40 °C) | −0.8609 | −0.8474 |
| r(capacity, respiration) | −0.6432 | −0.6201 |

Which files in `supp_data/` are byte-identical to the backup:

| byte-identical to backup | differs from backup |
|---|---|
| `calibration_measured.csv` | `bayes_posterior_samples.csv` |
| `capacity_bootstrap_draws.csv` | `bayes_predictive.csv` |
| `capacity_corr.csv` | `bayes_summary.csv` |
| `capacity_isolates.csv` | `calibration_curves.csv`, `calibration_r2.csv` |
| | `capacity_bootstrap.csv` |
| | `identifiability_optima.csv`, `identifiability_profiles.csv` |
| | `predicted_vs_observed.csv` |

`supp_data/` is a **mixture of two builds**. In particular
`capacity_bootstrap.csv` (new build; `capacity` = 0.6136…) and
`capacity_isolates.csv` (old build; isolates centred on 0.6703…) are compared
against each other in Fig 3 and in the supplement.

### 4.4 `figure3_data/capacity.csv` vs `calibration_r2.csv` — **CONFIRMED**

| clade | fig3 `capacity` | supp `kcat_scale` | Δ | backup `kcat_scale` | Δ vs backup |
|---|---:|---:|---:|---:|---:|
| I | 0.6702817373 | 0.6136443845 | 5.66e-02 | 0.6702817373 | **0** |
| II | 0.4428240312 | 0.4096579460 | 3.32e-02 | 0.4428240312 | **0** |
| III | 0.6788860349 | 0.6218047432 | 5.71e-02 | 0.6788860349 | **0** |
| IV | 0.7745872870 | 0.7104159095 | 6.42e-02 | 0.7745872870 | **0** |

`figure3_data/capacity.csv` is bit-for-bit the **backup** kcat_scale. Its `r2`
column likewise matches the backup `r2_calibrated` exactly and the current one
not at all. Same disagreement, same cause: two builds.

### 4.5 `apriori_Topt` differs across clades — **CONFIRMED**

The a-priori curve is clade-independent by construction (`curve(pm, [0,1,1], tf)`
takes no clade argument).

| | `supp_data` | `supp_data_original_backup` |
|---|---|---|
| max across-clade \|difference\| at any T | **0.000e+00** (bit-identical) | 9.07e-14 |
| a-priori peak | 0.7955200481 at 31.0 °C | 0.8518999712 at 30.0 °C |
| exact float ties at the peak | **6** (31.0 – 37.0 °C) | 1 |
| within 1e-6 of the peak | **30.0 – 40.0 °C (21 of 57 points)** | 30.0 – 46.0 °C (33 of 57) |
| within 5 % of the peak | 29.5 – 40.5 °C (23 of 57) | 29.5 – 46.0 °C (34 of 57) |
| `np.argmax` on that array must give | **31.0 °C for every clade** | 30.0 °C for every clade |
| `calibration_r2.csv` actually reports | **40 / 30 / 30 / 40** | 30 / 30 / 30 / 30 |
| consistent? | **NO** | yes |

Two findings.

1. **The plateau is enormous.** The a-priori curve is flat to within 1e-6 over
   **30–40 °C, 21 of its 57 sampled temperatures**. `argmax` on a plateau is
   meaningless: `apriori_Topt` is not an estimate of anything.
2. **The reported values are impossible from the file they sit beside.** Within
   `supp_data/calibration_curves.csv` the a-priori column is *bit-identical*
   across all four clades, so `np.argmax` must return 31.0 °C four times. It
   reports 40/30/30/40. The backup is internally consistent (30 everywhere);
   the current `supp_data/` is not.

The mechanism is directly demonstrable: in a **single fresh build** under the
current anchor, evaluating the same clade-independent a-priori curve four times
gave `apriori_Topt` = **32.5 / 32.5 / 30.0 / 32.5 °C**, and the C1 pipeline run
wrote **30 / 30 / 30 / 31**. The LP is degenerate on the plateau, so the solver
returns vertices that differ in the last bits and `argmax` lands wherever the
noise puts it.

### 4.6 `drop_excluded()` / `EXCLUDE_GROUPS = "glab"` — **CONFIRMED**

`drop_excluded()` is defined once in `config.R` and called in exactly **two**
places: `08_bayesian_models.R` (×1) and `09_bayesian_plots.R` (×1). Scripts
02, 03, 07, 10, 11, 12 and 13 never call it.

Consequence, by row count in the committed `results/tables/`:

| file | rows | glabrata rows | % |
|---|---:|---:|---:|
| `derived_N0_R_results_with_carbon.csv` | 991 | **168** | 17.0 |
| `fit_coefficients_long.csv` | 3240 | **540** | 16.7 |
| `fit_coefficients_wide.csv` | 1063 | **178** | 16.7 |
| `fit_metrics.csv` | 1080 | **180** | 16.7 |
| `Oxygen_Trimmed_Series_Metadata.csv` | 1080 | **180** | 16.7 |
| `Oxygen_Curve_Code_Key.csv` | 1080 | **180** | 16.7 |
| `group_lookup_with_delta.csv` | 1080 | **180** | 16.7 |
| `replication_summary_by_temperature.csv` | 204 | **34** | 16.7 |
| `activation_energy_summary.csv` | 18 | **3** | 16.7 |
| `arrhenius_growth_fgC_h_coefs.csv` | 18 | **3** | 16.7 |
| `arrhenius_respiration_fgC_h_coefs.csv` | 18 | **3** | 16.7 |
| `cue_quadratic_fit_coefs.csv` | 18 | **3** | 16.7 |
| `sharpe_schoolfield_*_coefs.csv` (×3) | 72 each | **12** each | 16.7 |
| `diagnostic_hot_temps_growth_check.csv` | 111 | **20** | 18.0 |
| `summary_hot_temps_growth_check.csv` | 26 | **4** | 15.4 |
| `failed_fits_at_bound.csv` | 17 | **2** | 11.8 |
| `removed_points.csv` | 72 | **10** | 13.9 |
| — Bayesian outputs — | | | |
| `bayes_clade_params.csv` | 5 | **0** | 0 |
| `bayes_cue_by_clade.csv` | 600 | **0** | 0 |
| `bayes_cue_by_isolate.csv` | 1800 | **0** | 0 |
| `bayes_E_resp_minus_growth.csv` | 5 | **0** | 0 |
| `carbon_tax_group.csv` / `_isolate.csv` / `_curves.csv` | 20 / 30 / 445 | **0** | 0 |
| `fig_values.csv` | 5 | **0** | 0 |

`08_bayesian_models.R` logged `EXCLUDE_GROUPS: dropped 168 of 991 derived rows
(glab; OTUs 13, 14, 15)` in the C1 run. So **every frequentist table and every
figure produced by 07 contains *C. glabrata*, and nothing downstream of 08
does.** That is the documented intent, but the two halves of `results/tables/`
describe different datasets, and nothing in the file names says so.

### 4.7 `capacity_bootstrap_draws.csv` pairs capacity and peak independently — **CONFIRMED (bug is real; the committed file predates it)**

Source (`stage_boot`):

```python
for c, pk in zip(*[caps[rng.choice(B, 1500, replace=False)],
                   peaks[rng.choice(B, 1500, replace=False)]]):
```

Two **independent** resamples of the same 4000 bootstrap replicates, so draw *i*'s
capacity and draw *i*'s peak come from different bootstrap samples. The pairing
is destroyed.

Within-clade `corr(capacity, peak)`, 1500 draws per clade:

| clade | **committed** file | **regenerated** run 1 | **regenerated** run 2 |
|---|---:|---:|---:|
| I | +0.1552 | −0.0250 | −0.0250 |
| II | **+0.9234** | −0.0471 | −0.0471 |
| III | −0.1600 | −0.0101 | −0.0101 |
| IV | −0.2666 | −0.0085 | −0.0085 |

The current code produces ≈ 0, as independent pairing must (SE ≈ 1/√1500 =
0.026). The committed file does not — Clade II's +0.923 is 35 SE away. So:

* the **bug is real and live** in the source, confirmed by the re-run;
* the **committed file was not produced by this code** — and indeed it is one of
  the four files byte-identical to `supp_data_original_backup/` (§4.3).

The pooled cross-clade correlation is +0.978 in the committed file and is driven
entirely by the between-clade spread; it is not evidence of within-clade pairing.

---

## 5. Blockers

### 5a. What had to change for the pipeline to run at all

Execution only. No analysis decision was touched. Nothing on the CONSTRAINTS
list was modified.

| # | change | why |
|---|---|---|
| 1 | Headless guards on 01, 04, 05, 06 | `run_all.R` source()d 06, whose bottom line was `if (interactive()) shinyApp(...) else runApp(...)`. Under `Rscript`, `interactive()` is FALSE, so the master runner started a Shiny server and **blocked forever**. This is the single reason the repository could not be run unattended. |
| 2 | `--file=` script-directory resolution in 13 files | `sys.frame(1)$ofile` is NULL under `Rscript`, so every script fell back to `getwd()` and died on `cannot open file .../config.R`. No script could be run on its own. |
| 3 | Un-mangling `~+~` in `--file=` | R encodes spaces in `--file=` as `~+~`. This project lives under `.../Candidas TPC/Candidas`, so absolute invocations produced `Candidas~+~TPC` and every `source()` failed. |
| 4 | `RESULTS_DIR` from `CANDIDAS_RESULTS` in `config.R` | Required by the brief, to run non-destructively. `models_dir` follows it, so brms genuinely refits. |
| 5 | `app_input()` in `config.R` | `results/tables/` mixes generated output with **five hand-made app decisions** that every script reads as data. Redirecting the results tree alone would have silently run the pipeline with no manual trim windows, no exclusions, the global fallback cell size and no isolate names — a different analysis. |
| 6 | `resolve_supp_data()` / `require_csv()` in `config.R` and 13 | The old `cand[which(vapply(cand, dir.exists, logical(1)))[1]]` returned `NA` when nothing matched and, worse, **accepted a directory that existed but was empty** — the unpopulated-submodule state this repo shipped in. Now named `ETCGEM_SUPP_DATA_MISSING` / `_EMPTY` / `ETCGEM_FILE_MISSING` errors at the point of resolution. |
| 7 | 11: removed the hard-coded `/Users/ilgazcakin/...` path | Unrunnable anywhere but one machine. |
| 8 | 11: removed the mid-script `setwd()` | It made every relative path depend on where R was started and left the caller's working directory changed after `run_all.R` sourced it. |
| 9 | 11/13: `CANDIDAS_EXPRESSION_OUT` | 11 wrote its outputs **into `data/expression/`**. A non-destructive run must not write into the data tree. |
| 10 | 12/13/14: figure output follows `CANDIDAS_RESULTS` | Both hard-coded `.../results/figures/manuscript`. |
| 11 | `run_all.R` runs 11, 13 and 14 | They were only *mentioned* in a closing message and never actually run. |
| 12 | `run_all.R` globals renamed `.RUNALL_*` | `08_bayesian_models.R` defines `.TIMINGS` and `fmt_dur` at top level and `source(local = FALSE)` let it clobber the runner's, so the timing table came out empty. |
| 13 | `generate_model_data.py`: `--out-suffix` / `--out-supp` / `--out-fig3` | Required by the brief. The DE fit cache `_percladefit.npy` follows the suffix, so a re-run genuinely re-calibrates. |
| 14 | `probe.py`: Kelvin → Celsius | It passed `303.15`/`313.15` to `compute_tpc`, which takes **Celsius** (`_supp_common.curve()` feeds it `temp_C` 22–44). It was probing the model at ~303 °C and ~313 °C. |
| 15 | `identifiability_stageA.py`: `DERIVED` path | Pointed at `{HERE}/../../tables/`, i.e. `cauris_etcgem/tables`, which has never existed. The file is in the parent project under `results/tables/`. |
| 16 | `identifiability_stageA.py`: `temp_C` → `T` | It grouped on `temp_C`; that CSV names the column `T`. The `groupby` raised `KeyError`, so `measured()` could never return. |

Also added, not fixes: `scripts/00_install.R`, `scripts/run_all.sh`,
`SETUP.md`, `env/versions.json`, `renv.lock`,
`cauris_etcgem/requirements.lock`, `reports/tools/*`, and per-stage timing +
emcee diagnostics (acceptance fraction, autocorrelation time, raw chains) in
`generate_model_data.py`.

### 5b. What still cannot be run

| item | status | reason |
|---|---|---|
| `01_convert_xlsx.R` | **intentionally not run** | Click-driven. Owns `data/*_Oxygen.csv` and `otu_names.csv`, both committed inputs. Re-running needs a human to retype 18 isolate names and 48 plate temperatures. |
| `04_trim_selector.R` | **intentionally not run** | Click-driven. Owns 1080 manual fit windows and 88 point exclusions — irreplaceable hand judgements. |
| `05_cell_sizes.R` | **intentionally not run** | Click-driven. Owns the 18 per-isolate cell volumes. |
| `identifiability_stageA.py` | **NOT RUNNABLE** | Its two genuine path/column bugs are fixed, so the growth side executes — but `respiration_hook()` is a `return np.full(len(temps), np.nan)` **placeholder**. Every `resp_heldout_logRMSE` is NaN, so the script's entire purpose (rank mechanisms M_E / M_M / M_U by held-out respiration) cannot be served. The stub was **deliberately left alone**: writing an O₂ read-out would be inventing an analysis. Two further mechanisms (M_C coupling, M_R respiratory cap) are `TODO` and not implemented at all. |
| `probe.py` | runs | Diagnostic only; writes nothing. |

No script was silently skipped. All 14 numbered scripts are accounted for: **11
ran** (02, 03, 06, 07, 08, 09, 10, 11, 12, 13, 14), **3 are documented
click-driven inputs** (01, 04, 05).

### 5c. Things that would have to change to make the numbers reproduce — **NOT done here**

Reporting only, per the constraints:

1. **The anchor.** `ANCHOR_MU`/`ANCHOR_T` is on the do-not-touch list, and it is
   also the entire explanation for §2d. Either the committed etc-GEM outputs must
   be regenerated under the current anchor, or the anchor must be reverted — but
   `supp_data/` cannot be left as a splice of both. **Blocker for C2+.**
2. **Four unseeded stochastic stages** (§3). Adding `set.seed()` to 10 and 12,
   and a `random_state` to `emcee.EnsembleSampler`, would make three of the four
   reproducible; the DE stage additionally needs a bit-stable objective.
3. **`figure3_data/`, `capacity_isolates.csv`, `capacity_bootstrap_draws.csv`
   and `capacity_corr.csv`** must be regenerated from the same build as the
   files they are plotted against.
4. **The manuscript's Spearman ρ** (§2c) disagrees with the committed results.

---

## 6. Verification

| # | requirement | result |
|---|---|---|
| 1 | `bash scripts/run_all.sh` completes on a clean checkout, zero interaction, per-stage timings, logs written | **PASS, with one honest qualification** — it completed in **48 m 50 s**, zero interaction, per-stage timings printed, four logs written (§7). `runs/C1_reproduction/`, `runs/C1/supp_data/`, `runs/C1/figure3_data/` and the DE fit cache were **deleted before the run**, so every output was regenerated from the committed raw data. It was run on this working tree, **not on a fresh `git clone`**: the R library and the Python venv were already installed. `scripts/00_install.R` was written and exercised, and every version it pins is recorded in `env/versions.json`, but a genuinely from-zero clone-and-install was not re-timed end to end. |
| 2 | Every one of the 14 scripts ran or is documented as intentionally not run | **PASS** — §5b |
| 3 | `runs/C1_reproduction/` and `runs/C1/supp_data/` fully populated; `results/` and `outputs/supp_data/` byte-identical to their pre-run state | **PASS** — `shasum -a 256 -c env/baseline_checksums_{results,etcgem,data}.txt` all OK; `git status` shows no modification to `data/`, `results/` or the committed etc-GEM outputs |
| 4 | Manuscript renders to .docx and .pdf with no missing figure and no unresolved reference | **PASS** — rendered twice (committed figures, and the regenerated C1 figures via a staging copy). 6 documents each time; 7/7 referenced figures resolve; 0 unresolved citations in any of the three `.tex`; no LaTeX undefined-reference warning; no quarto WARN/ERROR |
| 5 | `reports/REPRODUCTION.md` with all six sections and actual numbers | **PASS** — this file: §1 file-by-file comparison, §2 headline numbers, §3 seeds and determinism, §4 known inconsistencies, §5 blockers, plus `reports/RUNBOOK.md`. Every committed load-bearing number is classified. |
| 6 | `env/versions.json`, `SETUP.md`, `reports/RUNBOOK.md`, `requirements.lock`, `renv.lock` committed | **PASS** |
| 7 | Plain statement at the top | see **VERDICT** |

---

## 7. Wall time and logs

One clean end-to-end `bash scripts/run_all.sh` (through `run_c1.sh`, which sets
the non-destructive path variables), started 09:47:44, finished 10:36:34 on
2026-07-26. Zero user interaction. Logs: `logs/run_all_20260726_094744.log`,
`logs/etcgem_20260726_094744.log`, `logs/r_pipeline_20260726_094744.log`,
`logs/quarto_20260726_094744.log`.

### Total wall time: **48 min 50 s**

| stage | wall time | exit |
|---|---:|---:|
| 1/3 etc-GEM (python, all stages) | **35 m 29 s** | 0 |
| 2/3 R pipeline (02 → 14) | **12 m 37 s** | 0 |
| 3/3 quarto render (3 documents × docx + pdf) | **44 s** | 0 |
| **TOTAL** | **48 m 50 s** | 0 |

### etc-GEM, per stage

| stage | wall time |
|---|---:|
| model build (SBML → GECKO → sectors → budget bisection) | 2.5 s |
| `calibrate` (differential evolution × 4 clades) | **16 m 00 s** |
| `curves` | < 0.1 s |
| `fit` | < 0.1 s |
| `ident` (3 knobs × 25 values × 4 clades) | 61.7 s |
| `boot` (120-point grid, 4 000 bootstrap replicates) | 92.6 s |
| `apriori` | 8.6 s |
| `bayes` (15³ emulator grid, then 4 × 32 walkers × 4 000 steps) | **16 m 26 s** |

### R pipeline, per script

As reported by `run_all.R`'s own timing table (final clean run):

| script | wall time |
|---|---:|
| 02_longdata.R | 4.2 s |
| 03_trimming.R | 1 m 55 s |
| 06_inoculation.R | 0.3 s (app defined, not launched) |
| 07_oxygen_fits.R | 2 m 10 s |
| **08_bayesian_models.R** | **8 m 00 s** |
| — growth Sharpe-Schoolfield | 4 m 58 s |
| — respiration Arrhenius | 36 s |
| — respiration SS (for LOO) | 2 m 18 s |
| 09_bayesian_plots.R | 10.4 s |
| 10_carbon_tax.R | 1.2 s |
| 11_capacity_expression.R | 0.7 s |
| 12_main_figures.R | 6.0 s |
| 13_supplementary_figures.R | 1.8 s |
| 14_schematic.py | 2.6 s |
| **TOTAL (as printed by run_all.R)** | **12 m 33 s** |

(The 12 m 37 s in the stage table above is the same work plus R start-up and
the `Rscript` process wrapper.)

`SETUP.md` quotes 1–3 h for 08 as a safe upper bound for slower machines; on
this 16-core box it is 8 minutes.

### etc-GEM sampler configuration, as actually used

```
differential_evolution: seed=3 maxiter=16 popsize=10 polish=True
                        bounds=[(-12, 6), (0.2, 9), (0.2, 4)]
emcee EnsembleSampler : walkers=32 steps=4000 discard=1500 thin=8 ndim=3
                        seed(initial positions)=3, sampler NOT seeded
priors/bounds         : dTopt (-11.0, 1.0)  dCp_scale (0.5, 6.5)  kcat_scale (0.25, 1.05)
emulator grid         : 15 x 15 x 15 over 12 temperatures -> G.shape (15,15,15,12)
                        = 40 500 model evaluations
temperatures          : 22 24 26 28 30 32 34 36 38 40 42 44 degC
```

Convergence (from the final run):

| clade | mean acceptance | min / max acceptance | autocorrelation time τ (dTopt, dCp, kcat) |
|---|---:|---|---|
| I | 0.644 | 0.624 / 0.664 | 36.9, 37.6, 37.3 |
| II | 0.645 | 0.631 / 0.661 | 35.5, 37.1, 37.9 |
| III | 0.636 | 0.624 / 0.654 | 39.0, 38.9, 35.0 |
| IV | **0.588** | **0.304** / 0.623 | **53.6, 45.9, 40.4** |

`get_autocorr_time()` **returned without raising** for all four clades: 4 000
steps is 75–115 τ, comfortably past emcee's 50 τ threshold. Clade IV is the weak
one — one of its 32 walkers is effectively stuck (acceptance 0.304), its τ is
~45 % longer, its posterior median `kcat_scale` = 0.953 with an upper bound of
1.0397 against a **prior upper bound of 1.05**, and it is the clade whose DE fit
also degrades (R² 0.877 vs 0.94–0.96). Its capacity estimate is pressed against
the edge of its prior and should not be read as well-identified.

The full un-thinned chains, with walker and step indices, per-walker acceptance
fractions and τ, are in
`cauris_etcgem/runs/C1/supp_data/bayes_chain_raw.npz`
(`chain_<clade>`: shape `(step = 4000, walker = 32, param = 3)`), for the
convergence work in the C-series.

### brms, as actually used

`iter = 4000`, `warmup = 1000`, `chains = 4`, `seed = 1234`,
`adapt_delta = 0.99` (growth SS) / `0.95` (Arrhenius, resp SS),
`max_treedepth = 12`, `family = student()`, measurement error on growth
(`y | se(se_y, sigma = TRUE)`). Backend **rstan** (08 never sets `backend`, so
brms uses its default; cmdstanr is not installed).

| model | max R̂ | min ESS bulk | min ESS tail | divergences |
|---|---:|---:|---:|---|
| growth Sharpe-Schoolfield | 1.00297 | 2601 | 2823 | none reported |
| respiration Arrhenius | 1.00267 | 2933 | 5012 | none reported |

Data actually fitted: **823 growth points**, 15 isolates, 5 groups (no-growth
curves dropped; `EXCLUDE_GROUPS` dropped 168 of 991 derived rows — glabrata,
OTUs 13–15). Plot/SS point exclusions dropped 0 of 991 rows.

### Manuscript render

Rendered twice, deliberately:

1. **Against the committed figures**, by `run_all.sh` → `manuscript/draft/_output/`
   (44 s).
2. **Against the regenerated C1 figures**, by
   `reports/tools/render_with_C1_figures.sh` →
   `runs/C1_reproduction/manuscript_build/manuscript/draft/_output/`.

**Which method:** a **staging copy**. `_body_main.qmd` and `_body_supp.qmd`
reference figures as `../../results/figures/manuscript/…`, relative to
`manuscript/draft`. The script copies the Quarto project (identical `.qmd`
bytes — nothing in `manuscript/draft/` is edited) into
`runs/C1_reproduction/manuscript_build/manuscript/draft/`, whose sibling
`../../results/figures/manuscript` is a **symlink to
`runs/C1_reproduction/figures/manuscript`**. Same document, C1 figures. It also verifies
every referenced figure resolves *before* invoking quarto, because quarto
otherwise downgrades a missing image to a warning in docx and only fails in PDF.

Both renders produced all six outputs
(`manuscript`, `supplementary`, `manuscript_combined` × `.docx` + `.pdf`).

| check | result |
|---|---|
| figures referenced | 7 (FIG1_decoupling, FIG2_the_bill, FIG_model_schematic, FIG_MODEL, FIG_MODEL_SUPP, FIG_MODEL_SUPP_validation, FIG_MODEL_SUPP_consistency) |
| figures resolved | **7 / 7** |
| unresolved citations (`[?]` / `???` in the generated `.tex`) | **0** in all three documents |
| LaTeX undefined reference / citation warnings | **none** |
| quarto WARN / ERROR lines | **none** |

---

## 8. Run-to-run determinism, measured

Two complete, independent `run_all.sh` runs on the same machine with the same
code (09:47 and the earlier 08:49 run).

### R pipeline: bit-reproducible, with one exception

Of the 40 CSVs written to `runs/C1_reproduction/tables/`, **39 are byte-identical between
the two runs** — including every brms output (`bayes_growth_ss_summary.csv`,
`bayes_clade_params.csv`, `bayes_clade_contrasts.csv`, `bayes_cue_*`,
`fig_values.csv`, `fig_contrasts.csv`).

The exception is **`carbon_tax_curves.csv`**, which differed by up to
**0.205 absolute / 2.21 % relative**, from the unseeded 1500-draw subsample in
`10_carbon_tax.R:217`.

This is the decisive determinism result: **given identical inputs, the R half of
this pipeline is bit-reproducible.** The committed-vs-C1 differences in §2a are
therefore *version* differences, not sampling noise.

### etc-GEM: not reproducible

| file | max abs Δ between two runs | max rel Δ |
|---|---:|---:|
| `calibration_measured.csv` | 0 | 0 |
| `capacity_bootstrap.csv` | 1.17e-04 | 1.44e-04 |
| `capacity_isolates.csv` | 1.17e-04 | 1.44e-04 |
| `capacity_bootstrap_draws.csv` | 1.17e-04 | 1.44e-04 |
| `capacity_corr.csv` | 2.89e-05 | 3.16e-05 |
| `identifiability_profiles.csv` | 3.95e-04 | 4.32e-03 |
| `calibration_curves.csv` | 1.81e-03 | 5.18e-03 |
| `predicted_vs_observed.csv` | 1.81e-03 | 1.60e-01 |
| `identifiability_optima.csv` | 2.06e-02 | 5.45e-02 |
| `bayes_predictive.csv` | 2.04e-02 | 2.05e-01 |
| `bayes_summary.csv` | 3.32e-02 | 1.28e-01 |
| **`calibration_r2.csv`** | **2.00e+00** | 6.25e-02 |
| `bayes_posterior_samples.csv` | 2.32e+00 | 2.00e+00 |

Two things stand out.

* **`calibration_r2.csv` differs by 2.0 — in `apriori_Topt`.** Two runs of
  identical code on identical data gave

  | clade | run 1 | run 2 |
  |---|---:|---:|
  | I | 30.0 | 30.0 |
  | II | 30.0 | 30.0 |
  | III | **30.0** | **32.0** |
  | IV | **31.0** | **32.0** |

  This is §4.5 reproduced live: `argmax` on a flat plateau, decided by
  last-bit LP noise. `apriori_Topt` is not a measurement.

* **Differential evolution is seeded but not reproducible.** Clades II, III and
  IV came back bit-identical (`nfev` 550/550/546 both times); **Clade I did
  not** — `dTopt = −0.3785` (746 evaluations) versus `−0.3578` (778). The seed
  fixes the trial sequence, but the objective is an LP solve and degenerate LPs
  do not return bit-identical vertices; once one objective value differs in the
  last bits the search path diverges.

* **emcee is not seeded at all**, so `bayes_summary.csv` medians move by ~1e-3
  (e.g. Clade II `kcat_scale` 0.5354 → 0.5341) and individual posterior draws by
  up to 2.3 between runs.
