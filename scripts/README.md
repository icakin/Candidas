# Candida Oxygen Pipeline (temperature × isolate)

## Overview

This project analyses how temperature modulates the respiration of several
**Candida isolates**. Raw data are dissolved-oxygen time-series recorded on a
PreSens 24-well SDR plate at multiple temperatures. There is no pH or dose
treatment.

The isolates belong to six clade/species **groups**:

- `Clade1`–`Clade4` — *Candida auris*, four clades from four geographical regions
- `glab` — *Candida glabrata*
- `para` — *Candida parapsilosis*

- **Temperature**: 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44 °C (12 levels).
- **Isolate (OTU)**: a unique integer code `1`…`18` (3 isolates per group × 6
  groups). Each code carries a human-readable display name you set in Step 1.
- **Replicate**: `R1`, `R2`, …, `R5`.

### Plate layout — the key fact

Each raw `.xlsx` is **one plate = one group at one temperature**. On every plate:

- rows **A, B, C** are **3 different isolates** of that group,
- columns **1–5** are the **5 replicate wells** per isolate,
- row **D** and column **6** are empty (`No Sensor`) and are ignored.

So one file yields 3 isolates × 5 replicates = 15 oxygen series. The file name
(`<Group>_<Temp>_Oxygen.xlsx`, e.g. `Clade1_22_Oxygen.xlsx`) gives the group and
temperature; the row gives the isolate; the column gives the replicate.

Isolate codes are assigned in a fixed group order (see `STRAIN_GROUPS` in
`config.R`): `Clade1` rows A/B/C → 1/2/3, `Clade2` → 4/5/6, …, `para` → 16/17/18,
so every isolate across the whole dataset has a unique code.

## Step 1: convert the raw Excel files

`01_convert_xlsx.R` is an interactive Shiny app that turns the PreSens exports
into the tidy CSVs the pipeline expects. In the app you:

- tick the files to convert (all recognised groups are pre-ticked);
- type an **isolate name** for rows A/B/C of each group (free text — the
  strain/isolate ID). The CSV headers use safe codes `OTU1`…`OTU18`; your names
  are written to `tables/otu_names.csv` and shown in every result table and plot;
- set the **temperature** per file (auto-filled from the plate's measured `Tm`
  column; editable).

Run it in RStudio (open `01_convert_xlsx.R` → **Run App**) or with
`Rscript scripts/01_convert_xlsx.R`. It writes one `<name>_Oxygen.csv` per file
into `data/`, with columns `Time, T, OTU<c>_R1, …, OTU<c>_R5` for the plate's 3
isolates, plus `tables/otu_names.csv`.

## Script Map

| Script | What it does | Key inputs | Key outputs |
|---|---|---|---|
| `01_convert_xlsx.R` | Interactive Shiny app: raw PreSens `.xlsx` → tidy CSV; sets isolate display names & codes | `data/*.xlsx` | `data/*_Oxygen.csv`, `tables/otu_names.csv` |
| `05_cell_sizes.R` | OPTIONAL Shiny app: average cell volume per isolate → carbon per cell | (typed in) | `tables/otu_cell_sizes.csv` |
| `config.R` | Shared paths, constants, isolate registry, helpers | - | sourced by all |
| `02_longdata.R` | Wide → long; splits `OTU<c>_R<r>` into `OTU` + `Replicate` | `data/*_Oxygen.csv` | `tables/Oxygen_All_Long.csv` |
| `03_trimming.R` | Spline-based trimming with optional manual overrides | `tables/Oxygen_All_Long.csv` | `tables/Oxygen_Data_Filtered.csv`, trimming-diagnostics PDF |
| `04_trim_selector.R` | OPTIONAL Shiny app to review/override each curve's fit window | `03_trimming.R` outputs | `tables/manual_fit_windows.csv`, `tables/plot_exclude_points.csv` |
| `07_oxygen_fits.R` | nlsLM fits per series; carbon unit conversions; isolate-coloured plots; derives growth, respiration, CUE; **frequentist Sharpe-Schoolfield + Boltzmann-Arrhenius thermal fits per isolate** | filtered data + trim metadata | `tables/derived_N0_R_results_with_carbon.csv`, `tables/sharpe_schoolfield_*_coefs.csv`, `tables/arrhenius_*_coefs.csv`, `tables/activation_energy_summary.csv`, `tables/cue_quadratic_fit_coefs.csv`, overlay + CUE + `activation_energy_by_otu` PNGs |
| `08_auto_trim.R` | OPTIONAL automatic alternative to `04`: for each raw curve, fits the exponential over a family of candidate windows, picks the window with the cleanest/most-stable `r` (gated toward the isolate's trend), and discards curves no window can rescue. Writes the same files `04` does, so `06` uses them | `tables/Oxygen_Data_Filtered.csv` (from `03`) | `tables/manual_fit_windows.csv`, `tables/plot_exclude_points.csv`, `tables/auto_trim_log.csv` |
| `09_rmse_trim.R` | OPTIONAL, simpler alternative to `08`. **Windows are NOT changed.** The only lever is each curve's raw-fit RMSE: it scans RMSE cutoffs, fits the growth SS at each, and stops at the **knee** (diminishing returns), then discards any points still problematic. Guardrails cap how much can be dropped | `derived_N0_R_results_with_carbon.csv` + `fit_metrics.csv` (from `06`) | `tables/plot_exclude_points.csv`, `tables/rmse_trim_discarded.csv`, `tables/rmse_trim_log.csv`, `tables/rmse_trim_growth_ss.csv`, `figures/rmse_trim_diagnostics.pdf` |
| `07_tpc_refine.R` | OPTIONAL robust TPC refinement: growth = Sharpe-Schoolfield; respiration = best of Arrhenius / exponential / quadratic / SS by AICc. Robust (Tukey IRLS) fits **flag which points to trim to get the best model**, and mark them on per-isolate diagnostic plots | `tables/derived_N0_R_results_with_carbon.csv` | `tables/tpc_growth_ss_coefs.csv`, `tables/tpc_respiration_model_selection.csv`, `tables/tpc_respiration_aicc_table.csv`, `tables/tpc_flagged_points.csv`, `tables/plot_exclude_points_suggested.csv`, `figures/TPC_refined_growth.pdf` / `_respiration.pdf` (+ PNGs) |

## Schema

After `02_longdata.R`, the long-format table has columns:

```
File, Time, T, OTU, Replicate, Oxygen
```

`OTU` is the integer isolate code (1…18); its display name lives in
`tables/otu_names.csv`. Each series is uniquely identified by `(T, OTU, Replicate)`.

After `07_oxygen_fits.R`, `tables/derived_N0_R_results_with_carbon.csv` has one
row per `(T, OTU, Replicate)` combination, with the derived quantities
`growth_fgC_h`, `respiration_fgC_h`, `CUE`, `resp_over_growth`, plus their
biomass-corrected analogues `growth_C_per_C_h` and `respiration_C_per_C_h`.

## How to Run

Once the CSVs exist (Step 1), run the full pipeline:

```r
source("scripts/run_all.R")
```

Step by step:

```r
source("scripts/02_longdata.R")
source("scripts/03_trimming.R")
source("scripts/07_oxygen_fits.R")
```

## What to Edit Per Project

- `base_dir` in `config.R` is **auto-detected** as the parent of `scripts/`, so
  moving/renaming the project folder needs no edit. Override via `base_dir_manual`
  only if auto-detection fails.
- **Isolate registry** in `config.R`: `STRAIN_GROUPS` (the six group prefixes, in
  the order that fixes the isolate codes), `GROUP_SPECIES`, `PLATE_ROWS`,
  `N_REPLICATES`. Edit these only if the plate design changes.
- **Cell-carbon constants** in `config.R`: `CELL_WIDTH_UM`, `CELL_LENGTH_UM`,
  `CARBON_DENSITY_FG_PER_UM3`, `N_inoculation_cells_per_L`. These are Candida
  (yeast) values and set the carbon per cell, which scales growth (fg C/h) and
  CUE. Use `05_cell_sizes.R` for per-isolate sizes.
- `EXCLUDE_TEMPS` / `EXCLUDE_OTUS` / `EXCLUDE_T_OTU` in `config.R` to drop
  temperatures, isolates (by code 1…18), or specific (temperature, isolate) groups.
- `MANUAL_FIT_WINDOWS` / `PLOT_EXCLUDE_POINTS` in `config.R` (empty by default)
  after inspecting `figures/oxygen_trimming_diagnostics.pdf`.
  `USE_APP_TRIM_FILES` is `TRUE`, so `07_oxygen_fits.R` fits the exact
  windows/exclusions you set in `04_trim_selector.R` (auto-saved to `tables/`).
  Set it to `FALSE` to ignore the app and use automatic trimming for every curve.

> Thermal-performance fitting is **frequentist only** (nlsLM Sharpe-Schoolfield
> and lm-based Boltzmann-Arrhenius, one curve per isolate, in `07_oxygen_fits.R`).
> There is no Bayesian / MCMC stage in this build.

---

## Scripts NOT part of the current manuscript (etc-GEM model work)

The enzyme- and temperature-constrained genome-scale model (etc-GEM) analysis has
been **removed from the manuscript**. An audit found that the model was built once
from a single proteome and reused for every clade, so its "genome-only predictions
are identical across clades" result was a property of the build rather than a
finding; that the fitted `kcat_scale` parameter scales predicted growth exactly
linearly and correlates 0.996 with measured peak growth, so it restates the data it
was fitted to; that the proteome pool budget is calibrated rather than fixed, making
it non-identifiable against `kcat_scale`; that one clade's measured peak growth
exceeds the model's feasible maximum; and that a three-parameter empirical curve
reproduces the fit equally well. A flux-variability analysis had already shown that
oxygen consumption at fixed growth is not uniquely determined, so the model permits
rather than predicts respiration.

The code is retained here for a separate methods paper, where the per-clade build,
the `kcat_scale` x pool identifiability, the binding-constraint analysis and an
out-of-sample respiration prediction can be done properly.

Not part of the current manuscript pipeline:

- `13_capacity_expression.R` — RNA-seq test of the fitted capacity parameter
- `16_supplementary_figures.R` — etc-GEM diagnostic supplementaries
- `17_schematic.py` — the etc-GEM pipeline schematic
- the Figure 3 (etc-GEM) block inside `14_main_figures.R`; its Figure 1 and
  Figure 2 blocks remain part of the pipeline

`run_all.R` still calls `14_main_figures.R`, which will regenerate the unused
etc-GEM figure alongside Figures 1 and 2. That is harmless; the manuscript does not
reference it.

## Figures used by the current manuscript

| Manuscript figure | File | Produced by |
|---|---|---|
| Figure 1 | `results/figures/manuscript/FIG1_decoupling.png` | `14_main_figures.R` |
| Figure 2 | `results/figures/manuscript/FIG2_the_bill.png` | `14_main_figures.R` |
| Supplementary 1 | `results/figures/Fig_temperature_equilibration.png` | `08_temperature_equilibration_sensitivity.R` |
| Supplementary 2 | `results/figures/fig_bayes_resp_arrhenius.png` | `11_bayesian_plots.R` / `15_uncertainty_bands.R` |
| Supplementary 3 | `results/figures/fig_bayes_cue_by_clade.png` | `11_bayesian_plots.R` / `15_uncertainty_bands.R` |
| Supplementary 4 | `results/figures/Fig_n0_treatment_panel.png` | `18_n0_treatment_panel.R` |
