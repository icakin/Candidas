# scripts/ — the R pipeline (temperature × isolate respirometry)

## Overview

Raw data are dissolved-oxygen time series from a PreSens 24-well SDR plate reader, one
plate per group per temperature, at 22, 24, …, 44 °C (12 temperatures). Eight groups,
three isolates each, five replicate wells per isolate:

| group | species | isolates |
|---|---|---|
| `Clade1`–`Clade4` | *Candida auris*, clades I–IV | 2068–2070, 2071–2073, 2074–2076, 2077–2079 |
| `Hae` | *C. haemulonii* | 1724, 1768, 1769 |
| `Duo` | *C. duobushaemulonii* | 1770, 1771 (+ one failed) |
| `para` | *C. parapsilosis* | 2051–2053 |
| `glab` | *C. glabrata* | 2048–2050 |

Isolate codes 1–24 are assigned in the `STRAIN_GROUPS` order in `config.R`; display
names live in `results/tables/otu_names.csv`.

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
  strain/isolate ID). The CSV headers use safe codes `OTU1`…`OTU24`; your names
  are written to `results/tables/otu_names.csv` and shown in every result table and plot;
- set the **temperature** per file (auto-filled from the plate's measured `Tm`
  column; editable).

Run it in RStudio (open `01_convert_xlsx.R` → **Run App**) or with
`Rscript scripts/01_convert_xlsx.R`. It writes one `<name>_Oxygen.csv` per file
into `data/`, with columns `Time, T, OTU<c>_R1, …, OTU<c>_R5` for the plate's 3
isolates, plus `results/tables/otu_names.csv`.

## Script map (running order)

Every script sources `config.R` (shared paths, constants, isolate registry, the figure
whitelist). `run_all.R` runs the non-interactive ones in order; `run_all.sh` wraps that with
the C11 report and the manuscript render, logged and timed.

| # | script | what it does | writes |
|---|---|---|---|
| 00 | `00_install.R` | restore the pinned R library (`renv.lock`), check Stan | |
| 01 | `01_convert_xlsx.R` | **Shiny app**: PreSens `.xlsx` → tidy `.csv`, isolate names | `data/*_Oxygen.csv`, `results/tables/otu_names.csv` |
| 02 | `02_longdata.R` | wide → long | `results/tables/Oxygen_All_Long.csv` |
| 03 | `03_trimming.R` | spline-based trimming of each trace, diagnostics | `Oxygen_Data_Filtered.csv`, `Oxygen_Trimmed_Series_Metadata.csv` |
| 04 | `04_trim_selector.R` | **Shiny app**: review/override each curve's fit window | `manual_fit_windows.csv`, `plot_exclude_points.csv` |
| 05 | `05_cell_sizes.R` | **Shiny app**: cell volume per isolate → carbon per cell | `otu_cell_sizes.csv` |
| 06 | `06_inoculation.R` | **Shiny app**: inoculation density per group (sourced headless in a batch run) | `otu_inoc.csv` |
| 07 | `07_oxygen_fits.R` | per-well exponential fits; growth, respiration and CUE in carbon units; frequentist Sharpe-Schoolfield and Arrhenius per isolate | `derived_N0_R_results_with_carbon.csv`, `fit_coefficients_long.csv`, `fit_metrics.csv`, … |
| 08 | `08_temperature_equilibration_sensitivity.R` | how much each quantity moves with the onset-temperature assumption | `temperature_equilibration_*.csv`, **Supp Fig 1** |
| 09 | `09_bayesian_models.R` | hierarchical Bayesian TPCs (brms): growth Sharpe-Schoolfield, respiration Arrhenius | `results/rds/bayes_*.rds`, `bayes_*.csv` (slow) |
| 10 | `10_n0_term_test.R` | diagnostic: the r·δ back-projection term (refits 09; run deliberately) | |
| 11 | `11_bayesian_plots.R` | posterior curves, contrasts, CUE by clade/isolate | `fig_bayes_*.png` incl. **Supp Figs 2, 3** |
| 12 | `12_carbon_tax.R` | the carbon cost of 37 → 40 °C, per isolate and group | `carbon_tax_*.csv`, `FIG_carbon_tax.png` |
| 13 | `13_uncertainty_bands.R` | respiration and CUE with the equilibration uncertainty combined | `fig_*_uncertainty.png` |
| 14 | `14_rk_covariance_check.R` | diagnostic: does dropping cov(r, K) misstate the uncertainty? | `fig_rk_covariance_check.png` |
| 15 | `15_n0_treatment_panel.R` | the between-clade respiration E under three N₀ treatments (refits 09; run deliberately) | **Supp Fig 4** |
| — | `fig_common.R` | shared by 16 and 17: loads the 09 fits, draw-level quantities, theme, palette, `save_fig()`, Figure 1 panels | |
| 16 | `16_fig1.R` | **Figure 1** growth and respiration decouple | `results/figures/manuscript/FIG1_decoupling.png` |
| 17 | `17_fig2.R` | **Figure 2** the carbon cost of fever, + four supplementary panels | `FIG2_consequences.png`, `FIG_SUPP_*.png` |
| 18 | `18_fig3.R` | **Figure 3** high-temperature growth vs phylogeny (reads `phylo/trees/pg_rooted.nwk`) | `FIG3_discordance.{png,pdf}` |
| 19 | `19_fig4.R` | **Figure 4** the etcGEM counterfactual (reads `gem/tables/`; see `gem/README.md`) | `FIG4_etcgem_counterfactual.{png,pdf}` |
| 21 | `21_figures_all.R` | runs 16–19 in one command | |

The four Shiny apps are not run in a batch: their outputs are committed inputs. 10, 14
and 15 overwrite the published Bayesian respiration summary while they run; after any of
them, re-run 09 → 11 → 12 → 16 → 17.

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

## How to run

Once the CSVs exist (Step 1):

```
Rscript scripts/run_all.R          # 02 03 06 07 08 09 11 12 13, then Figures 1-4
bash scripts/run_all.sh            # the same, plus the C11 report and the manuscript render, logged
Rscript scripts/21_figures_all.R   # figures only
Rscript scripts/18_fig3.R          # any single script runs on its own
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

> `07_oxygen_fits.R` fits frequentist Sharpe-Schoolfield and Boltzmann-Arrhenius curves
> per isolate (fast, for diagnostics). The curves the manuscript reports come from the
> hierarchical Bayesian fits in `09_bayesian_models.R`.

---

## The etc-GEM that was cut, and the one that replaced it

An earlier single-proteome etc-GEM (one model reused for every clade, its capacity
parameter restating the growth data) was removed from the manuscript in August 2026; its
two R scripts, schematic, vendored outputs and the GEO expression test are in
`archive/old_etcgem/` for the record. The model in the current Figure 4 is different: one
reconstruction per species from its own genome, calibrated on *C. auris* alone and frozen,
then predicting the relatives. It lives in `gem/` (Python) and is documented in
`gem/README.md` and `gem/FIG4_LOCKED.md`.
