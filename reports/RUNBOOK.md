# RUNBOOK — Candidas, from scratch

> **2026-09-06.** The repository was reorganised after this runbook was last updated
> (scripts renumbered, the `cauris_etcgem` submodule and `outputs/supp_data/` removed, the
> etcGEM layer now in `gem/` with its own README). The current instructions are `README.md`,
> `SETUP.md`, `scripts/README.md` and `gem/README.md`; this file stands as the record of
> the August packaging and its stage names no longer match `run_all.sh`.


What to run, in what order, what it should take, and what to check at each stage.

**Updated 2026-08-11.** The original was written for C1, before the script
renumbering (14 → 19 scripts), before `run_all.sh` gained its current stage list,
and before the manuscript was flattened from `manuscript/draft/v3_quarto/` to
`manuscript/v3.qmd`. It has been corrected for all three rather than annotated:
this is instructions, not a finding, and its only value is being right. The C1-era
original is in the git history of this file.

For **installation and system prerequisites**, read
[`SETUP.md`](../SETUP.md) — it is the authority, and it is more complete than the
install section this runbook used to carry. What follows assumes the environment
is already in place.

---

## 0. Before you start

| | |
|---|---|
| Time | **~13 min** for a full run once installed, on 16 cores. The Bayesian stage (09) is ~8 min of it. |
| Disk | ~3 GB all in: repository, `renv/library`, TinyTeX. |
| Network | Not needed for a run. The GEO count matrix is committed; the etc-GEM outputs are vendored. |

Install first: `Rscript scripts/00_install.R`. It checks system prerequisites
before touching the library, restores `renv.lock`, and compiles and samples a toy
Stan model so a broken toolchain fails in 40 seconds rather than an hour into 09.

---

## 1. Clone — with submodules

```bash
git clone --recursive https://github.com/icakin/Candidas.git
cd Candidas
# already cloned without --recursive?
git submodule update --init --recursive
```

An unpopulated submodule is the single most common way this project fails. The
etc-GEM outputs are also vendored at `outputs/supp_data/`, so the model figures
resolve even without it, but the submodule is where the model itself lives.

---

## 2. Run everything

```bash
bash scripts/run_all.sh
```

Unattended: no prompts, no browser. Per-stage timing, tee'd logs in `logs/`,
`set -euo pipefail`, and a SUMMARY table at the end.

| Stage | What runs |
|---|---|
| 0 *(optional, off)* | etc-GEM Python pipeline — `RUN_ETCGEM=1` to enable |
| 1 | `Rscript scripts/run_all.R` → **02 → 03 → 06 → 07 → 08 → 09 → 11 → 12 → 14 → 15** |
| 2 | the C11 scale-free report (`reports/C11_scale_free/R/*.R`, then its render) |
| 3 *(optional)* | `13_capacity_expression.R`, `16_supplementary_figures.R`, `17_schematic.py` — the etc-GEM figure stages, gated off by default |
| 4 | `quarto render` on `manuscript/v3.qmd` |

Knobs, all optional: `RUN_ETCGEM`, `ETCGEM_STAGES`, `SKIP_R`, `SKIP_C11`,
`SKIP_RENDER`, `LOG_DIR`. Read the header of `scripts/run_all.sh`; it is the
authority on which script runs where and why.

### Which scripts do *not* run, and why

**The four click-driven Shiny apps — 01, 04, 05, 06.** Their outputs are
committed **inputs** and are treated as data. Sourcing any of them defines the
app and returns; none launches a browser in a batch run. 06 *is* sourced by
`run_all.R` (headless); 01, 04 and 05 are not run at all.

| App | Output it owns |
|---|---|
| `01_convert_xlsx.R` | `data/*_Oxygen.csv`, `results/tables/otu_names.csv` |
| `04_trim_selector.R` | `results/tables/manual_fit_windows.csv` (1 080 curves), `plot_exclude_points.csv` (88 points) |
| `05_cell_sizes.R` | `results/tables/otu_cell_sizes.csv` |
| `06_inoculation.R` | `results/tables/otu_inoc.csv` |

To open one deliberately:

```bash
Rscript scripts/04_trim_selector.R --app     # or CANDIDAS_RUN_APP=1 Rscript ...
```

**The two N₀ diagnostics — 10 and 18.** Both refit the Bayesian respiration model
(twice and three times) and write each alternative fit over
`results/tables/bayes_resp_arr_summary.csv` and `results/rds/`, restoring only
`derived_N0_R_results_with_carbon.csv`. A batch ending with either would leave the
*published* tables holding an alternative fit. Run them deliberately, then restore
the published set:

```bash
Rscript scripts/15_n0_treatment_panel.R       # Supplementary Figure 6, ~24 min
Rscript scripts/09_bayesian_models.R && Rscript scripts/11_bayesian_plots.R \
  && Rscript scripts/12_carbon_tax.R && Rscript scripts/16_fig1.R / 17_fig2.R \
  && Rscript scripts/13_uncertainty_bands.R
```

**`14_rk_covariance_check.R`** is a standalone diagnostic and is in no runner.

---

## 3. Stage by stage — timings and what to check

Measured on the machine in `env/versions.json` (Apple silicon, 16 cores).

### Stage 1 — the R pipeline · **~12 min**

| script | time | check |
|---|---|---|
| 02 longdata | ~1 min | 1 080 series; row count reported |
| 03 trimming | ~1 min | spline trimming; `manual_fit_windows.csv` loaded for 1 080 curves |
| 06 inoculation | seconds | prints "headless — app NOT launched" |
| 07 oxygen_fits | ~2 min | per-series fits; writes `derived_N0_R_results_with_carbon.csv` |
| 08 temperature_equilibration | ~30 s | writes `temperature_equilibration_percurve.csv`, an input to 11 and 15 |
| **09 bayesian_models** | **~8 min** | 3 models × 4 chains × 4 000 iter. **Check R̂ < 1.01 and no divergences** in the summary CSVs |
| 11 bayesian_plots | ~30 s | posterior/contrast/CUE, with the equilibration envelope from 08 |
| 12 carbon_tax | ~20 s | `carbon_tax_isolate.csv` |
| 14 main_figures | ~30 s | Fig 1 and Fig 2 into `results/figures/manuscript/` |
| 15 uncertainty_bands | ~10 s | `fig_respiration_uncertainty.png`, `fig_cue_uncertainty.png` |

### Stage 2 — the C11 scale-free report · **~1 min**

Four scripts under `reports/C11_scale_free/R/` then a Quarto render. It runs with
the pipeline on purpose: the manuscript's scale-free numbers are produced here
rather than under `scripts/`, so if this does not regenerate, those numbers drift.

**Check:** `reports/C11_scale_free/C11_scale_free.pdf`, 13 pages.

### Stage 3 — the etc-GEM figure stages · *off by default*

13, 16 and 17. The etc-GEM section has been cut from the manuscript, so these no
longer feed it; they are kept for the separate methods thread.

### Stage 4 — the manuscript · **~25 s**

`quarto render` on `manuscript/v3.qmd`, xelatex via quarto's TinyTeX.

**Check:** `manuscript/v3.pdf`, **15 pages, 6 figures**, and every figure
resolving through `../results/…` rather than an embedded bitmap. If a figure is
missing, the path is wrong, not the pipeline — see
[`MANUSCRIPT_FIGURE_LINKAGE.md`](MANUSCRIPT_FIGURE_LINKAGE.md).

---

## 4. What "success" looks like

* `run_all.sh` ends with a SUMMARY table, every stage `exit 0`.
* `results/tables/` regenerates with no table moving beyond floating-point noise.
  Two exceptions are known and documented, both traced to the environment rather
  than to the code: see `REPRODUCTION.md` and the C9 pull request.
* `results/tables/carbon_tax_curves.csv` **will** differ between runs — it uses an
  unseeded `sample.int`. So will several PNGs, from unseeded jitter. That is a
  known defect, recorded as residual item 2 in
  [`MANUSCRIPT_FIGURE_LINKAGE.md`](MANUSCRIPT_FIGURE_LINKAGE.md), not a failure.
* `manuscript/v3.pdf` is 15 pages with 6 figures.
* ~~`env/baseline_checksums_*.txt`~~ removed 2026-09-07: a checksum baseline cannot hold across machines (PDF timestamps, plot jitter, brms date stamps); git is the record, and the check is `git status` on the text tables after a re-run. Re-run and diff against git to prove nothing
  under `data/`, `results/` or the etc-GEM outputs moved.

---

## 5. Troubleshooting

Nearly everything is in [`SETUP.md`](../SETUP.md), which has a full
prerequisites section and a troubleshooting list. The three that bite most often:

| Symptom | Cause | Fix |
|---|---|---|
| `cannot open file '.../config.R'` | a script's directory resolution, on a path containing a space | every script un-mangles R's `~+~` now; if you see this, that script is unpatched |
| Stan will not compile | C++ toolchain | `Rscript scripts/00_install.R` and follow its output |
| A script hangs with no output | a Shiny app launched in a batch | check nothing sets `CANDIDAS_RUN_APP=1` |
| `renv::restore()` builds everything from source | 73 of 170 pins have no current CRAN binary | expected; see `SETUP.md` "System prerequisites" |
