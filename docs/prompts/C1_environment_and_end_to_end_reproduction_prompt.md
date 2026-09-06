# Claude Code prompt — C1: pin the environment, make the pipeline headless, reproduce every number end to end (autonomous)

Run from the project root (`.../Candidas TPC/Candidas`). This is the FIRST prompt in the Candidas
series and it is deliberately a REPRODUCTION task, not a science task: install and pin everything,
remove the things that block unattended execution, re-run the whole pipeline (R respirometry +
Python etc-GEM + Quarto manuscript) from the committed raw data, and write a reproduction report
saying which committed numbers came back and which did not. It changes NO analysis decision.

NOTE TO USER: launch in an auto-approving mode. Budget several hours of wall clock: the brms stage
(08) and the etc-GEM grid/emcee stage are both slow.

CONTEXT: the repo currently cannot be run unattended by anyone, including its author.
`run_all.R` sources `06_inoculation.R`, which is a Shiny app, so the master runner blocks forever;
four of the fourteen scripts (01, 04, 05, 06) are click-driven apps; `11_capacity_expression.R`
carries a hard-coded `/Users/ilgazcakin/...` path and calls `setwd()` mid-script; scripts 12 and 13
read `cauris_etcgem/strains/eci_cauris/outputs/supp_data/`, which was an empty submodule until now;
there is no renv lockfile, no Python pin, and the etc-GEM engine is an unversioned
`pip install git+https://github.com/GabYvonDurocher/etcGEMs` so the model changes under the results.
Separately, the committed outputs are known to be internally inconsistent (see the audit summary in
PART F.4) — the point of this prompt is to MEASURE that, not to fix it. The scientific fixes are
C2 onwards and must not be attempted here.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first: README.md,
GETTING_STARTED.md, scripts/{config.R,run_all.R,README.md}, every scripts/*.R header block,
cauris_etcgem/{README.md,SETUP.md,requirements.txt}, cauris_etcgem/strains/eci_cauris/
{MANIFEST.md,STATUS.md,strain.yaml,scripts/*.py}, and manuscript/draft/{_quarto.yml,manuscript.qmd,
_body_main.qmd,_body_supp.qmd}. Do NOT change any analysis decision (see CONSTRAINTS).

PART A - pin the environment; one command to install
- R: initialise renv, snapshot a lockfile covering every package the scripts actually load
  (tidyverse/dplyr/tidyr/readr/purrr/tibble/stringr/ggplot2/scales/ggrepel/ggdist/patchwork/rlang/
  grid, minpack.lm, zoo, shiny, rstudioapi, brms, posterior, and edgeR from Bioconductor). Record
  the R version. edgeR is NOT on CRAN — install via BiocManager and pin the Bioc release.
- Stan: brms needs a working backend. Detect which backend 08 uses (rstan vs cmdstanr), install and
  verify it by compiling and sampling a 10-line toy model. Record the backend, its version, and the
  C++ toolchain, and FAIL LOUDLY with instructions if it cannot compile.
- Python: create cauris_etcgem/.venv, install from requirements.txt, and PIN the engine — resolve
  github.com/GabYvonDurocher/etcGEMs to a concrete commit SHA and record it as
  `etcgem @ git+https://github.com/GabYvonDurocher/etcGEMs@<sha>` in a new
  cauris_etcgem/requirements.lock. Also pin an LP solver (glpk via cobra) and record its version.
  Record whether the pinned engine is the one that produced the committed outputs (PART F.4 tests it).
- Quarto: verify `quarto` and the PDF engine (`quarto install tinytex`); record versions.
- Deliverables: `scripts/00_install.R` (idempotent; installs R side, restores renv, checks Stan),
  `SETUP.md` at the project root (prerequisites, the two install commands, expected runtimes per
  stage, disk needed), and `env/versions.json` capturing every version and SHA above. SETUP.md must
  be sufficient for a new machine with nothing but R, Python and git.

PART B - make the pipeline headless (execution fixes ONLY)
- The four Shiny scripts (01_convert_xlsx, 04_trim_selector, 05_cell_sizes, 06_inoculation) must be
  SOURCEABLE without launching a browser. Guard each app launch behind an interactive check, e.g.
  `if (interactive() && !isTRUE(getOption("candidas.headless"))) shiny::runApp(...)`, so that under
  Rscript the file defines its objects, writes nothing, and returns. Their committed outputs
  (tables/manual_fit_windows.csv, plot_exclude_points.csv, otu_cell_sizes.csv, otu_inoc.csv) are
  INPUTS to this run and must be treated as data - do not regenerate or alter them.
- `run_all.R` must complete under `Rscript scripts/run_all.R` with no user interaction, running
  02 -> 03 -> 06 -> 07 -> 08 -> 09 -> 10 -> 11 -> 12 -> 13, then 14 (python3). Add 11 and 13 to the
  loop (currently they are only mentioned in a message). Keep the existing order otherwise.
- `11_capacity_expression.R`: replace the hard-coded `/Users/ilgazcakin/...` fallback with a path
  built from `base_dir`, and remove the `setwd()` (use explicit paths). Confirm it reads the
  committed data/expression/GSE165762/GSE165762_Raw_counts.txt.gz and does NOT need to download.
- `12_main_figures.R` and `13_supplementary_figures.R`: their `DD` candidate-path search must fail
  LOUDLY with a named error if the etc-GEM supp_data directory is missing or empty, instead of
  silently taking `NA` from `which(...)[1]`. Same for a missing capacity_isolates.csv.
- Fix `probe.py` (it passes Kelvin to `compute_tpc`, which expects Celsius) and the two broken paths
  in `identifiability_stageA.py` (`{HERE}/../../tables/` does not exist; the column is `T`, not
  `temp_C`). If identifiability_stageA cannot be made to run because its respiration hook is a
  `return np.full(len(temps), np.nan)` placeholder, leave the stub alone and record it in the report
  as NOT RUNNABLE - do not invent a respiration hook.
- Add `scripts/run_all.sh` at the project root that runs, in order: the Python etc-GEM pipeline, then
  `Rscript scripts/run_all.R`, then the Quarto render, with `set -euo pipefail`, per-stage timing,
  and a tee'd log to `logs/`.

PART C - run the Python etc-GEM pipeline (non-destructive)
- Run `python scripts/generate_model_data.py all` from cauris_etcgem/strains/eci_cauris/, with
  outputs redirected to a NEW directory `outputs/supp_data_C1/` (add a CLI/env override for the
  output dir; do not write into `outputs/supp_data/` or `outputs/supp_data_original_backup/`).
- This MUST include `stage_boot`, so capacity_isolates.csv / capacity_bootstrap*.csv / capacity_corr
  .csv are regenerated on the current calibration rather than inherited from an older run.
- Record wall time per stage, the DE seed, the emcee settings (walkers/steps/discard/thin) and the
  emulator grid shape actually used. Save the raw emcee chain WITH walker and step indices to
  `outputs/supp_data_C1/bayes_chain_raw.npz` so convergence can be assessed later (C-series). Report
  acceptance fraction and `get_autocorr_time()` (catch and report the exception if it raises).

PART D - run the R pipeline (non-destructive)
- Point the pipeline at a NEW results tree `results_C1/{tables,figures,rds}` via a single override in
  config.R (e.g. `RESULTS_DIR <- Sys.getenv("CANDIDAS_RESULTS", file.path(base_dir,"results"))`), so
  `results/` is left untouched for comparison. Everything else in config.R stays exactly as it is.
- Force recomputation: if 08/09 reuse cached `results/rds/*.rds`, make the cache path follow
  RESULTS_DIR so the brms fits are genuinely refitted into results_C1/rds. Record the brms seed,
  chains, iterations, adapt_delta actually used, and the resulting max Rhat / min ESS.
- Have 12/13 read the etc-GEM data from `outputs/supp_data_C1/` for this run.
- Record per-script wall time and capture ALL console output to `logs/NN_<script>.log`. Any `Error:`,
  any brms divergence, any "SS fit failed", any dropped-curve message must appear in the report.

PART E - render the manuscript
- `cd manuscript/draft && quarto render` -> `_output/` (manuscript, supplementary,
  manuscript_combined; .docx + .pdf). Confirm every figure referenced resolves and that the build
  uses the figures just regenerated (temporarily point the manuscript at the C1 figure directory, or
  copy results_C1 figures into a staging dir - state which you did).
- Report any missing figure, unresolved citation, or broken cross-reference.

PART F - the reproduction report (the actual deliverable)
1. `reports/REPRODUCTION.md` — for every committed output table, a row-by-row numeric comparison of
   `results/` vs `results_C1/` and `outputs/supp_data/` vs `outputs/supp_data_C1/`. Report max
   absolute and max relative difference per file, and classify each as IDENTICAL (<1e-10) /
   NUMERICALLY EQUIVALENT (<1e-6) / STOCHASTIC-BUT-CONSISTENT (differs, but every 95% interval
   overlaps) / DIVERGENT (does not overlap) / NOT REPRODUCED (could not be regenerated at all).
2. A headline table of the manuscript's load-bearing numbers, committed vs regenerated:
   E_G and E_R per taxon; dE and its CrI; T_opt(CUE) per taxon; tax37 and tax40 per taxon; the 10
   pairwise fever contrasts and how many resolve; Spearman rho; the etc-GEM R2 per clade; kcat_scale
   per clade; the 12 isolate capacities; ICC1, F(3,8), P; r = -0.86 and r = -0.64.
3. A SEEDS AND DETERMINISM section: for every stochastic stage (brms, differential_evolution, emcee,
   the two bootstraps in 12 and in generate_model_data.py), state whether a seed is set, where, and
   whether the stage is bit-reproducible on a re-run. Run the cheap stages TWICE to check.
4. A KNOWN INCONSISTENCIES section — verify or refute each of these by direct test, and report the
   numbers either way. Do NOT fix them here:
   * `outputs/supp_data/calibration_curves.csv` has a `calibrated` curve for Clade IV peaking at
     ~0.931 above its own `apriori` plateau of ~0.7955. The biosynthesis cap is a hard ceiling
     independent of all three knobs, so a single build should not be able to produce this. Test
     whether the two columns are consistent with one build.
   * Evaluating the committed parameters under the current `_supp_common.py` anchor
     (ANCHOR_MU=0.95 @ 307.15K) does not reproduce `calibration_r2.csv`; the older anchor
     (mu=1.1 @ 313.15K, still named in a comment) does. Report R2 under both.
   * `capacity_isolates.csv` is byte-identical between `supp_data/` and `supp_data_original_backup/`
     and its values are exact integer multiples of the 120-point grid step (0.0096639) away from the
     BACKUP clade kcat_scale, not the current one. Report the offsets in grid steps under both, and
     whether the PART C rerun changes ICC1, F, P, r=-0.86 and r=-0.64.
   * `figure3_data/capacity.csv` (0.670/0.443/0.679/0.775) disagrees with
     `supp_data/calibration_r2.csv` kcat_scale (0.614/0.410/0.622/0.710) for the same quantity.
   * `apriori_Topt` in calibration_r2.csv differs across clades (40/30/30/40) although the a-priori
     curve is clade-independent by construction. Report the across-clade spread of the a-priori
     curve and the size of the plateau.
   * `drop_excluded()` / EXCLUDE_GROUPS ("glab") is called only in 08 and 09, so `results/tables/`
     and every 07 figure still contain glabrata while the Bayesian results do not. Confirm and
     quantify (row counts per file).
   * `capacity_bootstrap_draws.csv` pairs capacity and peak via two independent `rng.choice` calls.
     Report the within-clade corr(capacity, peak) in the file.
5. A BLOCKERS section: everything that had to be changed to make the pipeline run, with a one-line
   justification each, and everything that still cannot be run (with the reason).
6. `reports/RUNBOOK.md` — the definitive from-scratch instructions: clone (with `--recursive`),
   install, run, expected wall time per stage, expected disk, and what to check at each stage.

VERIFY (report all)
1. `bash scripts/run_all.sh` completes on a clean checkout with zero user interaction; per-stage
   timings printed; logs written. State the total wall time.
2. Every one of the 14 scripts either ran, or is documented as intentionally not run (01/04/05 are
   click-driven and their outputs are committed inputs) — no silent skips.
3. results_C1/ and outputs/supp_data_C1/ are fully populated; results/ and outputs/supp_data/ are
   BYTE-IDENTICAL to their pre-run state (verify with git status / checksums).
4. The manuscript renders to .docx and .pdf with no missing figure and no unresolved reference.
5. reports/REPRODUCTION.md exists and contains all six sections, with ACTUAL numbers in the headline
   table; every committed load-bearing number is classified.
6. env/versions.json, SETUP.md, reports/RUNBOOK.md, cauris_etcgem/requirements.lock and renv.lock
   are committed.
7. State plainly, in one paragraph at the top of the report: can this repository be reproduced from
   raw data by a third party today, yes or no, and if no, exactly what is missing.

CONSTRAINTS
- CHANGE NO ANALYSIS DECISION. Specifically do not touch: N0_BACKPROJECT, USE_DRAWDOWN_WINDOW,
  FIT_DRAWDOWN_FRAC, USE_APP_TRIM_FILES, EXCLUDE_GROUPS/TEMPS/OTUS/T_OTU, N_inoculation_cells_per_L,
  the cell-carbon constants, RESPIRATORY_QUOTIENT, any brms formula/prior/seed, the DE bounds or
  seed, the emcee priors/bounds, THERMAL_MODE, ANCHOR_MU/ANCHOR_T, the sector fractions, or any
  committed trim window / exclusion file. If one of them MUST change for the code to execute, do not
  change it - report it as a blocker.
- Non-destructive: write only to new paths (results_C1/, outputs/supp_data_C1/, logs/, reports/,
  env/, SETUP.md, prompts/). Never overwrite results/, outputs/supp_data/, supp_data_original_backup/,
  data/, or manuscript/draft/*.qmd content (rendering to _output/ is fine).
- Fixes are limited to what blocks unattended execution: app guards, the hard-coded path, the
  setwd(), the fail-loud path resolution, the probe.py unit bug, the identifiability path bug, the
  output-dir override, and the runner. Nothing else.
- Report ACTUAL numbers everywhere; "matches" without a number is not acceptable. Where something
  does not reproduce, say so plainly rather than reconciling it.
- Autonomous; commit in parts: "C1: pin R/Python/Stan/Quarto environment + 00_install + SETUP",
  "C1: headless guards, path fixes, run_all.sh - pipeline runs unattended",
  "C1: rerun etc-GEM pipeline into supp_data_C1 (incl. stage_boot)",
  "C1: rerun R pipeline into results_C1 + render manuscript",
  "C1: reproduction report, runbook, known-inconsistency tests".
```
