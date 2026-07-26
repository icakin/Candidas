# RUNBOOK — Candidas, from scratch

The definitive from-nothing instructions. Every number below was measured on the
machine recorded in `env/versions.json` (macOS 26.2, Apple silicon, 16 cores,
R 4.5.2, Python 3.9.6) during the C1 reproduction run.

`SETUP.md` is the short version. This is the one that tells you what you should
see at each stage and what to do when you do not.

---

## 0. Before you start

You need, on `PATH`: **R ≥ 4.4**, **Python ≥ 3.9**, **git**, **quarto**, and a
working **C++17 compiler** (macOS: `xcode-select --install`; Linux:
`build-essential`; Windows: Rtools).

Disk: **~3 GB** total (repo ~200 MB, renv library ~1.5 GB, Python venv ~400 MB,
TinyTeX ~500 MB, one results tree ~70 MB, `_output/` ~20 MB).

Time: **49 minutes** for a full run on this hardware (16 cores) once installed —
etc-GEM 35 m, R pipeline 12.5 m, render 44 s. Add 15–45 min for the first-ever
install (rstan compiles). On a slower or lower-core machine, budget 2–4 h: the
two long poles are the etc-GEM `calibrate`/`bayes` stages and script 08.

---

## 1. Clone — with submodules

```bash
git clone --recursive <repo-url> Candidas
cd Candidas
```

If you already cloned without `--recursive`:

```bash
git submodule update --init --recursive
```

**Check:** `ls cauris_etcgem/strains/eci_cauris/scripts/` must list
`generate_model_data.py`. An empty `cauris_etcgem/` is the single most common
failure mode in this project — the etc-GEM figures used to fail far downstream
inside `read.csv()` with no clue why. They now stop immediately with
`ETCGEM_SUPP_DATA_MISSING` or `ETCGEM_SUPP_DATA_EMPTY`.

**Avoid a path containing spaces if you can.** This project lives under
`.../Candidas TPC/Candidas`, and R encodes spaces in `--file=` as `~+~`, which
broke every `source()` until C1 fixed it. The fix is in, but it is one less
thing to think about.

---

## 2. Install

### 2a. R (10–40 min cold — rstan and StanHeaders compile)

```bash
Rscript scripts/00_install.R
```

It restores `renv.lock` (170 packages, CRAN + Bioconductor 3.22), checks every
package the pipeline loads, works out which Stan backend brms will use, then
**compiles and samples a 10-line toy model**.

**Check:** the last lines read

```
  [ok]   all 22 packages present
  ==> 08_bayesian_models.R will use backend: rstan
  [ok]   toy model sampled: mu = 3.0837 (truth 3), sigma = 0.8454 (truth 1)
```

If Stan fails here, **stop**. 08 cannot run, and you will otherwise find that
out three hours in. `00_install.R` prints per-platform instructions.
`--no-stan` skips the compile check when you already know it works.

### 2b. Python (2–5 min)

```bash
cd cauris_etcgem
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.lock
cd ..
```

**Check:**

```bash
cauris_etcgem/.venv/bin/python -c "import etcgem, cobra; print(cobra.__version__)"
```

`requirements.lock` pins the etc-GEM **engine** to commit
`7d32383d902184c085662031b802406353e7b083`. Before C1 it was installed from a
bare branch reference, so the model could change under the results with no
record. If you deliberately move the pin, re-run the whole etc-GEM pipeline —
never mix engines across the outputs.

### 2c. PDF engine (once)

```bash
quarto install tinytex
```

---

## 3. Run everything

```bash
bash scripts/run_all.sh
```

Three stages, timed, with a tee'd log per stage in `logs/`:

1. the Python etc-GEM pipeline (`generate_model_data.py all`),
2. `Rscript scripts/run_all.R` (02 → 03 → 06 → 07 → 08 → 09 → 10 → 11 → 12 → 13, then 14),
3. `quarto render` in `manuscript/draft` → `_output/`.

To regenerate without touching the committed outputs, use `run_c1.sh` (or set
the same five variables yourself):

```bash
CANDIDAS_RESULTS="$PWD/results_C1" \
ETCGEM_OUT_SUFFIX=_C1 \
CANDIDAS_SUPP_DATA="$PWD/cauris_etcgem/strains/eci_cauris/outputs/supp_data_C1" \
CANDIDAS_EXPRESSION_OUT="$PWD/results_C1/expression" \
  bash scripts/run_all.sh
```

---

## 4. Stage by stage — timings, and what to check

### Stage 1 — Python etc-GEM · **35 min**

| step | wall time | writes |
|---|---:|---|
| model build (SBML → GECKO → sectors → budget bisection) | 2.4 s | — |
| `calibrate` — differential evolution × 4 clades | **16 min** | `_percladefit*.npy` |
| `curves` | < 1 s | `figure3_data*/` |
| `fit` | < 1 s | `predicted_vs_observed.csv` |
| `ident` — 3 knobs × 25 values × 4 clades | 62 s | `identifiability_*.csv` |
| `boot` — 120-point grid, 4 000 bootstrap replicates | 93 s | `capacity_*.csv` |
| `apriori` | 9 s | `calibration_*.csv` |
| `bayes` — 15³ emulator grid then 4 × (32 walkers × 4 000 steps) | **16 min** | `bayes_*.csv`, `bayes_chain_raw.npz` |

**Check:** the model build prints

```
[sectors] P_total=5.483  f=(metab 0.5, bio 0.20, maint 0.3)  translation_coeff=1.378 (mu*=0.7955)
```

`mu* = 0.7955` is the anchored growth rate and therefore the **ceiling of the
a-priori TPC**. If you see `mu* = 1.089`, you are on the OLD anchor
(`ANCHOR_MU = 1.1 @ 313.15 K`), and your numbers will match the committed
`calibration_r2.csv` instead of the current code. See
`REPRODUCTION.md` §4.2 — this is the single biggest provenance trap in the repo.

**Check:** each clade reports `R2 = 0.88–0.96`. Clade IV is the weak one
(0.88); its emcee posterior for `kcat_scale` presses the upper prior bound
(1.05) and one of its 32 walkers gets stuck at an acceptance fraction of 0.30.

**Check:** the emcee diagnostics print acceptance fractions of ≈ 0.59–0.65 and
`tau ≈ 35–54`. With 4 000 steps that is 75–115 autocorrelation times, so
`get_autocorr_time()` returns rather than raising.

### Stage 2 — R pipeline

| script | wall time | note |
|---|---:|---|
| 02 longdata | 4.2 s | 1 469 340 rows kept of 1 486 800 (Time ≤ 1455 min) |
| 03 trimming | 1 m 55 s | 1 080 curves; writes a large diagnostics PDF |
| 06 inoculation | 0.3 s | app NOT launched; `otu_inoc.csv` read as data |
| 07 oxygen_fits | 2 m 10 s | per-series `nlsLM` fits + descriptive plots |
| **08 bayesian_models** | **8 m 00 s** | 3 models × 4 chains × 4 000 iterations (growth SS 4 m 58 s, resp Arrhenius 36 s, resp SS 2 m 18 s) |
| 09 bayesian_plots | 10.4 s | |
| 10 carbon_tax | 1.2 s | |
| 11 capacity_expression | 0.7 s | reads the COMMITTED GEO matrix; no download |
| 12 main_figures | 6.0 s | needs 08 **and** the etc-GEM outputs |
| 13 supplementary_figures | 1.8 s | needs 11 **and** the etc-GEM outputs |
| 14 schematic.py | 2.6 s | |
| **total** | **12 m 33 s** | |

**Check after 02:** `Time cutoff applied: kept 1469340 of 1486800 rows`.

**Check after 03:** `Loaded MANUAL_FIT_WINDOWS from app file (1080 curves)` and
`Loaded PLOT_EXCLUDE_POINTS from app file (88 points)`. If you instead see
`0 curves`, the committed trim decisions are not being found and **you are
running a different analysis** — see §6.

**Check after 08:** no divergent transitions, `max Rhat ≤ 1.01`, `min ESS` in
the thousands. 08 prints its own timing table.

**Check after 12/13:** the console prints which `supp_data` directory it
resolved and how many CSVs it found.

### Stage 3 — quarto render · **~40 s**

**Check:** six files in `manuscript/draft/_output/` —
`manuscript`, `supplementary`, `manuscript_combined`, each `.docx` and `.pdf`.

---

## 5. The four click-driven apps (01, 04, 05, 06)

These are **not** part of an unattended run and **must not be** — they encode
hand judgements that are committed to the repository and read as data:

| app | what it owns |
|---|---|
| `01_convert_xlsx.R` | `data/*_Oxygen.csv`, `results/tables/otu_names.csv` |
| `04_trim_selector.R` | `results/tables/manual_fit_windows.csv` (1 080 curves), `plot_exclude_points.csv` (88 points) |
| `05_cell_sizes.R` | `results/tables/otu_cell_sizes.csv` |
| `06_inoculation.R` | `results/tables/otu_inoc.csv` |

Sourcing any of them defines the app and returns. To open one deliberately:

```bash
Rscript scripts/04_trim_selector.R --app       # or CANDIDAS_RUN_APP=1 Rscript ...
```

In RStudio, "Run App" behaves as before.

---

## 6. Path knobs

All are path-only. **None of them changes an analysis decision.**

| variable | effect | default |
|---|---|---|
| `CANDIDAS_RESULTS` | the whole results tree (`tables/`, `figures/`, `rds/`) | `results/` |
| `CANDIDAS_APP_INPUTS` | where the five committed app decisions are looked up **as a fallback** | `results/tables/` |
| `CANDIDAS_SUPP_DATA` | which etc-GEM outputs 12/13 read | submodule `outputs/supp_data` |
| `CANDIDAS_EXPRESSION_OUT` | where 11 writes | `data/expression/` |
| `ETCGEM_OUT_SUFFIX` | suffix for the etc-GEM output dirs | none |
| `ETCGEM_STAGES` | which etc-GEM stages to run | `all` |
| `CANDIDAS_SKIP` | script numbers to skip, e.g. `08,09` | none |
| `SKIP_PY` / `SKIP_R` / `SKIP_RENDER` | skip a stage of `run_all.sh` | off |

`results/tables/` mixes two kinds of file: generated output, and five hand-made
app decisions. Redirecting `CANDIDAS_RESULTS` alone would silently drop the
manual trim windows, the exclusions, the per-isolate cell sizes and the isolate
names — a *different analysis*, not a re-run. `app_input()` in `config.R`
therefore looks in the redirected tree first and falls back to the canonical
`results/tables/`. You will see, in the log:

```
app input 'manual_fit_windows.csv' not in the results tree; using committed copy: .../results/tables/manual_fit_windows.csv
```

That message is expected and correct for a redirected run.

---

## 7. Troubleshooting

| symptom | cause | fix |
|---|---|---|
| `ETCGEM_SUPP_DATA_MISSING` / `_EMPTY` | submodule not populated, or the Python stage never ran | `git submodule update --init --recursive`, then run the etc-GEM pipeline |
| `ETCGEM_FILE_MISSING: capacity_isolates.csv` | the `boot` stage was skipped | `generate_model_data.py boot` (needs `_percladefit*.npy`) |
| `cannot open file .../config.R` | script-directory resolution; a path containing spaces | fixed in C1; check you are on the current commit |
| a script hangs with no output | a Shiny app got launched | make sure nothing sets `CANDIDAS_RUN_APP=1` |
| Stan will not compile | C++ toolchain | `Rscript scripts/00_install.R` and follow its output |
| `renv::restore()` fails on edgeR | edgeR is Bioconductor | `BiocManager::install(version = "3.22")` |
| PDF render fails | no LaTeX | `quarto install tinytex` |
| `mu* = 1.089` in the build banner | you are on the OLD anchor | see `REPRODUCTION.md` §4.2 |

---

## 8. What "success" looks like

* `results_C1/` holds **37 tables** (40 with the three large oxygen
  intermediates, which `.gitignore` keeps out of the repo), **78 figures**,
  **5 `.rds`** model objects and **2** expression outputs.
* The committed `results/` holds 39 tables and 83 figures. Every difference is
  accounted for:
  * *only in `results/`* — the five committed app inputs
    (`manual_fit_windows.csv`, `plot_exclude_points.csv`, `otu_cell_sizes.csv`,
    `otu_inoc.csv`, `otu_names.csv`), which are inputs and are never
    regenerated; plus **7 stale top-level copies** of the manuscript figures
    (`results/figures/FIG1_decoupling.png` and friends) left over from an older
    layout. **No current script writes them** — 12, 13 and 14 all write to
    `figures/manuscript/`. They are duplicates and should be deleted; that is a
    housekeeping change, not a C1 one.
  * *only in `results_C1/`* — `Oxygen_All_Long.csv`,
    `Oxygen_Data_Filtered.csv`, `Oxygen_Data_Smoothed_Trimmed.csv` (gitignored
    by name, so absent from the committed tree) and the two large diagnostic
    PDFs `per_series_fits.pdf`, `oxygen_trimming_diagnostics.pdf`.
* `cauris_etcgem/strains/eci_cauris/outputs/supp_data_C1/` holds **13 CSVs**
  plus `bayes_chain_raw.npz`.
* `manuscript/draft/_output/` holds **6** documents.
* `git status` shows no modification to `data/`, `results/`, or the committed
  etc-GEM outputs — verify with:

```bash
shasum -a 256 -c env/baseline_checksums_results.txt
shasum -a 256 -c env/baseline_checksums_etcgem.txt
shasum -a 256 -c env/baseline_checksums_data.txt
```

**Getting the pipeline to run is not the same as reproducing the numbers.**
Read `reports/REPRODUCTION.md` before trusting any regenerated value: the
etc-GEM outputs do **not** reproduce from the current code, and the reason is
identified there.
