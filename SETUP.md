# Setup — Candidas

Everything you need to take a machine that has **only R, Python 3 and git** to a
state where `bash scripts/run_all.sh` runs the whole project unattended.

If you just want the commands, read [Install](#install) and [Run](#run).
If you want to know what to expect at each stage, read
[`reports/RUNBOOK.md`](reports/RUNBOOK.md).

---

## Prerequisites

| Thing | Minimum | What this project was built and verified with |
|---|---|---|
| **R** | 4.4 | **4.5.2** (`aarch64-apple-darwin20`) |
| **A C++17 toolchain** | any working one | Apple clang 17.0.0 (Xcode CLT); `R CMD config CXX17` → `clang++ -arch arm64` |
| **Python** | 3.9 | **3.9.6** (macOS system `/usr/bin/python3`, arm64) |
| **git** | 2.x | 2.50.1 |
| **Quarto** | 1.4 | **1.8.27** (bundles pandoc 3.6.3) |
| **A LaTeX engine** | xelatex | quarto's own TinyTeX **v2026.02** (XeTeX 3.141592653, TeX Live 2025) |

Platform notes:

* **macOS** — `xcode-select --install` gives you clang and `make`. That is the
  whole toolchain requirement; every CRAN package this project uses resolves to
  an arm64 **binary**, so no `gfortran` is needed. (If you force a source
  install of something Fortran-flavoured you will need it —
  see <https://mac.r-project.org/tools/>.)
* **Linux** — `sudo apt install build-essential` (needs `g++` and GNU make).
* **Windows** — install Rtools matching your R version.

Stan is the fussy dependency. `scripts/00_install.R` **compiles and samples a
toy model** to prove the toolchain works before you spend hours on the real fit;
if it cannot, it fails loudly with per-platform instructions.

---

## Install

Clone **with submodules** — the etc-GEM model lives in one, and an unpopulated
submodule is the single most common way this project fails:

```bash
git clone --recursive https://github.com/<you>/Candidas.git
cd Candidas
# already cloned without --recursive?
git submodule update --init --recursive
```

Then the two install commands:

```bash
# 1. R side  (restores renv.lock, then verifies Stan by compiling a toy model)
Rscript scripts/00_install.R

# 2. Python side  (the etc-GEM engine + its dependencies, fully pinned)
cd cauris_etcgem && python3 -m venv .venv \
  && .venv/bin/python -m pip install -r requirements.lock && cd ..
```

and, once, the PDF engine:

```bash
quarto install tinytex
```

### What gets pinned, and where

| Layer | Pin file | Notes |
|---|---|---|
| R packages | `renv.lock` | 170 packages, R 4.5.2, CRAN + **Bioconductor 3.22** (edgeR 4.8.2, limma 3.66.0). `renv::restore()` handles both. |
| Stan | via `renv.lock` | rstan 2.32.7 / StanHeaders 2.32.10 / Stan 2.32.2. **`08_bayesian_models.R` never sets `backend`, so brms uses its default: rstan.** cmdstanr is not used. |
| Python | `cauris_etcgem/requirements.lock` | Full `pip freeze`, plus the engine pinned to a commit (below). |
| etc-GEM **engine** | `cauris_etcgem/requirements.lock` | `etcgem @ git+https://github.com/GabYvonDurocher/etcGEMs@7d32383d902184c085662031b802406353e7b083` |
| LP solver | `cauris_etcgem/requirements.lock` | GLPK 5.0 via `swiglpk` 5.0.13. cobra sees only `glpk`, `glpk_exact`, `scipy` here, so the choice is unambiguous. |
| Everything, as measured | `env/versions.json` | Every version above, read from the live environment. |

**Why the engine pin matters.** `requirements.txt` asked for
`pip install git+https://github.com/GabYvonDurocher/etcGEMs` with no ref. That
tracks `main`, so the model could — and did — change under the committed
results with no version bump, no tag and no record. The lock now names a commit.
If you deliberately move it, change it in `requirements.lock`, note it in
`env/versions.json`, and re-run the whole etc-GEM pipeline; do not mix engines.

---

## Run

```bash
bash scripts/run_all.sh
```

That runs, in order, with per-stage timing and a tee'd log in `logs/`:

1. the Python etc-GEM pipeline (`generate_model_data.py all`),
2. `Rscript scripts/run_all.R` (02 → 03 → 06 → 07 → 08 → 09 → 10 → 11 → 12 → 13, then 14),
3. `quarto render` in `manuscript/draft` → `_output/`.

It is fully unattended: no prompts, no browser.

### Running non-destructively

Every path is an environment knob, so you can regenerate everything into a fresh
tree and leave the committed outputs untouched for comparison:

```bash
CANDIDAS_RESULTS="$PWD/results_C1" \
ETCGEM_OUT_SUFFIX=_C1 \
CANDIDAS_SUPP_DATA="$PWD/cauris_etcgem/strains/eci_cauris/outputs/supp_data_C1" \
CANDIDAS_EXPRESSION_OUT="$PWD/results_C1/expression" \
  bash scripts/run_all.sh
```

| Variable | Effect | Default |
|---|---|---|
| `CANDIDAS_RESULTS` | whole results tree (`tables/`, `figures/`, `rds/`) | `results/` |
| `CANDIDAS_SUPP_DATA` | which etc-GEM outputs 12/13 read | submodule `outputs/supp_data` |
| `CANDIDAS_EXPRESSION_OUT` | where 11 writes its expression CSV/PDF | `data/expression/` |
| `ETCGEM_OUT_SUFFIX` | suffix for the etc-GEM output dirs | none |
| `ETCGEM_STAGES` | which etc-GEM stages to run | `all` |
| `CANDIDAS_SKIP` | script numbers to skip, e.g. `08,09` | none |
| `SKIP_PY`, `SKIP_R`, `SKIP_RENDER` | skip a whole stage of `run_all.sh` | off |

The cached brms fits live under `CANDIDAS_RESULTS/rds/`, so pointing at an empty
tree genuinely **refits** rather than silently reusing `results/rds/*.rds`.

### The four click-driven apps (01, 04, 05, 06)

These are Shiny apps. Their outputs are **committed inputs** to the pipeline and
are treated as data — an unattended run never regenerates them:

| App | Output it owns |
|---|---|
| `01_convert_xlsx.R` | `data/*_Oxygen.csv`, `results/tables/otu_names.csv` |
| `04_trim_selector.R` | `results/tables/manual_fit_windows.csv`, `plot_exclude_points.csv` |
| `05_cell_sizes.R` | `results/tables/otu_cell_sizes.csv` |
| `06_inoculation.R` | `results/tables/otu_inoc.csv` |

Sourcing any of them defines the app and returns. To actually open one:

```bash
Rscript scripts/04_trim_selector.R --app     # or CANDIDAS_RUN_APP=1 Rscript ...
```

In RStudio, "Run App" works as it always did.

---

## Expected runtime and disk

Measured on the machine in `env/versions.json` (Apple silicon, 16 cores).
Per-stage numbers as actually observed are in
[`reports/RUNBOOK.md`](reports/RUNBOOK.md) and `reports/REPRODUCTION.md`.

| Stage | Wall time |
|---|---|
| `Rscript scripts/00_install.R` — first run, from nothing | 10–40 min (rstan/StanHeaders compile) |
| `00_install.R` — subsequent runs | < 1 min (`--no-stan`: seconds) |
| Python venv + `requirements.lock` | 2–5 min |
| etc-GEM: model build | ~4 s |
| etc-GEM: `calibrate` (differential evolution, 4 clades) | ~16 min |
| etc-GEM: `curves` / `fit` / `apriori` | seconds each |
| etc-GEM: `ident` (3 knobs × 25 values × 4 clades) | ~1 min |
| etc-GEM: `boot` (120-point grid + 4000 bootstrap draws) | ~1.5 min |
| etc-GEM: `bayes` (15³ emulator grid, then 4 × 32-walker × 4000-step emcee) | ~16 min |
| R 02 → 03 (reshape + spline trimming, 1 080 series) | ~2 min |
| R 07 (per-series fits + descriptive plots) | ~2 min |
| **R 08 (brms: 3 models × 4 chains × 4000 iter)** | **8 min measured on 16 cores; budget ~1 h on a slow machine** |
| R 09 / 10 / 11 / 12 / 13 / 14 | ~25 s total |
| `quarto render` (3 documents × docx + pdf) | ~1 min |

**Measured total on the reference machine: 48 min 50 s** (etc-GEM 35 m 29 s,
R pipeline 12 m 37 s, render 44 s) — see `reports/REPRODUCTION.md` §7. Budget
2–4 h on a slower or lower-core machine; the etc-GEM `calibrate`/`bayes` stages
and script 08 are the long poles.

Disk:

| Item | Size |
|---|---|
| repository (with submodule, no outputs) | ~200 MB |
| `renv/library` | ~1.5 GB |
| `cauris_etcgem/.venv` | ~400 MB |
| TinyTeX | ~500 MB |
| one full results tree (`results/` or `results_C1/`) | ~150 MB |
| one etc-GEM outputs tree | ~5 MB (+ ~15 MB if raw emcee chains are kept) |
| `manuscript/draft/_output/` | ~20 MB |
| **total, comfortably** | **~3 GB** |

---

## Troubleshooting

**`ETCGEM_SUPP_DATA_MISSING` / `ETCGEM_SUPP_DATA_EMPTY` from 12 or 13.**
The submodule is not populated, or the Python pipeline has not run.

```bash
git submodule update --init --recursive
cd cauris_etcgem/strains/eci_cauris/scripts && ../../../.venv/bin/python generate_model_data.py all
```

An *empty but existing* directory used to be accepted silently and fail much
later inside `read.csv()`; it is now a named error at the point of resolution.

**Stan will not compile.** Run `Rscript scripts/00_install.R` and read its
instructions. On macOS start with `xcode-select --install`, then
<https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Mac>.
To switch to cmdstanr instead:

```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
options(brms.backend = "cmdstanr")   # in ~/.Rprofile
```

**`renv::restore()` fails on edgeR.** edgeR is Bioconductor, not CRAN:

```r
install.packages("BiocManager"); BiocManager::install(version = "3.22")
```

**A script hangs with no output.** It should not any more — that was the four
Shiny apps, and they are guarded. If it happens, check nothing sets
`CANDIDAS_RUN_APP=1` in your environment.

**PDF render fails.** `quarto install tinytex`. If a system TeX is also on your
`PATH`, note that quarto prefers its own TinyTeX; both are recorded in
`env/versions.json` so a version difference is not mistaken for drift.
