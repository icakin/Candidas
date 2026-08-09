# Setup — Candidas

Everything you need to take a machine that has **only R, Python 3 and git** to a
state where `bash scripts/run_all.sh` runs the whole project unattended.

If you just want the commands, read [Install](#install) and [Run](#run).

---

## Prerequisites

| Thing | Minimum | What this project was built and verified with |
|---|---|---|
| **R** | 4.4 | **4.5.2** (`aarch64-apple-darwin20`) |
| **A C++17 toolchain** | any working one | Apple clang 17.0.0 (Xcode CLT); `R CMD config CXX17` → `clang++ -arch arm64` |
| **Python** | 3.9 | **3.9.6** (macOS system `/usr/bin/python3`, arm64) — with `numpy` and `matplotlib` |
| **git** | 2.x | 2.50.1 |
| **Quarto** | 1.4 | **1.8.27** (bundles pandoc 3.6.3) |
| **A LaTeX engine** | xelatex | quarto's own TinyTeX (XeTeX 3.141592653-2.6-0.999998, TeX Live 2026) |

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

Clone **with submodules** — the etc-GEM model lives in one:

```bash
git clone --recursive https://github.com/icakin/Candidas.git
cd Candidas
# already cloned without --recursive?
git submodule update --init --recursive
```

Then:

```bash
# 1. R side  (restores renv.lock, then verifies Stan by compiling a toy model)
Rscript scripts/00_install.R

# 2. Python side  (only 17_schematic.py needs Python in a default run)
python3 -m pip install numpy matplotlib
```

and, once, the PDF engine:

```bash
quarto install tinytex
```

### What gets pinned, and where

| Layer | Pin file | Notes |
|---|---|---|
| R packages | `renv.lock` | 170 packages, R 4.5.2, CRAN + **Bioconductor 3.22** (edgeR 4.8.2, limma 3.66.0). `renv::restore()` handles both. Verified to cover every package the eighteen numbered scripts load. |
| Stan | via `renv.lock` | rstan 2.32.7 / StanHeaders 2.32.10 / Stan 2.32.2. **`09_bayesian_models.R` never sets `backend`, so brms uses its default: rstan.** cmdstanr is not used. |
| Everything, as measured | `env/versions.json` | Every version above, read from the live environment. |
| Committed inputs and outputs | `env/baseline_checksums_*.txt` | SHA-256 of every tracked file under `data/`, `results/`, and the etc-GEM outputs. Re-run the pipeline and diff against these to prove nothing moved. |

**Not pinned: the Python side.** There is no Python lockfile on `main`.
`17_schematic.py` needs only numpy and matplotlib, so that is low risk — but the
etc-GEM engine is a different matter, and it is worth being explicit about it:

> `cauris_etcgem/requirements.txt` asks for
> `pip install git+https://github.com/GabYvonDurocher/etcGEMs` with **no ref**. That
> tracks `main`, so the model can change under the committed results with no version
> bump, no tag and no record. A `requirements.lock` that pins the engine to a commit
> exists, but on a **later submodule commit than the one this repository points at**,
> so it is not reachable from here. Until the submodule pointer is moved, regenerating
> the etc-GEM outputs is an unpinned operation.
>
> This is why `scripts/run_all.sh` leaves the etc-GEM stage **off by default**. The
> outputs it would produce are vendored at `outputs/supp_data/` and are what
> `14_main_figures.R` and `16_supplementary_figures.R` actually read.

---

## Run

```bash
bash scripts/run_all.sh
```

Fully unattended: no prompts, no browser. Per-stage timing, with a tee'd log in
`logs/`. The stages are:

| Stage | What it runs |
|---|---|
| 0 (optional, off) | etc-GEM Python pipeline — `RUN_ETCGEM=1` to enable, see above |
| 1 | `Rscript scripts/run_all.R` → **02 → 03 → 06 → 07 → 08 → 09 → 11 → 12 → 14 → 15** |
| 2 | `13_capacity_expression.R` |
| 3 | `16_supplementary_figures.R` |
| 4 | `17_schematic.py` |
| 5 | `quarto render` in `manuscript/draft` → `_output/` |

Knobs (all optional): `RUN_ETCGEM`, `ETCGEM_STAGES`, `SKIP_R`, `SKIP_RENDER`,
`LOG_DIR`.

### The four click-driven apps (01, 04, 05, 06)

These are Shiny apps. Their outputs are **committed inputs** to the pipeline and
are treated as data — an unattended run never regenerates them:

| App | Output it owns |
|---|---|
| `01_convert_xlsx.R` | `data/*_Oxygen.csv`, `results/tables/otu_names.csv` |
| `04_trim_selector.R` | `results/tables/manual_fit_windows.csv`, `plot_exclude_points.csv` |
| `05_cell_sizes.R` | `results/tables/otu_cell_sizes.csv` |
| `06_inoculation.R` | `results/tables/otu_inoc.csv` |

Sourcing any of them defines the app and returns. Before this was guarded,
`run_all.R` source()d 06 and blocked forever on its Shiny server, which is why
the master runner could never complete unattended.

To actually open one:

```bash
Rscript scripts/04_trim_selector.R --app     # or CANDIDAS_RUN_APP=1 Rscript ...
```

In RStudio, "Run App" works as it always did.

### The two diagnostics (10, 18) — run these deliberately, never in a batch

`10_n0_term_test.R` and `18_n0_treatment_panel.R` ask how much the between-clade
respiration ordering leans on the N0 back-projection. They answer it by **refitting**
the Bayesian respiration model two and three times respectively, and each alternative
fit is written over `results/tables/bayes_resp_arr_summary.csv` and `results/rds/`.
They restore `derived_N0_R_results_with_carbon.csv` but not those. So:

```bash
Rscript scripts/10_n0_term_test.R          # ~2 Bayesian fits
Rscript scripts/18_n0_treatment_panel.R    # ~3 Bayesian fits; Supplementary Fig. 6
# then put the published fit back:
Rscript scripts/09_bayesian_models.R && Rscript scripts/11_bayesian_plots.R \
  && Rscript scripts/12_carbon_tax.R && Rscript scripts/14_main_figures.R \
  && Rscript scripts/15_uncertainty_bands.R
```

`run_all.sh` deliberately leaves both out for that reason, and says so in its header.

---

## Expected runtime and disk

Measured on the machine in `env/versions.json` (Apple silicon, 16 cores),
`RUN_ETCGEM` off.

| Stage | Wall time |
|---|---|
| `Rscript scripts/00_install.R` — first run, from nothing | 10–40 min (rstan/StanHeaders compile) |
| `00_install.R` — subsequent runs | < 1 min (`--no-stan`: seconds) |
| R 02 → 03 (reshape + spline trimming, 1 080 series) | ~2 min |
| R 06 → 07 (inoculation, per-series fits, descriptive plots) | ~2 min |
| R 08 (temperature-equilibration sensitivity) | ~30 s |
| **R 09 (brms: 3 models × 4 chains × 4000 iter)** | **~8 min on 16 cores; budget ~1 h on a slow machine** |
| R 11 / 12 / 14 / 15 | ~1 min total |
| R 13 (capacity expression, edgeR) | ~30 s |
| R 16 (supplementary figures) | ~20 s |
| 17 (schematic.py) | ~5 s |
| `quarto render` (docx + pdf) | ~1 min |
| **Total (see `logs/run_all_*.log` for the measured figure)** | **see the SUMMARY block at the end of the log** |

If `RUN_ETCGEM=1`, add roughly 35 min: the `calibrate` (differential evolution,
4 clades) and `bayes` (15³ emulator grid, then 4 × 32-walker × 4000-step emcee)
stages dominate.

Disk:

| Item | Size |
|---|---|
| repository (with submodule) | ~350 MB |
| `renv/library` | ~1.5 GB (mostly symlinks into the renv cache) |
| TinyTeX | ~500 MB |
| `manuscript/draft/_output/` | ~20 MB |
| **total, comfortably** | **~3 GB** |

---

## Troubleshooting

**`supp_data/ not found` from 14 or 16.** The vendored copy at `outputs/supp_data/`
is missing, or the submodule is not populated:

```bash
git submodule update --init --recursive
```

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

**`cannot open file '.../config.R'`.** That was the `--file=` bug: run from a
project path containing a space and the script directory came back mangled
(R substitutes `~+~` for every space in `--file=`). Every script now un-mangles
it. If you see this again, the script you are running has an unpatched
script-directory block.

**PDF render fails.** `quarto install tinytex`. If a system TeX is also on your
`PATH`, note that quarto prefers its own TinyTeX; both are recorded in
`env/versions.json` so a version difference is not mistaken for drift.

**The manuscript did not pick up a regenerated figure.** It cannot: the v3
manuscript embeds extracted bitmaps rather than referencing
`results/figures/`. See `reports/MANUSCRIPT_FIGURE_LINKAGE.md` — that is a known
defect with a written remediation spec, not a configuration mistake.
