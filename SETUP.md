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

* **macOS** — see [System prerequisites](#system-prerequisites) below. `xcode-select
  --install` is necessary but **not sufficient**.
* **Linux** — `sudo apt install build-essential` (needs `g++` and GNU make). Not
  otherwise documented here; nobody has run this project on Linux.
* **Windows** — install Rtools matching your R version. Likewise undocumented.

> **An earlier version of this file said "every CRAN package this project uses
> resolves to an arm64 binary, so no `gfortran` is needed."** That was true on the
> machine that captured `renv.lock`, at the moment it was captured, and it is
> false for anyone restoring the lockfile later. See below for why, and what you
> actually need.

---

## System prerequisites (macOS) {#system-prerequisites}

**Read this before `renv::restore()`.** On a clean Mac the restore fails without
these. `scripts/00_install.R` now checks for each of them and stops with the
specific missing item rather than letting a build fail deep in a dependency tree.

### Why a *pinned* environment needs a compiler at all

`renv.lock` pins exact versions. **CRAN's macOS binary repository only ever
carries the *current* version of each package.** So any pin that has since been
superseded has no binary, and `renv` falls back to the source tarball — which
needs whatever that package needs to build.

Measured against this lockfile on 2026-08-10, R 4.5.2, `aarch64-apple-darwin20`:

| | packages |
|---|---|
| pinned versions with an exact CRAN **binary** today | **97** of 170 |
| pinned versions that would **build from source** | **73** of 170 |

That number grows as the lockfile ages. It is not a defect in the lockfile; it is
how CRAN binaries work.

### 1. Xcode command line tools — required, always

```bash
xcode-select --install
```

Gives `clang`, `clang++` and `make`. Every source build needs these, and so does
Stan.

### 2. The CRAN gfortran toolchain — required

Six pinned packages link the Fortran runtime **and** currently have no binary at
their pinned version, so they build from source:

| package | pin | what needs it |
|---|---|---|
| `Matrix` | 1.7-4 | base of almost everything |
| `mgcv` | 1.9-3 | dependency of `brms` |
| `nlme` | 3.1-168 | dependency of `mgcv` |
| `mvtnorm` | 1.3-3 | `brms` |
| `nleqslv` | 3.3.5 | `brms` |
| `edgeR` | 4.8.2 | `13_capacity_expression.R` |

(`minpack.lm`, `quadprog`, `RcppEigen` and `statmod` also link Fortran but do
currently have binaries at their pins. That can change at any time — install the
toolchain regardless.)

Download the **CRAN-built** toolchain, not Homebrew's:

<https://mac.r-project.org/tools/>

— specifically the *GNU Fortran* installer for your architecture
(`gfortran-…-universal.pkg`).

> **Homebrew's `gcc`/`gfortran` is not a substitute.** CRAN-built R is compiled
> against the CRAN toolchain and expects its `libgfortran`/`libquadmath` at the
> paths and ABI that toolchain provides. Mixing Homebrew's produces link errors,
> or worse, packages that load and misbehave.

### 3. Homebrew system libraries — required

Seven pinned packages need a named system library and currently build from source:

| Homebrew formula | R packages it serves |
|---|---|
| `openssl` | `openssl` (2.3.5), `curl` (7.0.0) |
| `freetype` | `ragg` (1.5.0), `textshaping` (1.0.4), `systemfonts` (1.3.1) |
| `harfbuzz`, `fribidi` | `textshaping` |
| `fontconfig` | `systemfonts` |
| `libpng`, `jpeg-turbo`, `libtiff` | `ragg` |
| `libxml2` | `xml2` (1.5.2) |
| `icu4c` | `stringi` (1.8.7) — bundles its own ICU if it must, slowly |

```bash
brew install openssl freetype harfbuzz fribidi fontconfig libpng jpeg-turbo libtiff libxml2 icu4c
```

`openssl` and `freetype` are the two that were actually reported as failing; the
rest are the remaining source builds' declared needs and are cheap to install at
the same time.

### 4. Stan — the heaviest single dependency {#stan-prereq}

**The backend is `rstan`.** `09_bayesian_models.R` calls `brms::brm()` without a
`backend` argument and nothing in the project sets `options(brms.backend)`, so
brms falls back to its default, which is rstan. **`cmdstanr` is not in
`renv.lock`, is not installed, and is not used** — it is not a CRAN package, so
`renv` could not restore it from the lockfile even if it were listed.

Good news: `rstan` **2.32.7** and `StanHeaders` **2.32.10** both have exact CRAN
binaries at their pinned versions today, so the notorious multi-hour rstan
compile is avoidable. What does still build from source in the Stan stack is
`Rcpp` (1.1.1), `RcppParallel` (5.1.11-1) and `QuickJSR` (1.9.0) — all C++ only,
needing nothing beyond the Xcode tools, but `RcppParallel` builds a bundled TBB
and takes a few minutes.

Verify Stan **before** you start a pipeline run:

```bash
Rscript scripts/00_install.R          # compiles and samples a toy model
```

Expected: about **40 seconds**, ending with a line reporting the recovered
parameters. Measured on the reference machine: 36.7 s, recovering
`mu = 3.0837` (truth 3) and `sigma = 0.8454` (truth 1), with rstan 2.32.7 /
StanHeaders 2.32.10 / Stan 2.32.2.

A restore that "succeeds" while Stan cannot compile is worse than one that fails,
because the failure then surfaces an hour into `09_bayesian_models.R`. That is
why the check is part of the installer and not an optional extra.

### 5. Bioconductor — a different failure mode {#bioconductor-prereq}

`edgeR` is Bioconductor, not CRAN, and `13_capacity_expression.R` needs it. This
fails differently from a compile problem: without Bioconductor repositories
configured, `renv::restore()` cannot **find** the package at all, whatever
toolchain you have.

What the lockfile records:

- the **release is pinned**: `renv.lock` carries `"Bioconductor": {"Version": "3.22"}`;
- `BiocManager` 1.30.27 is itself in the lockfile;
- three packages have `"Source": "Bioconductor"` — `BiocVersion` 3.22.0,
  `edgeR` 4.8.2, `limma` 3.66.0.

`renv` uses `BiocManager` to configure the Bioconductor repositories for the
pinned release, so **install `BiocManager` first** if the bootstrap has not
already:

```r
install.packages("BiocManager")
BiocManager::install(version = "3.22")   # configure, do not upgrade
```

`scripts/00_install.R` does this in the right order for you.

> **A finding worth knowing, not fixed here.** The lockfile records the
> Bioconductor packages' `Repository` as `https://bioc-release.r-universe.dev`.
> That mirror serves **none** of the three at their pinned versions as a macOS
> arm64 binary, so they build from source — which is why `edgeR` lands in the
> gfortran list above. The canonical mirror
> `https://bioconductor.org/packages/3.22/bioc` has **all three at exactly the
> pinned versions** as arm64 binaries. `00_install.R` therefore puts the
> canonical Bioconductor repository ahead of the recorded one when resolving.
> Correcting the `Repository` field in `renv.lock` would be the proper fix, but
> that is an edit to the pin file and belongs in its own change.

### 6. Known issue: conda breaks the `curl` build {#conda-libkrb5}

**Symptom.** `renv::restore()` fails while building `curl`, with a linker error
mentioning `libkrb5`, `gssapi`, or a Kerberos symbol.

**Cause.** An active conda environment puts its own `lib` directory ahead of the
system libraries on `PATH`/`LDFLAGS`. `curl` then compiles against conda's
`libkrb5` while linking against the system `libcurl`, and the two disagree.

**Workaround.** Take conda out of the way for the duration of the restore:

```bash
conda deactivate            # repeat until no (env) prefix remains
# if conda auto-activates 'base' from your shell profile:
conda config --set auto_activate_base false
# or, for one shell only:
env PATH="/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin" Rscript scripts/00_install.R
```

Re-activate conda afterwards; nothing in this project needs it.

### What has and has not been verified

**Verified on the reference machine** (the one in `env/versions.json`): the
package-by-package analysis above — which packages contain compiled code, which
link Fortran or a named system library, and which pinned versions do or do not
have a CRAN binary today — was measured directly from `renv.lock` and from the
restored library with `otool -L`. The Stan toy-model check was run and passed.

**NOT verified here: that `renv::restore()` now succeeds on a clean Mac.** No
clean machine was available. The prerequisite list is derived from a reported
clean-Mac failure plus the dependency analysis above; it is our best
reconstruction, not a tested procedure. If you are the person with the clean
machine, please work through this section and report anything missing or
unnecessary — that feedback is the only way this gets confirmed.

---

## Install

Clone **with submodules** — the etc-GEM model lives in one:

```bash
git clone --recursive https://github.com/icakin/Candidas.git
cd Candidas
# already cloned without --recursive?
git submodule update --init --recursive
```

Then — **after** working through
[System prerequisites](#system-prerequisites), which is not optional on macOS:

```bash
# 1. R side  (checks prerequisites, restores renv.lock, verifies Stan)
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

See [Bioconductor](#bioconductor-prereq) for why the recorded repository makes
edgeR build from source, and how `00_install.R` works around it.

**`renv::restore()` fails building a Fortran package** (`Matrix`, `mgcv`,
`nlme`, `mvtnorm`, `nleqslv`, `edgeR`) — install the CRAN gfortran toolchain,
not Homebrew's. See [System prerequisites](#system-prerequisites) §2.

**`renv::restore()` fails building `openssl`, `curl`, `ragg`, `textshaping`,
`systemfonts`, `xml2` or `stringi`** — a Homebrew system library is missing. See
[System prerequisites](#system-prerequisites) §3.

**`curl` fails to link, mentioning `libkrb5` or Kerberos** — that is conda on
`PATH`. See [Known issue](#conda-libkrb5).

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
