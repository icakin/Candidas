# Candidas

Thermal-performance and carbon-economy analysis of *Candida auris* (four clades)
and *C. parapsilosis*, from dissolved-oxygen respirometry through a
sequence-grounded, enzyme- and temperature-constrained genome-scale model (etc-GEM).

## Layout

| path | what it is |
|---|---|
| `data/` | raw instrument data (`original_data/` = exactly as received) and cleaned inputs |
| `scripts/` | the analysis pipeline (R, numbered 01–14, + `14_schematic.py`); `config.R` holds shared paths and is sourced by every script. `run_all.sh` runs everything; `run_c1.sh` and `run_c2_arms.sh` are the non-destructive re-runs behind the reports |
| **`results/`** | **CANONICAL.** The shipped outputs: `figures/`, `tables/`, `rds/` (cached model fits). Publication figures live in `results/figures/manuscript/`. This is what the manuscript renders against |
| **`runs/`** | **NON-CANONICAL.** One tree per analysis prompt, kept as evidence for a report. Each is manifested in [`runs/MANIFEST.md`](runs/MANIFEST.md) with its producer, command, environment, purpose and table checksums. Nothing here is read by the pipeline unless a `CANDIDAS_*` variable points at it |
| **`reports/`** | the written findings: `REPRODUCTION.md`, `RUNBOOK.md`, `N0_SENSITIVITY.md`, the rendered technical report in `n0/`, figures in `figures_N0/`, and the analysis tools in `tools/` |
| `manuscript/` | the write-up; `manuscript/draft/` is a Quarto project rendering to Word + PDF (main, supplement, and a combined document) |
| `env/`, `renv.lock` | the pinned environment and the pre-run checksum baselines |
| `logs/` | per-stage run logs (not tracked) |
| `genomes/` | clade reference assemblies (not tracked; re-downloadable from NCBI) |

Nothing in `results/` or `runs/` is hand-edited — it is all produced by
`scripts/`. **`results/` is the one the paper cites; `runs/` exists so that a
claim in a report can be checked against the tree that produced it.**

Getting the pipeline to run is not the same as reproducing the numbers — see
[`reports/REPRODUCTION.md`](reports/REPRODUCTION.md) before trusting any
regenerated value.

## Setup and running

```
Rscript scripts/00_install.R                     # R side, verifies Stan
cd cauris_etcgem && python3 -m venv .venv \
  && .venv/bin/python -m pip install -r requirements.lock && cd ..
quarto install tinytex                           # once
bash scripts/run_all.sh                          # ~49 min end to end
```

Full detail in [`SETUP.md`](SETUP.md); what to expect at each stage in
[`reports/RUNBOOK.md`](reports/RUNBOOK.md).

## Related repositories (model code)

The genome-scale model code is kept in separate repositories, not copied into this one:

- **etcGEMs** — base etc-GEM package by G. Yvon-Durocher. Clone alongside this repo:
  `git clone https://github.com/GabYvonDurocher/etcGEMs`
- **cauris_etcgem** — the *C. auris* application of the model (the actual model runs):
  <https://github.com/icakin/cauris_etcgem> — linked here as a git submodule, so
  `git clone --recursive` fetches it too.

## Rendering the manuscript

```
cd manuscript/draft
quarto install tinytex   # once, provides the PDF engine
quarto render            # -> _output/ : manuscript, supplementary, manuscript_combined (.docx + .pdf)
```
