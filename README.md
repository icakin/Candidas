# Candidas

Thermal-performance and carbon-economy analysis of *Candida auris* (four clades)
and *C. parapsilosis* from dissolved-oxygen respirometry: growth rate, per-cell
respiration and carbon-use efficiency across 22–44 °C.

## Layout

- `data/` — raw instrument data (`original_data/` = exactly as received) and cleaned inputs
- `scripts/` — analysis pipeline (R, numbered 00–18, plus `17_schematic.py`); `config.R` holds shared paths, sourced by every script. See `scripts/README.md` for what each does and which are out of scope
- `results/` — generated outputs: `figures/`, `tables/`, `rds/` (cached model fits). Publication figures live in `results/figures/manuscript/`
- `manuscript/` — the write-up: `v3.qmd` is the single canonical version, rendering to `v3.pdf`
- `outputs/supp_data/` — vendored etc-GEM outputs, retained for the separate methods thread
- `genomes/` — clade reference assemblies (not tracked; re-downloadable from NCBI)

Nothing in `results/` is hand-edited — it is all produced by `scripts/`.

## The etc-GEM model code

**The genome-scale model is not part of the current manuscript.** An audit found
that it was built once from a single proteome and reused for every clade, that its
fitted capacity parameter restates the growth data it was fitted to, and that a
three-parameter empirical curve reproduces the fit equally well; a flux-variability
analysis had already shown the model permits rather than predicts respiration.
`scripts/README.md` records the findings in full and what a proper treatment would
require. The code and its outputs are retained here for a separate methods paper.

The model code itself lives in separate repositories, not copied into this one:

- **etcGEMs** — base etc-GEM package by G. Yvon-Durocher. Clone alongside this repo:
  `git clone https://github.com/GabYvonDurocher/etcGEMs`
- **cauris_etcgem** — the *C. auris* application of the model (the actual model runs):
  <https://github.com/icakin/cauris_etcgem> — linked here as a git submodule, so
  `git clone --recursive` fetches it too.

## Rendering the manuscript

```
cd manuscript
quarto install tinytex   # once, provides the PDF engine
quarto render v3.qmd     # -> v3.pdf
```
