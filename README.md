# Candidas

Thermal-performance and carbon-economy analysis of *Candida auris* (four clades)
and *C. parapsilosis*, from dissolved-oxygen respirometry through a
sequence-grounded, enzyme- and temperature-constrained genome-scale model (etc-GEM).

## Layout

- `data/` — raw instrument data (`original_data/` = exactly as received) and cleaned inputs
- `scripts/` — analysis pipeline (R, numbered 01–14, + `14_schematic.py`); `config.R` holds shared paths, sourced by every script
- `results/` — generated outputs: `figures/`, `tables/`, `rds/` (cached model fits). Publication figures live in `results/figures/manuscript/`
- `manuscript/` — the write-up; `manuscript/draft/` is a Quarto project rendering to Word + PDF (main, supplement, and a combined document)
- `genomes/` — clade reference assemblies (not tracked; re-downloadable from NCBI)

Nothing in `results/` is hand-edited — it is all produced by `scripts/`.

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
