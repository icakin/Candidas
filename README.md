# Candidas

Thermal performance of *Candida auris* (four clades) and its relatives *C. haemulonii*,
*C. duobushaemulonii* and *C. parapsilosis* (plus *C. glabrata* in the respirometry set),
from dissolved-oxygen respirometry across 22–44 °C: growth rate, per-cell respiration and
carbon-use efficiency; the carbon cost of a febrile host; the discordance between
high-temperature growth and phylogeny; and an enzyme- and temperature-constrained
genome-scale model (etcGEM) that asks whether sequence-predicted enzyme thermal properties
can reproduce the divergence (they cannot, and Figure 4 documents by how much).

## Layout

```
Candidas/
├── data/          raw respirometry: <Group>_<T>_Oxygen.xlsx as exported by the PreSens SDR
│                  reader, and the tidy .csv 01_convert_xlsx.R makes from each (8 groups x 12 T)
├── scripts/       the R pipeline, 00-21 in running order (scripts/README.md is the script map)
├── gem/           the etcGEM layer, Python, 01-26 in running order (gem/README.md is the data flow)
├── phylo/         proteomes, the phylogenomic tree, thermal-machinery census (scripts 01-04)
├── results/       everything the pipeline writes: tables/, figures/ (manuscript/ holds Figs 1-4), rds/
├── manuscript/    v3.qmd (the single canonical text) -> v3.pdf; references, csl, plain summary, grant case
├── reports/       audit reports with their own code (reproduction, N0 sensitivity, scale-free CUE, clade contrasts)
├── docs/          getting started, handover notes, the pinned-environment record, prompt archive
└── archive/       the etc-GEM cut from the manuscript in August; the ITS phylogeny; superseded Python figures
```

Nothing under `results/` or `gem/tables/` is hand-edited: every file there is written by a
script, and every manuscript figure resolves to the script that draws it:

| figure | file | script | reads |
|---|---|---|---|
| Figure 1 | `results/figures/manuscript/FIG1_decoupling.png` | `scripts/16_fig1.R` | Bayesian fits from 09 (via `fig_common.R`) |
| Figure 2 | `results/figures/manuscript/FIG2_consequences.png` | `scripts/17_fig2.R` | same, plus 12_carbon_tax |
| Figure 3 | `results/figures/manuscript/FIG3_discordance.png` | `scripts/18_fig3.R` | per-well fits from 07, `phylo/trees/pg_rooted.nwk` |
| Figure 4 | `results/figures/manuscript/FIG4_etcgem_counterfactual.png` | `scripts/19_fig4.R` | `gem/tables/` (written by `gem/26_fig4_tables.py`) |
| Supp. Fig 1 | `results/figures/Fig_temperature_equilibration.png` | `scripts/08_…sensitivity.R` | |
| Supp. Figs 2, 3 | `results/figures/fig_bayes_resp_arrhenius.png`, `fig_bayes_cue_by_clade.png` | `scripts/11_bayesian_plots.R` | |
| Supp. Fig 4 | `results/figures/Fig_n0_treatment_panel.png` | `scripts/15_n0_treatment_panel.R` | |

## Reproducing

```
# R side: restores the pinned library (renv.lock), checks Stan
Rscript scripts/00_install.R
# whole project: data -> models -> figures -> C11 report -> manuscript PDF, logged to logs/
bash scripts/run_all.sh
# or only the figures
Rscript scripts/21_figures_all.R
```

`SETUP.md` has the system prerequisites (compilers, Stan, quarto). The Python etcGEM layer
has its own pinned environment (`gem/requirements.txt`) and is not part of `run_all.sh`:
its outputs are committed, and `gem/README.md` gives the run order and says which steps
are deterministic re-runs (minutes) and which are the two deep-learning predictors that
were run once (hours, and sensitive to input order). `gem/fetch_external.sh` fetches the
third-party code and weights those need.

## Where the numbers live

* Figure 4's every printed value, with the script that produced it: `gem/FIG4_LOCKED.md`.
  The full reasoning and every objection considered: `gem/FIG4_etcgem_caption.md`.
* The clade-level CUE optimum contrasts (why no within-*C. auris* spread is quoted):
  `phylo/notes/CLADE_SPREAD_VERDICT.md`, from `reports/C11_scale_free/R/05_clade_contrasts.R`.
* The scale-free CUE result the manuscript quotes: `reports/C11_scale_free/`.
* Whether the committed pipeline reproduces itself from a clean environment:
  `reports/REPRODUCTION.md`.

## Naming

Recent revisions place *C. auris* and the *haemulonii* complex in *Candidozyma* and
*C. parapsilosis* in *Lodderomyces*. The project uses *Candida* throughout, matching
`config.R` and the clinical literature.
