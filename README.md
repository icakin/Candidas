# Candidas

Thermal performance of *Candida auris* (four clades) and its relatives *C. haemulonii*,
*C. duobushaemulonii* and *C. parapsilosis*,
from dissolved-oxygen respirometry across 22–44 °C: growth rate, per-cell respiration and
carbon-use efficiency; the carbon cost of a febrile host; the discordance between
high-temperature growth and phylogeny; and the state of every well at fever
(growth, respiration without growth) from a per-well likelihood model. The etcGEM work was
removed from the manuscript in September 2026 and lives in `etcGEMs/` and `archive/`.

## Layout

```
Candidas/
├── data/          raw respirometry: <Group>_<T>_Oxygen.xlsx as exported by the PreSens SDR
│                  reader, and the tidy .csv 01_convert_xlsx.R makes from each (8 groups x 12 T)
├── scripts/       the R pipeline, 00-21 in running order (scripts/README.md is the script map)
├── gem/           the etcGEM layer, Python, 01-26 in running order (gem/README.md is the data flow)
├── phylo/         proteomes, the phylogenomic tree, thermal-machinery census (scripts 01-04)
├── results/       everything the pipeline writes: tables/, figures/ (manuscript/ holds Figs 1-4), rds/
├── manuscript/    v1.qmd (the manuscript) -> v1.pdf; earlier drafts in manuscript/archive/; references, csl, plain summary, grant case
├── reports/       audit reports with their own code (reproduction, N0 sensitivity, scale-free CUE, clade contrasts)
├── docs/          getting started, handover notes, the pinned-environment record, prompt archive
└── archive/       the etc-GEM cut from the manuscript in August; the ITS phylogeny; superseded Python figures
```

Nothing under `results/` or `gem/tables/` is hand-edited: every file there is written by a
script, and every manuscript figure resolves to the script that draws it:

| figure | file | script | reads |
|---|---|---|---|
| Figure 1 | `results/figures/manuscript/FIG4_growth_lik.png` | `scripts/57_fig4_growth_lik.py` | `reprocess/lik_results_labelled.csv` |
| Figure 2 | `results/figures/manuscript/FIG_thermal_limits_lik.png` | `scripts/56_fig3_limits_lik.py` | `results/tables/lik_transition.csv` (56_lik_transition.py), `phylo/trees/pg_rooted.nwk` |
| Figure 3 | `results/figures/manuscript/FIG5_respiration_lik.png` | `scripts/58_fig5_respiration_lik.py` | same, plus `results/tables/phase2b_consumption.csv` |
| Figure 4 | `results/figures/manuscript/FIG1_decoupling_py.png` | `scripts/79_fig4_tpc_py.py` | posterior draws (`isolate_params_draws.csv`, 68_...R), `fig_values.csv` |
| Figure 5 | `results/figures/manuscript/FIG2_consequences_py.png` | `scripts/80_fig5_cue_py.py` | same, plus `fig_contrasts.csv`, `wells_by_isolate.csv` |
| Supp. Figs 1-4 | `results/figures/manuscript/FIG_SUPP_{relative_RG,topt_vs_fever_cost,loo_correlation,pairwise_contrasts}_py.png` | `scripts/80_fig5_cue_py.py` | |
| Supp. Fig 5 | `results/figures/manuscript/FIG3_revised_displacement.png` | `scripts/23_thermal_headroom.R`, `scripts/51_fig3_revised.py` | |
| Supp. Fig 6 | `results/figures/manuscript/SUPP_threshold_sensitivity.png` | `scripts/25_supp_threshold_sensitivity.py` | `results/tables/thermal_headroom_draws.csv` |
| Supp. Fig 7 | `results/figures/manuscript/FIG_SUPP_equilibration.png` | `scripts/76_supp_equilibration.py` | `temperature_equilibration_sensitivity.csv` (08_...R) |
| Supp. Figs 8, 9 | `results/figures/manuscript/FIG_SUPP_resp_envelope.png`, `FIG_SUPP_cue_envelope.png` | `scripts/77_supp_envelopes.py` | `isolate_params_draws.csv`, `bayes_cue_by_clade.csv` |
| Supp. Fig 10 | `results/figures/manuscript/FIG_SUPP_n0_treatments.png` | `scripts/78_supp_n0_treatments.py` | `bayes_resp_arr_E_three_treatments.csv` (15_...R) |
| Supp. Fig 11 | `results/figures/manuscript/FIG_window_sensitivity.png` | `scripts/55_fig_window_sensitivity.py` | `reprocess/refit_growth_comparison.csv` |
| Supp. Fig 12 | `results/figures/manuscript/FIG_SUPP_safety_margin.png` | `scripts/67_safety_margin_lik.py` | `results/tables/isolate_crossings_summary.csv` (65_isolate_crossings.R), `lik_transition.csv` |
| Supp. Fig 13 | `results/figures/manuscript/FIG_SUPP_od_validation.png` | `scripts/81_od_validation.py` | `data/od_validation/od_{37,40,42}.csv`, `plate_map.csv`; writes `results/tables/od_validation.csv` |
| Supp. Note 4 | (text) | `reprocess/hiO2/run_hiO2.py`, `compare_hiO2.py` | `reprocess/blind_traces.csv.gz`; writes `reprocess/hiO2/lik_hiO2_0.50.csv`, `hiO2_comparison_0.50.csv` |

Every manuscript figure is drawn in Python (`bash scripts/make_figures.sh`); R scripts 16/17 still write `fig_values.csv` and `fig_contrasts.csv`, which the Python figures read so printed numbers are unchanged. All figures share `scripts/style_common.py` (palette identical to `fig_common.R`, same font stack, same greys); rerun them on a machine with Helvetica so the PDF and R figures match.

The per-well growth calls behind Figures 3-5 come from the likelihood model in
`reprocess/lik_model.py` (run by `reprocess/run_lik.py`); the pre-registered blind classifier
(`reprocess/classifier.py`, protocols v1-v4) is retained as the sensitivity arm described in
Methods. `reprocess/README.md` says which files are current and which are archived.

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

* The likelihood-model growth calls and limits: `reprocess/lik_results_labelled.csv`,
  `results/tables/lik_transition.csv`; agreement with the blind classifier is printed by
  `scripts/56_lik_transition.py` and summarised in the manuscript Methods.
* The clade-level CUE optimum contrasts (why no within-*C. auris* spread is quoted):
  `phylo/notes/CLADE_SPREAD_VERDICT.md`, from `reports/C11_scale_free/R/05_clade_contrasts.R`.
* The scale-free CUE result the manuscript quotes: `reports/C11_scale_free/`.
* Whether the committed pipeline reproduces itself from a clean environment:
  `reports/REPRODUCTION.md`.

## Naming

Recent revisions place *C. auris* and the *haemulonii* complex in *Candidozyma* and
*C. parapsilosis* in *Lodderomyces*. The project uses *Candida* throughout, matching
`config.R` and the clinical literature.
