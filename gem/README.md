# gem/ — the enzyme- and temperature-constrained metabolic models (etcGEM)

This directory builds a genome-scale metabolic model for each of the four species, gives
every enzyme a sequence-predicted molecular weight, turnover number, optimum temperature
and melting temperature, fits the three shape parameters of the temperature response to
*C. auris*'s measured growth curve alone, freezes them, and asks what the model then
predicts for the three relatives. Figure 4 of the manuscript is read off these outputs.
`FIG4_LOCKED.md` pins every number the figure prints to the script that produced it;
`FIG4_etcgem_caption.md` is the long-form caption and the record of every objection
considered.

Everything is Python (environment: `requirements.txt`). Every script imports its paths
from `gempaths.py`, so the layout is defined once:

```
gem/
├── 01_ … 26_*.py|sh     the pipeline, numbered in running order (below)
├── gempaths.py          directory layout + the species -> (model, medium) map
├── inputs/              hand-curated inputs (tracked)
├── models/              the four SBML models (tracked, 38 MB)
├── tables/              everything a script writes (tracked; fig4/ holds Figure 4's tables)
├── audits/              the sensitivity and falsification analyses reported in the caption
├── notes/               the written record: methods for review, audit verdicts, decisions
├── external/            third-party code, weights, databases (NOT tracked; fetch_external.sh)
├── fetch_external.sh    fetches external/
└── requirements.txt     pinned Python environment
```

## Data flow

```
 curated models            proteomes                     KEGG
 iRV973 (auris)            phylo/proteomes/*.faa         ko_reaction.list
 iDC1003 (parapsilosis)    inputs/parap_uniprot.tsv      kegg_smiles.tsv
      │                        │                              │
 01 fetch ──► models/     04 kofamscan ──► inputs/kofam/*_ko.txt
      │                        │                 │
      └────────────────────────┴─── 06 build_drafts ──► models/auris_iRV973_rekeyed.xml
                                                        models/{haemulonii,duobushaemulonii}_draft.xml
                                                        tables/*_evidence.csv
                                                 │
            ┌────────────────────────────────────┼──────────────────────────────┐
     08 enzyme_mw                        09 dlkcat_pairs                 13 seq2topt_input
     tables/enzyme_mw_<sp>.csv           tables/dlkcat_pairs_<sp>.tsv     tables/seq2topt_input.csv
                                         inputs/kegg_compounds_needed.txt          │
                                                 │  12 fetch_kegg_smiles      14 run Seq2Topt + Seq2Tm
                                         10 run DLKcat                        tables/thermal_topt.csv
                                         tables/kcat_<sp>.tsv                 tables/thermal_tm.csv
                                         11 aggregate_kcat
                                         tables/kcat_reaction_<sp>.csv
            └────────────────────────────────────┼──────────────────────────────┘
                                                 │
   results/tables (respirometry) ──► 17 build_measured_tpc ──► tables/measured_tpc_honest.csv
                                                 │
                                    18 build_etcgem_tpc  (fit sig, w, P, SCALE on C. auris; freeze)
                                                 │   tables/etcgem_calib.json, etcgem_tpc_pred*.csv
              ┌──────────────┬───────────────────┼─────────────────┬────────────────┐
     19 counterfactual   20 dyn_sparse     22 thermal_sens.   23 gap_robustness   25 loocv_nested
     (uniform shifts)    (greedy + beam)   (limit per °C Tm)  (3 defects fixed)   (species held out)
              │              │                                                     24 paired_dedup_audit
              └──────────────┴──────────► 26 fig4_tables ──► tables/fig4/*.csv ──► scripts/19_fig4.R
```

## Run order

| step | script | what it does | needs | writes |
|---|---|---|---|---|
| 01 | `01_fetch_curated_models.sh` | download iRV973 (*C. auris*) and iDC1003 (*C. parapsilosis*) SBML | network | `models/auris_iRV973.xml`, `models/parapsilosis_iDC1003.xml` |
| 02 | `02_inspect_models.py` | report namespaces, compartments, biomass, growth | cobra | `notes/model_report.md` |
| 03 | `03_apply_medium_fba.py` | map the YMS medium onto a model; min/proxy/rich scenarios | cobra | `tables/medium_fba_<model>.csv` |
| 04 | `04_kofam_annotate.sh` | KofamScan KO annotation of the three Candidozyma proteomes | kofamscan, `external/kofam` | `inputs/kofam/{auris,hae,duo}_ko.txt` |
| 06 | `06_build_drafts.py` | sequence-linked GPRs from each species' KOs on the iRV973 scaffold; the same procedure re-keys *C. auris* itself (its curated GPRs use locus tags that exist in no sequence database) | cobra | `models/auris_iRV973_rekeyed.xml`, `models/*_draft.xml`, `tables/*_evidence.csv` |
| 07 | `07_rbh_orthologs.py` | reciprocal-best-hit orthologs auris ↔ each relative, with identity (for the common-network control; not used to build models) | diamond | `tables/rbh/pid_<sp>.tsv` |
| 08 | `08_enzyme_mw.py` | molecular weight of every model enzyme from its sequence | cobra, biopython | `tables/enzyme_mw_<sp>.csv` |
| 09 | `09_dlkcat_pairs.py` | every (reaction, enzyme, substrate) triple for DLKcat; currency metabolites excluded | cobra | `tables/dlkcat_pairs_<sp>.tsv`, `inputs/kegg_compounds_needed.txt`, `tables/dlkcat_input_<sp>.tsv` |
| 10 | `10_run_dlkcat.sh` | DLKcat kcat prediction per triple | torch, rdkit, `external/DLKcat` | `tables/kcat_<sp>.tsv` |
| 11 | `11_aggregate_kcat.py` | one kcat per reaction (maximum over its enzymes and substrates) | | `tables/kcat_reaction_<sp>.csv` |
| 12 | `12_fetch_kegg_smiles.py` | SMILES for every substrate from KEGG (once) | network, rdkit | `inputs/kegg_smiles.tsv` |
| 13 | `13_seq2topt_input.py` | id, sequence for every model enzyme ≤ 2000 aa | cobra | `tables/seq2topt_input.csv` |
| 14 | `14_run_seq2topt_seq2tm.sh` | Seq2Topt and Seq2Tm over that table, in file order, batch 4 | torch, `external/Seq2Topt`, checkpoints, ESM2 | `tables/thermal_topt.csv`, `tables/thermal_tm.csv` (+ raw `topt_pred.csv`, `tm_pred.csv`) |
| 15 | `15_run_seq2tm.py` | Seq2Tm over a whole proteome or with batch = 1 (the padding audit) | torch | as named on the command line |
| 16 | `16_pool_binding_test.py` | does the protein pool bind at literature P? (precondition) | cobra | stdout (recorded in `notes/POOL_BINDING_RESULT.md`) |
| 17 | `17_build_measured_tpc.py` | measured growth curves with dead wells counted as zero | `results/tables/` | `tables/measured_tpc_honest.csv` |
| 18 | `18_build_etcgem_tpc.py` | calibrate σ, w, P, SCALE on *C. auris*; predict all four | cobra | `tables/etcgem_calib.json`, `tables/etcgem_tpc_pred.csv`, `etcgem_tpc_pred_fine.csv` |
| 19 | `19_etcgem_counterfactual.py` | uniform Tm / Topt shifts: how far must a relative's enzymes move to fail at 40 °C while keeping 34 °C growth? | cobra | `tables/counterfactual_results.json`, `counterfactual_sweep.csv` |
| 20 | `20_dyn_sparse.py Tm\|Topt` | dynamic greedy + beam search over the strongest bottleneck enzymes (hours) | cobra | `tables/dyn_sparse_{Tm,Topt}.json` |
| 21 | `21_counts_at_44.py` | isolate-level detection at 44 °C with Fig 3's criterion | `results/tables/` | `tables/counts_at_44.json` |
| 22 | `22_thermal_sensitivity.py` | thermal limit per °C of Tm, calibrated and physical width | cobra | stdout (in `FIG4_LOCKED.md`) |
| 23 | `23_gap_robustness.py` | the limit with RQ, energy cycle and unfolding width corrected together | cobra | `tables/gap_robustness_results.csv` |
| 24 | `24_paired_dedup_audit.py` | paired ortholog ΔTm/ΔTopt: reaction-level vs unique pairs vs cluster bootstrap | | stdout (in `FIG4_LOCKED.md`) |
| 25 | `25_loocv_nested.py` | species-held-out validation (~16 min) | cobra | `tables/loocv_nested_results.csv` |
| 26 | `26_fig4_tables.py` | flatten everything Figure 4 draws into CSV | | `tables/fig4/*.csv` |

`audits/` holds the analyses behind the mechanism table in the caption (maintenance
energetics, oxygen, proteome allocation, amino-acid pool, energy-generating cycle,
complex GPRs, medium and kcat aggregation, sparse counterfactuals, predictor-error
propagation, the Seq2Tm padding artefact, the proteostasis pre-check, the lineage-specific
gene screen). Each states its question in its docstring. They are run from the repository
root (`python3 gem/audits/<name>.py`) and write to `tables/`.

## Reproducing the committed tables

Steps 06, 07, 08, 09, 13, 17, 21, 24 and 26 were re-run from the committed inputs on
2026-09-06 and reproduce the committed tables exactly (models: 2510/2510 GPRs; MW: all
four tables; DLKcat triples: 4042/3974/3871/4667; predictor input: 3042 rows in order).
Steps 10 and 14 (the two deep-learning predictors) were run once, on 2026-09-02, and their
outputs are committed; they are deterministic for the committed checkpoints but slow, and
Seq2Topt/Seq2Tm's output depends on batch composition (see the header of 14 and
`FIG4_LOCKED.md`), so re-run them only on the unchanged input file. Steps 18 to 25 are
deterministic FBA (GLPK) and re-run in minutes to hours (20 is the slow one).

## What the model is, in one paragraph

Flux balance analysis with a single protein-pool constraint
Σ_r (MW_r / (kcat_r · 3600 · a_r(T))) · v_r ≤ P, where the activity factor
a(T) = exp(−(T − Topt)² / 2σ²) / (1 + exp((T − Tm) / w)) is a Gaussian peak at the enzyme's
predicted optimum times a sigmoid unfolding cut-off at its predicted melting temperature.
σ, w and P are shared by all enzymes and species and fitted, with one growth scale, to
*C. auris*'s measured curve (σ ≈ 10.2 °C, w ≈ 8.8 °C, P ≈ 0.36 g/gDW); every interspecies
difference the model then predicts comes only from the sequence-predicted per-enzyme
Topt and Tm. Full method and its limitations: `notes/ETCGEM_METHODS_FOR_REVIEW.md`.
