# reprocess/ — growth calls for every well

Two arms live here. Both start from `blind_traces.csv.gz`, the earliest immutable raw
oxygen exports keyed by a hashed WELLID (`blind_key.csv` maps WELLID to well; it was sealed
until the classifier arm was unblinded).

## Current: the per-well likelihood model (manuscript Figures 3-5)

| file | role |
|---|---|
| `lik_model.py` | the model: exponential-to-plateau vs linear-to-plateau, free plateau, AR(1) errors, exact likelihood, parametric-bootstrap null for the LR (499 draws) |
| `run_lik.py` | resumable 4-worker runner over every well (`python3 run_lik.py 150`, repeat until COMPLETE) |
| `lik_partial.csv`, `lik_results.csv` | raw output, one row per well |
| `lik_results_labelled.csv` | the same joined to the well labels and to the classifier's columns (suffix `_v4`) |

Growth rule used everywhere: bootstrap p < 0.01. No detectable respiration: drawdown < 1.0 mg/L.

## Archived sensitivity arm: the pre-registered blind classifier (protocols v1-v4)

`classifier.py` (frozen, hashed), `PROTOCOL_v1..v4.md` with their `.sha256`, `simulate.py`
and `sim_results*.csv` (calibration), `classify_real.py` / `run_par.py` / `run_chunk.py`
(runners), `blind_states*.csv` (blind output), `UNBLINDING_PLAN.md`, `unblinded_results.csv`
(states joined to labels). Files marked `_INVALID` or `_FAILED` are earlier protocol versions
that failed their own calibration and were archived rather than revised.

Agreement with the likelihood model: 930 of its 931 confirmed-growth wells are growth; of its
154 ambiguous wells, 90 failed only the blocked cross-validation criterion because the fit
window ran into the plateau (89 are growth), and 13 of the remaining 64 are growth; 31 of its
114 respiration-only wells show growth (window artefact, see manuscript Methods).

`40_refit_brms.R`, `refit_growth_*.rds/csv`: the manual-vs-automatic window sensitivity refit
(manuscript Supplementary Fig. 11).

`_amb_traces_check.png`: the nine ambiguous traces inspected during the diagnosis; scratch.
