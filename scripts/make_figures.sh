#!/usr/bin/env bash
# Regenerate every manuscript figure (v10) in one aesthetic, all in Python.
# Prerequisites written once by R: results/tables/isolate_params_draws.csv (68_isolate_params_dump.R),
# fig_values.csv + fig_contrasts.csv (17_fig2.R). Run on the Mac so Helvetica is used.
set -e; cd "$(dirname "$0")/.."
for s in 57_fig4_growth_lik 56_fig3_limits_lik 58_fig5_respiration_lik 79_fig4_tpc_py 80_fig5_cue_py \
         51_fig3_revised 25_supp_threshold_sensitivity 76_supp_equilibration 77_supp_envelopes \
         78_supp_n0_treatments 55_fig_window_sensitivity 67_safety_margin_lik; do
  python3 scripts/$s.py
done
echo "all figures regenerated"
python3 scripts/81_od_validation.py
python3 reprocess/hiO2/compare_hiO2.py 0.5
