#!/usr/bin/env bash
# =============================================================================
# 14_run_seq2topt_seq2tm.sh -- predict Topt and Tm for every model enzyme from sequence
# =============================================================================
#   bash gem/14_run_seq2topt_seq2tm.sh          (~15 min on a laptop CPU for 3042 proteins)
#
# Predictors: Seq2Topt and Seq2Tm (Qiu, S. et al. 2025, Brief Bioinform 26(2), bbaf114),
# https://github.com/SizheQiu/Seq2Topt release v1.0.0. Both are an attention head on
# ESM2-t6-8M layer-6 embeddings; checkpoints model_topt_window=3_r2=0.57.pth (reported
# test RMSE 12.3 C) and model_tm_window=3_r2=0.76.pth (RMSE 7.6 C). gem/fetch_external.sh
# clones the code, downloads the checkpoints and the ESM2 weights.
#
# Input   gem/tables/seq2topt_input.csv           from 13_seq2topt_input.py (id, sequence)
# Output  gem/tables/topt_pred.csv, tm_pred.csv   the predictors' raw output (id = row number)
#         gem/tables/thermal_topt.csv             id (species|gene), sequence, pred_topt
#         gem/tables/thermal_tm.csv               id (species|gene), sequence, pred_tm
#
# The upstream scripts run from Seq2Topt/code (they hard-code ../../large_model_pth/) and
# process the input in file order in fixed batches of four WITHOUT an attention mask, so a
# protein's prediction depends on which proteins share its batch (documented in
# gem/FIG4_LOCKED.md; measured in gem/audits/padding_sensitivity.py). Re-running this
# script on the same input file therefore reproduces the committed tables exactly; running
# it on a re-ordered or subsetted input does not. Do not sort or filter the input.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S2T="$HERE/external/Seq2Topt"
TAB="$HERE/tables"
PY="${PYTHON:-python3}"
export TORCH_HOME="${TORCH_HOME:-$HERE/external/torch_home}"   # esm2_t6_8M_UR50D.pt lives in $TORCH_HOME/hub/checkpoints

[ -d "$S2T/code" ] || { echo "!! $S2T not found. Run gem/fetch_external.sh first."; exit 1; }
[ -f "$HERE/external/large_model_pth/model_topt_window=3_r2=0.57.pth" ] || { echo "!! checkpoints missing; run gem/fetch_external.sh"; exit 1; }
[ -f "$TAB/seq2topt_input.csv" ] || { echo "!! run 13_seq2topt_input.py first"; exit 1; }
# upstream expects the checkpoints at Seq2Topt/large_model_pth relative to code/
[ -e "$S2T/large_model_pth" ] || ln -s "$HERE/external/large_model_pth" "$S2T/large_model_pth"

cd "$S2T/code"
echo "== Seq2Topt";  "$PY" seq2topt.py --input "$TAB/seq2topt_input.csv" --output "$TAB/topt_pred"
echo "== Seq2Tm";    "$PY" seq2tm.py   --input "$TAB/seq2topt_input.csv" --output "$TAB/tm_pred"

# put the species|gene ids back (the predictors return the row number as id)
"$PY" - "$TAB" <<'EOF'
import sys, pandas as pd
tab = sys.argv[1]
inp = pd.read_csv(f'{tab}/seq2topt_input.csv')
for raw, col, out in [('topt_pred.csv', 'pred_topt', 'thermal_topt.csv'), ('tm_pred.csv', 'pred_tm', 'thermal_tm.csv')]:
    p = pd.read_csv(f'{tab}/{raw}')
    assert len(p) == len(inp) and (p.sequence.values == inp.sequence.values).all(), f'{raw} does not line up with the input'
    pd.DataFrame({'id': inp.id, 'sequence': inp.sequence, col: p[col]}).to_csv(f'{tab}/{out}', index=False)
    print(f'{out}: {len(p)} rows, median {p[col].median():.1f} C')
EOF
echo "next: python3 gem/16_pool_binding_test.py, then 17 and 18"
