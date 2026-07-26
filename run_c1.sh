#!/usr/bin/env bash
# C1 non-destructive full run. Every output goes to a NEW path.
set -uo pipefail
cd "$(dirname "$0")"
export CANDIDAS_RESULTS="$PWD/runs/C1_reproduction"
export ETCGEM_OUT_SUFFIX=_C1
export CANDIDAS_SUPP_DATA="$PWD/cauris_etcgem/runs/C1/supp_data"
export CANDIDAS_EXPRESSION_OUT="$PWD/runs/C1_reproduction/expression"
export ETCGEM_STAGES=all
bash scripts/run_all.sh
