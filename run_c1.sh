#!/usr/bin/env bash
# C1 non-destructive full run. Every output goes to a NEW path.
set -uo pipefail
cd "$(dirname "$0")"
export CANDIDAS_RESULTS="$PWD/results_C1"
export ETCGEM_OUT_SUFFIX=_C1
export CANDIDAS_SUPP_DATA="$PWD/cauris_etcgem/strains/eci_cauris/outputs/supp_data_C1"
export CANDIDAS_EXPRESSION_OUT="$PWD/results_C1/expression"
export ETCGEM_STAGES=all
bash scripts/run_all.sh
