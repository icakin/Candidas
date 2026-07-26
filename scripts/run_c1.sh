#!/usr/bin/env bash
# C1 non-destructive full run. Every output goes to a NEW path.
set -uo pipefail
# lives in scripts/, but runs from the PROJECT ROOT
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export CANDIDAS_RESULTS="$PWD/runs/C1_reproduction"
# the etc-GEM outputs now live under the submodule's runs/C1/, so point at
# them explicitly rather than relying on the _C1 name suffix
export ETCGEM_OUT_SUPP="$PWD/cauris_etcgem/runs/C1/supp_data"
export ETCGEM_OUT_FIG3="$PWD/cauris_etcgem/runs/C1/figure3_data"
export ETCGEM_FITNPY="$PWD/cauris_etcgem/runs/C1/_percladefit.npy"
export CANDIDAS_SUPP_DATA="$PWD/cauris_etcgem/runs/C1/supp_data"
export CANDIDAS_EXPRESSION_OUT="$PWD/runs/C1_reproduction/expression"
export ETCGEM_STAGES=all
bash scripts/run_all.sh
