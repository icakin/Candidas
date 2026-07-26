#!/usr/bin/env bash
# =============================================================================
# run_c2_arms.sh -- C2 PART C. Three N0 treatments, end to end, side by side.
# =============================================================================
#   bash run_c2_arms.sh [current|ramp|nobp ...]     (default: all three)
#
#   ARM 1  current   N0 = N_inoc * exp(r*delta)              -> results_C2_current/
#   ARM 2  ramp      N0 = N_inoc * exp(r*delta) * f_ramp     -> results_C2_ramp/
#   ARM 3  nobp      N0 = N_inoc                             -> results_C2_nobp/
#                    LOWER BOUND on N0, not a candidate.
#
# Everything except the N0 treatment is held fixed:
#   * the same committed trim windows, exclusions, inoculum and cell constants
#     (they come from results/tables via config.R's app_input() fallback);
#   * the same brms formulae, priors and seed (1234);
#   * the same etc-GEM capacity values -- outputs/supp_data_C1/ for every arm,
#     because C1 established the etc-GEM is not reproducible, so letting
#     capacity vary would confound the comparison;
#   * ONE seed (CANDIDAS_SEED) for the four previously-unseeded R subsample /
#     bootstrap stages, identical across arms.
#
# 02/03/06 are N0-INDEPENDENT (N0 enters only at 07). Their outputs are copied
# from results_C1/ so all three arms literally share the same inputs, which is
# what makes the growth-side invariance check meaningful. 11 is likewise
# N0-independent; 13 reads its committed output from results_C1/expression/.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")"
ROOT="$PWD"

ARMS=("$@"); [ ${#ARMS[@]} -eq 0 ] && ARMS=(current ramp nobp)

export CANDIDAS_SEED=20260726          # recorded in reports/N0_SENSITIVITY.md
export CANDIDAS_SUPP_DATA="$ROOT/cauris_etcgem/strains/eci_cauris/outputs/supp_data_C1"
export CANDIDAS_EXPRESSION_OUT="$ROOT/results_C1/expression"
export CANDIDAS_N0_RAMP_CSV="$ROOT/reports/tools/c2_ramp_n0_factors.csv"
export CANDIDAS_SKIP="02,03,06,11,14"

# outputs of 02/03 that 07 consumes
CARRY=(Oxygen_All_Long.csv Oxygen_Data_Filtered.csv Oxygen_Data_Smoothed_Trimmed.csv
       Oxygen_Trimmed_Series_Metadata.csv Oxygen_Curve_Code_Key.csv
       removed_points.csv Skipped_Series_Log.csv)

mkdir -p logs
STAMP="$(date +%Y%m%d_%H%M%S)"

for ARM in "${ARMS[@]}"; do
  OUT="$ROOT/results_C2_${ARM}"
  LOG="$ROOT/logs/c2_arm_${ARM}_${STAMP}.log"
  echo ""
  echo "======================================================================"
  echo "ARM: $ARM   ->  $OUT"
  echo "  seed        : $CANDIDAS_SEED"
  echo "  etc-GEM     : $CANDIDAS_SUPP_DATA"
  echo "  log         : $LOG"
  echo "  started     : $(date '+%Y-%m-%d %H:%M:%S')"
  echo "======================================================================"

  rm -rf "$OUT"
  mkdir -p "$OUT/tables" "$OUT/figures" "$OUT/rds"
  for f in "${CARRY[@]}"; do
    [ -f "$ROOT/results_C1/tables/$f" ] && cp "$ROOT/results_C1/tables/$f" "$OUT/tables/"
  done

  t0=$(date +%s)
  CANDIDAS_RESULTS="$OUT" CANDIDAS_N0_MODE="$ARM" \
    Rscript "$ROOT/scripts/run_all.R" 2>&1 | tee "$LOG"
  rc=${PIPESTATUS[0]}
  dt=$(( $(date +%s) - t0 ))
  printf '  ARM %s finished in %dm %02ds (exit %d)\n' "$ARM" $((dt/60)) $((dt%60)) "$rc"
  [ "$rc" -ne 0 ] && { echo "!! ARM $ARM FAILED - see $LOG"; exit "$rc"; }
done

echo ""
echo "======================================================================"
echo "all arms done"
for ARM in "${ARMS[@]}"; do
  echo "  results_C2_${ARM}/tables : $(ls "$ROOT/results_C2_${ARM}/tables"/*.csv 2>/dev/null | wc -l | tr -d ' ') csv"
done
echo "======================================================================"
