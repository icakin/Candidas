#!/usr/bin/env bash
# =============================================================================
# run_all.sh - the whole project, unattended, in one command
# =============================================================================
#   bash scripts/run_all.sh
#
# Everything is tee'd to logs/ and every stage is timed. Three stages:
#
#   1. R pipeline (scripts/run_all.R): 02 03 06 07 08 09 11 12 13, then Figures 1-4
#      (16-19). The four Shiny apps (01, 04, 05, 06) are not launched: their outputs
#      (data/*_Oxygen.csv, results/tables/{otu_names,manual_fit_windows,
#      plot_exclude_points,otu_cell_sizes,otu_inoc}.csv) are committed INPUTS.
#      10, 14 and 15 are diagnostics that refit the respiration model under alternative
#      treatments and are deliberately NOT run here (see run_all.R's header).
#   2. C11 scale-free report (reports/C11_scale_free/R/01-05): the manuscript's
#      scale-free numbers are produced there, so it regenerates with the pipeline.
#   3. quarto render of manuscript/v3.qmd.
#
# THE PYTHON etcGEM LAYER (gem/) IS NOT RUN HERE. Figure 4 (19_fig4.R) reads the tables
# it writes, which are committed. Regenerating them is a deliberate act with its own
# environment (gem/requirements.txt) and run order (gem/README.md); the sequence
# predictors depend on batch composition, so a casual re-run could silently move the
# committed numbers. RUN_GEM=1 runs the cheap, deterministic part (17 -> 26) first.
#
# ENVIRONMENT (all optional):
#   RUN_GEM=1        regenerate gem/tables/ from the committed models and predictor
#                    outputs before the R pipeline (needs gem/requirements.txt tier 1)
#   SKIP_R=1         skip the R pipeline
#   SKIP_C11=1       skip the C11 scale-free report
#   SKIP_RENDER=1    skip the manuscript render (and the C11 report render)
#   LOG_DIR          log destination        (default: <root>/logs)
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$ROOT"

# Nothing in this run may block on a Shiny server. The four app scripts read
# this and define their app instead of launching it.
export CANDIDAS_HEADLESS=1

LOG_DIR="${LOG_DIR:-$ROOT/logs}"
mkdir -p "$LOG_DIR"

STAMP="$(date +%Y%m%d_%H%M%S)"
MAIN_LOG="$LOG_DIR/run_all_${STAMP}.log"

# --- pretty + timing helpers -------------------------------------------------
RUN_T0=$(date +%s)
declare -a STAGE_NAMES=() STAGE_SECS=() STAGE_RC=()

hms() {  # seconds -> 1h 02m 03s
  local s=$1
  printf '%dh %02dm %02ds' $((s/3600)) $(((s%3600)/60)) $((s%60))
}

banner() {
  echo ""
  echo "======================================================================"
  echo "$*"
  echo "======================================================================"
}

# stage <name> <logfile> <command...>
stage() {
  local name="$1"; shift
  local log="$1"; shift
  banner "STAGE: $name   ($(date '+%Y-%m-%d %H:%M:%S'))"
  echo "  log -> $log"
  local t0 rc
  t0=$(date +%s)
  set +e
  ( "$@" ) 2>&1 | tee "$log"
  rc=${PIPESTATUS[0]}
  set -e
  local dt=$(( $(date +%s) - t0 ))
  STAGE_NAMES+=("$name"); STAGE_SECS+=("$dt"); STAGE_RC+=("$rc")
  echo "--- $name finished in $(hms $dt)  (exit $rc)"
  if [ "$rc" -ne 0 ]; then
    echo ""
    echo "!! STAGE FAILED: $name (exit $rc). See $log"
    summary
    exit "$rc"
  fi
}

summary() {
  banner "SUMMARY"
  local i
  for i in "${!STAGE_NAMES[@]}"; do
    printf '  %-40s %-14s exit %s\n' \
      "${STAGE_NAMES[$i]}" "$(hms "${STAGE_SECS[$i]}")" "${STAGE_RC[$i]}"
  done
  printf '  %-40s %-14s\n' "TOTAL" "$(hms $(( $(date +%s) - RUN_T0 )))"
  echo "  full log: $MAIN_LOG"
}

# --- everything below is captured to MAIN_LOG as well ------------------------
exec > >(tee -a "$MAIN_LOG") 2>&1

banner "Candidas full pipeline   $(date '+%Y-%m-%d %H:%M:%S')"
echo "  root        : $ROOT"
echo "  RUN_GEM     : ${RUN_GEM:-0}"
echo "  logs        : $LOG_DIR"

# ---------------------------------------------------------------------------
# 0. etcGEM tables (python) - optional, off by default
# ---------------------------------------------------------------------------
if [ "${RUN_GEM:-0}" = "1" ]; then
  PY="${PYTHON:-python3}"
  for s in 17_build_measured_tpc 18_build_etcgem_tpc 19_etcgem_counterfactual 21_counts_at_44 26_fig4_tables; do
    stage "0/4 gem/$s.py" "$LOG_DIR/gem_${s}_${STAMP}.log" "$PY" "$ROOT/gem/$s.py"
  done
  echo "  (20_dyn_sparse.py, hours, and the predictors are not re-run here; see gem/README.md)"
else
  echo ""; echo "etcGEM stage skipped (RUN_GEM=1 to regenerate gem/tables/ from the committed models)."
fi

# ---------------------------------------------------------------------------
# 1. R pipeline: 02 03 06 07 08 09 11 12 13 16 17 18 19  (run_all.R)
# ---------------------------------------------------------------------------
if [ "${SKIP_R:-0}" != "1" ]; then
  stage "1/4 R pipeline + Figures 1-4 (run_all.R)" "$LOG_DIR/r_pipeline_${STAMP}.log" \
    Rscript "$HERE/run_all.R"
else
  echo ""; echo "SKIP_R=1 -> R pipeline skipped."
fi

# ---------------------------------------------------------------------------
# 2. C11 scale-free report  ->  reports/C11_scale_free/
# ---------------------------------------------------------------------------
# WHY THIS IS A PIPELINE STAGE AND NOT A ONE-OFF.
# The manuscript's scale-free results - E_r against E_K, the model-free CUE
# optimum and its margin below 37 C, the gap decomposition against the fitted
# optimum - are produced HERE, not by anything under scripts/. Those numbers are
# typed into v3.qmd by hand, so if this report does not regenerate alongside the
# pipeline they silently drift away from the tables they came from. That is the
# whole reason for the stage.
#
# Safe to run unattended: reports/C11_scale_free/R/00_common.R sources config.R
# for the carbon constants and the figure whitelist, reads the tables stage 1
# has just rewritten, and writes ONLY under reports/C11_scale_free/. It cannot
# disturb results/ or data/.
#
# Cheap: least squares and a bootstrap, no MCMC and no network. Seconds, not
# minutes. The four scripts are order-dependent (04 plots what 01-03 compute).
#
#   SKIP_C11=1   skip this stage
if [ "${SKIP_C11:-0}" != "1" ]; then
  C11_DIR="$ROOT/reports/C11_scale_free"
  if [ -d "$C11_DIR/R" ]; then
    for f in 01_balance_quantity 02_scale_free 03_gap_decomposition 04_figures 05_clade_contrasts; do
      stage "2/4 C11 $f" "$LOG_DIR/c11_${f}_${STAMP}.log" \
        Rscript "$C11_DIR/R/$f.R"
    done

    if [ "${SKIP_RENDER:-0}" != "1" ] && command -v quarto >/dev/null 2>&1; then
      stage "2/4 C11 report render" "$LOG_DIR/c11_render_${STAMP}.log" \
        bash -c "cd '$C11_DIR' && quarto render C11_scale_free.qmd"
    fi
  else
    echo ""; echo "!! reports/C11_scale_free/R not found - C11 stage skipped."
    echo "   The manuscript's scale-free numbers come from there; check the checkout."
  fi
else
  echo ""; echo "C11 scale-free report skipped (SKIP_C11=1)."
fi

# ---------------------------------------------------------------------------
# 3. Manuscript
# ---------------------------------------------------------------------------
if [ "${SKIP_RENDER:-0}" != "1" ]; then
  if command -v quarto >/dev/null 2>&1; then
    stage "4/4 quarto render (manuscript/v3.qmd)" "$LOG_DIR/quarto_${STAMP}.log" \
      bash -c "cd '$ROOT/manuscript' && quarto render v3.qmd"
  else
    echo ""; echo "!! quarto not on PATH - manuscript NOT rendered."
    echo "   Install quarto, then: quarto install tinytex"
    exit 1
  fi
else
  echo ""; echo "SKIP_RENDER=1 -> manuscript render skipped."
fi

summary
banner "ALL DONE"
