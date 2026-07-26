#!/usr/bin/env bash
# =============================================================================
# run_all.sh - the whole project, unattended, in one command
# =============================================================================
#   bash scripts/run_all.sh
#
# Order:
#   1. Python etc-GEM pipeline   (cauris_etcgem, incl. the slow emcee stage)
#   2. Rscript scripts/run_all.R (02 -> 03 -> 06 -> 07 -> ... -> 13, then 14)
#   3. quarto render             (manuscript, supplementary, combined; docx+pdf)
#
# Everything is tee'd to logs/ and every stage is timed.
#
# ENVIRONMENT (all optional; all are PATH knobs - no analysis decision):
#   CANDIDAS_RESULTS         results tree            (default: <root>/results)
#   CANDIDAS_SUPP_DATA       etc-GEM outputs 12/13 read
#   CANDIDAS_EXPRESSION_OUT  where 11 writes
#   ETCGEM_OUT_SUFFIX        suffix for the etc-GEM output dirs (e.g. _C1)
#   ETCGEM_STAGES            stages to run          (default: all)
#   SKIP_PY / SKIP_R / SKIP_RENDER   set to 1 to skip that stage
#   LOG_DIR                  log destination        (default: <root>/logs)
#
# Non-destructive C1-style run (writes nothing into results/ or supp_data/):
#   CANDIDAS_RESULTS="$PWD/runs/C1_reproduction" \
#   ETCGEM_OUT_SUFFIX=_C1 \
#   CANDIDAS_SUPP_DATA="$PWD/cauris_etcgem/runs/C1/supp_data" \
#   CANDIDAS_EXPRESSION_OUT="$PWD/runs/C1_reproduction/expression" \
#     bash scripts/run_all.sh
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$ROOT"

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
    printf '  %-34s %-14s exit %s\n' \
      "${STAGE_NAMES[$i]}" "$(hms "${STAGE_SECS[$i]}")" "${STAGE_RC[$i]}"
  done
  printf '  %-34s %-14s\n' "TOTAL" "$(hms $(( $(date +%s) - RUN_T0 )))"
  echo "  full log: $MAIN_LOG"
}

# --- everything below is captured to MAIN_LOG as well ------------------------
exec > >(tee -a "$MAIN_LOG") 2>&1

banner "Candidas full pipeline   $(date '+%Y-%m-%d %H:%M:%S')"
echo "  root                    : $ROOT"
echo "  CANDIDAS_RESULTS        : ${CANDIDAS_RESULTS:-<default: $ROOT/results>}"
echo "  CANDIDAS_SUPP_DATA      : ${CANDIDAS_SUPP_DATA:-<default: submodule outputs/supp_data>}"
echo "  CANDIDAS_EXPRESSION_OUT : ${CANDIDAS_EXPRESSION_OUT:-<default: data/expression>}"
echo "  ETCGEM_OUT_SUFFIX       : ${ETCGEM_OUT_SUFFIX:-<none>}"
echo "  logs                    : $LOG_DIR"

# ---------------------------------------------------------------------------
# 1. Python etc-GEM pipeline
# ---------------------------------------------------------------------------
PYBIN="$ROOT/cauris_etcgem/.venv/bin/python"
if [ ! -x "$PYBIN" ]; then
  echo ""
  echo "!! No virtualenv at $PYBIN"
  echo "   Create it first (see SETUP.md):"
  echo "     cd cauris_etcgem && python3 -m venv .venv \\"
  echo "       && .venv/bin/python -m pip install -r requirements.lock"
  exit 1
fi

if [ "${SKIP_PY:-0}" != "1" ]; then
  ETCGEM_STAGES="${ETCGEM_STAGES:-all}"
  stage "1/3 etc-GEM (python, $ETCGEM_STAGES)" "$LOG_DIR/etcgem_${STAMP}.log" \
    env PYTHONUNBUFFERED=1 "$PYBIN" \
        "$ROOT/cauris_etcgem/strains/eci_cauris/scripts/generate_model_data.py" \
        $ETCGEM_STAGES
else
  echo ""; echo "SKIP_PY=1 -> etc-GEM stage skipped."
fi

# ---------------------------------------------------------------------------
# 2. R pipeline
# ---------------------------------------------------------------------------
if [ "${SKIP_R:-0}" != "1" ]; then
  stage "2/3 R pipeline (02-14)" "$LOG_DIR/r_pipeline_${STAMP}.log" \
    Rscript "$HERE/run_all.R"
else
  echo ""; echo "SKIP_R=1 -> R pipeline skipped."
fi

# ---------------------------------------------------------------------------
# 3. Manuscript
# ---------------------------------------------------------------------------
if [ "${SKIP_RENDER:-0}" != "1" ]; then
  if command -v quarto >/dev/null 2>&1; then
    stage "3/3 quarto render" "$LOG_DIR/quarto_${STAMP}.log" \
      bash -c "cd '$ROOT/manuscript/draft' && quarto render"
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
