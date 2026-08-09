#!/usr/bin/env bash
# =============================================================================
# run_all.sh - the whole project, unattended, in one command
# =============================================================================
#   bash scripts/run_all.sh
#
# Everything is tee'd to logs/ and every stage is timed.
#
# ----------------------------------------------------------------------------
# WHERE EACH OF THE EIGHTEEN SCRIPTS RUNS, AND WHY
# ----------------------------------------------------------------------------
#   01 convert_xlsx        Shiny app.  NOT RUN. Its outputs (data/*_Oxygen.csv,
#                          results/tables/otu_names.csv) are committed INPUTS.
#   02 longdata            \
#   03 trimming             |
#   06 inoculation          |  run by `Rscript scripts/run_all.R`, in this order.
#   07 oxygen_fits          |  06 is a Shiny app: it is SOURCED headless, which
#   08 temperature_equil.   |  defines the app and returns without launching it
#   09 bayesian_models      |  (its results/tables/otu_inoc.csv is a committed
#   11 bayesian_plots       |  INPUT read by 07, not regenerated here).
#   12 carbon_tax           |
#   14 main_figures         |  08's per-curve table is an input to 15.
#   15 uncertainty_bands   /
#   04 trim_selector       Shiny app.  NOT RUN. results/tables/manual_fit_windows.csv
#                          and plot_exclude_points.csv are committed INPUTS that
#                          config.R reads when USE_APP_TRIM_FILES is TRUE.
#   05 cell_sizes          Shiny app.  NOT RUN. results/tables/otu_cell_sizes.csv
#                          is a committed INPUT.
#   13 capacity_expression run here, as its own stage. It is not in run_all.R
#                          because it used to download from GEO; the count matrix
#                          is now committed at data/expression/GSE165762/, so it
#                          runs offline and 16 needs its output.
#   16 supplementary_figs  run here, as its own stage (needs 13 + etc-GEM data).
#   17 schematic.py        run here, as its own stage (python3 + matplotlib).
#
#   10 n0_term_test        DELIBERATELY NOT RUN - diagnostic, not a pipeline stage.
#   18 n0_treatment_panel  DELIBERATELY NOT RUN - diagnostic, not a pipeline stage.
#       Both answer "how much does the between-clade ordering lean on the N0
#       back-projection?" by REFITTING the Bayesian respiration model two and
#       three times respectively, writing each alternative fit over
#       results/tables/bayes_resp_arr_summary.csv and results/rds/. They restore
#       derived_N0_R_results_with_carbon.csv but NOT those other outputs, so a
#       run that ended with either of them would leave the published bayes_*
#       tables holding an ALTERNATIVE fit. 10's own header says so: "re-run
#       09_bayesian_models.R normally afterwards to refresh the published ones."
#       They also cost tens of minutes each. Run them deliberately, then re-run
#       09 -> 11 -> 12 -> 14 -> 15 to restore the published set:
#           Rscript scripts/10_n0_term_test.R
#           Rscript scripts/18_n0_treatment_panel.R     # Supplementary Fig. 6
#
#   etc-GEM (python)       OPTIONAL, OFF BY DEFAULT - see RUN_ETCGEM below.
#
# ENVIRONMENT (all optional):
#   RUN_ETCGEM=1     also regenerate the etc-GEM outputs from the Python engine
#                    in the cauris_etcgem submodule. Off by default because the
#                    outputs it produces are VENDORED at outputs/supp_data/ (see
#                    commit "Vendor etc-GEM supp_data outputs so clean checkout
#                    reproduces model figures"), and because main does not pin
#                    the Python engine - the requirements.lock that pins it lives
#                    on a submodule commit this repository does not reference.
#                    Turning it on without that pin can silently change the model
#                    under the committed figures.
#   ETCGEM_STAGES    which etc-GEM stages to run when RUN_ETCGEM=1 (default: all)
#   SKIP_R=1         skip the R pipeline
#   SKIP_RENDER=1    skip the manuscript render
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
echo "  RUN_ETCGEM  : ${RUN_ETCGEM:-0}"
echo "  logs        : $LOG_DIR"

# ---------------------------------------------------------------------------
# 0. etc-GEM (python) - optional, off by default
# ---------------------------------------------------------------------------
if [ "${RUN_ETCGEM:-0}" = "1" ]; then
  PYBIN="$ROOT/cauris_etcgem/.venv/bin/python"
  if [ ! -x "$PYBIN" ]; then
    echo ""
    echo "!! RUN_ETCGEM=1 but there is no virtualenv at $PYBIN"
    echo "   cd cauris_etcgem && python3 -m venv .venv \\"
    echo "     && .venv/bin/python -m pip install -r requirements.txt"
    exit 1
  fi
  stage "0/5 etc-GEM (optional) (python, ${ETCGEM_STAGES:-all})" "$LOG_DIR/etcgem_${STAMP}.log" \
    env PYTHONUNBUFFERED=1 "$PYBIN" \
        "$ROOT/cauris_etcgem/strains/eci_cauris/scripts/generate_model_data.py" \
        ${ETCGEM_STAGES:-all}
else
  echo ""
  echo "etc-GEM stage SKIPPED (RUN_ETCGEM is not 1)."
  echo "  The outputs it would produce are vendored at outputs/supp_data/ and are"
  echo "  what 14 and 16 read. main does not pin the Python engine, so regenerating"
  echo "  them is a deliberate act, not part of an unattended run."
fi

# ---------------------------------------------------------------------------
# 1. R pipeline: 02, 03, 06, 07, 08, 09, 11, 12, 14, 15
# ---------------------------------------------------------------------------
if [ "${SKIP_R:-0}" != "1" ]; then
  stage "1/5 R pipeline 02-15 (run_all.R)" "$LOG_DIR/r_pipeline_${STAMP}.log" \
    Rscript "$HERE/run_all.R"

  stage "2/5 13_capacity_expression.R" "$LOG_DIR/r_13_${STAMP}.log" \
    Rscript "$HERE/13_capacity_expression.R"

  stage "3/5 16_supplementary_figures.R" "$LOG_DIR/r_16_${STAMP}.log" \
    Rscript "$HERE/16_supplementary_figures.R"

  if command -v python3 >/dev/null 2>&1; then
    stage "4/5 17_schematic.py" "$LOG_DIR/py_17_${STAMP}.log" \
      python3 "$HERE/17_schematic.py"
  else
    echo ""; echo "!! python3 not on PATH - 17_schematic.py NOT run."
    exit 1
  fi
else
  echo ""; echo "SKIP_R=1 -> R pipeline skipped."
fi

# ---------------------------------------------------------------------------
# 2. Manuscript
# ---------------------------------------------------------------------------
if [ "${SKIP_RENDER:-0}" != "1" ]; then
  if command -v quarto >/dev/null 2>&1; then
    stage "5/5 quarto render (manuscript/draft)" "$LOG_DIR/quarto_${STAMP}.log" \
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
