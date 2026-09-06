# =============================================================================
# 20_figures_all.R - every main manuscript figure, in one command
# =============================================================================
#   Rscript scripts/20_figures_all.R        or  source() it in RStudio
#
# A thin front end. It does not contain figure code: it runs each figure's own
# generator in order, so there is still exactly one script per figure and a
# failure in one does not silently take out the others. Every generator writes
# into results/figures/manuscript/.
#
#   Figure 1  growth/respiration decoupling  ->  16_fig1_fig2.R   (R)
#   Figure 2  the carbon cost of fever       ->  16_fig1_fig2.R   (R)
#   Figure 3  thermal/phylogeny discordance  ->  17_fig3.R           (R)
#   Figure 4  etcGEM counterfactual          ->  18_fig4.R   (python3)
#   (16_fig1_fig2.R also writes the capacity figure and its consistency supp.)
#
# WHY FIGURE 4 IS PYTHON. Its panels are read off etcGEM simulations. The
# simulation stages (gem/etcgem_counterfactual.py, gem/dyn_sparse.py) need cobra;
# the drawing stage does not, and reads only committed JSON/CSV from gem/, so
# this script can run it unattended. Regenerating its inputs is a deliberate act:
#
#   python3 gem/etcgem_counterfactual.py     # counterfactual_{results.json,sweep.csv}
#   python3 gem/dyn_sparse.py Tm             # dyn_sparse_Tm.json
#   python3 gem/dyn_sparse.py Topt           # dyn_sparse_Topt.json
#   python3 gem/counts_at_44.py              # counts_at_44.json (same rule as Fig 3)
#
# Each stage runs in its own environment so that one generator cannot clobber
# another's globals: config.R re-assigns FIG_KEEP on every source(), which is the
# bug documented at length in config.R's whitelist section.
#
# RELATION TO run_all. This is the manual entry point, for when you want every
# main figure rebuilt without remembering which script makes what. The full
# pipeline does NOT call this: run_all.R runs 16_fig1_fig2.R and 17_fig3.R as
# ordinary R stages and run_all.sh runs 18_fig4.R as its own stage, so
# each gets its own timed, separately-logged stage. Same generators either way,
# so the two paths cannot disagree about what a figure contains.
#
# ENVIRONMENT:
#   FIGS=1,3        run only these figures (default: all)
#   SKIP_PY=1       skip Figure 4 even if python3 is present
# =============================================================================

script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else
    if (file.exists("scripts/16_fig1_fig2.R")) normalizePath("scripts") else
      normalizePath(".")
})
root <- normalizePath(file.path(script_dir, ".."))

want <- Sys.getenv("FIGS", "")
want <- if (nzchar(want)) trimws(strsplit(want, ",")[[1]]) else NULL
skip_py <- identical(Sys.getenv("SKIP_PY", "0"), "1")

STAGES <- list(
  list(figs = c("1", "2"), what = "Figures 1-2 (+ capacity)", kind = "R",
       file = "16_fig1_fig2.R"),
  list(figs = "3",         what = "Figure 3 (discordance)",   kind = "R",
       file = "17_fig3.R"),
  list(figs = "4",         what = "Figure 4 (etcGEM)",        kind = "py",
       file = "18_fig4.R")
)

hr <- function() message(strrep("=", 70))
results <- data.frame(stage = character(), seconds = numeric(),
                      status = character(), stringsAsFactors = FALSE)

for (st in STAGES) {
  if (!is.null(want) && !any(st$figs %in% want)) next
  path <- file.path(script_dir, st$file)
  if (!file.exists(path)) {
    message("!! missing: ", st$file, " - ", st$what, " SKIPPED")
    results[nrow(results) + 1, ] <- list(st$what, 0, "missing")
    next
  }
  if (st$kind == "py" && skip_py) {
    message("Figure 4 skipped (SKIP_PY=1)."); next
  }
  if (st$kind == "py" && !nzchar(Sys.which("python3"))) {
    message("!! python3 not on PATH - ", st$what, " NOT run (Figure 4 stale).")
    results[nrow(results) + 1, ] <- list(st$what, 0, "no python3")
    next
  }

  hr(); message("Running: ", st$what, "   [", st$file, "]"); hr()
  t0 <- Sys.time()
  ok <- if (st$kind == "R") {
    # Own environment: config.R re-assigns FIG_KEEP on every source(), so two
    # generators sharing globals is exactly how figure writes got lost before.
    !inherits(try(sys.source(path, envir = new.env(parent = globalenv())),
                  silent = FALSE), "try-error")
  } else {
    system2("python3", shQuote(path)) == 0
  }
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  results[nrow(results) + 1, ] <- list(st$what, dt, if (ok) "ok" else "FAILED")
  message("--- ", st$what, ": ", if (ok) "ok" else "FAILED", " in ", dt, "s")
}

hr(); message("SUMMARY")
for (i in seq_len(nrow(results)))
  message(sprintf("  %-28s %8.1fs  %s",
                  results$stage[i], results$seconds[i], results$status[i]))
message("  figures -> ", file.path(root, "results/figures/manuscript"))
hr()
if (any(results$status == "FAILED")) quit(status = 1, save = "no")
