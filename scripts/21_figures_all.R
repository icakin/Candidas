# =============================================================================
# 21_figures_all.R - every manuscript figure, in one command
# =============================================================================
#   Rscript scripts/21_figures_all.R        or  source() it in RStudio
#
# A thin front end. It contains no figure code: it runs each figure's own script in
# order, in its own environment, so there is exactly one script per figure and a
# failure in one does not silently take out the others. Every generator writes into
# results/figures/manuscript/.
#
#   Figure 1  growth/respiration decoupling          ->  16_fig1.R
#   Figure 2  the carbon cost of fever (+ 4 supps)   ->  17_fig2.R
#   Figure 3  thermal/phylogeny discordance          ->  18_fig3.R
#   Figure 4  etcGEM counterfactual                  ->  19_fig4.R
#
# Figures 1 and 2 share their data and posterior draws through fig_common.R.
# Figure 4 reads tables written by the Python etcGEM layer (gem/tables/, see
# gem/README.md); regenerating THOSE is a deliberate act (python3 gem/26_fig4_tables.py
# after the simulation scripts), not part of drawing the figure.
#
# Supplementary Figures 1-4 are written by the analysis scripts that compute them
# (08, 11, 15), not here.
#
# ENVIRONMENT:
#   FIGS=1,3        run only these figures (default: all)
# =============================================================================

script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(gsub("~+~", " ", f[1], fixed = TRUE))) else
    if (file.exists("scripts/16_fig1.R")) normalizePath("scripts") else normalizePath(".")
})
root <- normalizePath(file.path(script_dir, ".."))

want <- Sys.getenv("FIGS", "")
want <- if (nzchar(want)) trimws(strsplit(want, ",")[[1]]) else NULL

STAGES <- list(
  list(fig = "1", what = "Figure 1 (decoupling)",     file = "16_fig1.R"),
  list(fig = "2", what = "Figure 2 (fever cost)",     file = "17_fig2.R"),
  list(fig = "3", what = "Figure 3 (discordance)",    file = "18_fig3.R"),
  list(fig = "4", what = "Figure 4 (etcGEM)",         file = "19_fig4.R")
)

hr <- function() message(strrep("=", 70))
results <- data.frame(stage = character(), seconds = numeric(),
                      status = character(), stringsAsFactors = FALSE)

for (st in STAGES) {
  if (!is.null(want) && !(st$fig %in% want)) next
  path <- file.path(script_dir, st$file)
  if (!file.exists(path)) {
    message("!! missing: ", st$file, " - ", st$what, " SKIPPED")
    results[nrow(results) + 1, ] <- list(st$what, 0, "missing"); next
  }
  hr(); message("Running: ", st$what, "   [", st$file, "]"); hr()
  t0 <- Sys.time()
  # Own environment: config.R re-assigns FIG_KEEP on every source(), so two generators
  # sharing globals is exactly how figure writes got lost before.
  ok <- !inherits(try(sys.source(path, envir = new.env(parent = globalenv())),
                      silent = FALSE), "try-error")
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  results[nrow(results) + 1, ] <- list(st$what, dt, if (ok) "ok" else "FAILED")
  message("--- ", st$what, ": ", if (ok) "ok" else "FAILED", " in ", dt, "s")
}

hr(); message("SUMMARY")
for (i in seq_len(nrow(results)))
  message(sprintf("  %-28s %8.1fs  %s", results$stage[i], results$seconds[i], results$status[i]))
message("  figures -> ", file.path(root, "results/figures/manuscript"))
hr()
if (any(results$status == "FAILED")) quit(status = 1, save = "no")
