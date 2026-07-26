# =============================================================================
# run_all.R - Master pipeline runner (Candida auris, temperature x isolate)
# =============================================================================
# UNATTENDED. Runs end to end under `Rscript scripts/run_all.R` with no user
# interaction and no browser. Scripts are numbered in running order; each
# sources config.R for shared paths.
#
#   01 convert_xlsx        (app)   PreSens .xlsx -> tidy CSV        [NOT RUN - see below]
#   02 longdata                    wide -> long oxygen
#   03 trimming                    spline trimming + diagnostics
#   04 trim_selector       (app)   review/override fit windows      [NOT RUN - see below]
#   05 cell_sizes          (app)   mean cell volume per isolate     [NOT RUN - see below]
#   06 inoculation         (app)   inoculation density per group    [sourced, app not launched]
#   07 oxygen_fits                 growth / respiration / CUE + SS/Arrhenius
#   08 bayesian_models             hierarchical Bayesian TPCs (brms)  [SLOW]
#   09 bayesian_plots              posterior / contrast / CUE diagnostics
#   10 carbon_tax                  fever-tax economics -> carbon_tax_isolate.csv
#   11 capacity_expression         GEO GSE165762 transcription test (data COMMITTED)
#   12 main_figures                Fig 1 (decoupling) + Fig 2 (the bill) + Fig 3 (etc-GEM)
#   13 supplementary_figures       etc-GEM diagnostics + validation  [needs 11 + etc-GEM]
#   14 schematic.py        (py)    pipeline schematic figure
#
# WHY 01 / 04 / 05 ARE NOT IN THE LOOP
#   They are click-driven Shiny apps whose outputs are COMMITTED INPUTS to this
#   pipeline and are treated as data:
#       01 -> data/*_Oxygen.csv, results/tables/otu_names.csv
#       04 -> results/tables/manual_fit_windows.csv, plot_exclude_points.csv
#       05 -> results/tables/otu_cell_sizes.csv
#   Re-running them would require a human to re-make hand judgements and would
#   overwrite committed decisions. All four apps are now guarded so that sourcing
#   them defines the app and returns; launch one deliberately with
#       Rscript scripts/04_trim_selector.R --app
#   06 IS in the loop - it is sourced (its output otu_inoc.csv is likewise a
#   committed input), it simply does not launch its app any more. Before this
#   change, run_all.R source()d 06 and blocked forever on its Shiny server, which
#   is why the master runner could never complete unattended.
#
# ETC-GEM DEPENDENCY
#   12 (Fig 3) and 13 read the Python etc-GEM outputs in
#   cauris_etcgem/strains/eci_cauris/outputs/supp_data/. Run that pipeline FIRST
#   (scripts/run_all.sh does). Without it, 12 and 13 now fail LOUDLY with a named
#   ETCGEM_SUPP_DATA_* error instead of silently mis-resolving the path.
#
# ENVIRONMENT KNOBS (paths only - no analysis decision depends on any of them)
#   CANDIDAS_RESULTS         redirect the whole results tree (default: results/)
#   CANDIDAS_SUPP_DATA       which etc-GEM outputs dir 12/13 read
#   CANDIDAS_EXPRESSION_OUT  where 11 writes its expression outputs
#   CANDIDAS_SKIP            comma-separated script numbers to skip, e.g. "08,09"
# =============================================================================

# Script directory: --file= FIRST, because under `Rscript scripts/run_all.R`
# neither sys.frame(1)$ofile nor rstudioapi resolves, and the old code fell back
# to getwd() - which made every source() below look in the project ROOT and fail.
script_dir <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) return(dirname(normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
})

# THE flag that stops every Shiny script from blocking this runner.
options(candidas.headless = TRUE)
Sys.setenv(CANDIDAS_HEADLESS = "1")

.SKIP <- trimws(unlist(strsplit(Sys.getenv("CANDIDAS_SKIP", ""), ",")))
.SKIP <- .SKIP[nzchar(.SKIP)]

fmt_dur <- function(secs) {
  secs <- as.numeric(secs)
  h <- floor(secs / 3600); m <- floor((secs %% 3600) / 60); s <- round(secs %% 60)
  if (h > 0) sprintf("%dh %02dm %02ds", h, m, s) else
  if (m > 0) sprintf("%dm %02ds", m, s) else sprintf("%.1fs", secs)
}

.T0      <- Sys.time()
.TIMINGS <- list()
.FAILED  <- character(0)

run_script <- function(name) {
  num <- sub("_.*$", "", name)
  if (num %in% .SKIP) {
    message("\n", strrep("-", 70), "\nSKIPPED (CANDIDAS_SKIP): ", name, "\n", strrep("-", 70))
    .TIMINGS[[name]] <<- NA_real_
    return(invisible(NULL))
  }
  message("\n", strrep("=", 70), "\nRunning: ", name,
          "   |  started ", format(Sys.time(), "%H:%M:%S"), "\n", strrep("=", 70), "\n")
  t0 <- Sys.time()
  source(file.path(script_dir, name), local = FALSE)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  .TIMINGS[[name]] <<- el
  message("\n--- DONE ", name, " in ", fmt_dur(el),
          "  (total elapsed ", fmt_dur(difftime(Sys.time(), .T0, units = "secs")), ")")
  invisible(NULL)
}

# The full unattended order. 11 and 13 used to be mentioned only in a closing
# message and were never actually run by this file; they are in the loop now.
SCRIPTS <- c("02_longdata.R",
             "03_trimming.R",
             "06_inoculation.R",
             "07_oxygen_fits.R",
             "08_bayesian_models.R",
             "09_bayesian_plots.R",
             "10_carbon_tax.R",
             "11_capacity_expression.R",
             "12_main_figures.R",
             "13_supplementary_figures.R")

for (s in SCRIPTS) run_script(s)

# ---- 14 is Python; run it here so `Rscript run_all.R` really is complete ------
if (!("14" %in% .SKIP)) {
  message("\n", strrep("=", 70), "\nRunning: 14_schematic.py (python3)\n", strrep("=", 70), "\n")
  t0 <- Sys.time()
  py  <- Sys.which("python3"); if (!nzchar(py)) py <- Sys.which("python")
  if (!nzchar(py)) {
    message("!! python3 not found on PATH - 14_schematic.py NOT run.")
    .FAILED <- c(.FAILED, "14_schematic.py (no python3)")
  } else {
    rc <- system2(py, shQuote(file.path(script_dir, "14_schematic.py")))
    if (!identical(as.integer(rc), 0L))
      .FAILED <- c(.FAILED, sprintf("14_schematic.py (exit %s)", rc))
  }
  .TIMINGS[["14_schematic.py"]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

# ---- summary -----------------------------------------------------------------
message("\n", strrep("=", 70))
message("PIPELINE TIMINGS")
message(strrep("=", 70))
for (nm in names(.TIMINGS)) {
  v <- .TIMINGS[[nm]]
  message(sprintf("  %-28s %s", nm, if (is.na(v)) "skipped" else fmt_dur(v)))
}
message(sprintf("  %-28s %s", "TOTAL",
                fmt_dur(difftime(Sys.time(), .T0, units = "secs"))))
message(strrep("=", 70))
message("NOT RUN (click-driven apps; their outputs are committed INPUTS):")
message("  01_convert_xlsx.R, 04_trim_selector.R, 05_cell_sizes.R")
message("  -> launch one deliberately with:  Rscript scripts/0N_name.R --app")
if (length(.FAILED)) {
  message(strrep("=", 70))
  message("FAILURES: ", paste(.FAILED, collapse = "; "))
  message(strrep("=", 70))
  quit(status = 1L)
}
message(strrep("=", 70))
