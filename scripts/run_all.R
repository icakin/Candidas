# =============================================================================
# run_all.R - Master pipeline runner (Candida auris, temperature x isolate)
# =============================================================================
# Scripts are numbered in running order; each sources config.R for shared paths.
#
#   01 convert_xlsx        (Run App)  PreSens .xlsx -> tidy CSV          [manual]
#   02 longdata                        wide -> long oxygen
#   03 trimming                        spline trimming + diagnostics
#   04 trim_selector       (Run App)   review/override fit windows       [optional]
#   05 cell_sizes          (Run App)   mean cell volume per isolate      [optional]
#   06 inoculation                     inoculation density per group
#   07 oxygen_fits                     growth / respiration / CUE + SS/Arrhenius
#   08 bayesian_models                 hierarchical Bayesian TPCs (brms)   [slow]
#   09 bayesian_plots                  posterior / contrast / CUE diagnostics
#   10 carbon_tax                      fever-tax economics -> carbon_tax_isolate.csv
#   11 capacity_expression             GEO GSE165762 transcription test  [downloads]
#   12 main_figures                    Fig 1 (decoupling) + Fig 2 (the bill) + Fig 3 (etc-GEM)
#   13 supplementary_figures           etc-GEM diagnostics + validation   [needs 11 + etc-GEM*]
#   14 schematic.py         (python3)  Fig 4 pipeline schematic
#
# * Fig 3 (inside 12) and the supplementary figures (13) read the etc-GEM outputs
#   produced by the Python pipeline in cauris_etcgem/strains/eci_cauris/scripts/.
#   Run that first. Without it, 12 still writes Fig 1 & 2 and stops before Fig 3.
# =============================================================================

script_dir <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
                  rstudioapi::isAvailable() &&
                  nzchar(rstudioapi::getActiveDocumentContext()$path)) {
  dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
}
run_script <- function(name) {
  message("\n", strrep("=", 70), "\nRunning: ", name, "\n", strrep("=", 70), "\n")
  source(file.path(script_dir, name), local = FALSE)
}
for (s in c("02_longdata.R","03_trimming.R","06_inoculation.R","07_oxygen_fits.R",
            "08_bayesian_models.R","09_bayesian_plots.R","10_carbon_tax.R","12_main_figures.R")) {
  run_script(s)
}
message("\n", strrep("=",70))
message("Done: data, models, diagnostics, and main figures (Fig 1-3).")
message("Run 11_capacity_expression.R (downloads GEO), 13_supplementary_figures.R")
message("(needs 11 + etc-GEM outputs), and 14_schematic.py (python3) separately.")
message(strrep("=",70))
