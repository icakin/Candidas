# =============================================================================
# run_all.R - Master pipeline runner (Candida auris, temperature x isolate)
# =============================================================================
# Scripts are numbered in running order; each sources config.R for shared paths.
#
#   01 convert_xlsx                 (Run App)  PreSens .xlsx -> tidy CSV        [manual]
#   02 longdata                      wide -> long oxygen
#   03 trimming                      spline trimming + diagnostics
#   04 trim_selector                (Run App)  review/override fit windows      [optional]
#   05 cell_sizes                   (Run App)  mean cell volume per isolate     [optional]
#   06 inoculation                   inoculation density per group
#   07 oxygen_fits                   growth / respiration / CUE + SS/Arrhenius
#   08 temperature_equilibration_sensitivity  onset-temperature uncertainty (T_internal)
#   09 bayesian_models               hierarchical Bayesian TPCs (brms)   [slow]
#   10 n0_term_test                  r*dt back-projection term diagnostic
#   11 bayesian_plots                posterior / contrast / CUE (+ temperature band from 08)
#   12 carbon_tax                    fever-tax economics -> carbon_tax_isolate.csv
#   13 capacity_expression           GEO GSE165762 transcription test    [downloads]
#   14 main_figures                  Fig 1 (decoupling) + Fig 2 (the bill) + Fig 3 (etc-GEM)
#   15 uncertainty_bands             respiration/CUE with combined uncertainty band
#   16 supplementary_figures         etc-GEM diagnostics + validation     [needs 13 + etc-GEM*]
#   17 schematic.py       (python3)  Fig 4 pipeline schematic
#
# * Fig 3 (inside 14) and the supplementary figures (16) read the etc-GEM outputs
#   produced by the Python pipeline in cauris_etcgem/strains/eci_cauris/scripts/.
#   Run that first. Without it, 14 still writes Fig 1 & 2 and stops before Fig 3.
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
            "08_temperature_equilibration_sensitivity.R","09_bayesian_models.R",
            "11_bayesian_plots.R","12_carbon_tax.R","14_main_figures.R",
            "15_uncertainty_bands.R")) {
  run_script(s)
}
message("\n", strrep("=",70))
message("Done: data, models, diagnostics, main figures, and uncertainty bands.")
message("Run 13_capacity_expression.R (downloads GEO), 16_supplementary_figures.R")
message("(needs 13 + etc-GEM outputs), and 17_schematic.py (python3) separately.")
message(strrep("=",70))
