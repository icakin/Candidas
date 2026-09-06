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

# Packaging: resolve the script directory from --file= FIRST. Under
# `Rscript scripts/<this>.R` neither rstudioapi nor sys.frame(1)$ofile resolves,
# so this fell back to getwd() and then died on "cannot open file .../config.R".
# Sourcing it from run_all.R was unaffected, which is why the bug stayed hidden.
script_dir <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) return(dirname(normalizePath(
    # R replaces every space in --file= with "~+~", so a project path that
    # contains a space comes back mangled and source() then fails on a path
    # that does not exist. Un-mangle it before normalising.
    gsub("~+~", " ", sub("^--file=", "", fa[1]), fixed = TRUE),
    mustWork = FALSE)))
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
})
# THE flag that stops 06 (a Shiny app) from blocking this runner. Sourcing 06
# under it defines the app and returns; results/tables/otu_inoc.csv is a
# committed INPUT and is not regenerated here.
options(candidas.headless = TRUE)
Sys.setenv(CANDIDAS_HEADLESS = "1")

run_script <- function(name) {
  message("\n", strrep("=", 70), "\nRunning: ", name, "\n", strrep("=", 70), "\n")
  source(file.path(script_dir, name), local = FALSE)
}
for (s in c("02_longdata.R","03_trimming.R","06_inoculation.R","07_oxygen_fits.R",
            "08_temperature_equilibration_sensitivity.R","09_bayesian_models.R",
            "11_bayesian_plots.R","12_carbon_tax.R","16_fig1_fig2.R",
            "13_uncertainty_bands.R","17_fig3.R")) {
  run_script(s)
}
message("\n", strrep("=",70))
message("Done: data, models, diagnostics, main figures (incl. Fig 3), bands.")
message("Run 13_capacity_expression.R (downloads GEO), 16_supplementary_figures.R")
message("(needs 13 + etc-GEM outputs), and 17_schematic.py (python3) separately.")
message("Fig 4 is python: run_all.sh runs scripts/18_fig4.R as its own stage.")
message(strrep("=",70))
