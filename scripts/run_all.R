# =============================================================================
# run_all.R - the R pipeline, start to finish (Candida, temperature x isolate)
# =============================================================================
# Scripts are numbered in running order; each sources config.R for shared paths.
#
#   00 install                        one-off: renv::restore() and checks         [manual]
#   01 convert_xlsx        (Run App)  PreSens .xlsx -> tidy CSV                   [manual]
#   02 longdata                       wide -> long oxygen
#   03 trimming                       spline trimming + diagnostics
#   04 trim_selector       (Run App)  review/override fit windows                 [optional]
#   05 cell_sizes          (Run App)  mean cell volume per isolate                [optional]
#   06 inoculation         (Run App)  inoculation density per group  (sourced headless)
#   07 oxygen_fits                    growth / respiration / CUE + SS/Arrhenius
#   08 temperature_equilibration_sensitivity   onset-temperature uncertainty  -> Supp Fig 1
#   09 bayesian_models                hierarchical Bayesian TPCs (brms)          [slow]
#   10 n0_term_test                   r*dt back-projection diagnostic            [not run here]
#   11 bayesian_plots                 posterior / contrast / CUE plots           -> Supp Figs 2, 3
#   12 carbon_tax                     fever-tax economics
#   13 uncertainty_bands              respiration/CUE with combined uncertainty band
#   14 rk_covariance_check            cov(r, K) diagnostic                       [not run here]
#   15 n0_treatment_panel             N0 back-projection, three treatments       -> Supp Fig 4 [not run here]
#   16 fig1                           Figure 1   (fig_common.R holds the shared draws)
#   17 fig2                           Figure 2 + four supplementary panels
#   18 fig3                           Figure 3
#   19 fig4                           Figure 4   (reads gem/tables/, written by the Python etcGEM layer)
#   20 (reserved)
#   21 figures_all                    16-19 in one command
#
# 10, 14 and 15 refit the Bayesian respiration model under alternative treatments and
# overwrite results/tables/bayes_resp_arr_summary.csv while they run; run them
# deliberately and then re-run 09 -> 11 -> 12 -> 16 -> 17 to restore the published set.
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
            "11_bayesian_plots.R","12_carbon_tax.R","13_uncertainty_bands.R",
            "16_fig1.R","17_fig2.R","18_fig3.R","19_fig4.R")) {
  run_script(s)
}
message("\n", strrep("=",70))
message("Done: data, models, diagnostics, bands, Figures 1-4 and Supplementary Figures 1-3.")
message("Supplementary Figure 4 is 15_n0_treatment_panel.R, run deliberately (see header).")
message("Figure 4's inputs come from the Python etcGEM layer: see gem/README.md.")
message(strrep("=",70))
