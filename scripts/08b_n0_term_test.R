# =============================================================================
# 08b_n0_term_test.R - Does the between-clade respiration ordering depend on the
#                      N0 back-projection (r * delta) term?
# =============================================================================
# Per-cell respiration is  R = K / N0  with  N0 = N_inoc * exp(r * delta).
# So R carries an explicit exp(-r * delta) factor: fast-growing curves get their
# respiration pulled DOWN the most. This script asks whether the manuscript's
# between-clade respiration-activation-energy ordering (e.g. "Clade IV is the
# cheapest at fever") is real biology or an artefact of that term.
#
# It re-fits the SAME Bayesian respiration model as 08_bayesian_models.R, twice:
#   (1) WITH the term     - the published derived table, untouched.
#   (2) WITHOUT the term  - respiration multiplied back up by exp(r * delta),
#                           i.e. N0 = N_inoc (no back-projection), CUE recomputed.
# then prints each clade's activation energy E and Clade4's rank each way.
#
# INPUT  : tables/derived_N0_R_results_with_carbon.csv  (run 02->03->06->07 first)
# OUTPUT : tables/bayes_resp_arr_WITH_term.csv
#          tables/bayes_resp_arr_WITHOUT_term.csv
#          console: the two rankings + the headline
#
# NOTE - it re-runs the full Bayesian fit (via 08) TWICE, so it is SLOW
#        (tens of minutes; longer on few cores). For a quick look you can lower
#        BAYES_ITER / BAYES_CHAINS in 08_bayesian_models.R first.
# NOTE - it restores the original derived table at the end, but the other bayes_*
#        / rds outputs will reflect the LAST (WITHOUT-term) fit; re-run
#        08_bayesian_models.R normally afterwards to refresh the published ones.
# RUN  : Rscript scripts/08b_n0_term_test.R   (or "Source" in RStudio)
# =============================================================================

# ---- locate scripts/ dir and shared config (same pattern as your other scripts)
.this_dir <- if (
  requireNamespace("rstudioapi", quietly = TRUE) &&
  rstudioapi::isAvailable() &&
  nzchar(rstudioapi::getActiveDocumentContext()$path)
) {
  dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1])))
  else tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
}
source(file.path(.this_dir, "config.R"))

suppressPackageStartupMessages({ library(dplyr); library(readr) })

stopifnot(file.exists(derived_csv))
orig <- readr::read_csv(derived_csv, show_col_types = FALSE)

summary_path <- file.path(tables_dir, "bayes_resp_arr_summary.csv")

# Re-run 08's Bayesian fit on whatever is currently in derived_csv, then pull the
# per-clade respiration activation energy (b_E_Group*) out of its summary.
run_and_get_E <- function(label, save_as) {
  message("\n>>> Fitting Bayesian respiration model: ", label, " ...")
  source(file.path(.this_dir, "08_bayesian_models.R"), local = new.env())
  s <- readr::read_csv(summary_path, show_col_types = FALSE)
  file.copy(summary_path, file.path(tables_dir, save_as), overwrite = TRUE)
  s <- s %>%
    dplyr::filter(grepl("^b_E_Group", variable)) %>%
    dplyr::transmute(clade = sub("^b_E_Group", "", variable),
                     E = mean, lwr = q2.5, upr = q97.5, rhat = rhat) %>%
    dplyr::arrange(E) %>%
    dplyr::mutate(rank = dplyr::row_number())
  message("=== ", label, " ===")
  print(as.data.frame(s), row.names = FALSE, digits = 3)
  s
}

# (1) WITH the term - published derived table, untouched
with_E <- run_and_get_E("WITH r*delta term (published)", "bayes_resp_arr_WITH_term.csv")

# (2) WITHOUT the term - undo the back-projection: R -> R * exp(r * delta)
noterm <- orig %>%
  dplyr::mutate(
    respiration_fgC_h = respiration_fgC_h * exp(r * delta_Ninoc_to_N0_min),
    CUE = growth_fgC_h / (growth_fgC_h + respiration_fgC_h)
  )
if ("respiration_C_per_C_h" %in% names(noterm))
  noterm <- dplyr::mutate(noterm,
    respiration_C_per_C_h = respiration_C_per_C_h * exp(r * delta_Ninoc_to_N0_min))
readr::write_csv(noterm, derived_csv)
without_E <- run_and_get_E("WITHOUT r*delta term", "bayes_resp_arr_WITHOUT_term.csv")

# ---- restore the original derived table so nothing is left modified ----------
readr::write_csv(orig, derived_csv)
message("\nRestored original ", basename(derived_csv),
        " (other bayes_*/rds outputs reflect the last fit; re-run 08 to refresh).")

# ---- headline ----------------------------------------------------------------
focus <- "Clade4"
rw <- with_E$rank[with_E$clade == focus]
ro <- without_E$rank[without_E$clade == focus]
n  <- nrow(with_E)
message(sprintf(
  "\n%s respiration-cost rank:  WITH term = %d,  WITHOUT term = %d  (of %d groups).",
  focus, rw, ro, n))
if (!is.na(rw) && !is.na(ro) && rw != ro)
  message("=> The ordering is NOT robust to the N0 back-projection term.")
