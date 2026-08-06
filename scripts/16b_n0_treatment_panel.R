# =============================================================================
# 16b_n0_treatment_panel.R  ->  Supplementary Figure 6
#   Between-clade respiration activation energy under THREE treatments of the
#   inoculum back-projection, as one figure.
# =============================================================================
# Per-cell respiration is  R = K / N0  with  N0 = N_inoc * exp(r * delta), so R
# carries an explicit exp(-r * delta) factor. This script shows how much the
# BETWEEN-CLADE respiration activation-energy ordering leans on that term, by
# refitting the SAME hierarchical Bayesian Arrhenius model (via 09) three ways
# and plotting each clade's E side by side:
#
#   (1) WITH term        published:  N0 = N_inoc * exp(r * delta)
#   (2) EQUILIBRATION    N0 recomputed from the logged internal vial temperature
#                        over the onset window (script 08's per-curve resp_corr).
#                        The physically-motivated middle case.
#   (3) TERM-FREE        N0 = N_inoc (term removed entirely). NOT the physical
#                        null, but bounds the correction from the opposite side.
#
# It reuses 10_n0_term_test.R's mechanism: write the chosen respiration into
# derived_csv, source 09 in a fresh env, harvest b_E_Group* from its summary,
# then restore the original table.
#
# INPUT  tables/derived_N0_R_results_with_carbon.csv
#        tables/temperature_equilibration_percurve.csv   (from 08; arm 2 only)
# OUTPUT tables/bayes_resp_arr_E_three_treatments.csv
#        figures/Fig_n0_treatment_panel.png              (Supplementary Fig. 6)
#
# NOTE  refits the Bayesian respiration model up to THREE times -> SLOW (tens of
#       minutes). Run 02->03->06->07 first (and 08 for arm 2). For a quick look,
#       lower BAYES_ITER / BAYES_CHAINS in 09_bayesian_models.R first.
# NOTE  leaves derived_csv restored to the original at the end, but the other
#       bayes_*/rds outputs reflect the LAST fit; re-run 09 to refresh them.
# RUN   Rscript scripts/16b_n0_treatment_panel.R
# =============================================================================

# ---- locate scripts/ dir + shared config (same pattern as the other scripts) -
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

# config.R installs a figure whitelist (FIG_KEEP) that silently drops any ggsave
# whose base name is not listed. Add this figure so it is actually written.
if (exists("FIG_KEEP")) FIG_KEEP <- unique(c(FIG_KEEP, "Fig_n0_treatment_panel"))

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2)
})
has_patch <- requireNamespace("patchwork", quietly = TRUE)

stopifnot(file.exists(derived_csv))
orig <- readr::read_csv(derived_csv, show_col_types = FALSE)
summary_path <- file.path(tables_dir, "bayes_resp_arr_summary.csv")

# Pretty clade labels for the axis (falls back to the raw group name).
clade_label <- function(x) {
  m <- c(Clade1 = "Clade I", Clade2 = "Clade II", Clade3 = "Clade III",
         Clade4 = "Clade IV", para = "C. parapsilosis", parapsilosis = "C. parapsilosis")
  ifelse(x %in% names(m), m[x], x)
}

# Refit 09 on whatever is currently in derived_csv, pull per-clade b_E_Group*.
run_and_get_E <- function(label, save_as) {
  message("\n>>> Fitting Bayesian respiration model: ", label, " ...")
  source(file.path(.this_dir, "09_bayesian_models.R"), local = new.env())
  s <- readr::read_csv(summary_path, show_col_types = FALSE)
  file.copy(summary_path, file.path(tables_dir, save_as), overwrite = TRUE)
  s %>%
    dplyr::filter(grepl("^b_E_Group", variable)) %>%
    dplyr::transmute(clade = sub("^b_E_Group", "", variable),
                     E = mean, lwr = q2.5, upr = q97.5, rhat = rhat,
                     treatment = label)
}

TREAT <- c(WITH  = "With term (published)",
           EQUIL = "Equilibration-corrected",
           FREE  = "Term-free")

# ---- (1) WITH term -- published derived table, untouched --------------------
E_with <- run_and_get_E(TREAT[["WITH"]], "bayes_resp_arr_E_WITH_term.csv")

# ---- (2) EQUILIBRATION -- swap in script 08's per-curve corrected respiration
percurve <- file.path(tables_dir, "temperature_equilibration_percurve.csv")
E_equil <- NULL
if (file.exists(percurve)) {
  pc <- readr::read_csv(percurve, show_col_types = FALSE) %>%
    dplyr::select(T, OTU, Replicate, resp_corr) %>%
    dplyr::mutate(T = as.numeric(T), OTU = as.integer(OTU),
                  Replicate = toupper(as.character(Replicate)))
  eq <- orig %>%
    dplyr::mutate(T = as.numeric(T), OTU = as.integer(OTU),
                  Replicate = toupper(as.character(Replicate))) %>%
    dplyr::left_join(pc, by = c("T", "OTU", "Replicate")) %>%
    dplyr::mutate(
      .scale = dplyr::if_else(is.finite(resp_corr) & respiration_fgC_h > 0,
                              resp_corr / respiration_fgC_h, 1),
      respiration_fgC_h = dplyr::if_else(is.finite(resp_corr), resp_corr, respiration_fgC_h),
      CUE = growth_fgC_h / (growth_fgC_h + respiration_fgC_h)) %>%
    dplyr::select(-resp_corr)
  if ("respiration_C_per_C_h" %in% names(eq))
    eq <- dplyr::mutate(eq, respiration_C_per_C_h = respiration_C_per_C_h * .scale)
  eq <- dplyr::select(eq, -.scale)
  readr::write_csv(eq, derived_csv)
  E_equil <- run_and_get_E(TREAT[["EQUIL"]], "bayes_resp_arr_E_EQUIL.csv")
} else {
  message("\n[skip] ", percurve, " not found -- run 08 first to include the ",
          "equilibration-corrected arm. Proceeding with WITH + TERM-FREE only.")
}

# ---- (3) TERM-FREE -- undo the back-projection: R -> R * exp(r * delta) ------
free <- orig %>%
  dplyr::mutate(
    respiration_fgC_h = respiration_fgC_h * exp(r * delta_Ninoc_to_N0_min),
    CUE = growth_fgC_h / (growth_fgC_h + respiration_fgC_h))
if ("respiration_C_per_C_h" %in% names(free))
  free <- dplyr::mutate(free,
    respiration_C_per_C_h = respiration_C_per_C_h * exp(r * delta_Ninoc_to_N0_min))
readr::write_csv(free, derived_csv)
E_free <- run_and_get_E(TREAT[["FREE"]], "bayes_resp_arr_E_FREE.csv")

# ---- restore the original derived table -------------------------------------
readr::write_csv(orig, derived_csv)
message("\nRestored original ", basename(derived_csv),
        " (other bayes_*/rds outputs reflect the last fit; re-run 09 to refresh).")

# ---- combine + save ---------------------------------------------------------
allE <- dplyr::bind_rows(E_with, E_equil, E_free) %>%
  dplyr::mutate(clade_lab = clade_label(clade),
                treatment = factor(treatment, levels = unname(TREAT)))
readr::write_csv(allE, file.path(tables_dir, "bayes_resp_arr_E_three_treatments.csv"))
cat("\n== respiration activation energy E per clade, by N0 treatment ==\n")
print(as.data.frame(allE %>%
  dplyr::group_by(treatment) %>% dplyr::arrange(E, .by_group = TRUE) %>%
  dplyr::mutate(rank = dplyr::row_number()) %>% dplyr::ungroup() %>%
  dplyr::select(treatment, clade, E, lwr, upr, rank)), row.names = FALSE, digits = 3)

# ---- figure (Supplementary Figure 6) ----------------------------------------
# Order clades by the published (WITH-term) E; highlight Clade IV.
ord <- E_with %>% dplyr::arrange(E) %>%
  dplyr::mutate(clade_lab = clade_label(clade)) %>% dplyr::pull(clade_lab)
allE$clade_lab <- factor(allE$clade_lab, levels = rev(ord))
allE$hi <- allE$clade == "Clade4"

p <- ggplot(allE, aes(E, clade_lab, colour = hi)) +
  geom_errorbar(aes(xmin = lwr, xmax = upr), orientation = "y", width = 0, linewidth = 0.9) +
  geom_point(size = 2.4) +
  facet_wrap(~ treatment, nrow = 1) +
  scale_colour_manual(values = c(`FALSE` = "grey45", `TRUE` = "#C0392B"), guide = "none") +
  labs(x = "respiration activation energy  E  (eV)", y = NULL,
       title = "Between-clade respiration activation energy under three N0 treatments",
       subtitle = paste0("C. auris (red = Clade IV) vs C. parapsilosis separation holds throughout; ",
                         "the fine within-auris order does not survive term removal.")) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        panel.grid.minor = element_blank())

ggsave(file.path(figures_dir, "Fig_n0_treatment_panel.png"),
       p, width = 12, height = 3.6, dpi = 200, bg = "white")
cat("\nwrote Fig_n0_treatment_panel.png (Supplementary Figure 6) and the E table.\n")
