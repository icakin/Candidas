# =============================================================================
# C11 common: paths, constants, data loading
# =============================================================================
# Shared by 01_balance_quantity.R, 02_scale_free.R and 03_gap_decomposition.R.
# READ-ONLY with respect to the pipeline: nothing here writes into results/ or
# data/, and no estimator, constant or committed number is touched. Every output
# lands under reports/C11_scale_free/.
# =============================================================================

.this_dir <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) return(dirname(normalizePath(
    # R replaces every space in --file= with "~+~"; this project path has one.
    gsub("~+~", " ", sub("^--file=", "", fa[1]), fixed = TRUE),
    mustWork = FALSE)))
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
})

C11_DIR  <- dirname(.this_dir)                  # reports/C11_scale_free
PROJ     <- dirname(dirname(C11_DIR))           # project root
SCRIPTS  <- file.path(PROJ, "scripts")

# config.R is the authority for the carbon constants and for EXCLUDE_GROUPS, so
# it is sourced rather than duplicated. It also installs the figure whitelist,
# which gates ggsave() on the figure's base name -- so this report declares its
# own figures through fig_keep_add(), the entry point added in C10.
.this_dir <- SCRIPTS                            # config.R reads this to find base_dir
source(file.path(SCRIPTS, "config.R"))

C11_TAB <- file.path(C11_DIR, "tables");  dir.create(C11_TAB, showWarnings = FALSE, recursive = TRUE)
C11_FIG <- file.path(C11_DIR, "figures"); dir.create(C11_FIG, showWarnings = FALSE, recursive = TRUE)

fig_keep_add(
  "c11_fig1_activation_energies", "c11_fig2_balance_vs_T",
  "c11_fig3_modelfree_vs_fitted", "c11_fig4_gap_decomposition")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

# ---- the conversion block, exactly as 07_oxygen_fits.R applies it -----------
#   respiration_fgC_h = R_O2_mg_cell_min * MG_TO_FG * O2_TO_C_MASS * RQ * MIN_TO_H
C11_C <- MG_TO_FG * O2_TO_C_MASS * RESPIRATORY_QUOTIENT * MIN_TO_H

# ---- inputs -----------------------------------------------------------------
C11_INPUTS <- c(
  derived  = file.path(tables_dir, "derived_N0_R_results_with_carbon.csv"),
  names    = file.path(tables_dir, "otu_names.csv"),
  percurve = file.path(tables_dir, "temperature_equilibration_percurve.csv"),
  figvals  = file.path(tables_dir, "fig_values.csv"),
  cuecoefs = file.path(tables_dir, "cue_quadratic_fit_coefs.csv"),
  fitwide  = file.path(tables_dir, "fit_coefficients_wide.csv")
)

c11_load <- function() {
  d  <- readr::read_csv(C11_INPUTS[["derived"]], show_col_types = FALSE)
  nm <- readr::read_csv(C11_INPUTS[["names"]],   show_col_types = FALSE) %>%
    dplyr::transmute(OTU = as.integer(OTU), group = as.character(group),
                     otu_label = as.character(otu_name))
  d %>%
    dplyr::mutate(OTU = as.integer(OTU), T = as.numeric(T),
                  Replicate = toupper(as.character(Replicate))) %>%
    dplyr::left_join(nm, by = "OTU") %>%
    # the pipeline's own exclusion (config.R EXCLUDE_GROUPS = "glab"), so this
    # analysis sees exactly the five taxa the manuscript reports
    drop_excluded("C11 derived") %>%
    dplyr::filter(is.finite(r), r > 0, is.finite(K), K > 0)
}

# Display order and labels, matching config.R's GROUP_LABEL_1L.
C11_TAXA  <- c("Clade1", "Clade2", "Clade3", "Clade4", "para")
c11_label <- function(g) unname(GROUP_LABEL_1L[g])

# Boltzmann-Arrhenius on the log scale: log X = a + E * (1/(k T_ref) - 1/(k TK))
c11_boltz <- function(TC) (1 / (k_B * T_ref)) - (1 / (k_B * (TC + 273.15)))

T_BODY <- 37
