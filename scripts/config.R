# =============================================================================
# config.R - Shared paths, constants, and helper functions
# =============================================================================
# Sourced by every script in the pipeline.
#
# Project: Candida oxygen-respiration as a function of temperature and isolate.
# The experimental variables are temperature and isolate; there is no pH or
# dose treatment.
#
# PLATE LAYOUT (important - this is what makes the data tidy):
#   Each raw Excel file is ONE 24-well PreSens plate recorded at ONE temperature
#   for ONE clade/species GROUP (the group is the file-name prefix, e.g.
#   "Clade1", "glab", "para"). On every plate:
#       - rows A, B, C   = 3 different ISOLATES of that group
#       - columns 1..5   = 5 replicate wells per isolate
#       - row D and column 6 are empty ("No Sensor") and are ignored.
#   So one file yields 3 isolates x 5 replicates = 15 oxygen series.
#
# The six groups are:
#   Clade1..Clade4 = Candida auris, four clades from four geographical regions
#   glab           = Candida glabrata
#   para           = Candida parapsilosis
# giving 6 groups x 3 isolates = 18 isolates in total. Each isolate is given a
# unique integer code (the OTU column, 1..18) by 01_convert_xlsx.R and carries a
# human-readable display name (tables/otu_names.csv) that you type in that app.
#
# Experimental-design columns expected throughout the pipeline:
#   T          numeric, temperature in degC (measured plate Tm: 22..44)
#   OTU        integer, unique ISOLATE code 1..18 (display name in otu_names.csv)
#   Replicate  character, "R1" / "R2" / ... / "R5"
#   Time       numeric, minutes from start of recording
#   Oxygen     numeric, dissolved O2 in mg/L
#
# Each series is uniquely identified by (T, OTU, Replicate). The pipeline
# reshapes, trims, and fits a respiration model per series, then derives
# growth, respiration and CUE, and fits per-isolate thermal-performance curves
# (frequentist Sharpe-Schoolfield + Boltzmann-Arrhenius). No Bayesian stage.
#
# Edit base_dir (and any settings below) in this one place only.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ===== Project root and output directories ====================================
# base_dir is normally AUTO-DETECTED as the parent of the scripts/ folder (this
# file lives in scripts/config.R), so you can move or rename the whole project
# folder - or run it on another computer - without editing anything here.
#
# If auto-detection ever fails, set base_dir_manual to the absolute path of the
# project folder (the one that contains data/, scripts/, tables/, figures/).
base_dir_manual <- ""   # e.g. "/Users/you/Desktop/Projects/Candidas"

.detect_base_dir <- function() {
  if (nzchar(base_dir_manual)) return(normalizePath(base_dir_manual, mustWork = FALSE))
  # Scripts define `.this_dir` (= the scripts/ folder) before sourcing this file.
  # base::exists / base::get are namespaced ON PURPOSE. A stray `get <- function(...)`
  # left in the global environment by any other script would otherwise mask base::get
  # and blow this up with a baffling "unused argument (inherits = TRUE)".
  if (base::exists(".this_dir", inherits = TRUE)) {
    d <- base::get(".this_dir", inherits = TRUE)
    if (length(d) == 1 && !is.na(d) && nzchar(d) && dir.exists(d)) {
      return(dirname(normalizePath(d, mustWork = FALSE)))
    }
  }
  # Fallbacks: this file's own folder (via source), else RStudio, else getwd().
  d <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  if ((is.null(d) || is.na(d) || !nzchar(d)) &&
      requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    d <- tryCatch(dirname(rstudioapi::getActiveDocumentContext()$path),
                  error = function(e) NA_character_)
  }
  if (is.null(d) || is.na(d) || !nzchar(d)) return(getwd())
  dirname(normalizePath(d, mustWork = FALSE))
}

base_dir    <- .detect_base_dir()
message("config.R: base_dir = ", base_dir)
data_dir    <- file.path(base_dir, "data")
tables_dir  <- file.path(base_dir, "tables")
figures_dir <- file.path(base_dir, "figures")
models_dir  <- file.path(base_dir, "models")

for (d in c(data_dir, tables_dir, figures_dir, models_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ===== Input files ============================================================

LONG_CSV      <- file.path(tables_dir, "Oxygen_All_Long.csv")
IN_CSV        <- file.path(tables_dir, "Oxygen_Data_Filtered.csv")
TRIM_META_CSV <- file.path(tables_dir, "Oxygen_Trimmed_Series_Metadata.csv")

# ===== Stage 1 output paths ==================================================

pdf_path         <- file.path(figures_dir, "per_series_fits.pdf")
coef_csv         <- file.path(tables_dir, "fit_coefficients_long.csv")
fit_metrics_csv  <- file.path(tables_dir, "fit_metrics.csv")
coef_wide_csv    <- file.path(tables_dir, "fit_coefficients_wide.csv")
derived_csv      <- file.path(tables_dir, "derived_N0_R_results_with_carbon.csv")
group_lookup_csv <- file.path(tables_dir, "group_lookup_with_delta.csv")

replication_summary_csv <- file.path(tables_dir,
                                     "replication_summary_by_temperature.csv")

# ===== Other plot outputs =====================================================

box_growth_png             <- file.path(figures_dir, "boxplot_growth_fgC_h.png")
box_resp_png               <- file.path(figures_dir, "boxplot_respiration_fgC_h.png")
scatter_growth_png         <- file.path(figures_dir, "scatter_growth_fgC_h_vs_T.png")
scatter_resp_png           <- file.path(figures_dir, "scatter_respiration_fgC_h_vs_T.png")
resp_over_growth_png       <- file.path(figures_dir, "resp_over_growth_vs_T.png")
cue_vs_T_png               <- file.path(figures_dir, "CUE_vs_T.png")
cue_shape_png              <- file.path(figures_dir, "CUE_temperature_response_shape.png")
box_growth_biomass_png     <- file.path(figures_dir, "boxplot_growth_C_per_C_h.png")
box_resp_biomass_png       <- file.path(figures_dir, "boxplot_respiration_C_per_C_h.png")
scatter_growth_biomass_png <- file.path(figures_dir, "scatter_growth_C_per_C_h_vs_T.png")
scatter_resp_biomass_png   <- file.path(figures_dir, "scatter_respiration_C_per_C_h_vs_T.png")
box_total_resp_png         <- file.path(figures_dir, "boxplot_total_respiration_O2_mg_per_L.png")
scatter_total_resp_png     <- file.path(figures_dir, "scatter_total_respiration_O2_mg_per_L_vs_T.png")
box_K_png                  <- file.path(figures_dir, "boxplot_K_O2_rate.png")
scatter_K_png              <- file.path(figures_dir, "scatter_K_O2_rate_vs_T.png")

# ===== Isolate / plate-group registry ========================================
# The six clade/species GROUPS, in a FIXED order. This order defines the unique
# isolate code (OTU) assignment used by 01_convert_xlsx.R: the first group's
# rows A/B/C become isolates 1,2,3; the second group's rows become 4,5,6; and so
# on, so every isolate across the whole dataset gets a unique code 1..18.
STRAIN_GROUPS <- c("Clade1", "Clade2", "Clade3", "Clade4", "glab", "para")

# Species behind each group (used only for labels/reference).
GROUP_SPECIES <- c(
  Clade1 = "Candida auris (Clade I)",
  Clade2 = "Candida auris (Clade II)",
  Clade3 = "Candida auris (Clade III)",
  Clade4 = "Candida auris (Clade IV)",
  glab   = "Candida glabrata",
  para   = "Candida parapsilosis"
)

# ---- Clade metadata (for figure labels and the discussion) -------------------
# Geographic origin of each C. auris clade. Established, stable designations:
#   Lancet Microbe 2024 (6th clade, Singapore) - reviews the clade I-V origins
#   https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(24)00101-0/fulltext
GROUP_ORIGIN <- c(
  Clade1 = "South Asia",
  Clade2 = "East Asia",
  Clade3 = "Africa",
  Clade4 = "South America",
  glab   = NA_character_,
  para   = NA_character_
)

# Published clade phenotypes, for the DISCUSSION only - never for a figure claim.
#   Clade I  : hyper-virulent, multi-drug resistant, major outbreak clade
#   Clade II : LEAST virulent (lowest mortality, murine bloodstream model), the
#              only generally drug-SUSCEPTIBLE clade, mostly ear infections rather
#              than invasive disease. The known phenotypic outlier.
#   Clade III: forms large cell aggregates
#   Clade IV : HIGHEST mortality in the murine model; hyper-virulent, MDR
#   Refs: https://journals.asm.org/doi/10.1128/spectrum.00498-23   (metabolic profiling)
#         https://journals.asm.org/doi/10.1128/jcm.00007-19        (East Asian clade)
#
# NOTE FOR THE WRITE-UP - READ THIS BEFORE COMPARING CLADES TO THE LITERATURE.
#
# From the REAL posterior (tables/fig_contrasts.csv), relative respiratory cost at
# 40 C. 8 of the 10 pairwise contrasts resolve:
#
#     Clade IV  vs para       0.42 [0.33-0.53]   credible
#     Clade III vs para       0.50 [0.39-0.64]   credible
#     Clade I   vs para       0.59 [0.45-0.76]   credible
#     Clade II  vs para       0.70 [0.53-0.91]   credible
#     Clade III vs Clade IV   1.20 [1.04-1.39]   credible
#     Clade II  vs Clade III  1.40 [1.14-1.71]   credible
#     Clade I   vs Clade IV   1.42 [1.20-1.67]   credible
#     Clade II  vs Clade IV   1.67 [1.39-2.03]   credible
#     Clade I   vs Clade II   0.85 [0.68-1.05]   NOT separable
#     Clade I   vs Clade III  1.18 [0.99-1.41]   NOT separable
#
# SO: all four C. auris clades credibly pay LESS at fever than C. parapsilosis, and
# Clade IV is credibly the cheapest of all. That is a clean species-level claim.
#
# (An earlier version of this note said "only 3 of 10 resolve". It was WRONG - it
# came from a normal-approximation cross-check that treated E, Eh, Th and E_R as
# independent. They are strongly correlated in the real posterior, which halves the
# interval on every derived quantity. Never rank or reject a contrast from a
# Gaussian approximation to this model; it is ~2x too conservative.)
#
# BUT STILL DO NOT CORRELATE WITH GEOGRAPHY. Our cost ranking is IV < III < I ~ II;
# published murine mortality is IV > I > III > II. The extremes agree (Clade IV:
# cheapest at fever AND most virulent; Clade II: dearest AND least virulent, the
# only drug-susceptible one). But the two clades in the MIDDLE - exactly where the
# rankings disagree - are the two we CANNOT separate (I vs III). Four points, with
# the disputed middle unresolvable, is not a correlation. Report the two extremes
# as a hypothesis; never draw a trend line through four clades.
GROUP_PHENOTYPE <- c(
  Clade1 = "hyper-virulent, MDR",
  Clade2 = "least virulent, drug-susceptible",
  Clade3 = "aggregating",
  Clade4 = "highest mortality, MDR",
  glab   = NA_character_,
  para   = NA_character_
)

# Display label: species + origin where there is one.
GROUP_LABEL <- c(
  Clade1 = "C. auris I\n(South Asia)",
  Clade2 = "C. auris II\n(East Asia)",
  Clade3 = "C. auris III\n(Africa)",
  Clade4 = "C. auris IV\n(South America)",
  glab   = "C. glabrata",
  para   = "C. parapsilosis"
)
GROUP_LABEL_1L <- c(
  Clade1 = "C. auris I (South Asia)",
  Clade2 = "C. auris II (East Asia)",
  Clade3 = "C. auris III (Africa)",
  Clade4 = "C. auris IV (South America)",
  glab   = "C. glabrata",
  para   = "C. parapsilosis"
)

# Plate geometry: isolates occupy rows A/B/C, replicates occupy columns 1..5.
PLATE_ROWS         <- c("A", "B", "C")   # one isolate per row
ISOLATES_PER_PLATE <- length(PLATE_ROWS) # 3
N_REPLICATES       <- 5                  # columns 1..5

# All isolate codes present in the dataset (1..18).
OTUS <- seq_len(length(STRAIN_GROUPS) * ISOLATES_PER_PLATE)

# Map a raw file-name prefix (any case, e.g. "clade3", "Para") to its canonical
# group name in STRAIN_GROUPS, or NA if it is not one of the known groups.
canonical_group <- function(prefix) {
  p <- tolower(trimws(as.character(prefix)))
  key <- tolower(STRAIN_GROUPS)
  idx <- match(p, key)
  ifelse(is.na(idx), NA_character_, STRAIN_GROUPS[idx])
}

# The globally-unique isolate code for (group, row). Rows are A/B/C -> 1/2/3
# within a group; groups are offset by 3 in STRAIN_GROUPS order.
isolate_code <- function(group, row) {
  g <- match(canonical_group(group), STRAIN_GROUPS)
  r <- match(toupper(row), PLATE_ROWS)
  as.integer((g - 1L) * ISOLATES_PER_PLATE + r)
}

# ===== Groups held out of the ANALYSIS =======================================
# C. glabrata is being re-run. Its per-cell respiration is ~2x noisier than every
# other taxon (replicate-level sd(log R) = 0.455 vs 0.221; growth 0.214 vs 0.090),
# because its low per-cell O2 consumption leaves the optode short of signal. That
# is an assay limit, not a fitting problem: the rt-cap fix in 04 already removed a
# 2.5x bias and 43% of the noise, and it is STILL 2x noisier than the rest.
#
# Nothing is deleted. The raw files, the 04 fit windows and the 06 per-curve fits
# all still include glab, so the moment the repeat lands you set this to
# character(0), re-run 06 onwards, and it is back. It is simply not carried into
# the ANALYSIS scripts (07, 10, and every figure downstream of them).
#
# >>> STRAIN_GROUPS IS DELIBERATELY NOT EDITED. <<<
# isolate_code() derives OTU numbers from a group's POSITION in STRAIN_GROUPS.
# Removing "glab" there would renumber para from 16,17,18 to 13,14,15 and silently
# invalidate every existing window, exclusion and table in tables/. Never do that.
EXCLUDE_GROUPS <- c("glab")

ACTIVE_GROUPS <- setdiff(STRAIN_GROUPS, EXCLUDE_GROUPS)

# The OTU codes belonging to the excluded groups (glab -> 13, 14, 15).
EXCLUDED_OTUS <- if (length(EXCLUDE_GROUPS)) {
  sort(as.integer(unlist(lapply(EXCLUDE_GROUPS, function(g)
    vapply(PLATE_ROWS, function(rw) isolate_code(g, rw), integer(1))))))
} else integer(0)

# Drop excluded groups from any data frame carrying Group / group / OTU.
# Returns the frame unchanged when EXCLUDE_GROUPS is empty.
drop_excluded <- function(df, what = "data") {
  if (!length(EXCLUDE_GROUPS) || is.null(df) || nrow(df) == 0) return(df)
  n0 <- nrow(df)
  if ("Group" %in% names(df)) {
    df <- df[!(as.character(df$Group) %in% EXCLUDE_GROUPS), , drop = FALSE]
  } else if ("group" %in% names(df)) {
    df <- df[!(as.character(df$group) %in% EXCLUDE_GROUPS), , drop = FALSE]
  } else if ("OTU" %in% names(df)) {
    df <- df[!(as.integer(df$OTU) %in% EXCLUDED_OTUS), , drop = FALSE]
  } else {
    return(df)
  }
  if (nrow(df) < n0) {
    message(sprintf("EXCLUDE_GROUPS: dropped %d of %d %s rows (%s; OTUs %s).",
                    n0 - nrow(df), n0, what,
                    paste(EXCLUDE_GROUPS, collapse = ", "),
                    paste(EXCLUDED_OTUS, collapse = ", ")))
  }
  df
}

if (length(EXCLUDE_GROUPS)) {
  message("EXCLUDE_GROUPS is set: ", paste(EXCLUDE_GROUPS, collapse = ", "),
          "  ->  analysis uses ", paste(ACTIVE_GROUPS, collapse = ", "))
}


# ===== User settings ==========================================================

# NOTE: there is no per-series RMSE quality filter. Every series that yields a
# usable oxygen-model fit is kept and appears in the results and on every plot.

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>> MANUAL OVERRIDES - paste the WHOLE block from the trim-selector app here >>>
# >>> (replace everything between this line and the matching END marker below)  >>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Two tables (optionally pasted from the trim-selector app):
#   PLOT_EXCLUDE_POINTS = samples dropped from the plots.
#   MANUAL_FIT_WINDOWS  = per-curve start/end (NA = automatic on that side).
# Start EMPTY: every curve is trimmed automatically by 03_trimming.R. Run
# 04_trim_selector.R to review curves and fill these in (or keep
# USE_APP_TRIM_FILES <- TRUE below to load the app's saved files directly).
PLOT_EXCLUDE_POINTS <- data.frame(
  T         = numeric(0),
  OTU       = numeric(0),
  Replicate = character(0),
  stringsAsFactors = FALSE
)

MANUAL_FIT_WINDOWS <- data.frame(
  T         = numeric(0),
  OTU       = numeric(0),
  Replicate = character(0),
  fit_start = numeric(0),
  fit_end   = numeric(0),
  stringsAsFactors = FALSE
)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# <<< END MANUAL OVERRIDES <<<
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Auto-connect the trim-selector app: if it has saved its files in tables/, load
# them DIRECTLY (no copy-paste). The app writes manual_fit_windows.csv and
# plot_exclude_points.csv; when present they REPLACE the two blocks above, so the
# pipeline just picks up whatever you set in the app. Set USE_APP_TRIM_FILES to
# FALSE to ignore the app files and use the hand-edited blocks above instead.
USE_APP_TRIM_FILES <- TRUE
if (isTRUE(USE_APP_TRIM_FILES)) {
  .mfw_csv <- file.path(tables_dir, "manual_fit_windows.csv")
  if (file.exists(.mfw_csv)) {
    .tmp <- tryCatch(read.csv(.mfw_csv, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(.tmp) &&
        all(c("T", "OTU", "Replicate", "fit_start", "fit_end") %in% names(.tmp))) {
      MANUAL_FIT_WINDOWS <- .tmp
      message("  Loaded MANUAL_FIT_WINDOWS from app file (", nrow(.tmp), " curves).")
    }
  }
  .pep_csv <- file.path(tables_dir, "plot_exclude_points.csv")
  if (file.exists(.pep_csv)) {
    .tmp <- tryCatch(read.csv(.pep_csv, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(.tmp) && all(c("T", "OTU", "Replicate") %in% names(.tmp))) {
      PLOT_EXCLUDE_POINTS <- .tmp
      message("  Loaded PLOT_EXCLUDE_POINTS from app file (", nrow(.tmp), " points).")
    }
  }
}


# Fit-window standardization (used in 07_oxygen_fits.R).
# Ending each fit window at the per-curve steepest-drop point made the window
# length vary and biased the fitted K (corr(K, window) ~ -0.55), a main driver
# of replicate-to-replicate respiration scatter. Instead, end the window when O2
# has fallen by a fixed FRACTION of its total drawdown from the green line, so
# every series is fit over the same slice of its curve and K is estimated
# consistently.
#   USE_DRAWDOWN_WINDOW  TRUE  = fractional-drawdown end (standardized).
#                        FALSE = legacy steepest-drop / chosen-end behaviour.
#   FIT_DRAWDOWN_FRAC    fraction (0-1) of total O2 drawdown to include. Lower =
#                        cleaner exponential phase but fewer points; ~0.45 stays
#                        below the inflection. Tune and re-check corr(K, window).
USE_DRAWDOWN_WINDOW <- FALSE
FIT_DRAWDOWN_FRAC   <- 0.45

# ===== USER INPUT =============================================================

# Inoculation density (cells/L). This anchors N0 and therefore the absolute
# growth / respiration (and CUE, which depends on cell carbon). Set it to your
# actual Candida inoculation density. Default: 1.5e7 cells/L (= 15,000 cells/mL).
N_inoculation_cells_per_L <- 1.5e7

# ---------------------------------------------------------------------------
# DATA-SCOPE FILTERS
# ---------------------------------------------------------------------------
# 02_longdata.R applies these filters immediately after loading, so every
# downstream script sees only the rows the user wants.
#
#   EXCLUDE_TEMPS   numeric vector of temperatures (degC) to drop. NULL = none.
#   EXCLUDE_OTUS    integer vector of ISOLATE codes (1..18) to drop. NULL = none.
EXCLUDE_TEMPS <- NULL
EXCLUDE_OTUS  <- NULL

# EXCLUDE_T_OTU drops specific (temperature, isolate) groups from the derived
# results, all downstream calculations, and every figure. It is applied in
# 07_oxygen_fits.R at the results stage, so the raw 02/03 trimming diagnostics
# still show every curve. Use a two-column data.frame (T in degC, OTU 1..18);
# set to NULL to drop none.
EXCLUDE_T_OTU <- NULL

# ===== Thermal constants ======================================================

k_B    <- 0.00008617
T_ref  <- 293.15
INV_KB <- 1 / k_B

# ===== Cell carbon conversion =================================================
# The values below describe a Candida (yeast) cell modelled as a prolate
# ellipsoid, width 3 um and length 4 um, at ~100 fg C/um3 (Bratbak 1985). They
# set the carbon per cell, which scales growth (fg C/h) and hence CUE. They do
# NOT change respiration_fgC_h or the temperature-response shapes.
#
# PER-ISOLATE sizes: run 05_cell_sizes.R to give each isolate its own average
# size. It writes tables/otu_cell_sizes.csv, which 07_oxygen_fits.R uses per
# isolate; the values below are the GLOBAL fallback for any isolate not listed.
CELL_WIDTH_UM  <- 3
CELL_LENGTH_UM <- 4
CARBON_DENSITY_FG_PER_UM3 <- 100
RESPIRATORY_QUOTIENT <- 1

cell_radius_um <- CELL_WIDTH_UM / 2
cell_cyl_length_um <- max(CELL_LENGTH_UM - CELL_WIDTH_UM, 0)

CELL_VOLUME_UM3 <-
  pi * (cell_radius_um^2) * cell_cyl_length_um +
  (4 / 3) * pi * (cell_radius_um^3)

CELL_CARBON_FG_PER_CELL <- CELL_VOLUME_UM3 * CARBON_DENSITY_FG_PER_UM3
MG_TO_FG <- 1e12
MIN_TO_H <- 60

# O2 mass -> C mass conversion. Carbon respired (as CO2) per unit O2 consumed:
#   mass_C = mass_O2 * (M_C / M_O2) * RQ
# Without this factor, respiration would be expressed in O2 mass, not carbon,
# and would not be comparable to growth_fgC_h (which is real cell carbon).
M_C_G_PER_MOL  <- 12.011
M_O2_G_PER_MOL <- 31.998
O2_TO_C_MASS   <- M_C_G_PER_MOL / M_O2_G_PER_MOL   # ~0.3754

# N0 anchoring for per-cell respiration (respiration = K / N0).
#   TRUE  = N0 = N_inoc * exp(r * delta), delta = time to the green line.
#           Back-projects biomass across the (largely lag-phase) interval before
#           the fit window; the exp(r * delta) term injects the inverted growth
#           curve and fit noise into respiration -> no thermal optimum, noisy.
#   FALSE = N0 = N_inoc (delta = 0): anchor at the consumption onset. Respiration
#           then scales with K (a shape-derived rate) divided by the same stock
#           constant for every curve, mirroring how growth = r * constant. This
#           removes the artifact so respiration shows a clean TPC comparable to
#           growth. fit_start_time is still recorded for reference.
N0_BACKPROJECT <- TRUE

# ===== Helper functions =======================================================

get_baseline <- function(x, head_max = 30, min_valid = 5) {
  x <- as.numeric(x)
  n <- length(x)
  if (n == 0) return(NA_real_)

  k <- min(head_max, n)
  hv <- x[seq_len(k)]
  hv <- hv[is.finite(hv)]
  if (length(hv) >= min_valid) {
    b <- median(hv, na.rm = TRUE)
    if (is.finite(b)) return(b)
  }
  suppressWarnings(mean(x, na.rm = TRUE))
}

resp_model <- function(r, K, t, O2_0) {
  O2_0 + (K / r) * (1 - exp(r * t))
}

# Sharpe-Schoolfield (log scale) prediction, used by the frequentist TPC fits.
pred_ss_log <- function(TK, lnB0, E, Eh, Th) {
  lnB0 -
    (E * 11604.51812) * ((1 / TK) - (1 / 293.15)) -
    log(1 + exp((Eh * 11604.51812) * ((1 / Th) - (1 / TK))))
}

# Boltzmann-Arrhenius (log scale) prediction.
pred_arr_log <- function(TK, alpha, E, k_B = 0.00008617, T_ref = 293.15) {
  boltz_shift <- (1 / (k_B * T_ref)) - (1 / (k_B * TK))
  alpha + E * boltz_shift
}

make_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}

message("config.R loaded: base_dir = ", base_dir)
message("  Project: Candida (temperature x isolate) | groups = {",
        paste(STRAIN_GROUPS, collapse = ", "), "} | isolates = ",
        length(OTUS), " | inoc = ", N_inoculation_cells_per_L, " cells/L")
message("  Cell volume: ", round(CELL_VOLUME_UM3, 2), " um^3 | C/cell: ",
        round(CELL_CARBON_FG_PER_CELL, 2), " fg")
