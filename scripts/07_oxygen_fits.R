# ===== Stage 1: Oxygen model fits, carbon unit conversions, and descriptive plots
#
# This script (Candida temperature x isolate design):
#  - Fits exponential-decay respiration models to dissolved oxygen time series
#  - Uses the trimmed time window from 03_trimming.R
#  - Fits the model to raw Oxygen values, not O2_fit
#  - Each series is one (T, OTU, Replicate) combination (OTU = isolate 1..18)
#  - Computes initial cell counts (N0) and carbon-based metrics from oxygen consumption
#  - Derives growth rate, respiration rate, and CUE (Carbon Use Efficiency)
#  - Filters outlier points and generates descriptive plots vs temperature,
#    coloured and facetted by isolate
#  - Fits per-isolate thermal-performance curves (Sharpe-Schoolfield + Arrhenius)
#
# Main output: tables/derived_N0_R_results_with_carbon.csv (one row per
# T x OTU x Replicate with growth_fgC_h, respiration_fgC_h, CUE, and their
# biomass-corrected analogues). This is the final stage of this build.


# =============================================================================
# 0) Source shared config
# =============================================================================

.this_dir <- if (
  requireNamespace("rstudioapi", quietly = TRUE) &&
  rstudioapi::isAvailable() &&
  nzchar(rstudioapi::getActiveDocumentContext()$path)
) {
  dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
}

source(file.path(.this_dir, "config.R"))


# =============================================================================
# 1) Load packages
# ===========================================x==================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(rlang)
  library(minpack.lm)
  library(grid)
})

options(mc.cores = max(1, parallel::detectCores() - 1))


# =============================================================================
# 1b) User knobs
# =============================================================================
# FIT_TO_SPLINE controls what the nls fitter sees.
#
#   FALSE -> fit nlsLM to the raw observed oxygen values (Oxygen column).
#            This is the default and matches the historical behaviour.
#   TRUE  -> fit nlsLM to the smoothed spline values (O2_fit column) that
#            03_trimming.R writes alongside Oxygen. Requires O2_fit to be
#            present in IN_CSV; the script stops with an error otherwise.
#
# RMSE, R2, AIC, AICc, and the "kept" decision are all computed against the
# same column the fit was run on, so the diagnostics stay self-consistent
# under either mode.
# =============================================================================

FIT_TO_SPLINE <- FALSE

# ---- Curvature ("is there any growth?") test --------------------------------
# The BEND in the O2 curve IS the growth: as r -> 0 the respiration model
# collapses to a straight line (O2 = O2_0 - K*t) = constant respiration with NO
# growth (maintenance metabolism, e.g. above the thermal maximum).
#
# A straight trace still fits the exponential model very well (low RMSE, high
# R2), so an RMSE / quality filter can NEVER detect this.
#
# THE TEST: the DIMENSIONLESS curvature  rt = r * window_length.
#   O2 = O2_0 + (K/r)(1 - exp(r*t)).  If r*t is small over the whole window then
#   exp(r*t) ~ 1 + r*t and the model collapses to O2 ~ O2_0 - K*t, a STRAIGHT
#   LINE. There is then no curvature to estimate r from, so r is noise.
#     has_curvature = (r * window_length) > CURV_MIN_RT
#   rt >~ 1 means the consumption rate at least e-folds across the window (real
#   exponential signal). rt < 0.5 means effectively linear.
#
# NOTE: do NOT use an AICc line-vs-curve test here. Each curve has hundreds of
# time points, so AICc differences run to 500-3000 even for visually straight
# traces - the test can never fire. delta_aicc_curv is still reported below as a
# diagnostic, but the FLAG is based on rt.
#
# With no curvature, 07/09 drop the curve from the GROWTH fit but KEEP it for
# respiration (K, the linear O2 slope, is still a perfectly valid rate).
CURV_MIN_RT <- 0.5

# ---- Minimum signal: how much O2 did the culture actually consume? -----------
# The curvature test above asks "is the shape curved?" - it never asks "is there
# any signal to HAVE a shape?". A culture that draws down only ~0.3 mg/L of O2 is
# essentially dead, and the model will happily fit exponential curvature to the
# noise, returning a confident, near-maximal growth rate. (Seen at 42-44 degC in
# the para isolates: 0.34-1.7 mg/L O2 consumed, yet "growth" of 1400-3900 - as
# high as at the optimum. Impossible: you cannot build biomass without paying for
# it in respiration.)
#
# So we require a minimum total O2 drawdown for a rate to be believable. In this
# dataset the median curve consumes ~4.7 mg/L, so 2.0 is a low bar; it removes
# ~3% of curves - the dead ones, plus a few very-low-signal cold curves whose r
# was never reliable either.
MIN_O2_DRAWDOWN <- 2.0   # mg/L. Set to 0 to disable this check.

# ---- Minimum fit-window length: reject tail-window selector misfires ---------
# For a subset of series (16 of them at 26 degC, plus a handful at 24/40/42/44)
# the trimming step latched onto a tiny late drawdown and returned a fit window
# of only ~9-22 min instead of the main ~200-400 min run. Fitting the 3-parameter
# exponential to a ~10 min sliver returns nonsense r/K (r ~0.06-0.09 vs a real
# ~0.007). These pass every other check, so they need their own guard. A real
# main-run window here is hundreds of minutes; 60 is a low, safe floor.
MIN_FIT_WINDOW_MIN <- 60   # min. Set to 0 to disable this check.

# Temperatures to drop from the descriptive plots. Set to numeric(0) to keep
# every temperature. (The dataset contains 22, 24, ..., 44 degC; none excluded
# by default.)
EXCLUDE_TEMPS_PLOT <- numeric(0)


# =============================================================================
# 2) Read oxygen data
# =============================================================================
# IMPORTANT:
# The input file should be the trimmed output from 03_trimming.R,
# for example Oxygen_Data_Filtered.csv or Oxygen_Data_Smoothed_Trimmed.csv.
#
# The trimming script may include both:
#   Oxygen = raw observed oxygen values
#   O2_fit = smoothed spline used for trimming diagnostics
#
# Which one is fed to the nls fitter is controlled by FIT_TO_SPLINE above.
# =============================================================================

tmp_cols <- try(
  readr::read_csv(IN_CSV, n_max = 1, show_col_types = FALSE),
  silent = TRUE
)

if (inherits(tmp_cols, "try-error")) {
  stop("Could not read input file: ", IN_CSV)
}

# 03_trimming.R writes one row per (T, OTU, Replicate, Time). Confirm the
# temperature and otu columns are present before reading with a typed spec.
if (!"T" %in% names(tmp_cols)) {
  stop("Input file is missing the temperature column `T`: ", IN_CSV)
}
if (!"OTU" %in% names(tmp_cols)) {
  stop("Input file is missing the OTU column `OTU`: ", IN_CSV)
}

colspec <- list(
  T         = readr::col_double(),
  OTU     = readr::col_integer(),
  Replicate = readr::col_character(),
  Time      = readr::col_double(),
  Oxygen    = readr::col_double()
)

has_o2fit <- "O2_fit" %in% names(tmp_cols)

if (has_o2fit) {
  colspec$O2_fit <- readr::col_double()
}

if (isTRUE(FIT_TO_SPLINE) && !has_o2fit) {
  stop(
    "FIT_TO_SPLINE = TRUE but the input file does not contain an O2_fit column: ",
    IN_CSV,
    "\nRe-run 03_trimming.R (which writes O2_fit) or set FIT_TO_SPLINE <- FALSE."
  )
}

o2f <- readr::read_csv(IN_CSV, col_types = do.call(readr::cols, colspec))

o2f <- o2f %>%
  dplyr::mutate(
    T         = as.numeric(T),
    OTU     = as.integer(OTU),
    Replicate = toupper(as.character(Replicate)),
    series_id = paste0("T=", T, " | OTU=", OTU, " | Rep=", Replicate),

    # Pick the fit target based on FIT_TO_SPLINE.
    #   FALSE -> raw Oxygen (default, historical behaviour)
    #   TRUE  -> smoothed O2_fit from the trimming step
    Oxygen_used = if (isTRUE(FIT_TO_SPLINE)) O2_fit else Oxygen
  ) %>%
  dplyr::arrange(T, OTU, Replicate, Time)

stopifnot(all(c("T", "OTU", "Replicate", "Time", "Oxygen_used") %in% names(o2f)))

message(sprintf(
  "Fit target: %s (FIT_TO_SPLINE = %s)",
  if (isTRUE(FIT_TO_SPLINE)) "spline (O2_fit)" else "raw (Oxygen)",
  isTRUE(FIT_TO_SPLINE)
))

if (nrow(o2f) == 0) {
  stop("No rows in input oxygen file.")
}


# =============================================================================
# 3) Plot bookkeeping
# =============================================================================

ymin_all <- suppressWarnings(min(o2f$Oxygen_used, na.rm = TRUE))
ymax_all <- suppressWarnings(max(o2f$Oxygen_used, na.rm = TRUE))

pad_all <- 0.02 * (ymax_all - ymin_all)

if (!is.finite(pad_all)) {
  pad_all <- 0.1
}

Y_LIMITS_SERIES <- c(ymin_all - pad_all, ymax_all + pad_all)


# =============================================================================
# 4) Read trimming metadata
# =============================================================================

trim_meta_raw <- readr::read_csv(TRIM_META_CSV, show_col_types = FALSE)

if (!("T" %in% names(trim_meta_raw))) {
  stop("Trimming metadata is missing the temperature column `T` in: ", TRIM_META_CSV)
}
if (!("OTU" %in% names(trim_meta_raw))) {
  stop("Trimming metadata is missing the OTU column `OTU` in: ", TRIM_META_CSV)
}

# Per-series fit window = the zone between the two DASHED diagnostic lines that
# 03_trimming.R draws (02 and its figure are left untouched; we only fit a
# sub-range of the window it produced):
#   green dashed  = main_run_start_time  (start of the main descending run)
#   purple dashed = steepest_drop_time   (inflection / steepest descent)
# Fitting only this exponential-growth phase keeps the data inside the
# respiration model's assumption (exponentially growing biomass).
needed_meta <- c("main_run_start_time", "steepest_drop_time", "chosen_end_time")
missing_meta <- setdiff(needed_meta, names(trim_meta_raw))
if (length(missing_meta) > 0) {
  stop("Trimming metadata is missing column(s) needed for the dashed-line fit window: ",
       paste(missing_meta, collapse = ", "))
}

fit_windows <- trim_meta_raw %>%
  dplyr::mutate(
    T = as.numeric(T),
    OTU = as.integer(OTU),
    Replicate = toupper(as.character(Replicate)),
    main_run_start_time = as.numeric(main_run_start_time),
    steepest_drop_time  = as.numeric(steepest_drop_time),
    chosen_end_time     = as.numeric(chosen_end_time)
  ) %>%
  dplyr::mutate(
    fit_start_time = main_run_start_time,
    # Fall back to the chosen end if the steepest-drop time is missing or does
    # not sit after the run start (keeps the window valid for odd curves).
    fit_end_time = dplyr::if_else(
      is.finite(steepest_drop_time) & steepest_drop_time > main_run_start_time,
      steepest_drop_time,
      chosen_end_time
    )
  ) %>%
  dplyr::select(T, OTU, Replicate, fit_start_time, fit_end_time) %>%
  dplyr::distinct()

stopifnot(all(c("T", "OTU", "Replicate", "fit_start_time", "fit_end_time") %in% names(fit_windows)))

# ---------------------------------------------------------------------------
# Standardize the fit-window END by fractional O2 drawdown (see config:
# USE_DRAWDOWN_WINDOW / FIT_DRAWDOWN_FRAC). For each series, end the window when
# O2 (the fit target) has fallen by FIT_DRAWDOWN_FRAC of its total drawdown from
# the green line, instead of at the per-curve steepest-drop point. This fits the
# same slice of every curve, removing the corr(K, window length) bias. The green
# line (fit_start_time) is unchanged; only the end moves. Auto-adapts to
# temperature (cold curves take longer to reach the same drawdown).
if (isTRUE(USE_DRAWDOWN_WINDOW)) {
  drawdown_end <- o2f %>%
    dplyr::inner_join(
      dplyr::select(fit_windows, T, OTU, Replicate, fit_start_time),
      by = c("T", "OTU", "Replicate")
    ) %>%
    dplyr::filter(is.finite(Oxygen_used), is.finite(Time), Time >= fit_start_time) %>%
    dplyr::arrange(T, OTU, Replicate, Time) %>%
    dplyr::group_by(T, OTU, Replicate) %>%
    dplyr::summarise(
      drawdown_end_time = {
        y  <- Oxygen_used
        tt <- Time
        o2_green <- y[1]                       # O2 at the green line (data is sorted)
        o2_min   <- min(y, na.rm = TRUE)       # lowest O2 reached
        target   <- o2_green - FIT_DRAWDOWN_FRAC * (o2_green - o2_min)
        hit <- which(y <= target)
        if (length(hit) > 0) tt[hit[1]] else tt[length(tt)]
      },
      .groups = "drop"
    )

  fit_windows <- fit_windows %>%
    dplyr::left_join(drawdown_end, by = c("T", "OTU", "Replicate")) %>%
    dplyr::mutate(
      fit_end_time = dplyr::if_else(
        is.finite(drawdown_end_time) & drawdown_end_time > fit_start_time,
        drawdown_end_time,
        fit_end_time   # fall back to steepest-drop/chosen-end if drawdown fails
      )
    ) %>%
    dplyr::select(T, OTU, Replicate, fit_start_time, fit_end_time)

  message(sprintf(
    "Fit-window END standardized by %.0f%% O2 drawdown (USE_DRAWDOWN_WINDOW = TRUE).",
    100 * FIT_DRAWDOWN_FRAC
  ))
}

# Manual fit-window overrides (see MANUAL_FIT_WINDOWS in config.R). For any
# listed curve, the given start/end replaces the automatic trimming; NA keeps the
# automatic value for that side. Applied last, so it wins over the automatic end.
if (exists("MANUAL_FIT_WINDOWS") && nrow(MANUAL_FIT_WINDOWS) > 0) {
  .mfw <- MANUAL_FIT_WINDOWS %>%
    dplyr::mutate(
      T = as.numeric(T),
      OTU = as.integer(OTU),
      Replicate = toupper(as.character(Replicate)),
      manual_start = as.numeric(fit_start),
      manual_end   = as.numeric(fit_end)
    ) %>%
    dplyr::select(T, OTU, Replicate, manual_start, manual_end)

  fit_windows <- fit_windows %>%
    dplyr::left_join(.mfw, by = c("T", "OTU", "Replicate")) %>%
    dplyr::mutate(
      fit_start_time = dplyr::if_else(!is.na(manual_start), manual_start, fit_start_time),
      fit_end_time   = dplyr::if_else(!is.na(manual_end),   manual_end,   fit_end_time)
    ) %>%
    dplyr::select(T, OTU, Replicate, fit_start_time, fit_end_time)

  message(sprintf(
    "Manual fit-window overrides applied to %d curve(s).", nrow(.mfw)
  ))
}

# N0 anchoring (see N0_BACKPROJECT in config.R):
#   TRUE  -> delta = fit-window start (green line); N0 = N_inoc * exp(r * delta).
#   FALSE -> delta = 0; N0 = N_inoc, anchored at the consumption onset. This
#            drops the exp(r * delta) term, so respiration (= K / N0) scales with
#            K and a fixed stock constant, giving a clean TPC comparable to
#            growth. fit_start_time is still kept as a reference column.
group_lookup <- o2f %>%
  dplyr::distinct(T, OTU, Replicate) %>%
  dplyr::left_join(fit_windows, by = c("T", "OTU", "Replicate")) %>%
  dplyr::mutate(
    delta_Ninoc_to_N0_min = if (isTRUE(N0_BACKPROJECT)) fit_start_time else 0,
    N_inoculation_cells_per_L = N_inoculation_cells_per_L
  ) %>%
  dplyr::arrange(T, OTU, Replicate)

# ---- PER-ISOLATE inoculation density (from 06_inoculation.R) ----------------
# Everything went in at the same OD, but OD measures BIOMASS, not cell number -
# so small-celled species (glabrata) were inoculated at a HIGHER cell density
# than large-celled ones (parapsilosis). A single global N_inoc is therefore
# wrong by construction. If tables/otu_inoc.csv exists, use each isolate's own
# value; otherwise fall back to the global constant in config.R.
#
# respiration = K / N0  and  N0 = N_inoc * exp(r * delta), so respiration scales
# as 1/N_inoc. GROWTH is unaffected (it depends on cell carbon, not cell count).
.inoc_csv <- file.path(tables_dir, "otu_inoc.csv")
if (file.exists(.inoc_csv)) {
  .inoc <- tryCatch(readr::read_csv(.inoc_csv, show_col_types = FALSE),
                    error = function(e) NULL)
  if (!is.null(.inoc) && all(c("OTU", "N_inoc_cells_per_L") %in% names(.inoc))) {
    .inoc <- .inoc %>%
      dplyr::transmute(OTU = as.integer(OTU),
                       .n_inoc = suppressWarnings(as.numeric(N_inoc_cells_per_L))) %>%
      dplyr::filter(is.finite(.n_inoc), .n_inoc > 0) %>%
      dplyr::distinct(OTU, .keep_all = TRUE)
    group_lookup <- group_lookup %>%
      dplyr::left_join(.inoc, by = "OTU") %>%
      dplyr::mutate(N_inoculation_cells_per_L =
                      dplyr::coalesce(.n_inoc, N_inoculation_cells_per_L)) %>%
      dplyr::select(-.n_inoc)
    message(sprintf("Per-isolate inoculation densities loaded from %s (%d isolates).",
                    basename(.inoc_csv), nrow(.inoc)))
  }
} else {
  message("No otu_inoc.csv - using the single global N_inoculation_cells_per_L = ",
          format(N_inoculation_cells_per_L, scientific = TRUE),
          " for every isolate. Run 06_inoculation.R to set it per species.")
}

readr::write_csv(group_lookup, group_lookup_csv)

if (any(is.na(group_lookup$delta_Ninoc_to_N0_min))) {
  warning("Missing fit-window start for some groups in trimming metadata. N0 and R will be NA for those groups.")
}


# =============================================================================
# 5) Fit one oxygen series
# =============================================================================

fit_one <- function(df, y_limits = NULL) {
  df0 <- df %>%
    dplyr::arrange(Time) %>%
    dplyr::mutate(Time0 = Time - min(Time, na.rm = TRUE))
  
  y <- df0$Oxygen_used
  
  base_plot <- function(subtitle_txt) {
    p <- ggplot2::ggplot(df0, ggplot2::aes(Time0, Oxygen_used)) +
      ggplot2::geom_point(size = 1.3) +
      ggplot2::labs(
        title = df$series_id[1],
        subtitle = subtitle_txt,
        x = "Time (min, rebased)",
        y = "O2 (mg/L)"
      ) +
      ggplot2::theme_classic(12)
    
    if (!is.null(y_limits) && all(is.finite(y_limits))) {
      p <- p + ggplot2::coord_cartesian(ylim = y_limits)
    }
    
    p
  }
  
  empty_metrics <- tibble::tibble(
    T = df$T[1],
    OTU = df$OTU[1],
    Replicate = df$Replicate[1],
    n = nrow(df0),
    r2 = NA_real_,
    rmse = NA_real_,
    rss = NA_real_,
    aic = NA_real_,
    aicc = NA_real_,
    aicc_lin = NA_real_,
    delta_aicc_curv = NA_real_,
    rt_curv = NA_real_,
    has_curvature = NA,
    r_at_bound = NA,
    K_at_bound = NA,
    fit_valid = FALSE,
    T_end_min = NA_real_,
    keep = FALSE
  )
  
  if (nrow(df0) < 6 || any(!is.finite(y))) {
    return(list(
      coefs = tibble::tibble(
        parameter = c("O2_0", "r", "K"),
        Estimate  = c(NA_real_, NA_real_, NA_real_),
        SE        = NA_real_,
        p_value   = NA_real_
      ),
      metrics = empty_metrics,
      keep = FALSE,
      plot = base_plot("Too few/invalid points")
    ))
  }
  
  O2_0_start <- get_baseline(y)
  
  if (!is.finite(O2_0_start)) {
    O2_0_start <- y[1]
  }
  
  seg_n <- max(6, floor(0.25 * nrow(df0)))
  seg <- head(df0, seg_n)
  
  slope0 <- suppressWarnings(
    median(diff(seg$Oxygen_used) / diff(seg$Time0), na.rm = TRUE)
  )
  
  if (!is.finite(slope0)) {
    slope0 <- -1e-5
  }
  
  r_start <- 1e-3
  
  K_est <- abs(slope0)
  
  if (!is.finite(K_est) || K_est <= 0) {
    K_est <- 1e-6
  }
  
  K_start <- pmin(pmax(K_est, 1e-8), 0.5)
  
  r_lower <- 1e-6
  r_upper <- 0.15
  
  K_lower <- max(1e-10, K_est / 50)
  K_upper <- max(1.0, K_est * 50)
  
  yr <- range(y[is.finite(y)], na.rm = TRUE)
  yspan <- diff(yr)
  
  if (!is.finite(yspan) || yspan <= 0) {
    yspan <- abs(y[1]) + 1
  }
  
  O2_lower <- yr[1] - 0.5 * yspan
  O2_upper <- yr[2] + 0.5 * yspan
  
  fit <- try(
    minpack.lm::nlsLM(
      Oxygen_used ~ resp_model(r, K, Time0, O2_0),
      data = df0,
      start = list(r = r_start, K = K_start, O2_0 = O2_0_start),
      lower = c(r = r_lower, K = K_lower, O2_0 = O2_lower),
      upper = c(r = r_upper, K = K_upper, O2_0 = O2_upper),
      control = minpack.lm::nls.lm.control(
        maxiter = 1500,
        ftol = 1e-12,
        ptol = 1e-12
      )
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(list(
      coefs = tibble::tibble(
        parameter = c("O2_0", "r", "K"),
        Estimate  = c(NA_real_, NA_real_, NA_real_),
        SE        = NA_real_,
        p_value   = NA_real_
      ),
      metrics = empty_metrics,
      keep = FALSE,
      plot = base_plot("Fit failed")
    ))
  }
  
  preds <- as.numeric(predict(fit, df0))
  n <- nrow(df0)
  
  rss <- sum((df0$Oxygen_used - preds)^2, na.rm = TRUE)
  rmse <- sqrt(rss / n)
  
  ss_tot <- sum(
    (df0$Oxygen_used - mean(df0$Oxygen_used, na.rm = TRUE))^2,
    na.rm = TRUE
  )
  
  r2 <- if (is.finite(ss_tot) && ss_tot > 1e-12) {
    1 - rss / ss_tot
  } else {
    NA_real_
  }
  
  k_param <- 3L
  
  aic <- n * log(rss / n) + 2 * k_param
  
  aicc <- if (n > k_param + 1) {
    aic + (2 * k_param * (k_param + 1)) / (n - k_param - 1)
  } else {
    NA_real_
  }
  
  T_end_min <- suppressWarnings(max(df0$Time0, na.rm = TRUE))

  co <- coef(fit)

  # ---- Curvature test: is there any evidence of GROWTH? ---------------------
  # Compare the exponential (curved) model against a STRAIGHT LINE on the same
  # window. If the line is just as good, there is no curvature -> no growth
  # (r unidentifiable), even though respiration (the slope) is real.
  lin_fit <- try(stats::lm(Oxygen_used ~ Time0, data = df0), silent = TRUE)
  aicc_lin <- NA_real_
  if (!inherits(lin_fit, "try-error")) {
    rss_lin <- sum(stats::residuals(lin_fit)^2, na.rm = TRUE)
    k_lin <- 2L
    if (is.finite(rss_lin) && rss_lin > 0 && n > k_lin + 1) {
      aic_lin  <- n * log(rss_lin / n) + 2 * k_lin
      aicc_lin <- aic_lin + (2 * k_lin * (k_lin + 1)) / (n - k_lin - 1)
    }
  }
  delta_aicc_curv <- aicc_lin - aicc          # reported only - NOT used as the flag

  # The real test: dimensionless curvature rt = r * window length.
  r_hat <- suppressWarnings(as.numeric(co[["r"]]))
  rt_curv <- if (is.finite(r_hat) && is.finite(T_end_min)) r_hat * T_end_min else NA_real_
  has_curvature <- isTRUE(is.finite(rt_curv) && rt_curv > CURV_MIN_RT)

  # ---- Boundary check: is this a FAILED fit rather than a bad one? ----------
  # A parameter sitting exactly on its constraint is not an estimate - the
  # optimiser ran out of room and railed against the bound. r = R_UPPER (0.15/min)
  # implies a doubling time of ~5 minutes, which is impossible for Candida.
  # This is NOT a quality filter (we removed those); it is a VALIDITY check: the
  # number is meaningless, so the curve is dropped from the results downstream.
  .tol <- 1e-8
  r_at_bound <- isTRUE(is.finite(co[["r"]]) &&
                         (co[["r"]] >= r_upper * (1 - .tol) ||
                            co[["r"]] <= r_lower * (1 + .tol)))
  K_at_bound <- isTRUE(is.finite(co[["K"]]) &&
                         (co[["K"]] >= K_upper * (1 - .tol) ||
                            co[["K"]] <= K_lower * (1 + .tol)))
  fit_valid <- !(r_at_bound || K_at_bound)

  # No RMSE quality filter: keep every series that produced a usable fit
  # (>= 6 points, needed for the 3-parameter model). Every point is retained.
  # Guard: drop tail-window selector misfires whose fit window is implausibly
  # short (see MIN_FIT_WINDOW_MIN) - these return spurious r/K.
  win_ok <- is.finite(T_end_min) && T_end_min >= MIN_FIT_WINDOW_MIN
  keep <- (n >= 6) && win_ok

  co_sum <- as.data.frame(summary(fit)$parameters) %>%
    tibble::rownames_to_column("parameter") %>%
    tibble::as_tibble()
  
  nm <- names(co_sum)
  
  if ("Std. Error" %in% nm) {
    nm[nm == "Std. Error"] <- "SE"
  }
  
  if ("Pr(>|t|)" %in% nm) {
    nm[nm == "Pr(>|t|)"] <- "p_value"
  }
  
  names(co_sum) <- nm
  
  if (!"SE" %in% names(co_sum)) {
    co_sum$SE <- NA_real_
  }
  
  if (!"p_value" %in% names(co_sum)) {
    co_sum$p_value <- NA_real_
  }
  
  metrics <- tibble::tibble(
    T = df$T[1],
    OTU = df$OTU[1],
    Replicate = df$Replicate[1],
    n = n,
    r2 = r2,
    rmse = rmse,
    rss = rss,
    aic = aic,
    aicc = aicc,
    aicc_lin = aicc_lin,
    delta_aicc_curv = delta_aicc_curv,
    rt_curv = rt_curv,
    has_curvature = has_curvature,
    r_at_bound = r_at_bound,
    K_at_bound = K_at_bound,
    fit_valid = fit_valid,
    T_end_min = T_end_min,
    keep = keep
  )
  
  # Say plainly whether the fitted line is CURVED (growth) or STRAIGHT (no growth).
  curv_txt <- if (!is.finite(rt_curv)) {
    " | curvature: ?"
  } else if (isTRUE(has_curvature)) {
    sprintf(" | CURVED -> growth (r*window=%.2f)", rt_curv)
  } else {
    sprintf(" | STRAIGHT -> NO GROWTH (r*window=%.2f)", rt_curv)
  }

  subtitle_txt <- paste0(
    "R2=", ifelse(is.na(r2), "NA", sprintf("%.3f", r2)),
    " | RMSE=", sprintf("%.3g", rmse),
    " | AIC=", sprintf("%.1f", aic),
    " | fit=", if (isTRUE(FIT_TO_SPLINE)) "spline" else "raw",
    curv_txt,
    if (!keep) " - FILTERED" else ""
  )
  
  # Grey dashed = the straight-line fit. If it lies on top of the red curve, the
  # trace has no curvature -> respiration but no growth.
  lin_preds <- if (!inherits(lin_fit, "try-error")) {
    as.numeric(stats::fitted(lin_fit))
  } else rep(NA_real_, nrow(df0))

  p <- ggplot2::ggplot(df0, ggplot2::aes(Time0, Oxygen_used)) +
    ggplot2::geom_point(size = 1.3) +
    ggplot2::geom_line(ggplot2::aes(y = lin_preds), linewidth = 0.7,
                       colour = "grey50", linetype = "dashed") +
    ggplot2::geom_line(ggplot2::aes(y = preds), linewidth = 0.9, color = "red") +
    ggplot2::labs(
      title = df$series_id[1],
      subtitle = subtitle_txt,
      x = "Time (min, rebased)",
      y = "O2 (mg/L)",
      caption = "red = exponential (growth) model · grey dashed = straight line (no growth)"
    ) +
    ggplot2::theme_classic(12)
  
  if (!is.null(y_limits) && all(is.finite(y_limits))) {
    p <- p + ggplot2::coord_cartesian(ylim = y_limits)
  }
  
  list(
    coefs = co_sum,
    metrics = metrics,
    keep = keep,
    plot = p
  )
}


# =============================================================================
# 6) Run oxygen fits
# =============================================================================

# Restrict each series to its dashed-line fit window (green -> purple) before
# splitting into per-series groups. The full trimmed series in IN_CSV already
# contains this sub-range, so we just filter to it.
o2f_fit <- o2f %>%
  dplyr::left_join(fit_windows, by = c("T", "OTU", "Replicate")) %>%
  dplyr::filter(
    is.finite(fit_start_time), is.finite(fit_end_time),
    Time >= fit_start_time, Time <= fit_end_time
  )

n_series_all <- nrow(dplyr::distinct(o2f, T, OTU, Replicate))
n_series_fit <- nrow(dplyr::distinct(o2f_fit, T, OTU, Replicate))
if (n_series_fit < n_series_all) {
  message(sprintf(
    "Note: %d of %d series had no points inside the dashed-line window and will be skipped.",
    n_series_all - n_series_fit, n_series_all
  ))
}

groups <- o2f_fit %>%
  dplyr::group_by(T, OTU, Replicate) %>%
  dplyr::group_split()

all_coef_rows <- list()
all_metrics <- list()

pdf(pdf_path, width = 6.8, height = 4.6)

for (g in groups) {
  res <- fit_one(
    g,
    y_limits = Y_LIMITS_SERIES
  )

  all_metrics[[length(all_metrics) + 1]] <- res$metrics

  # Plot EVERY curve to per_series_fits.pdf (its subtitle is marked "- FILTERED"
  # when it fails the RMSE / boundary check), so nothing is silently invisible.
  print(res$plot)

  # Collect coefficients for every fit that produced finite estimates. There is
  # no RMSE quality filter, so every usable fit is kept (`keep` is TRUE whenever
  # the series had enough points to fit) and appears in the results and plots.
  co <- res$coefs
  if (!is.null(co) && any(is.finite(co$Estimate))) {
    all_coef_rows[[length(all_coef_rows) + 1]] <-
      tibble::tibble(
        T = g$T[1],
        OTU = g$OTU[1],
        Replicate = g$Replicate[1]
      ) %>%
      dplyr::bind_cols(co) %>%
      dplyr::mutate(T_end_min = res$metrics$T_end_min[1],
                    keep = isTRUE(res$keep),
                    has_curvature = isTRUE(res$metrics$has_curvature[1]),
                    delta_aicc_curv = res$metrics$delta_aicc_curv[1],
                    fit_valid = isTRUE(res$metrics$fit_valid[1]))
  }
}

dev.off()

coef_out <- dplyr::bind_rows(all_coef_rows)
fit_metrics_out <- dplyr::bind_rows(all_metrics)

readr::write_csv(coef_out, coef_csv)
readr::write_csv(fit_metrics_out, fit_metrics_csv)

# Report the curvature (growth vs no-growth) verdict.
if ("has_curvature" %in% names(fit_metrics_out)) {
  .ng <- sum(!(fit_metrics_out$has_curvature %in% TRUE), na.rm = TRUE)
  message(sprintf(
    "Curvature test: %d of %d curves are STRAIGHT lines (respiration but NO growth); r is unidentifiable for those.",
    .ng, nrow(fit_metrics_out)))
  message("  They are dropped from the GROWTH fit in 07/09 but kept for respiration. ",
          "See the 'CURVED / STRAIGHT' verdict on each page of per_series_fits.pdf.")
}

if (nrow(coef_out) == 0) {
  stop("No oxygen model fits produced coefficients. Check the fit diagnostics PDF and the trimmed input.")
}


# =============================================================================
# 7) Wide coefficients
# =============================================================================

coef_wide <- coef_out %>%
  dplyr::select(T, OTU, Replicate, parameter, Estimate, T_end_min, keep,
                has_curvature, delta_aicc_curv, fit_valid) %>%
  tidyr::pivot_wider(names_from = parameter, values_from = Estimate) %>%
  dplyr::group_by(T, OTU, Replicate) %>%
  dplyr::summarise(
    r = dplyr::first(r),
    K = dplyr::first(K),
    O2_0 = dplyr::first(O2_0),
    T_end_min = dplyr::first(T_end_min),
    keep = dplyr::first(keep),
    has_curvature = dplyr::first(has_curvature),
    delta_aicc_curv = dplyr::first(delta_aicc_curv),
    fit_valid = dplyr::first(fit_valid),
    .groups = "drop"
  ) %>%
  dplyr::arrange(T, OTU, Replicate)

# Drop FAILED fits (a parameter railed against its bound => not an estimate).
.n_bad <- sum(!(coef_wide$fit_valid %in% TRUE))
if (.n_bad > 0) {
  message(sprintf(
    "Boundary check: %d of %d curves had a parameter pinned at its bound (failed fit) - REMOVED.",
    .n_bad, nrow(coef_wide)))
  readr::write_csv(coef_wide %>% dplyr::filter(!(fit_valid %in% TRUE)) %>%
                     dplyr::select(T, OTU, Replicate, r, K, O2_0),
                   file.path(tables_dir, "failed_fits_at_bound.csv"))
  coef_wide <- coef_wide %>% dplyr::filter(fit_valid %in% TRUE)
}

readr::write_csv(coef_wide, coef_wide_csv)


# =============================================================================
# 8) Compute sample-specific N0, respiration, carbon units, and CUE
# =============================================================================

# Per-OTU cell size / carbon. 05_cell_sizes.R writes tables/otu_cell_sizes.csv
# (one row per OTU with cell_volume_um3 and cell_carbon_fg). Any OTU not listed
# there falls back to the single global cell size from config.R. Cell carbon
# scales growth (fg C) and CUE; it does NOT change respiration_fgC_h or the
# temperature-response shapes.
otu_size_csv  <- file.path(tables_dir, "otu_cell_sizes.csv")   # from 05_cell_sizes.R
otu_names_csv <- file.path(tables_dir, "otu_names.csv")        # from 01_convert_xlsx.R
.otu_ids <- sort(unique(coef_wide$OTU))
otu_size_lookup <- tibble::tibble(
  OTU             = .otu_ids,
  otu_name        = paste0("OTU", .otu_ids),
  cell_volume_um3 = CELL_VOLUME_UM3,
  cell_carbon_fg  = CELL_CARBON_FG_PER_CELL
)

# (1) Strain display names typed in 01_convert_xlsx.R. Blank -> keep "OTU<n>".
if (file.exists(otu_names_csv)) {
  .nm <- tryCatch(readr::read_csv(otu_names_csv, show_col_types = FALSE),
                  error = function(e) NULL)
  if (!is.null(.nm) && all(c("OTU", "otu_name") %in% names(.nm))) {
    .nm2 <- tibble::tibble(
      OTU  = as.integer(.nm$OTU),
      .newname = trimws(as.character(.nm$otu_name))
    )
    .nm2$.newname <- ifelse(!is.na(.nm2$.newname) & nzchar(.nm2$.newname),
                            .nm2$.newname, NA_character_)
    otu_size_lookup <- otu_size_lookup %>%
      dplyr::left_join(.nm2, by = "OTU") %>%
      dplyr::mutate(otu_name = dplyr::coalesce(.newname, otu_name)) %>%
      dplyr::select(-.newname)
    message("Strain names loaded from ", otu_names_csv, ".")
  }
}

# (2) Per-OTU cell volume / carbon from 05_cell_sizes.R. Missing -> global size.
if (file.exists(otu_size_csv)) {
  .sz <- tryCatch(readr::read_csv(otu_size_csv, show_col_types = FALSE),
                  error = function(e) NULL)
  if (!is.null(.sz) && all(c("OTU", "cell_carbon_fg") %in% names(.sz))) {
    .has_vol <- "cell_volume_um3" %in% names(.sz)
    .sz2 <- tibble::tibble(
      OTU = as.integer(.sz$OTU),
      .cc = suppressWarnings(as.numeric(.sz$cell_carbon_fg)),
      .cv = if (.has_vol) suppressWarnings(as.numeric(.sz$cell_volume_um3)) else NA_real_
    )
    otu_size_lookup <- otu_size_lookup %>%
      dplyr::left_join(.sz2, by = "OTU") %>%
      dplyr::mutate(
        cell_carbon_fg  = dplyr::coalesce(.cc, cell_carbon_fg),
        cell_volume_um3 = dplyr::coalesce(.cv, cell_volume_um3)
      ) %>%
      dplyr::select(OTU, otu_name, cell_volume_um3, cell_carbon_fg)
    message("Per-OTU cell sizes loaded from ", otu_size_csv, ".")
  }
} else {
  message("No per-OTU cell sizes file; using the global cell size from config.R for every OTU.")
}
message("OTU lookup (name + cell size):")
print(otu_size_lookup)

# Name lookup used for plot facet strips and colour legends (falls back to
# "OTU<n>"). Internal grouping/keys still use the numeric OTU underneath.
otu_name_map <- stats::setNames(otu_size_lookup$otu_name,
                                as.character(otu_size_lookup$OTU))
otu_name_of  <- function(x) {
  nm <- unname(otu_name_map[as.character(x)])
  ifelse(is.na(nm), paste0("OTU", x), nm)
}

# Attach each OTU's cell size to a per-OTU table, so size travels with the TPC
# coefficient CSVs written later.
attach_otu_size <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  dplyr::left_join(df, otu_size_lookup, by = "OTU")
}

results <- coef_wide %>%
  dplyr::left_join(group_lookup, by = c("T", "OTU", "Replicate")) %>%
  dplyr::left_join(otu_size_lookup, by = "OTU") %>%
  dplyr::mutate(
    N0_cells_per_L = dplyr::if_else(
      is.finite(N_inoculation_cells_per_L) & N_inoculation_cells_per_L > 0 &
        is.finite(delta_Ninoc_to_N0_min) & delta_Ninoc_to_N0_min >= 0 &
        is.finite(r) & r > 0,
      N_inoculation_cells_per_L * exp(r * delta_Ninoc_to_N0_min),
      NA_real_
    ),
    C_tot_O2_mg_per_L = dplyr::if_else(
      is.finite(K) & is.finite(r) & r > 0 &
        is.finite(T_end_min) & T_end_min > 0,
      (K / r) * (exp(r * T_end_min) - 1),
      NA_real_
    ),
    biomass_integral_cells_min_per_L = dplyr::if_else(
      is.finite(N0_cells_per_L) & N0_cells_per_L > 0 &
        is.finite(r) & r > 0 &
        is.finite(T_end_min) & T_end_min > 0,
      N0_cells_per_L * (exp(r * T_end_min) - 1) / r,
      NA_real_
    ),
    R_O2_mg_cell_min = dplyr::if_else(
      is.finite(C_tot_O2_mg_per_L) & C_tot_O2_mg_per_L > 0 &
        is.finite(biomass_integral_cells_min_per_L) & biomass_integral_cells_min_per_L > 0,
      C_tot_O2_mg_per_L / biomass_integral_cells_min_per_L,
      NA_real_
    ),
    # cell_volume_um3 and cell_carbon_fg now come from otu_size_lookup (per isolate).
    growth_fgC_h = dplyr::if_else(
      is.finite(r) & r > 0,
      r * cell_carbon_fg * MIN_TO_H,
      NA_real_
    ),
    respiration_fgC_h = dplyr::if_else(
      is.finite(R_O2_mg_cell_min) & R_O2_mg_cell_min > 0,
      # mg O2/cell/min -> fg C/cell/h: mg->fg, O2 mass -> C mass (molar ratio
      # x RQ), min->h. The O2_TO_C_MASS factor was previously missing, so this
      # was reporting fg O2 (mislabeled as carbon, ~2.7x too high).
      R_O2_mg_cell_min * MG_TO_FG * O2_TO_C_MASS * RESPIRATORY_QUOTIENT * MIN_TO_H,
      NA_real_
    ),
    growth_C_per_C_h = dplyr::if_else(
      is.finite(growth_fgC_h) & growth_fgC_h > 0 &
        is.finite(cell_carbon_fg) & cell_carbon_fg > 0,
      growth_fgC_h / cell_carbon_fg,
      NA_real_
    ),
    respiration_C_per_C_h = dplyr::if_else(
      is.finite(respiration_fgC_h) & respiration_fgC_h > 0 &
        is.finite(cell_carbon_fg) & cell_carbon_fg > 0,
      respiration_fgC_h / cell_carbon_fg,
      NA_real_
    ),
    CUE = dplyr::if_else(
      is.finite(growth_fgC_h) & growth_fgC_h > 0 &
        is.finite(respiration_fgC_h) & respiration_fgC_h > 0,
      growth_fgC_h / (growth_fgC_h + respiration_fgC_h),
      NA_real_
    ),
    resp_over_growth = dplyr::if_else(
      is.finite(respiration_fgC_h) & respiration_fgC_h > 0 &
        is.finite(growth_fgC_h) & growth_fgC_h > 0,
      respiration_fgC_h / growth_fgC_h,
      NA_real_
    )
  ) %>%
  dplyr::arrange(T, OTU, Replicate)

message(sprintf(
  "Per-series fits: %d series with a usable fit (no RMSE filtering; every point is retained).",
  nrow(results)
))


# =============================================================================
# 9) Exclude specific outlier points
# =============================================================================

EXCLUDE_POINTS <- tibble::tibble(
  T         = numeric(),
  OTU     = integer(),
  Replicate = character()
)

results <- results %>%
  dplyr::anti_join(EXCLUDE_POINTS, by = c("T", "OTU", "Replicate"))

message(sprintf(
  "Exclusion filter applied: %d point(s) removed. %d rows retained.",
  nrow(EXCLUDE_POINTS),
  nrow(results)
))

# Drop whole (temperature, otu) groups listed in EXCLUDE_T_OTU (config.R).
if (exists("EXCLUDE_T_OTU") && !is.null(EXCLUDE_T_OTU) && nrow(EXCLUDE_T_OTU) > 0) {
  excl_tc <- EXCLUDE_T_OTU %>%
    dplyr::mutate(T = as.numeric(T), OTU = as.integer(OTU)) %>%
    dplyr::distinct()
  before_tc <- nrow(results)
  results <- results %>% dplyr::anti_join(excl_tc, by = c("T", "OTU"))
  message(sprintf(
    "EXCLUDE_T_OTU: removed %d row(s) for %s. %d rows retained.",
    before_tc - nrow(results),
    paste(sprintf("(T=%g, OTU=%g)", excl_tc$T, excl_tc$OTU), collapse = ", "),
    nrow(results)
  ))
}

# =============================================================================
# 9c) THE ONE FILTER: remove points that should not appear anywhere
# =============================================================================
# Removed from BOTH growth and respiration, and therefore from every table and
# every plot downstream (07, 09, 10 all read this file):
#   (1) curves YOU excluded in 04_trim_selector.R  (plot_exclude_points.csv)
#   (2) ZERO-GROWTH curves - straight O2 line, no curvature, r not measurable
#   (3) failed fits (a parameter pinned at its bound) - already dropped above
# Nothing else is filtered. Everything that survives is real data.
# =============================================================================

.n_before <- nrow(results)
.removed <- list()

# (1) your 04 exclusions
if (exists("PLOT_EXCLUDE_POINTS") && is.data.frame(PLOT_EXCLUDE_POINTS) &&
    nrow(PLOT_EXCLUDE_POINTS) > 0 &&
    all(c("T", "OTU", "Replicate") %in% names(PLOT_EXCLUDE_POINTS))) {
  .ex <- PLOT_EXCLUDE_POINTS %>%
    dplyr::mutate(T = as.numeric(T), OTU = as.integer(OTU),
                  Replicate = toupper(as.character(Replicate))) %>%
    dplyr::distinct(T, OTU, Replicate)
  .removed$excluded_in_04 <- results %>% dplyr::semi_join(.ex, by = c("T", "OTU", "Replicate")) %>%
    dplyr::transmute(T, OTU, Replicate, reason = "excluded_in_04")
  results <- results %>% dplyr::anti_join(.ex, by = c("T", "OTU", "Replicate"))
}

# NOTE: zero-growth curves are NOT discarded automatically here. YOU decide, in
# 04_trim_selector.R ("Exclude ALL zero-growth curves"), and 06 simply obeys
# whatever plot_exclude_points.csv says. The has_curvature / rt_curv / O2-drawdown
# flags are still computed and written out - they are what 04 uses to FIND the
# zero-growth curves for you - but nothing acts on them behind your back.

.rm_tbl <- dplyr::bind_rows(.removed)
if (nrow(.rm_tbl) > 0) {
  readr::write_csv(.rm_tbl %>% dplyr::arrange(OTU, T, Replicate),
                   file.path(tables_dir, "removed_points.csv"))
}
message(sprintf(
  "Removed %d of %d curves - all of them excluded by YOU in 04. %d curves remain.",
  .n_before - nrow(results), .n_before, nrow(results)))
message("  -> gone from BOTH growth and respiration, and from every plot. See tables/removed_points.csv")
message("  -> nothing else is discarded automatically. Use 04 to mark what you want gone.")

readr::write_csv(results, derived_csv)


# =============================================================================
# 9b) Diagnostics for 43°C and 44°C
# =============================================================================
# This block helps check whether unexpected 43°C vs 44°C differences come from
# r, K, O2_0, fit duration, or N0 calculations.
# =============================================================================

# Generalised from the legacy 43/44 check: inspect the two HOTTEST temperatures
# present in the data (per otu), where thermal stress is most likely to make
# r, K, O2_0, fit duration, or N0 behave unexpectedly.
hot_temps <- sort(unique(results$T), decreasing = TRUE)
hot_temps <- utils::head(hot_temps, 2)

diagnostic_hot <- results %>%
  dplyr::filter(T %in% hot_temps) %>%
  dplyr::select(
    T, OTU, Replicate,
    r, K, O2_0, T_end_min,
    delta_Ninoc_to_N0_min,
    N0_cells_per_L,
    growth_fgC_h,
    respiration_fgC_h,
    CUE,
    resp_over_growth
  ) %>%
  dplyr::arrange(T, OTU, Replicate)

readr::write_csv(
  diagnostic_hot,
  file.path(tables_dir, "diagnostic_hot_temps_growth_check.csv")
)

summary_hot <- diagnostic_hot %>%
  dplyr::group_by(T, OTU) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_r = mean(r, na.rm = TRUE),
    median_r = median(r, na.rm = TRUE),
    mean_growth_fgC_h = mean(growth_fgC_h, na.rm = TRUE),
    median_growth_fgC_h = median(growth_fgC_h, na.rm = TRUE),
    mean_T_end_min = mean(T_end_min, na.rm = TRUE),
    median_T_end_min = median(T_end_min, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  summary_hot,
  file.path(tables_dir, "summary_hot_temps_growth_check.csv")
)

if (nrow(diagnostic_hot) > 0) {
  p_check_hot_r <- ggplot2::ggplot(
    diagnostic_hot %>% dplyr::filter(is.finite(r)),
    ggplot2::aes(factor(T), r, colour = factor(OTU))
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    ggplot2::geom_jitter(
      position = ggplot2::position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
      size = 2
    ) +
    ggplot2::labs(
      title = sprintf("Fitted oxygen-model r at %s degC", paste(hot_temps, collapse = " and ")),
      x = "Temperature (°C)",
      y = "r fitted from oxygen model",
      colour = "Isolate"
    ) +
    ggplot2::theme_classic(12)

  ggplot2::ggsave(
    file.path(figures_dir, "check_hot_temps_fitted_r.png"),
    p_check_hot_r,
    width = 6,
    height = 4,
    dpi = 300
  )
}


# =============================================================================
# 10) Descriptive plots (temperature x isolate)
# =============================================================================
# Every figure distinguishes isolates by colour and gives each isolate its own
# facet panel so the thermal responses can be compared side by side.
# =============================================================================

results_plot <- results %>%
  dplyr::mutate(OTU = as.integer(OTU))

if (length(EXCLUDE_TEMPS_PLOT) > 0) {
  n_before <- nrow(results_plot)
  results_plot <- dplyr::filter(results_plot, !(T %in% EXCLUDE_TEMPS_PLOT))
  message(sprintf(
    "Plot/SS temperature filter: dropped %s degC -> %d of %d rows kept.",
    paste(EXCLUDE_TEMPS_PLOT, collapse = ", "),
    nrow(results_plot),
    n_before
  ))
}

# Replicate counts per temperature AND otu.
replication_summary <- results_plot %>%
  dplyr::count(T, OTU, name = "n_replicates") %>%
  dplyr::arrange(T, OTU)

readr::write_csv(replication_summary, replication_summary_csv)

# results_plot_full keeps the full data and feeds the Bayesian frames in
# section 11 (04 unaffected). results_plot is the plot/SS subset, with only the
# manual point exclusions below applied (Bayesian keeps those points).
results_plot_full <- results_plot

# Manual point exclusions for the plots + SS fits only (Bayesian data unaffected).
if (exists("PLOT_EXCLUDE_POINTS") && nrow(PLOT_EXCLUDE_POINTS) > 0) {
  .excl <- PLOT_EXCLUDE_POINTS %>%
    dplyr::mutate(
      T = as.numeric(T),
      OTU = as.integer(OTU),
      Replicate = toupper(as.character(Replicate))
    )
  .n_pre_excl <- nrow(results_plot)
  results_plot <- results_plot %>%
    dplyr::anti_join(.excl, by = c("T", "OTU", "Replicate"))
  message(sprintf(
    "Plot/SS point exclusions: dropped %d of %d rows.",
    .n_pre_excl - nrow(results_plot), .n_pre_excl
  ))
}

# Shared isolate colour scale + facet helper used by every descriptive figure.
# Legend and facet strips show each isolate's display name (otu_name_map, from
# 01_convert_xlsx.R); grouping still uses the numeric OTU code underneath.
#
# There are up to 18 isolates, which exceeds the 8-colour cap of the old
# RColorBrewer "Dark2" palette (that would error). Build a manual colour scale
# sized to the isolates actually present, using grDevices::hcl.colors(), which
# generates any number of well-separated hues.
.otu_all <- sort(unique(otu_size_lookup$OTU))
.otu_pal <- stats::setNames(
  grDevices::hcl.colors(max(length(.otu_all), 2L), palette = "Dark 3"),
  as.character(.otu_all)
)
otu_colour_scale <- ggplot2::scale_colour_manual(values = .otu_pal, name = "Isolate",
                                                 labels = otu_name_of)
otu_fill_scale   <- ggplot2::scale_fill_manual(values = .otu_pal, name = "Isolate",
                                               labels = otu_name_of)
otu_facet        <- ggplot2::facet_wrap(
  ~ OTU, labeller = ggplot2::as_labeller(otu_name_of))


box_by_otu <- function(df, value_col, title_txt, ylab) {
  v <- rlang::sym(value_col)
  d <- df %>% dplyr::filter(is.finite(!!v))

  ggplot2::ggplot(d, ggplot2::aes(factor(T), !!v, colour = factor(OTU))) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    ggplot2::geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
    otu_colour_scale +
    otu_facet +
    ggplot2::labs(title = title_txt, x = "Temperature (°C)", y = ylab) +
    ggplot2::theme_classic(12)
}

scatter_by_otu <- function(df, value_col, title_txt, ylab, positive_only = TRUE) {
  v <- rlang::sym(value_col)
  d <- df %>% dplyr::filter(is.finite(!!v))
  if (positive_only) d <- d %>% dplyr::filter(!!v > 0)

  ggplot2::ggplot(d, ggplot2::aes(T, !!v, colour = factor(OTU))) +
    ggplot2::geom_point(size = 2.1, alpha = 0.9) +
    otu_colour_scale +
    otu_facet +
    ggplot2::labs(title = title_txt, x = "Temperature (°C)", y = ylab) +
    ggplot2::theme_classic(12)
}


p_box_growth <- box_by_otu(
  results_plot, "growth_fgC_h",
  "growth in carbon units (by isolate)",
  expression("Growth (fg C h"^{-1}*")")
)

p_box_resp <- box_by_otu(
  results_plot, "respiration_fgC_h",
  "respiration in carbon units (by isolate)",
  expression("Respiration (fg C h"^{-1}*")")
)


ratio_dat <- results_plot %>%
  dplyr::filter(
    is.finite(growth_fgC_h), growth_fgC_h > 0,
    is.finite(respiration_fgC_h), respiration_fgC_h > 0
  ) %>%
  dplyr::mutate(
    log_resp_over_growth = log(resp_over_growth),
    TK = as.numeric(T) + 273.15
  )


p_rg_t <- scatter_by_otu(
  ratio_dat, "resp_over_growth",
  "Respiration / Growth vs Temperature (by isolate)",
  "Respiration / Growth"
)

p_growth_vs_T <- scatter_by_otu(
  results_plot, "growth_fgC_h",
  "growth vs Temperature (by isolate)",
  expression("Growth (fg C h"^{-1}*")")
)

p_resp_vs_T <- scatter_by_otu(
  results_plot, "respiration_fgC_h",
  "respiration vs Temperature (by isolate)",
  expression("Respiration (fg C h"^{-1}*")")
)

# Total (whole-sample) respiration = total O2 consumed over each curve's
# green->purple window (C_tot_O2_mg_per_L). Same facet-by-otu style as the
# per-cell figures, but NOT divided by cell count.
p_box_total_resp <- box_by_otu(
  results_plot, "C_tot_O2_mg_per_L",
  "total respiration (by isolate)",
  expression("Total O"[2]*" consumed (mg L"^{-1}*")")
)

p_total_resp_vs_T <- scatter_by_otu(
  results_plot, "C_tot_O2_mg_per_L",
  "total respiration vs Temperature (by isolate)",
  expression("Total O"[2]*" consumed (mg L"^{-1}*")")
)

p_box_growth_biomass <- box_by_otu(
  results_plot %>% dplyr::filter(growth_C_per_C_h > 0), "growth_C_per_C_h",
  "Biomass-corrected growth (by isolate)",
  expression("Growth (C C"^{-1}*" h"^{-1}*")")
)

p_box_resp_biomass <- box_by_otu(
  results_plot %>% dplyr::filter(respiration_C_per_C_h > 0), "respiration_C_per_C_h",
  "Biomass-corrected respiration (by isolate)",
  expression("Respiration (C C"^{-1}*" h"^{-1}*")")
)

p_growth_biomass_vs_T <- scatter_by_otu(
  results_plot, "growth_C_per_C_h",
  "Biomass-corrected growth vs Temperature (by isolate)",
  expression("Growth (C C"^{-1}*" h"^{-1}*")")
)

p_resp_biomass_vs_T <- scatter_by_otu(
  results_plot, "respiration_C_per_C_h",
  "Biomass-corrected respiration vs Temperature (by isolate)",
  expression("Respiration (C C"^{-1}*" h"^{-1}*")")
)


# =============================================================================
# 10b) Sharpe-Schoolfield fits for growth and respiration (frequentist)
# =============================================================================
# Frequentist nlsLM fit of log(rate) ~ pred_ss_log(TK, lnB0, E, Eh, Th), one
# thermal performance curve per OTU. No Bayesian / MCMC. pred_ss_log is defined
# in config.R. Outputs a fitted-curve overlay PNG per rate plus coef CSVs.
# =============================================================================

fit_ss_nls <- function(df, value_col,
                       lnB0_start = NULL, E_start = 0.65,
                       Eh_start = 3.5, Th_start = 315) {
  d <- df %>%
    dplyr::transmute(
      Replicate = as.character(Replicate),
      TK = as.numeric(T) + 273.15,
      y  = as.numeric(.data[[value_col]])
    ) %>%
    dplyr::filter(is.finite(TK), is.finite(y), y > 0)

  if (nrow(d) < 5) return(list(ok = FALSE, reason = "Too few finite positive points"))
  if (is.null(lnB0_start) || !is.finite(lnB0_start)) lnB0_start <- log(stats::median(d$y))

  lower <- c(lnB0 = -50, E = 0.01, Eh = 0.05, Th = 285)
  # Cap Th (temperature of half high-T deactivation) just above the hottest
  # MEASURED temperature. Without this cap the optimiser can push Th to a far-off
  # boundary (e.g. 340 K / 67 degC) and return a near-monotonic curve that
  # ignores the high-T decline - the failure that made OTU 2 look terrible.
  Tmax  <- max(d$TK, na.rm = TRUE)
  upper <- c(lnB0 = 50, E = 5, Eh = 30, Th = min(340, Tmax + 8))

  fit_with <- function(s) {
    th0 <- min(max(s$Th, lower["Th"] + 1), upper["Th"] - 1)
    try(minpack.lm::nlsLM(
      log(y) ~ pred_ss_log(TK, lnB0, E, Eh, Th), data = d,
      start = list(lnB0 = lnB0_start, E = s$E, Eh = s$Eh, Th = th0),
      lower = lower, upper = upper,
      control = minpack.lm::nls.lm.control(maxiter = 1000, ftol = 1e-10, ptol = 1e-10)
    ), silent = TRUE)
  }

  # Try several start sets ALWAYS and keep the lowest-RSS fit. A single default
  # start can "converge" to a poor boundary solution without erroring, so we
  # never rely on just one.
  start_sets <- list(
    list(E = E_start, Eh = Eh_start, Th = Th_start),
    list(E = 0.6, Eh = 0.5, Th = Tmax - 2),
    list(E = 0.8, Eh = 3,   Th = Tmax - 5),
    list(E = 1.0, Eh = 6,   Th = Tmax - 8),
    list(E = 0.4, Eh = 10,  Th = max(Tmax - 12, 300)),
    list(E = 1.5, Eh = 15,  Th = Tmax)
  )
  best <- NULL; best_rss <- Inf
  for (s in start_sets) {
    ft <- fit_with(s)
    if (!inherits(ft, "try-error")) {
      rss <- sum(stats::residuals(ft)^2, na.rm = TRUE)
      if (is.finite(rss) && rss < best_rss) { best <- ft; best_rss <- rss }
    }
  }
  if (is.null(best)) return(list(ok = FALSE, reason = "SS failed to converge (all starts)"))
  fit <- best

  co <- as.list(coef(fit))
  # Fit quality on the (log) scale the model is fit on.
  res_log <- as.numeric(stats::residuals(fit))
  ss_res  <- sum(res_log^2, na.rm = TRUE)
  ss_tot  <- sum((log(d$y) - mean(log(d$y)))^2, na.rm = TRUE)
  r2      <- if (is.finite(ss_tot) && ss_tot > 1e-12) 1 - ss_res / ss_tot else NA_real_
  rmse    <- sqrt(mean(res_log^2, na.rm = TRUE))

  T_grid  <- seq(min(d$TK) - 273.15, max(d$TK) - 273.15, length.out = 200)
  TK_grid <- T_grid + 273.15
  list(
    ok = TRUE, reason = NA_character_,
    coefs = tibble::tibble(parameter = c("lnB0", "E", "Eh", "Th"),
                           Estimate  = c(co$lnB0, co$E, co$Eh, co$Th),
                           r_squared = r2,
                           rmse_log  = rmse,
                           n         = nrow(d)),
    preds = tibble::tibble(T = T_grid,
                           y_fit = exp(pred_ss_log(TK_grid, co$lnB0, co$E, co$Eh, co$Th))),
    data  = d
  )
}

fit_ss_by_otu <- function(df, value_col) {
  pts <- list(); prd <- list(); cof <- list()
  for (cl in sort(unique(df$OTU))) {
    out <- fit_ss_nls(df %>% dplyr::filter(OTU == cl), value_col)
    if (!isTRUE(out$ok)) {
      message(sprintf("SS fit (%s) failed for OTU %s: %s", value_col, cl, out$reason)); next
    }
    pts[[length(pts) + 1]] <- out$data  %>% dplyr::mutate(OTU = cl)
    prd[[length(prd) + 1]] <- out$preds %>% dplyr::mutate(OTU = cl)
    cof[[length(cof) + 1]] <- out$coefs %>% dplyr::mutate(OTU = cl)
  }
  list(points = dplyr::bind_rows(pts), preds = dplyr::bind_rows(prd),
       coefs = dplyr::bind_rows(cof))
}

ss_growth_out <- fit_ss_by_otu(results_plot, "growth_fgC_h")
ss_resp_out   <- fit_ss_by_otu(results_plot, "respiration_fgC_h")
ss_K_out      <- fit_ss_by_otu(results_plot, "K")

ss_plot_by_otu <- function(ss_out, y_lab, title_txt) {
  if (nrow(ss_out$points) == 0) {
    message(title_txt, ": no OTU-level SS fits succeeded - plotting nothing.")
    return(ggplot2::ggplot() +
             ggplot2::labs(title = paste0(title_txt, " (no SS fit)"),
                           x = "Temperature (°C)", y = y_lab) +
             ggplot2::theme_classic(12))
  }
  label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(
      data = ss_out$points,
      ggplot2::aes(TK - 273.15, y, label = Replicate, colour = factor(OTU)),
      size = 2.5, show.legend = FALSE, max.overlaps = Inf,
      min.segment.length = 0, segment.size = 0.2,
      box.padding = 0.3, point.padding = 0.2, seed = 1)
  } else {
    ggplot2::geom_text(
      data = ss_out$points,
      ggplot2::aes(TK - 273.15, y, label = Replicate, colour = factor(OTU)),
      size = 2.5, vjust = -0.6, show.legend = FALSE, check_overlap = TRUE)
  }
  ggplot2::ggplot() +
    ggplot2::geom_point(data = ss_out$points,
      ggplot2::aes(TK - 273.15, y, colour = factor(OTU)), size = 2.1, alpha = 0.85) +
    label_layer +
    ggplot2::geom_line(data = ss_out$preds,
      ggplot2::aes(T, y_fit, colour = factor(OTU)), linewidth = 1.0) +
    otu_colour_scale + otu_facet +
    ggplot2::labs(title = title_txt, x = "Temperature (°C)", y = y_lab) +
    ggplot2::theme_classic(12)
}

p_ss_growth <- ss_plot_by_otu(ss_growth_out,
  expression("Growth (fg C h"^{-1}*")"), "Sharpe-Schoolfield fit - growth (per isolate)")
p_ss_resp   <- ss_plot_by_otu(ss_resp_out,
  expression("Respiration (fg C h"^{-1}*")"), "Sharpe-Schoolfield fit - respiration (per isolate)")

K_lab    <- expression("O"[2]*" rate K (mg L"^{-1}*" min"^{-1}*")")
p_box_K  <- box_by_otu(results_plot, "K", "O2 consumption rate K (by isolate)", K_lab)
p_K_vs_T <- scatter_by_otu(results_plot, "K", "O2 consumption rate K vs Temperature (by isolate)", K_lab)
p_ss_K   <- ss_plot_by_otu(ss_K_out, K_lab, "Sharpe-Schoolfield fit - K (per isolate)")

ss_growth_png <- file.path(figures_dir, "sharpe_schoolfield_growth_fgC_h.png")
ss_resp_png   <- file.path(figures_dir, "sharpe_schoolfield_respiration_fgC_h.png")
ss_K_png      <- file.path(figures_dir, "sharpe_schoolfield_K_O2_rate.png")
ss_growth_csv <- file.path(tables_dir, "sharpe_schoolfield_growth_fgC_h_coefs.csv")
ss_resp_csv   <- file.path(tables_dir, "sharpe_schoolfield_respiration_fgC_h_coefs.csv")
ss_K_csv      <- file.path(tables_dir, "sharpe_schoolfield_K_O2_rate_coefs.csv")

if (nrow(ss_growth_out$coefs) > 0) readr::write_csv(attach_otu_size(ss_growth_out$coefs), ss_growth_csv)
if (nrow(ss_resp_out$coefs)   > 0) readr::write_csv(attach_otu_size(ss_resp_out$coefs),   ss_resp_csv)
if (nrow(ss_K_out$coefs)      > 0) readr::write_csv(attach_otu_size(ss_K_out$coefs),      ss_K_csv)


# =============================================================================
# 10d) Boltzmann-Arrhenius fits for growth and respiration (frequentist)
# =============================================================================
# log(rate) is linear in the Boltzmann term, so a simple lm() per OTU gives the
# activation energy E (eV) and the normalization alpha. No Bayesian / MCMC.
# pred_arr_log, k_B and T_ref are defined in config.R.
# =============================================================================

fit_arr_lm <- function(df, value_col) {
  d <- df %>%
    dplyr::transmute(
      Replicate = as.character(Replicate),
      TK = as.numeric(T) + 273.15,
      y  = as.numeric(.data[[value_col]])
    ) %>%
    dplyr::filter(is.finite(TK), is.finite(y), y > 0) %>%
    dplyr::mutate(boltz = (1 / (k_B * T_ref)) - (1 / (k_B * TK)), ln_y = log(y))

  if (nrow(d) < 3 || length(unique(d$TK)) < 2) {
    return(list(ok = FALSE, reason = "Too few temperatures"))
  }
  fit <- try(stats::lm(ln_y ~ boltz, data = d), silent = TRUE)
  if (inherits(fit, "try-error")) return(list(ok = FALSE, reason = "lm failed"))
  sm    <- summary(fit)
  alpha <- unname(coef(fit)[1]); E <- unname(coef(fit)[2])
  se_E  <- tryCatch(sm$coefficients["boltz", "Std. Error"], error = function(e) NA_real_)
  T_grid  <- seq(min(d$TK) - 273.15, max(d$TK) - 273.15, length.out = 200)
  list(
    ok = TRUE, reason = NA_character_,
    coefs = tibble::tibble(alpha = alpha, E_eV = E, se_E_eV = se_E,
                           r_squared = sm$r.squared, n = nrow(d)),
    preds = tibble::tibble(T = T_grid, y_fit = exp(pred_arr_log(T_grid + 273.15, alpha, E))),
    data  = d %>% dplyr::mutate(T = TK - 273.15)
  )
}

fit_arr_by_otu <- function(df, value_col) {
  pts <- list(); prd <- list(); cof <- list()
  for (cl in sort(unique(df$OTU))) {
    out <- fit_arr_lm(df %>% dplyr::filter(OTU == cl), value_col)
    if (!isTRUE(out$ok)) {
      message(sprintf("Arrhenius fit (%s) failed for OTU %s: %s", value_col, cl, out$reason)); next
    }
    pts[[length(pts) + 1]] <- out$data  %>% dplyr::mutate(OTU = cl)
    prd[[length(prd) + 1]] <- out$preds %>% dplyr::mutate(OTU = cl)
    cof[[length(cof) + 1]] <- out$coefs %>% dplyr::mutate(OTU = cl)
  }
  list(points = dplyr::bind_rows(pts), preds = dplyr::bind_rows(prd),
       coefs = dplyr::bind_rows(cof))
}

arr_growth_out <- fit_arr_by_otu(results_plot, "growth_fgC_h")
arr_resp_out   <- fit_arr_by_otu(results_plot, "respiration_fgC_h")

arr_plot_by_otu <- function(arr_out, y_lab, title_txt) {
  if (nrow(arr_out$points) == 0) {
    message(title_txt, ": no OTU-level Arrhenius fits succeeded - plotting nothing.")
    return(ggplot2::ggplot() +
             ggplot2::labs(title = paste0(title_txt, " (no fit)"),
                           x = "Temperature (°C)", y = y_lab) +
             ggplot2::theme_classic(12))
  }
  ggplot2::ggplot() +
    ggplot2::geom_point(data = arr_out$points,
      ggplot2::aes(T, y, colour = factor(OTU)), size = 2.1, alpha = 0.85) +
    ggplot2::geom_line(data = arr_out$preds,
      ggplot2::aes(T, y_fit, colour = factor(OTU)), linewidth = 1.0) +
    otu_colour_scale + otu_facet +
    ggplot2::scale_y_log10() +
    ggplot2::labs(title = title_txt, x = "Temperature (°C)", y = y_lab) +
    ggplot2::theme_classic(12)
}

p_arr_growth <- arr_plot_by_otu(arr_growth_out,
  expression("Growth (fg C h"^{-1}*")"), "Arrhenius fit - growth (per isolate)")
p_arr_resp   <- arr_plot_by_otu(arr_resp_out,
  expression("Respiration (fg C h"^{-1}*")"), "Arrhenius fit - respiration (per isolate)")

arr_growth_png <- file.path(figures_dir, "arrhenius_growth_fgC_h.png")
arr_resp_png   <- file.path(figures_dir, "arrhenius_respiration_fgC_h.png")
arr_growth_csv <- file.path(tables_dir, "arrhenius_growth_fgC_h_coefs.csv")
arr_resp_csv   <- file.path(tables_dir, "arrhenius_respiration_fgC_h_coefs.csv")

if (nrow(arr_growth_out$coefs) > 0) readr::write_csv(attach_otu_size(arr_growth_out$coefs), arr_growth_csv)
if (nrow(arr_resp_out$coefs)   > 0) readr::write_csv(attach_otu_size(arr_resp_out$coefs),   arr_resp_csv)


# =============================================================================
# 10e) Combined activation-energy (E) summary across OTUs
# =============================================================================
# Side-by-side Boltzmann-Arrhenius activation energy E (eV) for growth vs
# respiration, one row per OTU. Written to tables/activation_energy_summary.csv,
# with a forest-style E +/- SE figure.
# =============================================================================

arr_summary_prep <- function(out, resp) {
  if (nrow(out$coefs) == 0) return(NULL)
  out$coefs %>%
    dplyr::transmute(OTU, response = resp, E_eV, se_E_eV, r_squared, n)
}

arr_e_long <- dplyr::bind_rows(
  arr_summary_prep(arr_growth_out, "growth"),
  arr_summary_prep(arr_resp_out,   "respiration")
)

activation_energy_csv <- file.path(tables_dir, "activation_energy_summary.csv")
activation_energy_png <- file.path(figures_dir, "activation_energy_by_otu.png")
p_activation_E <- NULL

if (!is.null(arr_e_long) && nrow(arr_e_long) > 0) {
  activation_energy_summary <- arr_e_long %>%
    tidyr::pivot_wider(
      id_cols     = OTU,
      names_from  = response,
      values_from = c(E_eV, se_E_eV, r_squared, n),
      names_glue  = "{response}_{.value}"
    ) %>%
    dplyr::arrange(OTU) %>%
    attach_otu_size()

  readr::write_csv(activation_energy_summary, activation_energy_csv)
  message("Wrote activation-energy summary: ", activation_energy_csv)

  p_activation_E <- ggplot2::ggplot(
      arr_e_long,
      ggplot2::aes(x = factor(OTU), y = E_eV, colour = response)
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.4), size = 2.6) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = E_eV - se_E_eV, ymax = E_eV + se_E_eV),
      position = ggplot2::position_dodge(width = 0.4), width = 0.2, linewidth = 0.7
    ) +
    ggplot2::labs(title = "Arrhenius activation energy by isolate",
                  x = "Isolate", y = "Activation energy E (eV)", colour = NULL) +
    ggplot2::theme_classic(12)
} else {
  message("No Arrhenius coefficients available - skipping activation-energy summary.")
}


# =============================================================================
# 10f) CUE plots + per-OTU temperature fit
# =============================================================================
# CUE is a bounded ratio (0-1), not a rate, so Arrhenius / Sharpe-Schoolfield do
# not apply. Instead fit a per-OTU QUADRATIC in temperature (peak-capable: it
# exposes a thermal optimum when the response is humped, or a monotonic trend
# otherwise). Outputs a boxplot, a CUE-vs-T scatter with the fitted curve, and a
# coefficient CSV (with the optimum temperature). To use a smoother instead,
# swap the geom_line below for ggplot2::geom_smooth(method = "gam"/"loess").
# =============================================================================

fit_cue_lm <- function(df) {
  d <- df %>%
    dplyr::transmute(T = as.numeric(T), y = as.numeric(CUE)) %>%
    dplyr::filter(is.finite(T), is.finite(y))
  if (nrow(d) < 4 || length(unique(d$T)) < 3) {
    return(list(ok = FALSE, reason = "Too few temperatures"))
  }
  fit <- try(stats::lm(y ~ T + I(T^2), data = d), silent = TRUE)
  if (inherits(fit, "try-error")) return(list(ok = FALSE, reason = "lm failed"))
  co <- coef(fit); b0 <- co[[1]]; b1 <- co[[2]]; b2 <- co[[3]]
  sm <- summary(fit)
  T_opt   <- if (is.finite(b2) && b2 != 0) -b1 / (2 * b2) else NA_real_
  cue_opt <- if (is.finite(T_opt)) b0 + b1 * T_opt + b2 * T_opt^2 else NA_real_
  Tg <- seq(min(d$T), max(d$T), length.out = 200)
  list(
    ok = TRUE, reason = NA_character_,
    coefs = tibble::tibble(
      b0 = b0, b1_T = b1, b2_T2 = b2,
      T_opt_C = T_opt, CUE_at_opt = cue_opt,
      concave = is.finite(b2) && b2 < 0,
      r_squared = sm$r.squared, n = nrow(d)),
    preds = tibble::tibble(T = Tg, y_fit = pmin(pmax(b0 + b1 * Tg + b2 * Tg^2, 0), 1)),
    data  = d
  )
}

fit_cue_by_otu <- function(df) {
  pts <- list(); prd <- list(); cof <- list()
  for (cl in sort(unique(df$OTU))) {
    out <- fit_cue_lm(df %>% dplyr::filter(OTU == cl))
    if (!isTRUE(out$ok)) {
      message(sprintf("CUE fit failed for OTU %s: %s", cl, out$reason)); next
    }
    pts[[length(pts) + 1]] <- out$data  %>% dplyr::mutate(OTU = cl)
    prd[[length(prd) + 1]] <- out$preds %>% dplyr::mutate(OTU = cl)
    cof[[length(cof) + 1]] <- out$coefs %>% dplyr::mutate(OTU = cl)
  }
  list(points = dplyr::bind_rows(pts), preds = dplyr::bind_rows(prd),
       coefs = dplyr::bind_rows(cof))
}

cue_out  <- fit_cue_by_otu(results_plot)
cue_lab  <- "CUE (carbon-use efficiency)"

p_box_cue <- box_by_otu(results_plot, "CUE", "Carbon-use efficiency (by isolate)", cue_lab)

p_cue_vs_T <- if (nrow(cue_out$points) > 0) {
  ggplot2::ggplot() +
    ggplot2::geom_point(data = cue_out$points,
      ggplot2::aes(T, y, colour = factor(OTU)), size = 2.1, alpha = 0.85) +
    ggplot2::geom_line(data = cue_out$preds,
      ggplot2::aes(T, y_fit, colour = factor(OTU)), linewidth = 1.0) +
    otu_colour_scale + otu_facet +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = "CUE vs Temperature (quadratic fit, by OTU)",
                  x = "Temperature (°C)", y = cue_lab) +
    ggplot2::theme_classic(12)
} else {
  scatter_by_otu(results_plot, "CUE", "CUE vs Temperature (by isolate)", cue_lab)
}

box_cue_png  <- file.path(figures_dir, "boxplot_CUE.png")
cue_coefs_csv <- file.path(tables_dir, "cue_quadratic_fit_coefs.csv")
if (nrow(cue_out$coefs) > 0) {
  readr::write_csv(attach_otu_size(cue_out$coefs), cue_coefs_csv)
}


# =============================================================================
# 10g) Individual per-isolate thermal-performance curves (growth & respiration)
# =============================================================================
# One clean TPC per isolate: replicate points (each labelled with its replicate
# number) plus the Sharpe-Schoolfield fitted curve, titled with the isolate name.
# Saved as (a) a multi-page PDF, one isolate per page, and (b) an individual PNG
# per isolate under figures/TPC_<rate>/. This is the readable, per-isolate view
# (the faceted overview PNGs pack all isolates into one busy grid).
# =============================================================================

.repel_ok <- requireNamespace("ggrepel", quietly = TRUE)

make_one_tpc <- function(iso, value_col, ss_out, arr_out, rate_lab) {
  v <- rlang::sym(value_col)
  # Every series with a usable fit is shown (there is no RMSE filtering).
  pts <- results_plot %>%
    dplyr::filter(OTU == iso, is.finite(!!v), !!v > 0) %>%
    dplyr::mutate(rep_lab = sub("^R", "", as.character(Replicate)))

  ln <- if (!is.null(ss_out$preds) && nrow(ss_out$preds) > 0) {
    dplyr::filter(ss_out$preds, OTU == iso)
  } else ss_out$preds[0, , drop = FALSE]

  ss_r2 <- NA_real_
  if (!is.null(ss_out$coefs) && nrow(ss_out$coefs) > 0) {
    .r <- ss_out$coefs %>% dplyr::filter(OTU == iso) %>% dplyr::pull(r_squared)
    if (length(.r)) ss_r2 <- .r[1]
  }
  E_eV <- NA_real_
  if (!is.null(arr_out$coefs) && nrow(arr_out$coefs) > 0) {
    .e <- arr_out$coefs %>% dplyr::filter(OTU == iso) %>% dplyr::pull(E_eV)
    if (length(.e)) E_eV <- .e[1]
  }

  sub <- paste0(
    if (is.finite(E_eV)) sprintf("E = %.2f eV", E_eV) else "E = NA",
    if (is.finite(ss_r2)) sprintf("   |   SS R² = %.2f", ss_r2) else "",
    sprintf("   |   n = %d", nrow(pts))
  )

  lab_geom <- function(d) {
    if (nrow(d) == 0) return(NULL)
    if (.repel_ok) {
      ggrepel::geom_text_repel(
        data = d, ggplot2::aes(T, !!v, label = rep_lab),
        size = 3, colour = "grey30", max.overlaps = Inf, seed = 1,
        min.segment.length = 0, segment.size = 0.2, box.padding = 0.25)
    } else {
      ggplot2::geom_text(
        data = d, ggplot2::aes(T, !!v, label = rep_lab),
        size = 3, colour = "grey30", vjust = -0.6, check_overlap = TRUE)
    }
  }

  p <- ggplot2::ggplot()
  if (nrow(ln) > 0) {
    p <- p + ggplot2::geom_line(data = ln, ggplot2::aes(T, y_fit),
                                colour = "#d95f0e", linewidth = 1.1)
  }
  # Every point is shown (no RMSE filtering): solid blue, labelled by replicate.
  p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(T, !!v),
                               size = 2.6, alpha = 0.85, colour = "#2c7fb8") +
    lab_geom(pts)
  p +
    ggplot2::labs(title = otu_name_of(iso), subtitle = sub,
                  x = "Temperature (°C)", y = rate_lab) +
    ggplot2::theme_bw(13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   panel.grid.minor = ggplot2::element_blank())
}

save_tpc_set <- function(value_col, ss_out, arr_out, rate_lab, tag) {
  isos <- sort(unique(results_plot$OTU))
  if (length(isos) == 0) return(invisible(NULL))

  pdf_out <- file.path(figures_dir, sprintf("TPC_%s_by_isolate.pdf", tag))
  grDevices::pdf(pdf_out, width = 7, height = 5)
  for (iso in isos) print(make_one_tpc(iso, value_col, ss_out, arr_out, rate_lab))
  grDevices::dev.off()

  png_dir <- make_dir(file.path(figures_dir, sprintf("TPC_%s", tag)))
  for (iso in isos) {
    fn <- sprintf("%02d_%s.png", iso, gsub("[^A-Za-z0-9]+", "_", otu_name_of(iso)))
    ggplot2::ggsave(file.path(png_dir, fn),
                    make_one_tpc(iso, value_col, ss_out, arr_out, rate_lab),
                    width = 6, height = 4.2, dpi = 300)
  }
  message("Saved individual ", tag, " TPCs -> ", pdf_out, " and ", png_dir, "/")
}

save_tpc_set("growth_fgC_h", ss_growth_out, arr_growth_out,
             expression("Growth (fg C h"^{-1}*")"), "growth")
save_tpc_set("respiration_fgC_h", ss_resp_out, arr_resp_out,
             expression("Respiration (fg C h"^{-1}*")"), "respiration")


# =============================================================================
# 10c) Save descriptive + TPC plots to figures_dir
# =============================================================================
# The faceted overview PNGs cram all 18 isolates into one grid, so they are
# saved large. For the readable per-isolate view use the TPC_*_by_isolate.pdf /
# figures/TPC_growth, figures/TPC_respiration outputs from section 10g above.

save_plot <- function(p, path, width = 6, height = 4.5, dpi = 300) {
  ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi)
}

# Larger canvas for the 18-panel faceted overviews so each panel is legible.
save_facet <- function(p, path) save_plot(p, path, width = 14, height = 9)

save_facet(p_box_growth,           box_growth_png)
save_facet(p_box_resp,             box_resp_png)
save_facet(p_box_total_resp,       box_total_resp_png)
save_facet(p_growth_vs_T,          scatter_growth_png)
save_facet(p_resp_vs_T,            scatter_resp_png)
save_facet(p_total_resp_vs_T,      scatter_total_resp_png)
save_facet(p_rg_t,                 resp_over_growth_png)
save_facet(p_box_growth_biomass,   box_growth_biomass_png)
save_facet(p_box_resp_biomass,     box_resp_biomass_png)
save_facet(p_growth_biomass_vs_T,  scatter_growth_biomass_png)
save_facet(p_resp_biomass_vs_T,    scatter_resp_biomass_png)
save_facet(p_box_K,                box_K_png)
save_facet(p_K_vs_T,               scatter_K_png)
save_facet(p_ss_growth,            ss_growth_png)
save_facet(p_ss_resp,              ss_resp_png)
save_facet(p_ss_K,                 ss_K_png)
save_facet(p_arr_growth,           arr_growth_png)
save_facet(p_arr_resp,             arr_resp_png)
save_facet(p_box_cue,              box_cue_png)
save_facet(p_cue_vs_T,             cue_vs_T_png)
if (!is.null(p_activation_E)) save_plot(p_activation_E, activation_energy_png)

message("Saved descriptive + TPC plots to: ", figures_dir)




message("\nDone: 07_oxygen_fits.R completed successfully.")
message("Key output: ", derived_csv)
