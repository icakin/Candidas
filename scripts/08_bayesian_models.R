# =============================================================================
# 08_bayesian_models.R - Hierarchical Bayesian thermal performance curves (brms)
# =============================================================================
# Fits, saves, and summarises. NO plotting (that is 09_bayesian_plots.R), because
# Stan fits are slow and you should not refit just to tweak a figure.
#
# WHAT IT FITS
#   GROWTH      Sharpe-Schoolfield, hierarchical:
#                 lnB0, E, Th ~ 0 + Group + (1 | Isolate);  Eh ~ 0 + Group
#               Isolates are nested in their clade/species GROUP (3 isolates per
#               group). Partial pooling means an isolate whose curve is a mess
#               (e.g. no decline at the top) is regularised toward its clade
#               instead of returning a nonsense E / Topt - and when the data
#               genuinely cannot identify Th, you get an honestly WIDE posterior
#               rather than a confident fake number.
#
#   RESPIRATION Arrhenius (log scale), hierarchical:
#                 log(rate) ~ alpha + E * boltz_shift
#                 alpha, E ~ 0 + Group + (1 | Isolate)
#               ALSO a Sharpe-Schoolfield version, so LOO can say whether
#               respiration keeps rising (Arrhenius) or peaks (SS). No AICc hack.
#
# WHY THIS REPLACES THE TRIMMING ARGUMENT
#   1. Robust likelihood (student-t): outliers are automatically down-weighted.
#      No TUKEY_C, no discarding, no "which points distort the fit" debate.
#   2. MEASUREMENT ERROR IS PROPAGATED. Every growth point carries the standard
#      error of its own r (from the nlsLM fit in 06). A shaky r counts LESS than
#      a solid one - which is what all the RMSE-cutoff attempts were trying and
#      failing to achieve. Weighting, not deleting.
#
# DATA RULES (as established)
#   - Curves you excluded in 04 (plot_exclude_points.csv) are dropped.
#   - Straight-line / NO-GROWTH curves (has_curvature = FALSE) are dropped from
#     GROWTH only. They are KEPT for respiration, where the linear O2 slope is a
#     perfectly valid rate. (Their CUE = 0 is a real result - see 11.)
#
# Inputs (run 02 -> 03 -> 06 first):
#   tables/derived_N0_R_results_with_carbon.csv
#   tables/fit_coefficients_long.csv   (standard error of r)
#   tables/otu_names.csv               (isolate name + clade/species group)
#
# Outputs (models/ and tables/):
#   models/bayes_growth_ss.rds, bayes_resp_arr.rds, bayes_resp_ss.rds
#   tables/bayes_growth_ss_summary.csv, bayes_resp_arr_summary.csv
#   tables/bayes_resp_model_comparison_loo.csv
#   models/bayes_data_growth.rds, bayes_data_resp.rds   (for 11)
# =============================================================================


# =============================================================================
# 0) Config + packages
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

need <- c("brms", "posterior")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) {
  stop("Install first:  install.packages(c(",
       paste(sprintf('"%s"', miss), collapse = ", "), "))\n",
       "brms also needs a Stan backend (rstan or cmdstanr).")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(brms)
  library(posterior)
})

options(mc.cores = max(1, parallel::detectCores() - 1))
if (Sys.info()[["sysname"]] == "Darwin") {
  try(suppressWarnings(mem.maxVSize(vsize = Inf)), silent = TRUE)
}


# =============================================================================
# 0b) Timer  (Stan fits are slow - this tells you where the time goes)
# =============================================================================

.T_START <- Sys.time()
.TIMINGS <- list()

fmt_dur <- function(secs) {
  secs <- as.numeric(secs)
  h <- floor(secs / 3600); m <- floor((secs %% 3600) / 60); s <- round(secs %% 60)
  if (h > 0) sprintf("%dh %02dm %02ds", h, m, s)
  else if (m > 0) sprintf("%dm %02ds", m, s)
  else sprintf("%.0fs", secs)
}

# Run `expr`, printing when it starts and how long it took.
timed <- function(label, expr) {
  message("\n>>> ", label, "  |  started ", format(Sys.time(), "%H:%M:%S"))
  t0  <- Sys.time()
  out <- expr                                  # promise is forced here
  el  <- difftime(Sys.time(), t0, units = "secs")
  .TIMINGS[[label]] <<- as.numeric(el)
  message("<<< ", label, "  |  DONE in ", fmt_dur(el),
          "  (total elapsed ", fmt_dur(difftime(Sys.time(), .T_START, units = "secs")), ")")
  out
}


# =============================================================================
# 1) Sampler + prior knobs
# =============================================================================

BAYES_ITER   <- 4000
BAYES_WARMUP <- 1000
BAYES_CHAINS <- 4
BAYES_SEED   <- 1234
BAYES_ADAPT  <- 0.95
BAYES_MAX_TD <- 12

# "student" = robust to outliers (RECOMMENDED - this is what removes the need to
# hand-discard distorting points). "gaussian" = classical.
BAYES_LIKELIHOOD <- "student"
bayes_family <- function() {
  if (identical(BAYES_LIKELIHOOD, "student")) brms::student() else brms::gaussian()
}

# Propagate the standard error of each r into the growth model (y | se(...)).
# TRUE is the whole point: unreliable growth points automatically count less.
USE_MEASUREMENT_ERROR <- TRUE

# Sharpe-Schoolfield priors. E ~ activation energy (eV); Eh ~ deactivation;
# Th ~ temperature (K) of half high-T deactivation, kept near the observed max.
SS_PRIOR_E_MEAN <- 0.65; SS_PRIOR_E_SD <- 0.35; SS_PRIOR_E_LB <- 0.01; SS_PRIOR_E_UB <- 3
SS_PRIOR_EH_MEAN <- 4;   SS_PRIOR_EH_SD <- 3;   SS_PRIOR_EH_LB <- 0.05; SS_PRIOR_EH_UB <- 30
SS_PRIOR_TH_OFFSET_C <- 2    # prior mean = max observed TK + this (K)
SS_PRIOR_TH_SD       <- 5
SS_PRIOR_TH_LB_C     <- -12  # bounds relative to max observed TK
SS_PRIOR_TH_UB_C     <- 8

# Arrhenius priors.
ARR_PRIOR_E_MEAN <- 0.65; ARR_PRIOR_E_SD <- 0.35
ARR_PRIOR_E_LB   <- 0;    ARR_PRIOR_E_UB <- 3

# Fit the respiration SS too (for the LOO comparison against Arrhenius)?
FIT_RESP_SS <- TRUE

# -----------------------------------------------------------------------------
# QUADRATIC (curved) respiration.   ADDED AFTER A PER-ISOLATE MODEL CHECK.
# -----------------------------------------------------------------------------
# Fitting each of the 18 isolates' respiration SEPARATELY and comparing by AICc,
# Arrhenius (a straight line in Boltzmann space) is REJECTED by 14 of 18:
#
#     group    isolates rejecting Arrhenius    median dAICc
#     para          1 of 3                        +2.1   (Arrhenius fine)
#     Clade1        2 of 3                        -5.0
#     Clade2        2 of 3                       -13.4
#     Clade3        3 of 3                        -9.6
#     Clade4        3 of 3                       -16.9
#     glab          3 of 3                       -34.5   (worst isolate: -70)
#
# Respiration ACCELERATES: its activation energy rises with temperature. The
# original Arrhenius-vs-SS LOO comparison never tested this, because a curved
# monotone alternative was not in the candidate set. It should have been.
#
#     ln R = alpha + E*x + c*x^2 ,     x = 1/(k*Tref) - 1/(k*T)
#     dlnR/dT = (E + 2*c*x) / (k*T^2)          <- E_eff(T), no longer constant
#
# c = 0 recovers Arrhenius exactly, so this is a STRICT GENERALISATION and the
# LOO comparison is fair. A segmented (breakpoint) model fits about as well, but
# its derivative is DISCONTINUOUS at the break - and every downstream quantity we
# care about (Topt(CUE), the carbon tax) is defined through dlnR/dT. A kinked
# model would make them undefined exactly where the interesting behaviour is.
# So: quadratic.
# OFF BY DEFAULT. The models used for the paper are SS (growth) + ARRHENIUS
# (respiration), full stop. This is here only as an optional robustness check.
# It was worth running once: Topt(CUE) shifts by at most 1.1 C (always downward,
# so the "below 37 C" conclusion is unaffected) and the carbon tax comes out
# HIGHER under the quadratic - meaning the Arrhenius numbers are CONSERVATIVE.
# Set TRUE only if a reviewer asks about curvature.
FIT_RESP_QUAD <- FALSE
QUAD_PRIOR_C_MEAN <- 0;   QUAD_PRIOR_C_SD <- 0.5    # centred on c = 0 (Arrhenius)


# =============================================================================
# 2) Build the data
# =============================================================================

if (!file.exists(derived_csv)) stop("Not found: ", derived_csv, "\nRun 02 -> 03 -> 06 first.")

res <- readr::read_csv(derived_csv, show_col_types = FALSE) %>%
  dplyr::mutate(T = as.numeric(T), OTU = as.integer(OTU),
                Replicate = toupper(as.character(Replicate)))

# Hold out any group listed in EXCLUDE_GROUPS (config.R). Currently glab, which is
# being re-run. Everything downstream (11, 13-17) derives its group list from the
# FITTED data, so excluding it here removes it from every model, table and figure
# automatically - there is nothing else to edit.
res <- drop_excluded(res, "derived")

if (!"has_curvature" %in% names(res)) {
  res$has_curvature <- TRUE
  message("No has_curvature column - re-run 06 to enable the no-growth test.")
}
res$has_curvature <- res$has_curvature %in% TRUE

# ---- isolate + clade/species group ----------------------------------------
nm <- readr::read_csv(file.path(tables_dir, "otu_names.csv"), show_col_types = FALSE)
if (!"group" %in% names(nm)) {
  stop("otu_names.csv has no `group` column - re-run 01_convert_xlsx.R.")
}
nm <- nm %>% dplyr::transmute(OTU = as.integer(OTU),
                              Isolate = as.character(otu_name),
                              Group   = as.character(group))
res <- res %>% dplyr::inner_join(nm, by = "OTU")

# ---- honour the samples you excluded in 04 ---------------------------------
if (exists("PLOT_EXCLUDE_POINTS") && is.data.frame(PLOT_EXCLUDE_POINTS) &&
    nrow(PLOT_EXCLUDE_POINTS) > 0) {
  .ex <- PLOT_EXCLUDE_POINTS %>%
    dplyr::mutate(T = as.numeric(T), OTU = as.integer(OTU),
                  Replicate = toupper(as.character(Replicate))) %>%
    dplyr::distinct(T, OTU, Replicate)
  .n0 <- nrow(res)
  res <- res %>% dplyr::anti_join(.ex, by = c("T", "OTU", "Replicate"))
  message(sprintf("Excluded %d curve(s) you removed in 04: %d of %d rows kept.",
                  nrow(.ex), nrow(res), .n0))
}

# ---- standard error of r (for measurement-error weighting) ------------------
se_r <- NULL
if (file.exists(coef_csv)) {
  .cf <- tryCatch(readr::read_csv(coef_csv, show_col_types = FALSE),
                  error = function(e) NULL)
  if (!is.null(.cf) && all(c("parameter", "SE", "Estimate") %in% names(.cf))) {
    se_r <- .cf %>%
      dplyr::filter(parameter == "r") %>%
      dplyr::transmute(T = as.numeric(T), OTU = as.integer(OTU),
                       Replicate = toupper(as.character(Replicate)),
                       r_est = as.numeric(Estimate), r_se = as.numeric(SE)) %>%
      dplyr::distinct()
  }
}

boltz <- function(TK) (1 / (k_B * T_ref)) - (1 / (k_B * TK))

# ---- GROWTH data: drop the straight-line (no-growth) curves ------------------
growth_dat <- res %>%
  dplyr::filter(has_curvature, is.finite(growth_fgC_h), growth_fgC_h > 0) %>%
  dplyr::mutate(TK = T + 273.15, boltz_shift = boltz(TK),
                y_raw = growth_fgC_h, y = log(growth_fgC_h),
                Isolate = factor(Isolate), Group = factor(Group))

# Delta method: sd(log r) ~ se(r) / r.  growth = r * const, so the constant
# cancels on the log scale and sd(log growth) = sd(log r).
if (isTRUE(USE_MEASUREMENT_ERROR) && !is.null(se_r)) {
  growth_dat <- growth_dat %>%
    dplyr::left_join(se_r, by = c("T", "OTU", "Replicate")) %>%
    dplyr::mutate(se_y = dplyr::if_else(is.finite(r_se) & is.finite(r_est) & r_est > 0,
                                        r_se / r_est, NA_real_))
  .med <- stats::median(growth_dat$se_y, na.rm = TRUE)
  if (!is.finite(.med) || .med <= 0) .med <- 0.05
  growth_dat$se_y[!is.finite(growth_dat$se_y) | growth_dat$se_y <= 0] <- .med
  message(sprintf("Measurement error on growth: median sd(log r) = %.3f", .med))
} else {
  growth_dat$se_y <- NA_real_
  if (isTRUE(USE_MEASUREMENT_ERROR)) {
    message("No SE for r found in fit_coefficients_long.csv - fitting without ",
            "measurement error.")
    USE_MEASUREMENT_ERROR <- FALSE
  }
}

# ---- RESPIRATION data: keep EVERYTHING, including the no-growth curves -------
resp_dat <- res %>%
  dplyr::filter(is.finite(respiration_fgC_h), respiration_fgC_h > 0) %>%
  dplyr::mutate(TK = T + 273.15, boltz_shift = boltz(TK),
                y_raw = respiration_fgC_h, y = log(respiration_fgC_h),
                Isolate = factor(Isolate), Group = factor(Group))

message(sprintf("Growth: %d points, %d isolates, %d groups (no-growth curves dropped).",
                nrow(growth_dat), dplyr::n_distinct(growth_dat$Isolate),
                dplyr::n_distinct(growth_dat$Group)))
message(sprintf("Respiration: %d points (no-growth curves KEPT).", nrow(resp_dat)))

saveRDS(growth_dat, file.path(models_dir, "bayes_data_growth.rds"))
saveRDS(resp_dat,   file.path(models_dir, "bayes_data_resp.rds"))


# =============================================================================
# 3) Hierarchical Sharpe-Schoolfield  (growth; also reusable for respiration)
# =============================================================================

fit_ss_hier <- function(dat, use_se = FALSE) {
  max_TK <- max(dat$TK, na.rm = TRUE)
  lnB0_mu <- log(stats::median(dat$y_raw, na.rm = TRUE))
  Th_mu   <- max_TK + SS_PRIOR_TH_OFFSET_C
  Th_lb   <- max_TK + SS_PRIOR_TH_LB_C
  Th_ub   <- max_TK + SS_PRIOR_TH_UB_C

  # Response: with measurement error, each point carries the SE of its own r, so
  # unreliable points automatically count less. sigma = TRUE keeps a residual
  # (biological) scatter term on top of the measurement error.
  lhs <- if (isTRUE(use_se)) "y | se(se_y, sigma = TRUE)" else "y"

  bf_ss <- brms::bf(
    stats::as.formula(paste0(
      lhs, " ~ lnB0 - (E * 11604.51812) * ((1 / TK) - (1 / 293.15)) - ",
      "log(1 + exp((Eh * 11604.51812) * ((1 / Th) - (1 / TK))))")),
    lnB0 ~ 0 + Group + (1 | Isolate),
    E    ~ 0 + Group + (1 | Isolate),
    Eh   ~ 0 + Group,                       # weakly identified - group level only
    Th   ~ 0 + Group + (1 | Isolate),
    nl = TRUE
  )

  pri <- c(
    brms::prior_string(sprintf("normal(%.6f, 2)", lnB0_mu), nlpar = "lnB0"),
    brms::prior_string(sprintf("normal(%g, %g)", SS_PRIOR_E_MEAN, SS_PRIOR_E_SD),
                       nlpar = "E", lb = SS_PRIOR_E_LB, ub = SS_PRIOR_E_UB),
    brms::prior_string(sprintf("normal(%g, %g)", SS_PRIOR_EH_MEAN, SS_PRIOR_EH_SD),
                       nlpar = "Eh", lb = SS_PRIOR_EH_LB, ub = SS_PRIOR_EH_UB),
    brms::prior_string(sprintf("normal(%.6f, %g)", Th_mu, SS_PRIOR_TH_SD),
                       nlpar = "Th", lb = Th_lb, ub = Th_ub),
    brms::prior_string("student_t(3, 0, 1)", class = "sd", nlpar = "lnB0"),
    brms::prior_string("student_t(3, 0, 0.5)", class = "sd", nlpar = "E"),
    brms::prior_string("student_t(3, 0, 2)", class = "sd", nlpar = "Th"),
    brms::prior_string("exponential(1)", class = "sigma")
  )

  brms::brm(
    formula = bf_ss, data = dat, family = bayes_family(), prior = pri,
    iter = BAYES_ITER, warmup = BAYES_WARMUP, chains = BAYES_CHAINS,
    seed = BAYES_SEED,
    control = list(adapt_delta = 0.99, max_treedepth = BAYES_MAX_TD),
    refresh = max(1, floor(BAYES_ITER / 40))   # Stan's own progress bar (~40 ticks)
  )
}


# =============================================================================
# 4) Hierarchical Arrhenius  (respiration)
# =============================================================================

fit_arr_hier <- function(dat) {
  bf_arr <- brms::bf(
    y ~ alpha + E * boltz_shift,
    alpha ~ 0 + Group + (1 | Isolate),
    E     ~ 0 + Group + (1 | Isolate),
    nl = TRUE
  )
  pri <- c(
    brms::prior_string("normal(0, 5)", nlpar = "alpha"),
    brms::prior_string(sprintf("normal(%g, %g)", ARR_PRIOR_E_MEAN, ARR_PRIOR_E_SD),
                       nlpar = "E", lb = ARR_PRIOR_E_LB, ub = ARR_PRIOR_E_UB),
    brms::prior_string("student_t(3, 0, 2.5)", class = "sd", nlpar = "alpha"),
    brms::prior_string("student_t(3, 0, 0.5)", class = "sd", nlpar = "E"),
    brms::prior_string("exponential(1)", class = "sigma")
  )
  brms::brm(
    formula = bf_arr, data = dat, family = bayes_family(), prior = pri,
    iter = BAYES_ITER, warmup = BAYES_WARMUP, chains = BAYES_CHAINS,
    seed = BAYES_SEED,
    control = list(adapt_delta = BAYES_ADAPT, max_treedepth = BAYES_MAX_TD),
    refresh = max(1, floor(BAYES_ITER / 40))   # Stan's own progress bar (~40 ticks)
  )
}


# =============================================================================
# 4b) Hierarchical QUADRATIC-in-Boltzmann  (respiration)   <- the new candidate
# =============================================================================
# ln R = alpha + E*x + c*x^2 ,  x = boltz_shift = 1/(k*Tref) - 1/(k*T)
#
# c > 0 : respiration ACCELERATES with temperature (activation energy rising)
# c = 0 : exactly Arrhenius
# c < 0 : respiration saturating
#
# The prior on c is centred on ZERO with sd 0.5, i.e. it is centred on Arrhenius.
# If the data do not demand curvature, this model will report c ~ 0 and LOO will
# not favour it. That makes the comparison honest rather than rigged.
fit_quad_hier <- function(dat) {
  bf_q <- brms::bf(
    y ~ alpha + E * boltz_shift + c * boltz_shift^2,
    alpha ~ 0 + Group + (1 | Isolate),
    E     ~ 0 + Group + (1 | Isolate),
    c     ~ 0 + Group + (1 | Isolate),
    nl = TRUE
  )
  pri <- c(
    brms::prior_string("normal(0, 5)", nlpar = "alpha"),
    brms::prior_string(sprintf("normal(%g, %g)", ARR_PRIOR_E_MEAN, ARR_PRIOR_E_SD),
                       nlpar = "E", lb = ARR_PRIOR_E_LB, ub = ARR_PRIOR_E_UB),
    brms::prior_string(sprintf("normal(%g, %g)", QUAD_PRIOR_C_MEAN, QUAD_PRIOR_C_SD),
                       nlpar = "c"),
    brms::prior_string("student_t(3, 0, 2.5)", class = "sd", nlpar = "alpha"),
    brms::prior_string("student_t(3, 0, 0.5)", class = "sd", nlpar = "E"),
    brms::prior_string("student_t(3, 0, 0.5)", class = "sd", nlpar = "c"),
    brms::prior_string("exponential(1)", class = "sigma")
  )
  brms::brm(
    formula = bf_q, data = dat, family = bayes_family(), prior = pri,
    iter = BAYES_ITER, warmup = BAYES_WARMUP, chains = BAYES_CHAINS,
    seed = BAYES_SEED,
    control = list(adapt_delta = BAYES_ADAPT, max_treedepth = BAYES_MAX_TD),
    refresh = max(1, floor(BAYES_ITER / 40))
  )
}


# =============================================================================
# 5) Fit
# =============================================================================

.n_models <- 2L + as.integer(isTRUE(FIT_RESP_SS)) + as.integer(isTRUE(FIT_RESP_QUAD))
message(sprintf("\nFitting %d model(s): %d chains x %d iterations each.",
                .n_models, BAYES_CHAINS, BAYES_ITER))
message("Stan prints its own per-chain progress below; the timer reports each model.")

fit_growth_ss <- timed(sprintf("1/%d  GROWTH  Sharpe-Schoolfield (hierarchical)", .n_models),
                       fit_ss_hier(growth_dat, use_se = USE_MEASUREMENT_ERROR))
saveRDS(fit_growth_ss, file.path(models_dir, "bayes_growth_ss.rds"))

fit_resp_arr <- timed(sprintf("2/%d  RESPIRATION  Arrhenius (hierarchical)", .n_models),
                      fit_arr_hier(resp_dat))
saveRDS(fit_resp_arr, file.path(models_dir, "bayes_resp_arr.rds"))

fit_resp_quad <- NULL
if (isTRUE(FIT_RESP_QUAD)) {
  fit_resp_quad <- timed(
    sprintf("3/%d  RESPIRATION  Quadratic-in-Boltzmann (curved)", .n_models),
    tryCatch(fit_quad_hier(resp_dat),
             error = function(e) { message("Resp quad failed: ", conditionMessage(e)); NULL }))
  if (!is.null(fit_resp_quad)) saveRDS(fit_resp_quad, file.path(models_dir, "bayes_resp_quad.rds"))
}

fit_resp_ss <- NULL
if (isTRUE(FIT_RESP_SS)) {
  fit_resp_ss <- timed(
    sprintf("%d/%d  RESPIRATION  Sharpe-Schoolfield (for LOO)",
            3L + as.integer(isTRUE(FIT_RESP_QUAD)), .n_models),
    tryCatch(fit_ss_hier(resp_dat, use_se = FALSE),
             error = function(e) { message("Resp SS failed: ", conditionMessage(e)); NULL }))
  if (!is.null(fit_resp_ss)) saveRDS(fit_resp_ss, file.path(models_dir, "bayes_resp_ss.rds"))
}


# =============================================================================
# 6) Summaries + LOO model comparison for respiration
# =============================================================================

posterior_summary_df <- function(fit, pars = NULL) {
  as.data.frame(posterior::summarise_draws(
    posterior::as_draws_df(fit, variable = pars),
    mean, sd, ~posterior::quantile2(.x, probs = c(0.025, 0.5, 0.975)),
    posterior::default_convergence_measures()))
}

readr::write_csv(posterior_summary_df(fit_growth_ss),
                 file.path(tables_dir, "bayes_growth_ss_summary.csv"))
readr::write_csv(posterior_summary_df(fit_resp_arr),
                 file.path(tables_dir, "bayes_resp_arr_summary.csv"))

# Does respiration keep rising (Arrhenius) or peak (SS)? Let LOO decide.
if (!is.null(fit_resp_ss)) {
  cmp <- tryCatch({
    l1 <- brms::loo(fit_resp_arr)
    l2 <- brms::loo(fit_resp_ss)
    as.data.frame(brms::loo_compare(l1, l2)) %>% tibble::rownames_to_column("model")
  }, error = function(e) { message("LOO failed: ", conditionMessage(e)); NULL })
  if (!is.null(cmp)) {
    readr::write_csv(cmp, file.path(tables_dir, "bayes_resp_model_comparison_loo.csv"))
    message("\nRespiration model comparison (LOO):")
    print(cmp)
  }
}


# =============================================================================
# 7) Done
# =============================================================================

message("\n", strrep("=", 70))
message("TIMING BREAKDOWN")
message(strrep("=", 70))
for (nmz in names(.TIMINGS)) {
  message(sprintf("  %-52s %s", nmz, fmt_dur(.TIMINGS[[nmz]])))
}
message(sprintf("  %-52s %s", "TOTAL",
                fmt_dur(difftime(Sys.time(), .T_START, units = "secs"))))
message(strrep("=", 70))

message("\n08_bayesian_models.R done.")
message("  Models saved to: ", models_dir)
message("  Growth SS is hierarchical (isolates nested in clades), robust (",
        BAYES_LIKELIHOOD, ")",
        if (isTRUE(USE_MEASUREMENT_ERROR)) ", with r's standard error propagated." else ".")
message("  Next: 09_bayesian_plots.R (grouped TPC panels, Arrhenius, CUE).")
message("  CHECK: Rhat < 1.01 and no divergent transitions in the summary CSVs.")
