#!/usr/bin/env Rscript
# =============================================================================
# 19_rk_covariance_check.R - does discarding cov(r, K) misstate the uncertainty?
# =============================================================================
#   Rscript scripts/19_rk_covariance_check.R
#
# READ-ONLY with respect to the pipeline. Reads the committed trimmed traces and
# the derived table, writes only its own tables and one figure. Touches no
# estimator, no constant and no published number.
#
# -----------------------------------------------------------------------------
# WHY THIS EXISTS
# -----------------------------------------------------------------------------
# Each oxygen trace is fitted ONCE and yields TWO parameters:
#
#     O2(t) = O2_0 + (K/r) (1 - e^{r t})
#
# so r and K are not independent measurements. They are two coordinates of a
# single solution, and the fit trades one against the other along a ridge. The
# pipeline records only the marginal standard errors (fit_coefficients_long.csv
# has Estimate and SE, no covariance), so everything downstream treats them as
# independent. This script asks what that costs.
#
# It matters because respiration is a COMBINATION of the two. With
# N0 = N_inoc e^{r delta} and R = c K / N0,
#
#     log R = log c - log N_inoc + log K - r * delta                     (1)
#
# The variance of a difference carries a cross term:
#
#     Var(log R) = Var(log K) + delta^2 Var(r) - 2 delta Cov(log K, r)   (2)
#
# and Cov(log K, r) is NEGATIVE here (the ridge), so the -2 delta Cov term is
# POSITIVE and the true variance EXCEEDS the independent one. Dropping the
# covariance understates the uncertainty.
#
# -----------------------------------------------------------------------------
# WHY THE SISTER REPOSITORY GETS THE OPPOSITE ANSWER
# -----------------------------------------------------------------------------
# OxygenModel anchors N0 at the END and back-projects:
#
#     N0 = FC_Final * C * e^{-r Delta}   =>   log R = log K - log(FC_Final C) + r * Delta
#
# a PLUS on the r term, so its cross term is +2 Delta Cov, which is NEGATIVE and
# SHRINKS the variance. Its D5 study found exactly that: sd 0.0035 against 0.0092
# if treated as independent, and concluded the r-K ridge is "protective".
#
# Candidas anchors at the START and projects FORWARD, giving equation (1) with a
# MINUS. Same ridge, same correlation, opposite consequence. The D5 result does
# NOT transfer, and assuming it does would be a mistake in the reassuring
# direction. That asymmetry is the reason this script exists as its own check
# rather than a citation of the sister repository.
#
# -----------------------------------------------------------------------------
# WHAT IT DOES
# -----------------------------------------------------------------------------
#   1. Refits every analysed trace with minpack.lm::nlsLM, STARTING FROM THE
#      COMMITTED SOLUTION so it converges immediately to the same optimum, and
#      recovers the full 3x3 covariance matrix the pipeline throws away.
#   2. Reports the within-curve correlation between r-hat and K-hat.
#   3. Evaluates equation (2) three ways: covariance ignored, this project's
#      minus sign, and the sister project's plus sign.
#   4. Compares the within-curve term against the BETWEEN-REPLICATE scatter,
#      which is what actually decides whether any of it matters.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2)
})

.this_dir <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  d  <- if (length(fa)) dirname(normalizePath(
    # R replaces every space in --file= with "~+~"; this project path has one.
    gsub("~+~", " ", sub("^--file=", "", fa[1]), fixed = TRUE), mustWork = FALSE))
        else tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
  if (!length(d) || is.na(d) || !nzchar(d)) d <- getwd()
  d
})
source(file.path(.this_dir, "config.R"))

fig_keep_add("fig_rk_covariance_check")

hdr <- function(x) message("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74))


# =============================================================================
# 1) The analysed curve set, built exactly as 09_bayesian_models.R builds it
# =============================================================================

if (!file.exists(derived_csv)) stop("Not found: ", derived_csv, "\nRun 02 -> 03 -> 06 -> 07 first.")

d <- readr::read_csv(derived_csv, show_col_types = FALSE) %>%
  dplyr::mutate(T = as.numeric(T), OTU = as.integer(OTU),
                Replicate = toupper(as.character(Replicate)))

nm <- readr::read_csv(file.path(tables_dir, "otu_names.csv"), show_col_types = FALSE) %>%
  dplyr::transmute(OTU = as.integer(OTU), Group = as.character(group),
                   Isolate = as.character(otu_name))
d <- dplyr::inner_join(d, nm, by = "OTU")
d <- drop_excluded(d, "derived")

if (exists("PLOT_EXCLUDE_POINTS") && is.data.frame(PLOT_EXCLUDE_POINTS) &&
    nrow(PLOT_EXCLUDE_POINTS) > 0) {
  .ex <- PLOT_EXCLUDE_POINTS %>%
    dplyr::mutate(T = as.numeric(T), OTU = as.integer(OTU),
                  Replicate = toupper(as.character(Replicate))) %>%
    dplyr::distinct(T, OTU, Replicate)
  d <- dplyr::anti_join(d, .ex, by = c("T", "OTU", "Replicate"))
}

d <- d %>% dplyr::filter(is.finite(r), r > 0, is.finite(K), K > 0,
                         is.finite(delta_Ninoc_to_N0_min))
hdr("1  Curve set")
message(sprintf("  %d analysed series", nrow(d)))


# =============================================================================
# 2) Recover the covariance the pipeline discards
# =============================================================================
# The trimmed traces are what 07 fitted (FIT_TO_SPLINE = FALSE, so the raw
# Oxygen column). Starting nlsLM from the committed (r, K, O2_0) means it lands
# on the same optimum in a step or two; the point is vcov(), not re-estimation.

trim_csv <- file.path(tables_dir, "Oxygen_Data_Smoothed_Trimmed.csv")
if (!file.exists(trim_csv)) stop("Not found: ", trim_csv, "\nRun 03_trimming.R first.")

ox <- readr::read_csv(trim_csv, show_col_types = FALSE,
                      col_select = c(Time, T, OTU, Replicate, Oxygen)) %>%
  dplyr::mutate(T = as.numeric(T), OTU = as.integer(OTU),
                Replicate = toupper(as.character(Replicate)))

oxs  <- split(ox, list(ox$T, ox$OTU, ox$Replicate), drop = TRUE)
kfun <- function(T, O, R) paste(T, O, R, sep = ".")

hdr("2  Refitting to recover cov(r, K)")
message("  (starts from the committed solution; this is not a re-estimation)")

out <- vector("list", nrow(d))
n_fail <- 0L
for (i in seq_len(nrow(d))) {
  row <- d[i, ]
  g   <- oxs[[kfun(row$T, row$OTU, row$Replicate)]]
  if (is.null(g) || nrow(g) < 10) { n_fail <- n_fail + 1L; next }
  g   <- g[order(g$Time), ]
  df0 <- data.frame(Oxygen_used = g$Oxygen, Time0 = g$Time - min(g$Time))

  fit <- try(minpack.lm::nlsLM(
    Oxygen_used ~ resp_model(r, K, Time0, O2_0), data = df0,
    start   = list(r = row$r, K = row$K, O2_0 = row$O2_0),
    control = minpack.lm::nls.lm.control(maxiter = 200, ftol = 1e-12, ptol = 1e-12)),
    silent = TRUE)
  if (inherits(fit, "try-error")) { n_fail <- n_fail + 1L; next }
  V <- try(vcov(fit), silent = TRUE)
  if (inherits(V, "try-error") || any(!is.finite(V))) { n_fail <- n_fail + 1L; next }

  cf <- coef(fit)
  out[[i]] <- tibble::tibble(
    T = row$T, OTU = row$OTU, Replicate = row$Replicate,
    Group = row$Group, Isolate = row$Isolate,
    r = cf[["r"]], K = cf[["K"]], delta = row$delta_Ninoc_to_N0_min,
    var_r = V["r", "r"], var_K = V["K", "K"], cov_rK = V["r", "K"],
    corr_rK = V["r", "K"] / sqrt(V["r", "r"] * V["K", "K"]))
}
res <- dplyr::bind_rows(out)
message(sprintf("  refitted %d of %d  (%d could not be refitted)",
                nrow(res), nrow(d), n_fail))
stopifnot(nrow(res) > 0)

hdr("3  Within-curve correlation of r-hat and K-hat")
qc <- stats::quantile(res$corr_rK, c(0, .05, .25, .5, .75, .95, 1), na.rm = TRUE)
for (i in seq_along(qc)) message(sprintf("  %-5s %8.4f", names(qc)[i], qc[[i]]))
message(sprintf("\n  negative in %.1f%% of curves", 100 * mean(res$corr_rK < 0, na.rm = TRUE)))
message(sprintf("  median delta (min): %.1f", stats::median(res$delta, na.rm = TRUE)))


# =============================================================================
# 3) Equation (2), three ways
# =============================================================================
# Delta method on the log scale:
#     Var(log K)     = Var(K) / K^2
#     Cov(log K, r)  = Cov(K, r) / K

res <- res %>% dplyr::mutate(
  v_logK    = var_K / K^2,
  v_rdelta  = delta^2 * var_r,
  cov_logKr = cov_rK / K,
  var_indep = v_logK + v_rdelta,                        # what the code assumes
  var_minus = v_logK + v_rdelta - 2 * delta * cov_logKr, # this project
  var_plus  = v_logK + v_rdelta + 2 * delta * cov_logKr  # sister project
) %>%
  dplyr::filter(is.finite(var_indep), var_indep > 0,
                is.finite(var_minus), var_minus > 0,
                is.finite(var_plus),  var_plus  > 0) %>%
  dplyr::mutate(sd_indep = sqrt(var_indep),
                sd_minus = sqrt(var_minus),
                sd_plus  = sqrt(var_plus),
                ratio_minus = sd_minus / sd_indep,
                ratio_plus  = sd_plus  / sd_indep)

hdr("4  sd of log(respiration) per curve")
message(sprintf("  covariance ignored (current code)      %.5f", stats::median(res$sd_indep)))
message(sprintf("  this project   log K - r*delta         %.5f   x%.2f",
                stats::median(res$sd_minus), stats::median(res$ratio_minus)))
message(sprintf("  sister project log K + r*delta         %.5f   x%.2f",
                stats::median(res$sd_plus),  stats::median(res$ratio_plus)))
message("\n  So the covariance INFLATES the uncertainty here and DEFLATES it there.")
message("  The sister repository's 'protective ridge' does not transfer.")


# =============================================================================
# 4) Does it matter? Against the between-replicate scatter
# =============================================================================
# The hierarchical models estimate replicate scatter from the data. A
# within-curve term far below it cannot move a posterior appreciably, however
# wrong it is. This is the same argument the sister repository's D5 made, and it
# is what decides the question.

bt <- res %>%
  dplyr::mutate(logR_rel = log(K) - r * delta) %>%
  dplyr::group_by(Group, T, OTU) %>%
  dplyr::filter(dplyr::n() >= 3) %>%
  dplyr::summarise(within  = stats::median(sd_minus),
                   between = stats::sd(logR_rel), .groups = "drop") %>%
  dplyr::filter(is.finite(between), between > 0)

between_med <- stats::median(bt$between)
w_now  <- stats::median(res$sd_indep)
w_true <- stats::median(res$sd_minus)
tot_now  <- sqrt(between_med^2 + w_now^2)
tot_true <- sqrt(between_med^2 + w_true^2)

hdr("5  Does it matter?")
message(sprintf("  isolate x temperature cells with >= 3 replicates : %d", nrow(bt)))
message(sprintf("  between-replicate sd                             %.4f   <- dominant", between_med))
message(sprintf("  within-curve sd, as coded                        %.4f", w_now))
message(sprintf("  within-curve sd, covariance included             %.4f   (+%.0f%%)",
                w_true, 100 * (w_true / w_now - 1)))
message(sprintf("  between / within                                 %.1fx",
                stats::median(bt$between / bt$within)))
message(sprintf("\n  total sd, as coded                               %.5f", tot_now))
message(sprintf("  total sd, covariance included                    %.5f", tot_true))
message(sprintf("  => reported intervals widen by                   %.2f%%",
                100 * (tot_true / tot_now - 1)))

w95 <- stats::quantile(res$sd_minus, .95); i95 <- stats::quantile(res$sd_indep, .95)
message(sprintf("  worst 5%% of curves                              %.2f%%",
                100 * (sqrt(between_med^2 + w95^2) / sqrt(between_med^2 + i95^2) - 1)))

verdict <- 100 * (tot_true / tot_now - 1)
message("\n  VERDICT: ",
        if (verdict < 1) {
          "real but immaterial. The correction is smaller than the rounding\n           in every reported number; no conclusion is affected."
        } else if (verdict < 5) {
          "small but worth reporting as a supplementary robustness note."
        } else {
          "MATERIAL. Propagate the joint (r, K) uncertainty before submission."
        })


# =============================================================================
# 5) Outputs
# =============================================================================

readr::write_csv(res %>% dplyr::select(T, OTU, Replicate, Group, Isolate, r, K, delta,
                                       var_r, var_K, cov_rK, corr_rK,
                                       sd_indep, sd_minus, sd_plus,
                                       ratio_minus, ratio_plus),
                 file.path(tables_dir, "rk_covariance_percurve.csv"))

summ <- tibble::tibble(
  quantity = c("n curves", "median corr(r,K)", "fraction negative",
               "median sd_indep", "median sd_minus", "median sd_plus",
               "ratio minus/indep", "ratio plus/indep",
               "between-replicate sd", "between/within",
               "interval widening %", "worst 5% widening %"),
  value = c(nrow(res), stats::median(res$corr_rK), mean(res$corr_rK < 0),
            w_now, w_true, stats::median(res$sd_plus),
            stats::median(res$ratio_minus), stats::median(res$ratio_plus),
            between_med, stats::median(bt$between / bt$within),
            verdict, 100 * (sqrt(between_med^2 + w95^2) / sqrt(between_med^2 + i95^2) - 1)))
readr::write_csv(summ, file.path(tables_dir, "rk_covariance_summary.csv"))

pdat <- dplyr::bind_rows(
  res %>% dplyr::transmute(Group, sd = sd_indep, which = "covariance ignored"),
  res %>% dplyr::transmute(Group, sd = sd_minus, which = "this project  (log K - r*delta)"),
  res %>% dplyr::transmute(Group, sd = sd_plus,  which = "sister project (log K + r*delta)"))
pdat$which <- factor(pdat$which, levels = unique(pdat$which))

p <- ggplot2::ggplot(pdat, ggplot2::aes(x = which, y = sd, fill = which)) +
  ggplot2::geom_boxplot(outlier.size = 0.4, alpha = 0.75, show.legend = FALSE) +
  ggplot2::geom_hline(yintercept = between_med, linetype = 2, colour = "grey30") +
  ggplot2::annotate("text", x = 0.6, y = between_med * 1.12, hjust = 0, size = 3,
                    label = sprintf("between-replicate sd = %.3f", between_med)) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    title = "Discarding cov(r, K) misstates the within-curve uncertainty",
    subtitle = sprintf(
      "%d series | median corr(r,K) = %.3f | intervals widen %.2f%% once corrected",
      nrow(res), stats::median(res$corr_rK), verdict),
    x = NULL, y = "sd of log(respiration), log scale") +
  ggplot2::theme_classic(11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 12, hjust = 1))

ggplot2::ggsave(file.path(figures_dir, "fig_rk_covariance_check.png"), p,
                width = 7.5, height = 5, dpi = 300)

hdr("19_rk_covariance_check.R done")
message("  tables/rk_covariance_percurve.csv")
message("  tables/rk_covariance_summary.csv")
message("  figures/fig_rk_covariance_check.png")
