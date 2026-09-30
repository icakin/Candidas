#!/usr/bin/env Rscript
# =============================================================================
# 05_clade_contrasts.R - does the within-C. auris clade spread in the MODEL-FREE
# CUE optimum survive pairwise contrasts?
# =============================================================================
# WHY THIS EXISTS
# The phylogenomics write-up quotes a 2.74 C spread in CUE optimum across the
# four C. auris clades, taken from the FITTED (Sharpe-Schoolfield, per-cell)
# optima. That is within-auris fine structure, which is the class of result the
# paper cut for N0 fragility: the fitted per-cell quantities inherit the
# N0 = N_inoc exp(r delta) back-projection, and the fine within-auris ordering
# does not survive its removal (the between-taxon spread in E_K collapses from
# 0.22 to 0.07 eV term-free).
#
# The model-free optimum does NOT inherit that problem. CUE = 1/(1 + a K/(r e^{r delta}))
# with a temperature-independent, so the optimum is argmin of K/r and every
# conversion constant cancels. So the honest test is: in the model-free optima,
# do the pairwise clade differences exclude zero?
#
# Marginal 95% intervals are NOT that test - two overlapping intervals can still
# have a difference excluding zero, and two non-overlapping ones always do. This
# bootstraps the DIFFERENCE directly, resampling series within temperature within
# each clade independently, exactly as 02_scale_free.R bootstraps each optimum.
#
# Run: Rscript reports/C11_scale_free/R/05_clade_contrasts.R
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tibble); library(readr)})

HERE  <- "reports/C11_scale_free"
d     <- readr::read_csv(file.path(HERE, "tables/partA_series.csv"), show_col_types = FALSE)
NBOOT <- 2000
set.seed(11)                              # same seed/scheme as 02_scale_free.R
AURIS <- c("Clade1", "Clade2", "Clade3", "Clade4")

grid_argmin <- function(TT, logB) {       # verbatim from 02_scale_free.R
  i <- which.min(logB); gridT <- TT[i]; refined <- gridT
  if (i > 1 && i < length(TT)) {
    x <- TT[(i-1):(i+1)]; y <- logB[(i-1):(i+1)]
    cf <- stats::coef(stats::lm(y ~ poly(x, 2, raw = TRUE)))
    if (!is.na(cf[3]) && cf[3] > 0) {
      v <- -cf[2] / (2 * cf[3])
      if (is.finite(v) && v >= min(x) && v <= max(x)) refined <- unname(v)
    }
  }
  c(grid = gridT, refined = refined)
}
taxon_curve <- function(df, col = "B_free") {
  df %>% dplyr::group_by(T) %>%
    dplyr::summarise(logB = mean(log(.data[[col]])), n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(T)
}
# one bootstrap replicate of a clade's optimum, resampling series within temperature
boot_draws <- function(df, col = "B_free") {
  keys <- split(seq_len(nrow(df)), df$T)
  vapply(seq_len(NBOOT), function(b) {
    idx <- unlist(lapply(keys, function(ii) sample(ii, length(ii), replace = TRUE)))
    cv  <- taxon_curve(df[idx, ], col)
    grid_argmin(cv$T, cv$logB)[["refined"]]
  }, numeric(1))
}

cat("Model-free CUE optimum, per C. auris clade (bootstrap, NBOOT =", NBOOT, ")\n\n")
draws <- list(); pt <- c()
for (g in AURIS) {
  df <- dplyr::filter(d, group == g)
  cv <- taxon_curve(df)
  pt[g]    <- grid_argmin(cv$T, cv$logB)[["refined"]]
  draws[[g]] <- boot_draws(df)
  cat(sprintf("  %-8s T_opt = %5.2f C   [%5.2f, %5.2f]   n_series = %d\n",
              g, pt[g], quantile(draws[[g]], .025), quantile(draws[[g]], .975), nrow(df)))
}
cat(sprintf("\n  point-estimate spread across the four clades: %.2f C\n\n",
            max(pt) - min(pt)))

cat("Pairwise contrasts (difference of bootstrap draws; excludes zero = real)\n\n")
cat(sprintf("  %-20s %8s %18s %10s\n", "contrast", "diff", "95% CI", "P(same sign)"))
cat("  ", strrep("-", 60), "\n", sep = "")
res <- list()
for (i in 1:3) for (j in (i+1):4) {
  a <- AURIS[i]; b <- AURIS[j]
  dd <- draws[[a]] - draws[[b]]
  lo <- unname(quantile(dd, .025)); hi <- unname(quantile(dd, .975))
  psign <- max(mean(dd > 0), mean(dd < 0))
  excl  <- (lo > 0) || (hi < 0)
  cat(sprintf("  %-20s %8.2f  [%6.2f, %6.2f] %8.3f %s\n",
              paste(a, "-", b), pt[a] - pt[b], lo, hi, psign,
              if (excl) "  <== excludes 0" else ""))
  res[[length(res)+1]] <- tibble(contrast = paste(a, "-", b), diff = pt[a] - pt[b],
                                 lo = lo, hi = hi, P_same_sign = psign, excludes_zero = excl)
}
out <- bind_rows(res)
readr::write_csv(out, file.path(HERE, "tables/partB_clade_contrasts.csv"))

n_excl <- sum(out$excludes_zero)
cat(sprintf("\n  %d of 6 pairwise contrasts exclude zero (nominal, uncorrected).\n", n_excl))

# multiplicity: 6 contrasts tested at once. Bonferroni on the two-sided 5% level
# means a contrast must have P(same sign) >= 1 - 0.025/6 = 0.9958 to count.
thr <- 1 - 0.025 / 6
out$survives_bonferroni <- out$P_same_sign >= thr
n_bonf <- sum(out$survives_bonferroni)
cat(sprintf("  %d of 6 survive Bonferroni (P(same sign) >= %.4f).\n\n", n_bonf, thr))
readr::write_csv(out, file.path(HERE, "tables/partB_clade_contrasts.csv"))

if (n_bonf == 0) {
  cat("VERDICT: no clade contrast survives multiplicity. Do not quote the spread.\n")
} else {
  surv <- out[out$survives_bonferroni, ]
  cat("VERDICT: the bare spread must not be quoted. The only contrast(s) that survive:\n")
  for (k in seq_len(nrow(surv)))
    cat(sprintf("    %s  %.2f C  [%.2f, %.2f]\n", surv$contrast[k], surv$diff[k], surv$lo[k], surv$hi[k]))
  cat("  An effect below 1 C inside an interval several degrees wide is not a quotable result.\n")
}
