# =============================================================================
# C11 PART B - scale-free E_r, E_K and the model-free CUE optimum, per taxon
# =============================================================================
# Works from the fitted r and K ONLY. Nothing here converts to per-cell carbon,
# and nothing here uses N0, the carbon quota, the O2-to-C factor, RQ or the
# inoculation density. Anything that needs those belongs in PART C.
# =============================================================================

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
       value = TRUE)[1])), "00_common.R"))

set.seed(11)          # bootstrap only; no estimator here is stochastic
NBOOT <- 2000

d <- readr::read_csv(file.path(C11_TAB, "partA_series.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(group = factor(group, levels = C11_TAXA),
                x = c11_boltz(T),
                log_r = log(r), log_K = log(K), log_B = log(B_free))

hdr <- function(x) message("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74))

# =============================================================================
# B.1  Which functional form for which quantity
# =============================================================================
# r turns over inside the measured range, so Sharpe-Schoolfield is identifiable
# for it. K does NOT turn over, so SS is not identifiable for K: the pipeline's
# own per-isolate SS fits to K run to the E = 5.0 bound for 2 of the 15 analysed
# isolates (Clade3_2076 and para_2053) and reach 1.38 for a third. The pipeline's
# own model comparison agrees - LOO prefers Arrhenius over Sharpe-Schoolfield for
# respiration (elpd_diff -3.04 +/- 0.36). So:
#
#     E_r  <- Sharpe-Schoolfield        (the unimodal alternative, where supported)
#     E_K  <- Boltzmann-Arrhenius, full measured range
#
# Both are scale-free. The pipeline already fits SS per isolate to growth_fgC_h
# = r * q * 60; q and 60 are constant within an isolate, so they shift lnB0 and
# leave E untouched - that E IS the activation energy of r, reused rather than
# refitted.
hdr("B.1  Functional form: SS is not identifiable for K")

ss_K <- readr::read_csv(file.path(tables_dir, "sharpe_schoolfield_K_O2_rate_coefs.csv"),
                        show_col_types = FALSE) %>%
  dplyr::filter(parameter == "E") %>%
  dplyr::transmute(OTU = as.integer(OTU), otu_name, E_ss_K = Estimate, r2 = r_squared) %>%
  dplyr::inner_join(dplyr::distinct(d, OTU, group), by = "OTU")
message(sprintf("  SS-on-K isolates at or above E = 1.3 eV : %d of %d",
                sum(ss_K$E_ss_K >= 1.3), nrow(ss_K)))
message(sprintf("  SS-on-K isolates pinned at the E = 5 bound: %s",
                paste(ss_K$otu_name[ss_K$E_ss_K >= 4.999], collapse = ", ")))
loo <- readr::read_csv(file.path(tables_dir, "bayes_resp_model_comparison_loo.csv"),
                       show_col_types = FALSE)
message(sprintf("  pipeline LOO on respiration: best = %s, SS elpd_diff = %.2f +/- %.2f",
                loo$model[1], loo$elpd_diff[2], loo$se_diff[2]))
readr::write_csv(ss_K, file.path(C11_TAB, "partB_ss_on_K_unstable.csv"))

# =============================================================================
# B.2  E_r, E_K and E_K - E_r, per taxon
# =============================================================================
hdr("B.2  Activation energies (scale-free)")

# E_r per isolate: the pipeline's own SS fit to growth (q cancels in E)
E_r_iso <- readr::read_csv(file.path(tables_dir, "sharpe_schoolfield_growth_fgC_h_coefs.csv"),
                           show_col_types = FALSE) %>%
  dplyr::filter(parameter == "E") %>%
  dplyr::transmute(OTU = as.integer(OTU), E_r = Estimate, r2_r = r_squared)

# E_K per isolate: Boltzmann-Arrhenius on K over the full measured range
E_K_iso <- d %>% dplyr::group_by(OTU) %>%
  dplyr::group_modify(function(df, key) {
    m <- stats::lm(log_K ~ x, data = df)
    tibble::tibble(E_K = unname(coef(m)[2]),
                   E_K_se = unname(summary(m)$coefficients[2, 2]),
                   r2_K = summary(m)$r.squared, n = nrow(df))
  }) %>% dplyr::ungroup()

iso <- dplyr::distinct(d, OTU, otu_label, group) %>%
  dplyr::inner_join(E_r_iso, by = "OTU") %>%
  dplyr::inner_join(E_K_iso, by = "OTU") %>%
  dplyr::mutate(dE = E_K - E_r)          # negative => growth more T-sensitive
readr::write_csv(iso, file.path(C11_TAB, "partB_E_isolate.csv"))

E_tax <- iso %>% dplyr::group_by(group) %>%
  dplyr::group_modify(function(df, key) {
    bs <- replicate(NBOOT, {
      i <- sample.int(nrow(df), nrow(df), replace = TRUE)
      c(mean(df$E_r[i]), mean(df$E_K[i]), mean(df$dE[i]))
    })
    q <- function(v, p) unname(stats::quantile(v, p))
    tibble::tibble(
      n_isolates = nrow(df),
      E_r = mean(df$E_r), E_r_lo = q(bs[1, ], .025), E_r_hi = q(bs[1, ], .975),
      E_K = mean(df$E_K), E_K_lo = q(bs[2, ], .025), E_K_hi = q(bs[2, ], .975),
      dE  = mean(df$dE),  dE_lo  = q(bs[3, ], .025), dE_hi  = q(bs[3, ], .975),
      P_growth_more_sensitive = mean(bs[3, ] < 0))
  }) %>% dplyr::ungroup()
print(as.data.frame(E_tax %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
readr::write_csv(E_tax, file.path(C11_TAB, "partB_activation_energies.csv"))

# Cross-check E_K against Ilgaz's own Term-free arm (18_n0_treatment_panel.R).
# Term-free sets N0 = N_inoc, so respiration is proportional to K and that arm's
# hierarchical Bayesian E IS the scale-free E_K. Reused, not reimplemented.
free_arm <- readr::read_csv(file.path(tables_dir, "bayes_resp_arr_E_three_treatments.csv"),
                            show_col_types = FALSE) %>%
  dplyr::filter(treatment == "Term-free") %>%
  dplyr::transmute(group = clade, E_K_termfree = E,
                   E_K_termfree_lo = lwr, E_K_termfree_hi = upr)
xchk <- E_tax %>% dplyr::select(group, E_K) %>%
  dplyr::inner_join(free_arm, by = "group") %>%
  dplyr::mutate(diff = E_K - E_K_termfree)
message("\n  cross-check of E_K against 18's Term-free Bayesian arm:")
print(as.data.frame(xchk %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
message(sprintf("  max |difference| = %.4f eV", max(abs(xchk$diff))))
readr::write_csv(xchk, file.path(C11_TAB, "partB_EK_crosscheck.csv"))

# =============================================================================
# B.3  The model-free CUE optimum
# =============================================================================
# argmax CUE == argmin B_free, and B_free = K/r needs no conversion constant.
#
# MODEL-FREE means exactly that: the taxon-level B at each measured temperature
# is the geometric mean over replicate series (B is a ratio of positive rates and
# is close to log-normal), and the optimum is the argmin over the measured grid.
# No functional form is imposed on either r or K. A three-point parabola in
# log B through the grid minimum and its neighbours gives a sub-grid refinement;
# where the minimum sits on an edge the grid value stands.
hdr("B.3  Model-free T_opt(CUE) from argmin K/r")

grid_argmin <- function(TT, logB) {
  i <- which.min(logB)
  gridT <- TT[i]; refined <- gridT
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
boot_opt <- function(df, col = "B_free") {
  keys <- split(seq_len(nrow(df)), df$T)
  vapply(seq_len(NBOOT), function(b) {
    idx <- unlist(lapply(keys, function(ii) sample(ii, length(ii), replace = TRUE)))
    cv <- taxon_curve(df[idx, ], col)
    grid_argmin(cv$T, cv$logB)[["refined"]]
  }, numeric(1))
}

opt_rows <- list(); curves <- list()
for (g in C11_TAXA) {
  df <- dplyr::filter(d, group == g)
  cv <- taxon_curve(df); a <- grid_argmin(cv$T, cv$logB); bb <- boot_opt(df)
  opt_rows[[g]] <- tibble::tibble(
    group = g, label = c11_label(g), n_series = nrow(df), n_T = nrow(cv),
    T_opt_grid = a[["grid"]], T_opt = a[["refined"]],
    lo = unname(stats::quantile(bb, .025)), hi = unname(stats::quantile(bb, .975)),
    margin_to_37 = T_BODY - a[["refined"]], P_below_37 = mean(bb < T_BODY),
    fold_change_B = exp(max(cv$logB) - min(cv$logB)),
    T_min_meas = min(cv$T), T_max_meas = max(cv$T))
  curves[[g]] <- cv %>% dplyr::mutate(group = g, B = exp(logB))
}
opt <- dplyr::bind_rows(opt_rows)
print(as.data.frame(opt %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
readr::write_csv(opt, file.path(C11_TAB, "partB_modelfree_optimum.csv"))
readr::write_csv(dplyr::bind_rows(curves), file.path(C11_TAB, "partB_balance_curves.csv"))

hdr("B.4  The headline test")
for (i in seq_len(nrow(opt))) with(opt[i, ], message(sprintf(
  "  %-28s T_opt = %5.2f C  [%5.2f, %5.2f]   margin to 37 C = %5.2f C   P(<37) = %.4f",
  label, T_opt, lo, hi, margin_to_37, P_below_37)))
message(sprintf("\n  all five below 37 C: %s   (min margin %.2f C, min P = %.4f)",
                all(opt$T_opt < T_BODY), min(opt$margin_to_37), min(opt$P_below_37)))
message(sprintf("  fold-change in K/r across the measured range: %.2f - %.2f",
                min(opt$fold_change_B), max(opt$fold_change_B)))

hdr("PART B done")
