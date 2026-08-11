# =============================================================================
# C11 PART C - compare against the published optimum, and decompose the gap
# =============================================================================
# Quantities here DO depend on N0, so by the PART B/C split they belong in this
# file and are labelled as such. The three N0 treatments are Ilgaz's own, reused
# from 08_temperature_equilibration_sensitivity.R and 18_n0_treatment_panel.R
# rather than reimplemented.
# =============================================================================

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
       value = TRUE)[1])), "00_common.R"))

set.seed(11)
NBOOT <- 2000

d <- readr::read_csv(file.path(C11_TAB, "partA_series.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(group = factor(group, levels = C11_TAXA))
hdr <- function(x) message("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74))

# ---- the published values ---------------------------------------------------
pub <- readr::read_csv(file.path(tables_dir, "fig_values.csv"), show_col_types = FALSE) %>%
  dplyr::transmute(group = Group, T_pub = Tcue, T_pub_lo = Tcue_lo, T_pub_hi = Tcue_hi,
                   P_pub_below_37 = P_below_body, Topt_growth = Topt) %>%
  dplyr::filter(group %in% C11_TAXA)

# ---- ARM 2: the equilibration-corrected CUE, straight from script 08 ---------
# 08 writes a per-curve corrected respiration (resp_corr) and the corrected CUE
# (cue_corr) that follows from it. 18 uses resp_corr for exactly this treatment.
# Nothing is recomputed here; cue_corr is read and its argmax taken.
pc <- readr::read_csv(C11_INPUTS[["percurve"]], show_col_types = FALSE) %>%
  dplyr::transmute(T = as.numeric(T), OTU = as.integer(OTU),
                   Replicate = toupper(as.character(Replicate)),
                   cue_corr = as.numeric(cue_corr), cue_pipe = as.numeric(cue_pipe))
d <- d %>% dplyr::left_join(pc, by = c("T", "OTU", "Replicate"))

# sanity: 08's cue_pipe must be the pipeline CUE we already have
.chk <- d %>% dplyr::filter(is.finite(cue_pipe), is.finite(CUE))
message(sprintf("  08's cue_pipe vs the derived table's CUE: max rel. diff %.3e over %d series",
                max(abs(.chk$cue_pipe - .chk$CUE) / .chk$CUE), nrow(.chk)))

# =============================================================================
# The four estimates of T_opt(CUE)
# =============================================================================
# All three model-free estimates use the SAME argmin machinery; they differ only
# in which N0 treatment defines the quantity being optimised.
#
#   TERM-FREE       argmin K/r                    no N0 at all      (PART B)
#   WITH TERM       argmin K/(r e^{r delta})      published N0
#   EQUILIBRATION   argmax cue_corr               08's corrected N0
#   PUBLISHED       fitted Bayesian CUE curve     published N0 + functional form
grid_arg <- function(TT, y, minimise = TRUE) {
  z <- if (minimise) y else -y
  i <- which.min(z); out <- TT[i]
  if (i > 1 && i < length(TT)) {
    x <- TT[(i-1):(i+1)]; v <- z[(i-1):(i+1)]
    cf <- stats::coef(stats::lm(v ~ poly(x, 2, raw = TRUE)))
    if (!is.na(cf[3]) && cf[3] > 0) {
      p <- -cf[2] / (2 * cf[3])
      if (is.finite(p) && p >= min(x) && p <= max(x)) out <- unname(p)
    }
  }
  out
}
curve_of <- function(df, col) df %>%
  dplyr::filter(is.finite(.data[[col]]), .data[[col]] > 0) %>%
  dplyr::group_by(T) %>%
  dplyr::summarise(y = mean(log(.data[[col]])), .groups = "drop") %>% dplyr::arrange(T)

opt_of <- function(df, col, minimise) {
  cv <- curve_of(df, col)
  if (nrow(cv) < 3) return(c(est = NA_real_, lo = NA_real_, hi = NA_real_))
  est <- grid_arg(cv$T, cv$y, minimise)
  sub <- dplyr::filter(df, is.finite(.data[[col]]), .data[[col]] > 0)
  keys <- split(seq_len(nrow(sub)), sub$T)
  bb <- vapply(seq_len(NBOOT), function(b) {
    idx <- unlist(lapply(keys, function(ii) sample(ii, length(ii), replace = TRUE)))
    c2 <- curve_of(sub[idx, ], col)
    if (nrow(c2) < 3) return(NA_real_)
    grid_arg(c2$T, c2$y, minimise)
  }, numeric(1))
  c(est = est, lo = unname(stats::quantile(bb, .025, na.rm = TRUE)),
    hi = unname(stats::quantile(bb, .975, na.rm = TRUE)),
    p37 = mean(bb < T_BODY, na.rm = TRUE))
}

hdr("C.1  Four estimates of T_opt(CUE) per taxon")
rows <- list()
for (g in C11_TAXA) {
  df <- dplyr::filter(d, group == g)
  a_free  <- opt_of(df, "B_free", TRUE)
  a_bp    <- opt_of(df, "B_bp",   TRUE)
  a_eq    <- opt_of(df, "cue_corr", FALSE)
  rows[[g]] <- tibble::tibble(
    group = g, label = c11_label(g),
    T_free = a_free[["est"]], T_free_lo = a_free[["lo"]], T_free_hi = a_free[["hi"]],
    P37_free = a_free[["p37"]],
    T_bp   = a_bp[["est"]],   T_bp_lo   = a_bp[["lo"]],   T_bp_hi   = a_bp[["hi"]],
    P37_bp = a_bp[["p37"]],
    T_eq   = a_eq[["est"]],   T_eq_lo   = a_eq[["lo"]],   T_eq_hi   = a_eq[["hi"]],
    P37_eq = a_eq[["p37"]])
}
est <- dplyr::bind_rows(rows) %>% dplyr::inner_join(pub, by = "group")

# =============================================================================
# The decomposition. Additive and exhaustive by construction:
#     T_free - T_pub  =  (T_free - T_bp)  +  (T_bp - T_pub)
#                         \___ N0 ___/       \_ functional form _/
# The solubility/unit term that D8 needed is ZERO here: PART A established that
# K is already volumetric because the fit is to raw oxygen.
# =============================================================================
est <- est %>% dplyr::mutate(
  gap_total       = T_free - T_pub,
  gap_N0          = T_free - T_bp,
  gap_form        = T_bp - T_pub,
  gap_solubility  = 0,
  gap_N0_equil    = T_free - T_eq,     # the physically-motivated middle case
  margin_free     = T_BODY - T_free,
  margin_pub      = T_BODY - T_pub,
  margin_widens   = margin_free - margin_pub)

print(as.data.frame(est %>% dplyr::select(label, T_free, T_free_lo, T_free_hi, P37_free,
                                          T_bp, P37_bp, T_eq, T_eq_lo, T_eq_hi, P37_eq,
                                          T_pub, P_pub_below_37) %>%
                      dplyr::mutate(across(where(is.numeric), ~round(.x, 3)))), row.names = FALSE)
message("\n  the tightest case, under EVERY treatment:")
for (i in seq_len(nrow(est))) with(est[i, ], {
  m <- c(free = T_BODY - T_free, with_term = T_BODY - T_bp, equil = T_BODY - T_eq,
         published = T_BODY - T_pub)
  message(sprintf("    %-28s min margin %+6.2f C  (%s)", label, min(m), names(m)[which.min(m)]))
})
hdr("C.2  Gap decomposition (degrees C)")
print(as.data.frame(est %>% dplyr::select(label, gap_total, gap_N0, gap_form, gap_solubility) %>%
                      dplyr::mutate(across(where(is.numeric), ~round(.x, 2)))), row.names = FALSE)
message(sprintf("\n  mean total gap        %+.2f C   (model-free minus published)", mean(est$gap_total)))
message(sprintf("  mean N0 component     %+.2f C   (%.0f%% of the total)",
                mean(est$gap_N0), 100 * mean(est$gap_N0) / mean(est$gap_total)))
message(sprintf("  mean form component   %+.2f C   (%.0f%% of the total)",
                mean(est$gap_form), 100 * mean(est$gap_form) / mean(est$gap_total)))
message(sprintf("  solubility component  %+.2f C   (K is volumetric; see PART A)", 0))

hdr("C.3  Does the claim strengthen or weaken?")
for (i in seq_len(nrow(est))) with(est[i, ], message(sprintf(
  "  %-28s margin: published %4.2f C -> model-free %5.2f C   (%+.2f C)",
  label, margin_pub, margin_free, margin_widens)))
message(sprintf("\n  margin widens for %d of %d taxa; minimum model-free margin %.2f C",
                sum(est$margin_widens > 0), nrow(est), min(est$margin_free)))

readr::write_csv(est, file.path(C11_TAB, "partC_gap_decomposition.csv"))

# =============================================================================
# C.4  Where the fitted curve puts the optimum if the DATA are model-free
# =============================================================================
# Splits "functional form" one level further: how much of gap_form is the
# Sharpe-Schoolfield + Arrhenius pair, as opposed to the Bayesian hierarchy?
# Uses the pipeline's own per-isolate frequentist SS(growth) and BA(K) fits to
# build a parametric K/r curve and take its analytic argmin, with the SAME N0
# treatment as term-free.
hdr("C.4  Parametric vs model-free argmin, under the term-free treatment")
ssg <- readr::read_csv(file.path(tables_dir, "sharpe_schoolfield_growth_fgC_h_coefs.csv"),
                       show_col_types = FALSE) %>%
  dplyr::select(parameter, Estimate, OTU) %>%
  tidyr::pivot_wider(names_from = parameter, values_from = Estimate) %>%
  dplyr::mutate(OTU = as.integer(OTU))
bak <- readr::read_csv(file.path(C11_TAB, "partB_E_isolate.csv"), show_col_types = FALSE)

TG <- seq(22, 44, by = 0.05)
par_rows <- list()
for (g in C11_TAXA) {
  ids <- dplyr::filter(bak, group == g)$OTU
  s <- dplyr::filter(ssg, OTU %in% ids); b <- dplyr::filter(bak, OTU %in% ids)
  if (!all(c("lnB0", "E", "Eh", "Th") %in% names(s))) { next }
  # mean over isolates of log r(T) from SS, and of log K(T) from BA
  logr <- rowMeans(vapply(seq_len(nrow(s)), function(i)
    pred_ss_log(TG + 273.15, s$lnB0[i], s$E[i], s$Eh[i], s$Th[i]), numeric(length(TG))))
  logK <- rowMeans(vapply(seq_len(nrow(b)), function(i)
    b$E_K[i] * c11_boltz(TG), numeric(length(TG))))
  lb <- logK - logr
  par_rows[[g]] <- tibble::tibble(group = g, T_param = TG[which.min(lb)])
}
par_est <- dplyr::bind_rows(par_rows)
est2 <- est %>% dplyr::inner_join(par_est, by = "group") %>%
  dplyr::mutate(form_parametric = T_param - T_free,
                form_bayes_hier = T_bp - T_param - gap_N0)
print(as.data.frame(est2 %>% dplyr::select(label, T_free, T_param, T_bp, T_pub,
                                           form_parametric) %>%
                      dplyr::mutate(across(where(is.numeric), ~round(.x, 2)))), row.names = FALSE)
message(sprintf("\n  mean shift from imposing SS(growth) + BA(K) alone: %+.2f C",
                mean(est2$form_parametric)))
readr::write_csv(est2, file.path(C11_TAB, "partC_form_split.csv"))

hdr("PART C done")
