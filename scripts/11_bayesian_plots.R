# =============================================================================
# 11_bayesian_plots.R - Figures + contrasts from the Bayesian fits
# =============================================================================
# Reads the models fitted by 09_bayesian_models.R. No refitting.
#
# FIGURES
#   fig_bayes_growth_tpc_by_clade.png   6 panels (Clade1-4, glab, para). Within
#                                       each panel: the 3 isolates' data points +
#                                       their posterior curves, plus the CLADE
#                                       curve in black with a credible ribbon.
#   fig_bayes_resp_arrhenius.png        Respiration on a LOG scale vs temperature,
#                                       same 6-panel layout. (LOO says Arrhenius
#                                       beats Sharpe-Schoolfield: respiration keeps
#                                       rising and never peaks.)
#   fig_bayes_cue_by_clade.png          CUE = G/(G+R), DERIVED from the growth and
#                                       respiration posteriors - so it is bounded
#                                       in [0,1] by construction and carries proper
#                                       credible intervals. No separate CUE model.
#   fig_bayes_contrasts.png             Clade contrasts with 95% credible intervals.
#
# TABLES
#   bayes_clade_params.csv        E, Th, Topt per clade (growth) + E (respiration)
#   bayes_clade_contrasts.csv     pairwise differences + posterior probabilities
#   bayes_E_resp_minus_growth.csv E_R - E_G per clade: the single number that
#                                 explains WHY CUE collapses with warming.
#
# NOTE on combining the two posteriors: growth and respiration were fitted as
# separate models, so their draws are independent. Pairing draw i with draw i to
# build CUE and E_R - E_G is standard posterior propagation and is fine - it just
# assumes no residual correlation between the two fits.
# =============================================================================


# =============================================================================
# 0) Config + packages
# =============================================================================

# Packaging: resolve the script directory from --file= FIRST. Under
# `Rscript scripts/<this>.R` neither rstudioapi nor sys.frame(1)$ofile resolves,
# so this fell back to getwd() and then died on "cannot open file .../config.R".
# Sourcing it from run_all.R was unaffected, which is why the bug stayed hidden.
.this_dir <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) return(dirname(normalizePath(
    # R replaces every space in --file= with "~+~", so a project path that
    # contains a space comes back mangled and source() then fails on a path
    # that does not exist. Un-mangle it before normalising.
    gsub("~+~", " ", sub("^--file=", "", fa[1]), fixed = TRUE),
    mustWork = FALSE)))
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
})

source(file.path(.this_dir, "config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(brms)
  library(posterior)
})

TGRID_N <- 120          # temperature resolution of the fitted curves
CI      <- c(0.025, 0.975)


# =============================================================================
# 1) Load
# =============================================================================

req <- function(p) if (!file.exists(p)) stop("Not found: ", p, "\nRun 09_bayesian_models.R first.") else p

growth_dat <- readRDS(req(file.path(models_dir, "bayes_data_growth.rds")))
resp_dat   <- readRDS(req(file.path(models_dir, "bayes_data_resp.rds")))
fit_g      <- readRDS(req(file.path(models_dir, "bayes_growth_ss.rds")))
fit_r      <- readRDS(req(file.path(models_dir, "bayes_resp_arr.rds")))

groups   <- levels(factor(growth_dat$Group))
iso_tbl  <- growth_dat %>% dplyr::distinct(Isolate, Group)
Trange   <- range(c(growth_dat$T, resp_dat$T))
Tseq     <- seq(Trange[1], Trange[2], length.out = TGRID_N)

# Consistent colours for the 3 isolates within a panel.
iso_cols <- c("#1b9e77", "#d95f02", "#7570b3")

# Palette sized to the isolates ACTUALLY PRESENT, not to an assumed 3-per-group.
# `rep(iso_cols, length(groups))` breaks the moment the isolate count changes -
# e.g. excluding glab left 15 colours for a layer that still carried 18 isolates:
#   "Insufficient values in manual scale. 18 needed but only 15 provided."
# Recycling to the real number of levels cannot go out of sync.
iso_palette <- function(...) {
  lv <- unique(unlist(lapply(list(...), function(d) {
    if (is.null(d) || !"Isolate" %in% names(d)) return(NULL)
    as.character(d$Isolate)
  })))
  lv <- lv[!is.na(lv)]
  stats::setNames(rep(iso_cols, length.out = max(length(lv), 1L)), sort(lv))
}


# =============================================================================
# 2) Prediction grids + posterior draws
# =============================================================================

# (a) per-ISOLATE curves (uses each isolate's random effects)
grid_iso <- tidyr::expand_grid(iso_tbl, T = Tseq) %>%
  dplyr::mutate(TK = T + 273.15,
                boltz_shift = (1 / (k_B * T_ref)) - (1 / (k_B * TK)),
                se_y = 0)          # se() term must exist; 0 = predict the mean

# (b) CLADE-level curves (population level: re_formula = NA, isolate effects off)
grid_grp <- tidyr::expand_grid(Group = groups, T = Tseq) %>%
  dplyr::mutate(Isolate = iso_tbl$Isolate[1],   # placeholder, ignored when re_formula = NA
                TK = T + 273.15,
                boltz_shift = (1 / (k_B * T_ref)) - (1 / (k_B * TK)),
                se_y = 0)

# posterior_epred returns draws x rows, on the LOG-rate scale (the response is
# log(rate)), so exp() converts to natural units.
ep <- function(fit, nd, re) brms::posterior_epred(fit, newdata = nd, re_formula = re)

epred_g_iso <- ep(fit_g, grid_iso, NULL)   # growth, per isolate
epred_g_grp <- ep(fit_g, grid_grp, NA)     # growth, per clade
epred_r_iso <- ep(fit_r, grid_iso, NULL)   # respiration, per isolate
epred_r_grp <- ep(fit_r, grid_grp, NA)     # respiration, per clade

summ <- function(draws, nd, val = "y") {
  q <- t(apply(exp(draws), 2, function(x)
    c(med = stats::median(x), lo = stats::quantile(x, CI[1]), hi = stats::quantile(x, CI[2]))))
  dplyr::bind_cols(nd, tibble::as_tibble(q, .name_repair = ~c("med", "lo", "hi")))
}

g_iso <- summ(epred_g_iso, grid_iso)
g_grp <- summ(epred_g_grp, grid_grp)
r_iso <- summ(epred_r_iso, grid_iso)
r_grp <- summ(epred_r_grp, grid_grp)


# =============================================================================
# 3) FIGURE: growth Sharpe-Schoolfield, 6 panels, 3 isolates each
# =============================================================================

p_growth <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(data = g_grp, ggplot2::aes(T, ymin = lo, ymax = hi),
                       fill = "grey70", alpha = 0.35) +
  ggplot2::geom_point(data = growth_dat,
                      ggplot2::aes(T, y_raw, colour = Isolate),
                      size = 1.5, alpha = 0.55) +
  ggplot2::geom_line(data = g_iso,
                     ggplot2::aes(T, med, colour = Isolate), linewidth = 0.8) +
  ggplot2::geom_line(data = g_grp, ggplot2::aes(T, med),
                     colour = "black", linewidth = 1.3) +
  ggplot2::scale_colour_manual(values = iso_palette(growth_dat, resp_dat),
                               na.value = "grey60") +
  ggplot2::facet_wrap(~ Group, scales = "free_y", nrow = 2) +
  ggplot2::labs(
    title = "Growth thermal performance (Bayesian Sharpe-Schoolfield)",
    subtitle = "black = clade curve (95% ribbon) | coloured = the 3 isolates in that clade",
    x = "Temperature (°C)", y = expression("Growth (fg C h"^{-1}*")")) +
  ggplot2::theme_bw(12) +
  ggplot2::theme(legend.position = "bottom",
                 plot.title = ggplot2::element_text(face = "bold"),
                 panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave(file.path(figures_dir, "fig_bayes_growth_tpc_by_clade.png"),
                p_growth, width = 13, height = 7.5, dpi = 300)


# =============================================================================
# 4) FIGURE: respiration Arrhenius, LOG scale
# =============================================================================

# ---- temperature-equilibration band (from 08_..sensitivity.R, runs before this) --
# The credible ribbon reflects replicate + fit uncertainty but treats N0 as
# exact. The onset-temperature effect adds a further, temperature-structured
# uncertainty. Here we widen the ribbon by that component (quadrature) and draw
# it as a lighter outer band UNDER the credible interval. If the sensitivity
# table is absent, the figure falls back to the plain credible ribbon.
.teq_file <- file.path(tables_dir, "temperature_equilibration_sensitivity.csv")
have_teq  <- file.exists(.teq_file)
if (have_teq) {
  .teq <- readr::read_csv(.teq_file, show_col_types = FALSE)
  rf <- stats::approx(.teq$T, abs(.teq$resp_swing_med) / 100, xout = r_grp$T, rule = 2)$y
  hw <- sqrt(((log(r_grp$hi) - log(r_grp$lo)) / 2)^2 + log1p(rf)^2)
  r_grp$lo_wide <- exp(log(r_grp$med) - hw)
  r_grp$hi_wide <- exp(log(r_grp$med) + hw)
} else {
  message("sensitivity table not found - drawing plain credible ribbon (run 08 first for the temperature band).")
}
resp_wide <- if (have_teq) ggplot2::geom_ribbon(
  data = r_grp, ggplot2::aes(T, ymin = lo_wide, ymax = hi_wide),
  fill = "grey80", alpha = 0.10, colour = "grey30", linetype = "22", linewidth = 0.3) else NULL

p_resp <- ggplot2::ggplot() +
  resp_wide +
  ggplot2::geom_ribbon(data = r_grp, ggplot2::aes(T, ymin = lo, ymax = hi),
                       fill = "grey70", alpha = 0.35) +
  ggplot2::geom_point(data = resp_dat,
                      ggplot2::aes(T, y_raw, colour = Isolate),
                      size = 1.5, alpha = 0.55) +
  ggplot2::geom_line(data = r_iso,
                     ggplot2::aes(T, med, colour = Isolate), linewidth = 0.8) +
  ggplot2::geom_line(data = r_grp, ggplot2::aes(T, med),
                     colour = "black", linewidth = 1.3) +
  ggplot2::scale_colour_manual(values = iso_palette(growth_dat, resp_dat),
                               na.value = "grey60") +
  ggplot2::scale_y_log10() +
  ggplot2::facet_wrap(~ Group, nrow = 2) +
  ggplot2::labs(
    title = "Respiration thermal response (Bayesian Arrhenius, log scale)",
    subtitle = "Inner grey = credible interval; dashed outer envelope adds temperature-equilibration uncertainty (widens at the extremes)",
    x = "Temperature (°C)",
    y = expression("Respiration (fg C h"^{-1}*", log scale)")) +
  ggplot2::theme_bw(12) +
  ggplot2::theme(legend.position = "bottom",
                 plot.title = ggplot2::element_text(face = "bold"),
                 panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave(file.path(figures_dir, "fig_bayes_resp_arrhenius.png"),
                p_resp, width = 13, height = 7.5, dpi = 300)


# =============================================================================
# 5) FIGURE: CUE, DERIVED from the two posteriors
# =============================================================================
# CUE = G / (G + R), computed draw-by-draw, so it is bounded in [0,1] and its
# credible interval inherits the uncertainty of BOTH fits. No separate CUE model,
# no Arrhenius-on-a-ratio (CUE is not a rate - Arrhenius does not apply to it).

G <- exp(epred_g_grp); R <- exp(epred_r_grp)
CUE_draws <- G / (G + R)

cue_grp <- dplyr::bind_cols(
  grid_grp,
  tibble::as_tibble(t(apply(CUE_draws, 2, function(x)
    c(stats::median(x), stats::quantile(x, CI[1]), stats::quantile(x, CI[2])))),
    .name_repair = ~c("med", "lo", "hi")))

# widen the CUE ribbon by the temperature-equilibration component (quadrature)
if (have_teq) {
  cf <- stats::approx(.teq$T, abs(.teq$cue_swing_med) / 100, xout = cue_grp$T, rule = 2)$y
  hw <- sqrt(((cue_grp$hi - cue_grp$lo) / 2)^2 + (cf * cue_grp$med)^2)
  cue_grp$lo_wide <- pmax(0, cue_grp$med - hw)
  cue_grp$hi_wide <- pmin(1, cue_grp$med + hw)
}
cue_wide <- if (have_teq) ggplot2::geom_ribbon(
  data = cue_grp, ggplot2::aes(T, ymin = lo_wide, ymax = hi_wide),
  fill = "#2c7fb8", alpha = 0.07, colour = "#08306b", linetype = "22", linewidth = 0.3) else NULL

p_cue <- ggplot2::ggplot(cue_grp, ggplot2::aes(T, med)) +
  cue_wide +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "#2c7fb8", alpha = 0.25) +
  ggplot2::geom_line(colour = "#2c7fb8", linewidth = 1.2) +
  ggplot2::facet_wrap(~ Group, nrow = 2) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Carbon-use efficiency, derived from the growth and respiration posteriors",
    subtitle = "Inner band = credible interval; dashed outer envelope adds temperature-equilibration uncertainty (widens toward the extremes)",
    x = "Temperature (°C)", y = "CUE") +
  ggplot2::theme_bw(12) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                 panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave(file.path(figures_dir, "fig_bayes_cue_by_clade.png"),
                p_cue, width = 13, height = 7.5, dpi = 300)

readr::write_csv(cue_grp %>% dplyr::select(Group, T, CUE = med, lo, hi),
                 file.path(tables_dir, "bayes_cue_by_clade.csv"))


# =============================================================================
# 6) Clade parameters + contrasts
# =============================================================================

dg <- posterior::as_draws_df(fit_g)
dr <- posterior::as_draws_df(fit_r)

pull_grp <- function(dd, stem) {
  cols <- grep(paste0("^b_", stem, "_Group"), names(dd), value = TRUE)
  out <- as.matrix(dd[, cols, drop = FALSE])
  colnames(out) <- sub(paste0("^b_", stem, "_Group"), "", cols)
  out
}

E_g  <- pull_grp(dg, "E")     # growth activation energy, draws x clade
Th_g <- pull_grp(dg, "Th")
E_r  <- pull_grp(dr, "E")     # respiration activation energy

# Topt per draw: peak of the fitted SS curve, found on the temperature grid.
Topt_draws <- matrix(NA_real_, nrow = nrow(G), ncol = length(groups),
                     dimnames = list(NULL, groups))
for (j in seq_along(groups)) {
  idx <- which(grid_grp$Group == groups[j])
  sub <- G[, idx, drop = FALSE]
  Topt_draws[, j] <- grid_grp$T[idx][apply(sub, 1, which.max)]
}

qsum <- function(x) c(mean = mean(x), median = stats::median(x),
                      lo = unname(stats::quantile(x, CI[1])),
                      hi = unname(stats::quantile(x, CI[2])))

clade_params <- dplyr::bind_rows(
  lapply(groups, function(gname) {
    tibble::tibble(
      Group = gname,
      growth_E_eV   = qsum(E_g[, gname])[["median"]],
      growth_E_lo   = qsum(E_g[, gname])[["lo"]],
      growth_E_hi   = qsum(E_g[, gname])[["hi"]],
      growth_Th_K   = qsum(Th_g[, gname])[["median"]],
      growth_Topt_C = qsum(Topt_draws[, gname])[["median"]],
      growth_Topt_lo = qsum(Topt_draws[, gname])[["lo"]],
      growth_Topt_hi = qsum(Topt_draws[, gname])[["hi"]],
      resp_E_eV     = qsum(E_r[, gname])[["median"]],
      resp_E_lo     = qsum(E_r[, gname])[["lo"]],
      resp_E_hi     = qsum(E_r[, gname])[["hi"]])
  }))
readr::write_csv(clade_params, file.path(tables_dir, "bayes_clade_params.csv"))

# ---- E_R - E_G : why CUE collapses with warming -----------------------------
dE <- E_r[seq_len(min(nrow(E_r), nrow(E_g))), , drop = FALSE] -
  E_g[seq_len(min(nrow(E_r), nrow(E_g))), , drop = FALSE]

e_diff <- dplyr::bind_rows(lapply(groups, function(gname) {
  x <- dE[, gname]
  tibble::tibble(Group = gname,
                 E_resp_minus_growth = stats::median(x),
                 lo = unname(stats::quantile(x, CI[1])),
                 hi = unname(stats::quantile(x, CI[2])),
                 P_resp_more_sensitive = mean(x > 0))
}))
readr::write_csv(e_diff, file.path(tables_dir, "bayes_E_resp_minus_growth.csv"))

# ---- pairwise clade contrasts ------------------------------------------------
contrast <- function(mat, label) {
  out <- list()
  for (i in seq_along(groups)) for (j in seq_along(groups)) {
    if (i >= j) next
    a <- groups[i]; b <- groups[j]
    d <- mat[, a] - mat[, b]
    out[[length(out) + 1]] <- tibble::tibble(
      parameter = label, group_a = a, group_b = b,
      diff = stats::median(d),
      lo = unname(stats::quantile(d, CI[1])),
      hi = unname(stats::quantile(d, CI[2])),
      P_a_greater = mean(d > 0))
  }
  dplyr::bind_rows(out)
}

contrasts_tbl <- dplyr::bind_rows(
  contrast(E_g, "growth_E_eV"),
  contrast(Th_g, "growth_Th_K"),
  contrast(Topt_draws, "growth_Topt_C"),
  contrast(E_r, "resp_E_eV"))
readr::write_csv(contrasts_tbl, file.path(tables_dir, "bayes_clade_contrasts.csv"))


# =============================================================================
# 7) FIGURE: contrasts (E growth vs respiration, per clade)
# =============================================================================

e_long <- dplyr::bind_rows(
  clade_params %>% dplyr::transmute(Group, rate = "growth",
                                    E = growth_E_eV, lo = growth_E_lo, hi = growth_E_hi),
  clade_params %>% dplyr::transmute(Group, rate = "respiration",
                                    E = resp_E_eV, lo = resp_E_lo, hi = resp_E_hi))

p_e <- ggplot2::ggplot(e_long, ggplot2::aes(E, Group, colour = rate)) +
  ggplot2::geom_pointrange(ggplot2::aes(xmin = lo, xmax = hi),
                           position = ggplot2::position_dodge(width = 0.5), size = 0.6) +
  ggplot2::scale_colour_manual(values = c(growth = "#1b9e77", respiration = "#d95f02")) +
  ggplot2::labs(
    title = "Activation energy: growth vs respiration",
    subtitle = "Growth is ~2x more temperature-sensitive than respiration - this is why CUE collapses on warming",
    x = "Activation energy E (eV), 95% credible interval", y = NULL, colour = NULL) +
  ggplot2::theme_bw(12) +
  ggplot2::theme(legend.position = "bottom",
                 plot.title = ggplot2::element_text(face = "bold"))

ggplot2::ggsave(file.path(figures_dir, "fig_bayes_contrasts.png"),
                p_e, width = 8, height = 5, dpi = 300)


# =============================================================================
# 8) CUE per INDIVIDUAL ISOLATE  (3 isolates per clade, + observed points)
# =============================================================================
# The clade-level CUE hides the isolate-to-isolate spread. Here each isolate gets
# its own curve, derived draw-by-draw from ITS OWN growth and respiration
# posteriors (random effects included), with the observed per-curve CUE overlaid.

G_i <- exp(epred_g_iso); R_i <- exp(epred_r_iso)
CUE_iso <- G_i / (G_i + R_i)

cue_iso <- dplyr::bind_cols(
  grid_iso,
  tibble::as_tibble(t(apply(CUE_iso, 2, function(x)
    c(stats::median(x), stats::quantile(x, CI[1]), stats::quantile(x, CI[2])))),
    .name_repair = ~c("med", "lo", "hi")))

# Observed CUE per curve (straight from 06), so you can see the raw scatter.
cue_obs <- NULL
if (file.exists(derived_csv)) {
  .d <- tryCatch(readr::read_csv(derived_csv, show_col_types = FALSE),
                 error = function(e) NULL)
  if (!is.null(.d) && "CUE" %in% names(.d)) {
    .nm <- readr::read_csv(file.path(tables_dir, "otu_names.csv"), show_col_types = FALSE) %>%
      dplyr::transmute(OTU = as.integer(OTU), Isolate = as.character(otu_name),
                       Group = as.character(group))
    cue_obs <- .d %>%
      dplyr::mutate(OTU = as.integer(OTU), T = as.numeric(T)) %>%
      dplyr::inner_join(.nm, by = "OTU") %>%
      dplyr::filter(is.finite(CUE), CUE >= 0, CUE <= 1) %>%
      dplyr::select(T, Isolate, Group, CUE)
    # This layer reads derived_csv DIRECTLY, so it bypasses the EXCLUDE_GROUPS
    # filter applied in 10. Without this it still carries glab, and the plot then
    # asks for 18 isolate colours from a 15-colour palette:
    #   "Insufficient values in manual scale. 18 needed but only 15 provided."
    # Any script that reads derived_csv directly must filter it directly too.
    cue_obs <- drop_excluded(cue_obs, "observed CUE")
  }
}

p_cue_iso <- ggplot2::ggplot()
if (!is.null(cue_obs) && nrow(cue_obs) > 0) {
  p_cue_iso <- p_cue_iso +
    ggplot2::geom_point(data = cue_obs,
                        ggplot2::aes(T, CUE, colour = Isolate),
                        size = 1.2, alpha = 0.35)
}
p_cue_iso <- p_cue_iso +
  ggplot2::geom_ribbon(data = cue_grp, ggplot2::aes(T, ymin = lo, ymax = hi),
                       fill = "grey70", alpha = 0.30) +
  ggplot2::geom_line(data = cue_iso,
                     ggplot2::aes(T, med, colour = Isolate), linewidth = 0.9) +
  ggplot2::geom_line(data = cue_grp, ggplot2::aes(T, med),
                     colour = "black", linewidth = 1.3) +
  ggplot2::scale_colour_manual(values = iso_palette(growth_dat, resp_dat),
                               na.value = "grey60") +
  ggplot2::facet_wrap(~ Group, nrow = 2) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Carbon-use efficiency per isolate",
    subtitle = "coloured = each isolate's posterior CUE | black = clade (95% ribbon) | points = observed per-curve CUE",
    x = "Temperature (°C)", y = "CUE", colour = NULL) +
  ggplot2::theme_bw(12) +
  ggplot2::theme(legend.position = "bottom",
                 plot.title = ggplot2::element_text(face = "bold"),
                 panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave(file.path(figures_dir, "fig_bayes_cue_by_isolate.png"),
                p_cue_iso, width = 13, height = 7.5, dpi = 300)

readr::write_csv(cue_iso %>% dplyr::select(Group, Isolate, T, CUE = med, lo, hi),
                 file.path(tables_dir, "bayes_cue_by_isolate.csv"))


# =============================================================================
# 9) POSTERIOR DISTRIBUTIONS  (half-eye) - the point of going Bayesian
# =============================================================================
# A pointrange hides the SHAPE of a posterior. If a clade's Th posterior is wide
# and skewed, that means the data cannot pin down its high-T decline - and you
# want to see that, not a tidy interval that implies false precision.

.gg <- requireNamespace("ggdist", quietly = TRUE)

draws_long <- dplyr::bind_rows(
  tibble::as_tibble(E_g)        %>% tidyr::pivot_longer(dplyr::everything(),
                    names_to = "Group", values_to = "value") %>%
    dplyr::mutate(parameter = "Growth: E (eV)"),
  tibble::as_tibble(Th_g)       %>% tidyr::pivot_longer(dplyr::everything(),
                    names_to = "Group", values_to = "value") %>%
    dplyr::mutate(parameter = "Growth: Th (K)"),
  tibble::as_tibble(Topt_draws) %>% tidyr::pivot_longer(dplyr::everything(),
                    names_to = "Group", values_to = "value") %>%
    dplyr::mutate(parameter = "Growth: Topt (°C)"),
  tibble::as_tibble(E_r)        %>% tidyr::pivot_longer(dplyr::everything(),
                    names_to = "Group", values_to = "value") %>%
    dplyr::mutate(parameter = "Respiration: E (eV)"))

p_post <- ggplot2::ggplot(draws_long, ggplot2::aes(value, Group, fill = Group))
p_post <- if (.gg) {
  p_post + ggdist::stat_halfeye(.width = c(0.5, 0.95), point_size = 1.6,
                                slab_alpha = 0.75)
} else {
  p_post + ggplot2::geom_violin(alpha = 0.75, colour = NA) +
    ggplot2::stat_summary(fun = stats::median, geom = "point", size = 1.8)
}
p_post <- p_post +
  ggplot2::facet_wrap(~ parameter, scales = "free_x", nrow = 2) +
  ggplot2::labs(
    title = "Posterior distributions by clade",
    subtitle = "the full posterior, not just an interval - width and skew tell you what the data could NOT pin down",
    x = NULL, y = NULL) +
  ggplot2::theme_bw(12) +
  ggplot2::theme(legend.position = "none",
                 plot.title = ggplot2::element_text(face = "bold"))

ggplot2::ggsave(file.path(figures_dir, "fig_bayes_posteriors.png"),
                p_post, width = 12, height = 7, dpi = 300)


# =============================================================================
# 10) CONTRAST POSTERIORS  - do the clades actually DIFFER?
# =============================================================================
# For every clade pair, the posterior of the DIFFERENCE, with a line at zero.
# If the mass sits clearly off zero, the difference is real. If it straddles
# zero, it is not - however far apart the point estimates looked.

contr_draws <- function(mat, label) {
  out <- list()
  for (i in seq_along(groups)) for (j in seq_along(groups)) {
    if (i >= j) next
    a <- groups[i]; b <- groups[j]
    d <- mat[, a] - mat[, b]
    out[[length(out) + 1]] <- tibble::tibble(
      parameter = label, pair = paste(a, "-", b), value = as.numeric(d))
  }
  dplyr::bind_rows(out)
}

cd_long <- dplyr::bind_rows(
  contr_draws(E_g, "Growth: E (eV)"),
  contr_draws(Topt_draws, "Growth: Topt (°C)"),
  contr_draws(E_r, "Respiration: E (eV)"))

# Label each pair with P(difference > 0).
plab <- cd_long %>% dplyr::group_by(parameter, pair) %>%
  dplyr::summarise(P = mean(value > 0), .groups = "drop") %>%
  dplyr::mutate(lab = sprintf("P=%.2f", P),
                strong = P > 0.95 | P < 0.05)

p_contr <- ggplot2::ggplot(cd_long, ggplot2::aes(value, pair))
p_contr <- if (.gg) {
  p_contr + ggdist::stat_halfeye(.width = 0.95, point_size = 1.3,
                                 slab_alpha = 0.7, fill = "#2c7fb8")
} else {
  p_contr + ggplot2::geom_violin(fill = "#2c7fb8", alpha = 0.7, colour = NA)
}
p_contr <- p_contr +
  ggplot2::geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
  ggplot2::geom_text(data = plab,
                     ggplot2::aes(x = Inf, y = pair, label = lab,
                                  fontface = ifelse(strong, "bold", "plain")),
                     hjust = 1.05, size = 3, inherit.aes = FALSE) +
  ggplot2::facet_wrap(~ parameter, scales = "free_x") +
  ggplot2::labs(
    title = "Clade contrasts: posterior of the difference",
    subtitle = "red line = no difference. P = posterior probability the difference is positive (bold = >0.95 or <0.05)",
    x = "difference (clade A - clade B)", y = NULL) +
  ggplot2::theme_bw(11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

ggplot2::ggsave(file.path(figures_dir, "fig_bayes_contrast_posteriors.png"),
                p_contr, width = 13, height = 8, dpi = 300)


# =============================================================================
# 11) E_R - E_G : WHY CUE collapses with warming
# =============================================================================
# If respiration is more temperature-sensitive than growth (E_R > E_G), then
# warming costs more than it returns and CUE must fall. This single posterior is
# the mechanism behind the CUE curve.

dE_long <- tibble::as_tibble(dE) %>%
  tidyr::pivot_longer(dplyr::everything(), names_to = "Group", values_to = "value")

p_dE <- ggplot2::ggplot(dE_long, ggplot2::aes(value, Group, fill = Group))
p_dE <- if (.gg) {
  p_dE + ggdist::stat_halfeye(.width = c(0.5, 0.95), point_size = 1.8, slab_alpha = 0.8)
} else {
  p_dE + ggplot2::geom_violin(alpha = 0.8, colour = NA)
}
p_dE <- p_dE +
  ggplot2::geom_vline(xintercept = 0, colour = "red", linetype = "dashed", linewidth = 0.8) +
  ggplot2::labs(
    title = expression(bold("Differential activation energy:  "*E[respiration] - E[growth])),
    subtitle = "left of the red line = GROWTH is the more temperature-sensitive process. This is why CUE falls as it warms.",
    x = "E(respiration) - E(growth)   [eV]", y = NULL) +
  ggplot2::theme_bw(12) +
  ggplot2::theme(legend.position = "none",
                 plot.title = ggplot2::element_text(face = "bold"))

ggplot2::ggsave(file.path(figures_dir, "fig_bayes_E_resp_minus_growth.png"),
                p_dE, width = 8, height = 5, dpi = 300)


# =============================================================================
# 12) Report
# =============================================================================

message("11_bayesian_plots.R done.\n")
message("Clade parameters:")
print(as.data.frame(clade_params %>%
                      dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 2)))))
message("\nE(respiration) - E(growth)  [negative => growth is MORE temperature-sensitive]:")
print(as.data.frame(e_diff %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 3)))))
message("\nStrongest clade contrasts (posterior probability > 0.95 or < 0.05):")
print(as.data.frame(contrasts_tbl %>%
                      dplyr::filter(P_a_greater > 0.95 | P_a_greater < 0.05) %>%
                      dplyr::arrange(parameter, dplyr::desc(abs(P_a_greater - 0.5))) %>%
                      dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 3)))))
message("\nFigures written to: ", figures_dir)
