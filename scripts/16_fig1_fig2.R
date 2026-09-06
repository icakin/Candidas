# =============================================================================
# 16_fig1_fig2.R  --  THE TWO RESPIROMETRY MAIN FIGURES (Fig 1 decoupling, Fig 2 the carbon bill)
# =============================================================================
#
# One script. Three figures. Everything computed draw-by-draw from the posteriors
# fitted in 09_bayesian_models.R (growth = Sharpe-Schoolfield, respiration =
# Boltzmann-Arrhenius, both hierarchical with isolates nested in clades).
#
#   FIG 1  DECOUPLING          Growth turns over; respiration does not. Growth is
#                              ~2x more temperature-sensitive. Therefore the
#                              carbon-economy optimum of every taxon lies BELOW
#                              human body temperature.
#
#   FIG 2  THE BILL            What the host actually costs: growth lost and
#                              carbon burned at 37 C (body) and 40 C (fever) -
#                              plus, honestly, which clade differences resolve.
#
#   (Fig 3 RETIRED.) The Casadevall correlation - warmer environmental optimum ->
#   lower fever cost, posterior Spearman rho ~ -0.9 - is now a PANEL of Fig 2, not
#   its own figure: at n = 5 taxa it is suggestive, not decisive, and does not yet
#   deserve to anchor a standalone figure. When more species land it becomes the
#   headline. The mechanism panel (morphing curves) moves to Supplementary.
#
# ---------------------------------------------------------------------------
# THE ONE THING THAT MAKES THIS PUBLISHABLE
# ---------------------------------------------------------------------------
# Absolute CUE is not identifiable here: it scales with cell volume (unknown to a
# factor of several) and with inoculum density. But every quantity in Figs 1c-1d,
# 2 and 3 is a RATIO, and the unknown constants cancel EXACTLY:
#
#   T_opt(CUE)    : depends only on E, E_R, Eh, Th   (shape, not level)
#   carbon tax    : [R/G](T) / min[R/G]              (alpha, lnB0 cancel)
#   fever cost    : [R/G](40) / [R/G](37)            (even the reference cancels)
#   growth % peak : G(T) / G(Topt)                   (lnB0 cancels)
#   activation E  : a SLOPE in log space             (a constant shifts intercept)
#
# Figs 1a/1b are absolute rates and DO depend on those constants. Everything you
# CLAIM lives in the ratio panels. Say so in the legend. Proof: header of
# 13_thermal_economics.R.
#
# ---------------------------------------------------------------------------
# WHAT THE DATA DOES NOT SUPPORT - do not let a figure imply otherwise
# ---------------------------------------------------------------------------
#   * (CORRECTED) The clades ARE mostly rankable. 8 of the 10 pairwise contrasts in
#     relative respiratory cost at 40 C resolve. ALL FOUR C. auris clades credibly
#     pay less than C. parapsilosis; Clade IV is credibly cheapest of all. Only
#     I-vs-II [0.68-1.05] and I-vs-III [0.99-1.41] cannot be separated.
#     An earlier note here said "only 3 of 10 resolve". That was WRONG. It came from
#     a Python sanity check that drew E, Eh, Th and E_R as INDEPENDENT normals. In
#     the real posterior they are strongly CORRELATED, and those correlations HALVE
#     the interval on every derived quantity. Trust the R output; a normal-
#     approximation cross-check on this model is ~2x too conservative, every time.
#   * Correlating thermal traits with clade geography. n = 4 clades. Even though the
#     ordering is now mostly resolved, four points is not a correlation, and the two
#     that DON'T resolve (I vs II, I vs III) are the middle of the putative gradient.
#     Origins stay LABELS. Say "consistent with", never "correlated".
#   * The n = 5 Casadevall correlation (Fig 2, Topt vs fever cost) is SUGGESTIVE,
#     not decisive: rho = -0.90 [-1.00, -0.30]. Report it as such. And note Topt and
#     fever cost are PARTLY geometrically linked (a low optimum puts 40 C further
#     down the decline) - but respiration E_R is independent and could have broken
#     it; it did not, so the bundling is a real result, not a tautology.
#   * C. glabrata is excluded (EXCLUDE_GROUPS in config.R) pending a repeat; its
#     respiration was ~2x noisier than any other taxon.
#
# INPUTS   models/bayes_growth_ss.rds, bayes_resp_arr.rds,
#          models/bayes_data_growth.rds, bayes_data_resp.rds
# OUTPUTS  figures/FIG1_decoupling.png/.pdf
#          figures/FIG2_consequences.png/.pdf   (name follows save_fig() below, not this comment)
#          tables/fig_values.csv, tables/fig_contrasts.csv
#
# RUN AFTER 09_bayesian_models.R.
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

need <- c("posterior", "ggplot2", "dplyr", "tidyr", "tibble", "readr", "patchwork")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) stop("Install first:\n  install.packages(c(",
                           paste(sprintf('"%s"', miss), collapse = ", "), "))")
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble); library(patchwork)
})

CINV   <- 11604.51812     # 1 / k_B  [K/eV] - exactly as used in 10
TREF   <- 293.15
T_BODY <- 37; T_FEVER <- 40; T_MIN <- 22; T_MAX <- 44
E_MTE  <- 0.65            # metabolic theory's "universal" activation energy

message("\n", strrep("=", 66))
message("16_fig1_fig2.R   Growth = Sharpe-Schoolfield.  Resp = Arrhenius.")
message(strrep("=", 66))

req <- function(p) { if (!file.exists(p)) stop("Not found: ", p,
                     "\nRun 09_bayesian_models.R first."); p }
fit_g <- readRDS(req(file.path(models_dir, "bayes_growth_ss.rds")))
fit_r <- readRDS(req(file.path(models_dir, "bayes_resp_arr.rds")))
gdat  <- readRDS(req(file.path(models_dir, "bayes_data_growth.rds")))
rdat  <- readRDS(req(file.path(models_dir, "bayes_data_resp.rds")))

dg <- posterior::as_draws_df(fit_g)
dr <- posterior::as_draws_df(fit_r)
nd <- min(nrow(dg), nrow(dr)); dg <- dg[seq_len(nd), ]; dr <- dr[seq_len(nd), ]

GRPS <- STRAIN_GROUPS[STRAIN_GROUPS %in% unique(as.character(gdat$Group))]
iso_map <- gdat %>% dplyr::distinct(Isolate, Group) %>%
  dplyr::mutate(Isolate = as.character(Isolate), Group = as.character(Group))
message(sprintf("Groups: %s | isolates: %d | draws: %d",
                paste(GRPS, collapse = ", "), nrow(iso_map), nd))


# =============================================================================
# 1) MODEL MATHS  (all closed forms; see 13_thermal_economics.R for derivations)
# =============================================================================
dcol <- function(d, nm) if (nm %in% names(d)) as.numeric(d[[nm]]) else rep(NA_real_, nrow(d))
qq   <- function(x, p) stats::quantile(x[is.finite(x)], p, names = FALSE)
med  <- function(x) stats::median(x[is.finite(x)])

# Sharpe-Schoolfield log-growth. lnB0 is OMITTED: it cancels in every ratio we
# report, and including it would only invite a cell-volume argument.
lnG <- function(TC, E, Eh, Th) {
  TK <- TC + 273.15; u <- Eh * CINV * (1 / Th - 1 / TK)
  E * CINV * (1 / TREF - 1 / TK) - log1p(exp(pmin(u, 700)))
}
# Arrhenius log-respiration, alpha omitted for the same reason.
lnR <- function(TC, ER) ER * CINV * (1 / TREF - 1 / (TC + 273.15))

# Topt of an SS curve whose Boltzmann slope is `Eslope`.
#   Eslope = E      -> T_opt(growth)
#   Eslope = E - ER -> T_opt(CUE)      (the temperature where R/G is cheapest)
topt <- function(Eslope, Eh, Th) {
  s <- Eslope / Eh; o <- rep(NA_real_, length(s))
  k <- is.finite(s) & s > 0 & s < 1
  o[k] <- 1 / (1 / Th[k] - stats::qlogis(s[k]) / (Eh[k] * CINV)) - 273.15
  o
}
ln_cost <- function(TC, E, ER, Eh, Th) {          # log(R/G), up to a constant
  TK <- TC + 273.15; u <- Eh * CINV * (1 / Th - 1 / TK)
  (ER - E) * CINV * (1 / TREF - 1 / TK) + ifelse(u > 30, u, log1p(exp(pmin(u, 30))))
}
tax <- function(TC, E, ER, Eh, Th) {              # carbon tax vs OWN cheapest temp
  t0 <- topt(E - ER, Eh, Th); t0[!is.finite(t0)] <- T_MIN
  t0 <- pmin(pmax(t0, T_MIN), T_MAX)
  exp(ln_cost(TC, E, ER, Eh, Th) - ln_cost(t0, E, ER, Eh, Th))
}
gpct <- function(TC, E, Eh, Th) {                 # growth as % of OWN peak
  tp <- topt(E, Eh, Th); out <- rep(NA_real_, length(E)); k <- is.finite(tp)
  out[k] <- 100 * exp(lnG(TC, E[k], Eh[k], Th[k]) - lnG(tp[k], E[k], Eh[k], Th[k]))
  out
}
# Th that puts the optimum of (E, Eh) at target_C. Topt is monotone in Th at fixed
# E, Eh, so bisection is exact.
th_for_topt <- function(E, Eh, target_C, iter = 60) {
  L <- rep(285, length(E)); H <- rep(340, length(E))
  for (i in seq_len(iter)) {
    M <- (L + H) / 2; t <- topt(E, Eh, M)
    up <- is.finite(t) & t < target_C
    L[up] <- M[up]; H[!up] <- M[!up]
  }
  (L + H) / 2
}

P <- lapply(GRPS, function(g) tibble::tibble(
  Group = g, .draw = seq_len(nd),
  lnB0 = dcol(dg, paste0("b_lnB0_Group",  g)),
  E    = dcol(dg, paste0("b_E_Group",     g)),
  Eh   = dcol(dg, paste0("b_Eh_Group",    g)),
  Th   = dcol(dg, paste0("b_Th_Group",    g)),
  alpha = dcol(dr, paste0("b_alpha_Group", g)),
  ER    = dcol(dr, paste0("b_E_Group",     g)))) %>% dplyr::bind_rows()

Piso <- lapply(seq_len(nrow(iso_map)), function(i) {
  iso <- iso_map$Isolate[i]; g <- iso_map$Group[i]
  re <- function(d, par) {
    nm <- sprintf("r_Isolate__%s[%s,Intercept]", par, iso)
    if (nm %in% names(d)) as.numeric(d[[nm]]) else rep(0, nrow(d))
  }
  tibble::tibble(Isolate = iso, Group = g,
    E  = dcol(dg, paste0("b_E_Group",  g)) + re(dg, "E"),
    Eh = dcol(dg, paste0("b_Eh_Group", g)),
    Th = dcol(dg, paste0("b_Th_Group", g)) + re(dg, "Th"),
    ER = dcol(dr, paste0("b_E_Group",  g)) + re(dr, "E"))
}) %>% dplyr::bind_rows()


# =============================================================================
# 2) EVERY NUMBER IN THE PAPER, ONCE
# =============================================================================
S <- lapply(GRPS, function(g) {
  p   <- dplyr::filter(P, Group == g)
  tg  <- topt(p$E, p$Eh, p$Th)
  tc  <- topt(p$E - p$ER, p$Eh, p$Th)
  x37 <- tax(T_BODY,  p$E, p$ER, p$Eh, p$Th)
  x40 <- tax(T_FEVER, p$E, p$ER, p$Eh, p$Th)
  g37 <- gpct(T_BODY,  p$E, p$Eh, p$Th)
  g40 <- gpct(T_FEVER, p$E, p$Eh, p$Th)
  # fever cost = [R/G](40) / [R/G](37). The reference temperature cancels, so this
  # is the single most assumption-free number in the paper.
  fev <- exp((lnR(T_FEVER, p$ER) - lnG(T_FEVER, p$E, p$Eh, p$Th)) -
             (lnR(T_BODY,  p$ER) - lnG(T_BODY,  p$E, p$Eh, p$Th)))
  tibble::tibble(Group = g,
    E = med(p$E),  E_lo = qq(p$E, .025),  E_hi = qq(p$E, .975),
    ER = med(p$ER), ER_lo = qq(p$ER, .025), ER_hi = qq(p$ER, .975),
    dE = med(p$E - p$ER), dE_lo = qq(p$E - p$ER, .025), dE_hi = qq(p$E - p$ER, .975),
    dE_credible = qq(p$E - p$ER, .025) > 0,
    Eh = med(p$Eh), Eh_lo = qq(p$Eh, .025), Eh_hi = qq(p$Eh, .975),
    Topt = med(tg), Topt_lo = qq(tg, .025), Topt_hi = qq(tg, .975),
    Tcue = med(tc), Tcue_lo = qq(tc, .025), Tcue_hi = qq(tc, .975),
    P_below_body = mean(tc < T_BODY, na.rm = TRUE),
    tax37 = med(x37), tax37_lo = qq(x37, .025), tax37_hi = qq(x37, .975),
    tax40 = med(x40), tax40_lo = qq(x40, .025), tax40_hi = qq(x40, .975),
    growth37 = med(g37), growth40 = med(g40),
    g40_lo = qq(g40, .025), g40_hi = qq(g40, .975),
    # DRAW-WISE, then summarise. `med(g37) - med(g40)` would be a difference of
    # medians, which is NOT the median of the difference and carries no interval.
    # g37 and g40 are strongly correlated within a draw (same curve), so the
    # draw-wise version is far more precise than the two marginals imply.
    growth_lost    = med(g37 - g40),
    growth_lost_lo = qq(g37 - g40, .025),
    growth_lost_hi = qq(g37 - g40, .975),
    # GROWTH KEPT THROUGH THE FEVER, as a plain ratio: growth at 40 / growth at 37.
    # This replaces "percentage points", which is correct but a trip hazard - and it
    # makes the two halves of the bill the numerator and denominator of one fraction:
    #     fever cost = (respiration 40/37) / (growth 40/37)
    # so panel b now shows what panel b's own arithmetic is made of.
    growth_kept    = med(exp(lnG(T_FEVER, p$E, p$Eh, p$Th) -
                             lnG(T_BODY,  p$E, p$Eh, p$Th))),
    gk_lo = qq(exp(lnG(T_FEVER, p$E, p$Eh, p$Th) - lnG(T_BODY, p$E, p$Eh, p$Th)), .025),
    gk_hi = qq(exp(lnG(T_FEVER, p$E, p$Eh, p$Th) - lnG(T_BODY, p$E, p$Eh, p$Th)), .975),
    fever_cost = med(fev), fc_lo = qq(fev, .025), fc_hi = qq(fev, .975))
}) %>% dplyr::bind_rows() %>% dplyr::arrange(match(Group, STRAIN_GROUPS))
readr::write_csv(S, file.path(tables_dir, "fig_values.csv"))

Siso <- Piso %>% dplyr::group_by(Isolate, Group) %>%
  dplyr::group_modify(~ tibble::tibble(
    E = med(.x$E), ER = med(.x$ER),
    Topt = med(topt(.x$E, .x$Eh, .x$Th)),
    Tcue = med(topt(.x$E - .x$ER, .x$Eh, .x$Th)),
    growth40 = med(gpct(T_FEVER, .x$E, .x$Eh, .x$Th)),
    tax40 = med(tax(T_FEVER, .x$E, .x$ER, .x$Eh, .x$Th)))) %>%
  dplyr::ungroup()

# ---- pairwise fever-tax contrasts (the honesty panel) ------------------------
wide <- P %>% dplyr::mutate(x = tax(T_FEVER, E, ER, Eh, Th)) %>%
  dplyr::select(.draw, Group, x) %>%
  tidyr::pivot_wider(names_from = Group, values_from = x)

CON <- utils::combn(GRPS, 2, simplify = FALSE) %>%
  lapply(function(p) {
    d <- wide[[p[1]]] / wide[[p[2]]]; d <- d[is.finite(d)]
    tibble::tibble(a = p[1], b = p[2],
      ratio = stats::median(d), lo = qq(d, .025), hi = qq(d, .975),
      # `credible` MUST agree with the interval drawn beside it. A one-sided
      # P > 0.95 would flag contrasts whose 95% interval still spans 1.
      credible = (qq(d, .975) < 1) | (qq(d, .025) > 1))
  }) %>% dplyr::bind_rows() %>% dplyr::arrange(ratio)
readr::write_csv(CON, file.path(tables_dir, "fig_contrasts.csv"))

message("\n--- headline numbers ---------------------------------------------")
for (i in seq_len(nrow(S))) {
  r <- S[i, ]
  message(sprintf("%-8s Topt %.1f | Tcue %.1f (P<37 = %.3f) | growth@40 %2.0f%% | tax@40 %.2fx | fever cost %.2fx",
    r$Group, r$Topt, r$Tcue, r$P_below_body, r$growth40, r$tax40, r$fever_cost))
}
message(sprintf("\nfever-tax contrasts: %d of %d resolve", sum(CON$credible), nrow(CON)))


# =============================================================================
# 3) STYLE  (Okabe-Ito; colour-blind safe and greyscale-separable)
# =============================================================================
# Clade1 and Clade2 used to be two blues (#0072B2 / #56B4E9) - indistinguishable at
# a glance, and Clade2 is the clade that behaves differently, so it needs to READ
# as different. With glab excluded its pink is free, so Clade2 takes it.
OI <- c(Clade1 = "#0072B2",   # blue
        Clade2 = "#CC79A7",   # pink   (was light blue - too close to Clade1)
        Clade3 = "#009E73",   # green
        Clade4 = "#E69F00",   # orange
        glab   = "#56B4E9",   # light blue (unused while excluded)
        para   = "#D55E00",   # vermillion
        Hae    = "#984EA3",   # violet (C. haemulonii; was grey)
        Duo    = "#00A6A6")   # teal   (C. duobushaemulonii; was black)
SHORT <- c(Clade1 = "C. auris I", Clade2 = "C. auris II", Clade3 = "C. auris III",
           Clade4 = "C. auris IV", glab = "C. glabrata", para = "C. parapsilosis",
           Hae = "C. haemulonii", Duo = "C. duobushaemulonii")
FULL  <- GROUP_LABEL_1L
RED   <- "#B22222"

th <- theme_classic(base_size = 8) +
  theme(
    axis.text  = element_text(size = 6.8, colour = "grey25"),
    axis.title = element_text(size = 7.5, colour = "grey10"),
    axis.line  = element_line(linewidth = .3, colour = "grey30"),
    axis.ticks = element_line(linewidth = .3, colour = "grey30"),
    legend.text = element_text(size = 6.4), legend.title = element_blank(),
    legend.key.height = unit(7.5, "pt"), legend.key.width = unit(13, "pt"),
    legend.background = element_rect(fill = alpha("white", .85), colour = NA),
    legend.key        = element_rect(fill = alpha("white", .85), colour = NA),
    legend.margin = margin(0, 0, 0, 0),
    plot.title    = element_text(size = 8, face = "bold", hjust = 0, colour = "grey5"),
    plot.subtitle = element_text(size = 6.8, colour = "grey35", lineheight = 1.05),
    plot.tag      = element_text(size = 9, face = "bold"),
    plot.title.position = "plot",
    plot.margin = margin(4, 6, 3, 4))

S$Group    <- factor(S$Group,    levels = GRPS)
Siso$Group <- factor(Siso$Group, levels = GRPS)
TG <- seq(T_MIN, T_MAX, by = 0.25)

save_fig <- function(f, name, h_mm, w_mm = 183) {
  W <- w_mm / 25.4; H <- h_mm / 25.4
  figMSdir <- file.path(figures_dir, "manuscript"); dir.create(figMSdir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(figMSdir, paste0(name, ".png")), f, width = W, height = H,
         dpi = 600, bg = "white")
  ok <- tryCatch({
    dv <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
    ggsave(file.path(figMSdir, paste0(name, ".pdf")), f, width = W, height = H,
           device = dv, bg = "white"); TRUE
  }, error = function(e) { message("  PDF failed (", name, "): ",
                                   conditionMessage(e)); FALSE })
  message("  ", name, ".png", if (ok) " + .pdf" else "")
}


# =============================================================================
# FIGURE 1 - THE PHENOMENON  (full page: the complete argument, in order)
# =============================================================================
#   a  growth turns over            (data + Sharpe-Schoolfield)
#   b  respiration does not         (data + Arrhenius; LOO-favoured over SS)
#   c  WHY: growth is ~2x more temperature-sensitive than respiration
#   d  SO: carbon-use efficiency collapses - and peaks early
#   e  THEREFORE: every economic optimum lies below human body temperature
#
# a, b and d are ABSOLUTE (fg C/h, and CUE's LEVEL), so they depend on cell volume
# and inoculum density - which we cannot pin down. c and e are SLOPES/RATIOS and
# are invariant to both. The CLAIM lives in c and e; a, b and d are the picture.
# Panel d's subtitle says so outright: the peak POSITION is invariant even though
# the height is not (proof: header of 13_thermal_economics.R).
# =============================================================================
cur <- lapply(GRPS, function(g) {
  p <- dplyr::filter(P, Group == g)
  i <- if (nrow(p) > 1200) sample.int(nrow(p), 1200) else seq_len(nrow(p))
  Gm <- vapply(TG, function(t)
    exp(p$lnB0[i] + lnG(t, p$E[i], p$Eh[i], p$Th[i])), numeric(length(i)))
  Rm <- vapply(TG, function(t)
    exp(p$alpha[i] + lnR(t, p$ER[i])), numeric(length(i)))
  Cm <- Gm / (Gm + Rm)                       # CUE, derived draw-by-draw
  mk <- function(M, w) tibble::tibble(Group = g, rate = w, T_C = TG,
    med = apply(M, 2, stats::median),
    lo  = apply(M, 2, stats::quantile, .025, names = FALSE),
    hi  = apply(M, 2, stats::quantile, .975, names = FALSE))
  dplyr::bind_rows(mk(Gm, "Growth"), mk(Rm, "Respiration"), mk(Cm, "CUE"))
}) %>% dplyr::bind_rows() %>% dplyr::mutate(Group = factor(Group, levels = GRPS))

pts <- dplyr::bind_rows(
  gdat %>% dplyr::transmute(Group = as.character(Group), rate = "Growth",
                            T_C = as.numeric(T), y = y_raw),
  rdat %>% dplyr::transmute(Group = as.character(Group), rate = "Respiration",
                            T_C = as.numeric(T), y = y_raw)) %>%
  dplyr::mutate(Group = factor(Group, levels = GRPS))

# ---- CLIP FITTED CURVES AT EACH TAXON'S LAST OBSERVED TEMPERATURE ------------
# No extrapolated predictions are drawn: e.g. Duo has no kept data above 38 C
# (its 40-44 C curves showed no growth and are excluded), so its curve stops
# there instead of being extrapolated to 44 C. Purely cosmetic - the posterior,
# Topt, activation energies and all contrasts are untouched.
TMAX_G <- pts %>% dplyr::group_by(Group) %>%
  dplyr::summarise(Tmax_obs = max(T_C, na.rm = TRUE), .groups = "drop")
cur <- cur %>% dplyr::left_join(TMAX_G, by = "Group") %>%
  dplyr::filter(T_C <= Tmax_obs) %>% dplyr::select(-Tmax_obs)

# ---- DISPLAY GROWTH AS THE RATE CONSTANT r (h^-1), NOT A CARBON FLUX ----------
# The oxygen-dynamics model recovers r; the growth carbon flux is G = r x q, where
# q is a temperature-INDEPENDENT per-cell carbon quota (constant per taxon). So
# growth_fgC_h / q = r (h^-1), a pure rescale of the growth panel that leaves E_G,
# T_opt and CUE completely unchanged (q shifts the intercept, not the slope).
# Panel a is therefore shown in r (h^-1) so that E_G is visibly the activation
# energy of a rate constant; respiration (b) and CUE stay as carbon fluxes, and
# the r->carbon conversion enters only where CUE is formed (Cm above, computed
# from the carbon fluxes BEFORE this rescale).
.der <- readr::read_csv(derived_csv, show_col_types = FALSE)
# quota per group, straight from the derived table (growth_fgC_h / r-in-h^-1)
.q <- .der %>%
  dplyr::filter(is.finite(growth_fgC_h), is.finite(growth_C_per_C_h),
                growth_fgC_h > 0, growth_C_per_C_h > 0) %>%
  dplyr::mutate(OTU = as.integer(OTU), q = growth_fgC_h / growth_C_per_C_h) %>%
  dplyr::left_join(
    readr::read_csv(file.path(tables_dir, "otu_names.csv"), show_col_types = FALSE) %>%
      dplyr::transmute(OTU = as.integer(OTU), Group = as.character(group)),
    by = "OTU") %>%
  dplyr::group_by(Group) %>% dplyr::summarise(q = stats::median(q), .groups = "drop")
QUOTA <- stats::setNames(.q$q, .q$Group)

.rescale_growth <- function(df, col) {
  df[[col]] <- ifelse(df$rate == "Growth",
                      df[[col]] / QUOTA[as.character(df$Group)], df[[col]])
  df
}
cur <- .rescale_growth(.rescale_growth(.rescale_growth(cur, "med"), "lo"), "hi")
pts <- .rescale_growth(pts, "y")

# ---- per-isolate fitted curves (boss: "show all isolates") -------------------
# One posterior-median TPC per isolate, from the SAME hierarchical fit: clade
# fixed effect + that isolate's random effects (shape via E/Th/ER, and height via
# the lnB0/alpha random intercepts if the model carries them; re() returns 0 if not).
cur_iso <- lapply(seq_len(nrow(iso_map)), function(k) {
  iso <- iso_map$Isolate[k]; g <- iso_map$Group[k]
  re <- function(d, par) { nm <- sprintf("r_Isolate__%s[%s,Intercept]", par, iso)
    if (nm %in% names(d)) as.numeric(d[[nm]]) else rep(0, nrow(d)) }
  lnB0 <- dcol(dg, paste0("b_lnB0_Group", g)) + re(dg, "lnB0")
  E    <- dcol(dg, paste0("b_E_Group",    g)) + re(dg, "E")
  Eh   <- dcol(dg, paste0("b_Eh_Group",   g))
  Th   <- dcol(dg, paste0("b_Th_Group",   g)) + re(dg, "Th")
  alph <- dcol(dr, paste0("b_alpha_Group", g)) + re(dr, "alpha")
  ER   <- dcol(dr, paste0("b_E_Group",     g)) + re(dr, "E")
  Gm <- vapply(TG, function(t) exp(lnB0 + lnG(t, E, Eh, Th)), numeric(nd))
  Rm <- vapply(TG, function(t) exp(alph + lnR(t, ER)),        numeric(nd))
  mk <- function(M, w) tibble::tibble(Isolate = iso, Group = g, rate = w, T_C = TG,
                                      med = apply(M, 2, stats::median))
  dplyr::bind_rows(mk(Gm, "Growth"), mk(Rm, "Respiration"))
}) %>% dplyr::bind_rows() %>% dplyr::mutate(Group = factor(Group, levels = GRPS))
cur_iso <- .rescale_growth(cur_iso, "med")   # same r (h^-1) rescale as the clade curve
cur_iso <- cur_iso %>% dplyr::left_join(TMAX_G, by = "Group") %>%
  dplyr::filter(T_C <= Tmax_obs) %>% dplyr::select(-Tmax_obs)  # same clip as cur

# ---- a, b: the two rates -----------------------------------------------------
# ND = assayed-but-no-quantifiable-growth: temperatures a group was run at (every
# group was assayed at all 12 temperatures) with NO kept growth curve. Drawn as
# crosses at the panel's lower margin so non-detections are visible rather than
# silently absent (e.g. Duo at 40-44 C).
ND_TEMPS <- seq(T_MIN, T_MAX, by = 2)
nd_growth <- pts %>% dplyr::filter(rate == "Growth") %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(kept = list(unique(round(T_C))), .groups = "drop") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(nd = list(setdiff(ND_TEMPS, kept))) %>%
  dplyr::ungroup() %>%
  dplyr::select(Group, nd) %>% tidyr::unnest(nd) %>%
  dplyr::rename(T_C = nd) %>%
  dplyr::mutate(Group = factor(Group, levels = GRPS))

tpc_panel <- function(w, ylab, ttl, sub, legend) {
  cv <- dplyr::filter(cur, rate == w); pp <- dplyr::filter(pts, rate == w)
  cvi <- dplyr::filter(cur_iso, rate == w)
  g <- ggplot(cv, aes(T_C, med, colour = Group, fill = Group)) +
    annotate("rect", xmin = T_BODY, xmax = T_FEVER, ymin = 0, ymax = Inf,
             fill = RED, alpha = .06) +
    # size .5 / alpha .30: these are the only MEASUREMENTS in the whole figure and
    # the first version rendered them as a ghost (size .3, alpha .16).
    geom_point(data = pp, aes(y = y), shape = 16, size = .5, alpha = .30) +
    geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = .10) +
    geom_line(linewidth = .8) +
    { if (w == "Growth" && nrow(nd_growth) > 0)
        geom_point(data = nd_growth %>%
                     dplyr::mutate(y_nd = min(pp$y, na.rm = TRUE) * .55),
                   aes(T_C, y_nd, colour = Group), shape = 4, size = 1.5,
                   stroke = .6, alpha = .9, inherit.aes = FALSE,
                   show.legend = FALSE)
      else NULL } +
    scale_colour_manual(values = OI, labels = SHORT) +
    scale_fill_manual(values = OI, guide = "none") +
    scale_x_continuous(breaks = seq(24, 44, 4)) +
    scale_y_log10() +
    coord_cartesian(xlim = c(T_MIN, T_MAX)) +
    labs(x = "Temperature (°C)", y = ylab, title = ttl, subtitle = sub) + th
  if (legend) g + guides(colour = guide_legend(nrow = 1))  # collected at assembly
  else g + guides(colour = "none")
}
# Panel a: SPECIFIC GROWTH RATE r (h^-1) - the rate constant the model recovers,
# so E_G is unambiguously the activation energy of a rate. Panel b: per-cell
# respiration carbon flux (fg C cell^-1 h^-1), the quantity that pairs with growth
# in CUE. They are deliberately NOT in the same units - that is the point.
f1a <- tpc_panel("Growth", expression("Specific growth rate, "*italic(r)*" (h"^-1*")"),
                 "Growth rates turn over with temperature",
                 "Sharpe-Schoolfield.  × = assayed, no detectable growth", TRUE)
f1b <- tpc_panel("Respiration", expression("O"[2]*"-derived per-cell respiration (fg C cell"^-1*" h"^-1*")"),
                 "No respiratory downturn detected over 22–44 °C",
                 "Arrhenius posterior fits", FALSE)

# ---- c: the mechanism --------------------------------------------------------
El <- dplyr::bind_rows(
  S %>% dplyr::transmute(Group, what = "Growth",      e = E,  lo = E_lo,  hi = E_hi),
  S %>% dplyr::transmute(Group, what = "Respiration", e = ER, lo = ER_lo, hi = ER_hi)) %>%
  dplyr::mutate(what  = factor(what, levels = c("Growth", "Respiration")),
                Group = factor(Group, levels = rev(GRPS)))

f1c <- ggplot(El, aes(e, Group, colour = Group, shape = what)) +
  geom_vline(xintercept = E_MTE, colour = "grey78", linetype = "22", linewidth = .3) +
  geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = .62),
                 linewidth = .9) +
  geom_point(aes(fill = Group), position = position_dodge(width = .62),
             size = 2.4, stroke = .5, colour = "white") +
  scale_shape_manual(values = c(Growth = 21, Respiration = 24),
                     guide = guide_legend(override.aes = list(fill = "grey35"))) +
  scale_colour_manual(values = OI, guide = "none") +
  scale_fill_manual(values = OI, guide = "none") +
  scale_y_discrete(labels = SHORT, expand = expansion(add = c(.6, 1.2))) +
  # parked ABOVE the top row - it used to sit ON the line it annotates
  annotate("text", x = E_MTE, y = length(GRPS) + 1.05, size = 2.1, colour = "grey30",
           label = "0.65 eV — MTE reference", hjust = .5) +
  scale_x_continuous(breaks = seq(0.2, 1.4, 0.2)) +   # was only 0.5 and 1.0
  labs(x = "Activation energy (eV)", y = NULL,
       title = "Growth has higher ascending-limb activation energy than respiration",
       subtitle = sprintf("Ascending-limb activation energies; difference credible in %d of %d taxa.",
                          sum(S$dE_credible), nrow(S))) +
  th + theme(axis.text.y = element_text(size = 7, face = "italic"),
             axis.line.y = element_blank(), axis.ticks.y = element_blank(),
             panel.grid.major.y = element_line(linewidth = .2, colour = "grey94"),
             legend.position = c(.87, .13))

# ---- d: CUE collapses --------------------------------------------------------
cue <- dplyr::filter(cur, rate == "CUE")
cue_pk <- cue %>% dplyr::group_by(Group) %>%
  dplyr::slice_max(med, n = 1, with_ties = FALSE) %>% dplyr::ungroup()

f1d <- ggplot(cue, aes(T_C, med, colour = Group, fill = Group)) +
  annotate("rect", xmin = T_BODY, xmax = T_FEVER, ymin = -Inf, ymax = Inf,
           fill = RED, alpha = .07) +
  geom_vline(xintercept = c(T_BODY, T_FEVER), colour = RED,
             linetype = "22", linewidth = .32) +
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = .10) +
  geom_line(linewidth = .8) +
  geom_point(data = cue_pk, aes(fill = Group), shape = 21, size = 2.4,
             colour = "white", stroke = .6) +
  scale_colour_manual(values = OI, guide = "none") +
  scale_fill_manual(values = OI, guide = "none") +
  scale_x_continuous(breaks = seq(24, 44, 4)) +
  coord_cartesian(xlim = c(T_MIN, T_MAX), ylim = c(0, 1)) +
  annotate("text", x = (T_BODY + T_FEVER)/2, y = .985, vjust = 1, size = 2.1,
           colour = RED, label = "body 37 \u2192 fever 40 \u00b0C") +
  labs(x = "Temperature (°C)", y = "Apparent carbon-use efficiency,  G / (G + R)",
       title = "Apparent carbon-use efficiency declines through the febrile range",
       subtitle = "● peak of each curve") + th

# ---- e: the result -----------------------------------------------------------
ord1 <- rev(as.character(GRPS))   # fixed taxonomic order, matches Fig 1c
S1  <- S    %>% dplyr::mutate(Group = factor(Group, levels = ord1))
Si1 <- Siso %>% dplyr::mutate(Group = factor(Group, levels = ord1))

f1e <- ggplot(S1, aes(y = Group, colour = Group)) +
  annotate("rect", xmin = T_BODY, xmax = T_FEVER, ymin = -Inf, ymax = Inf,
           fill = RED, alpha = .07) +
  geom_vline(xintercept = c(T_BODY, T_FEVER), colour = RED,
             linetype = "22", linewidth = .32) +
  geom_linerange(aes(xmin = Topt_lo, xmax = Topt_hi), colour = "grey65",
                 linewidth = .4) +
  geom_segment(aes(x = Tcue, xend = Topt, yend = Group), colour = "grey40",
               linewidth = .35, linetype = "13") +
  geom_linerange(aes(xmin = Tcue_lo, xmax = Tcue_hi), linewidth = 1.15,
                 lineend = "round") +
  # faint isolate-level CUE optima behind the group summary (nudged up a touch)
  geom_point(data = Si1, aes(x = Tcue), shape = 16, size = 1.1, alpha = .45,
             position = position_nudge(y = .30)) +
  geom_point(aes(x = Topt), shape = 22, size = 2.2, fill = "white",
             colour = "grey40", stroke = .5) +
  geom_point(aes(x = Tcue, fill = Group), shape = 21, size = 3.1,
             colour = "white", stroke = .6) +
  geom_text(aes(x = Tcue, label = sprintf("%.1f", Tcue)), vjust = 2.3, size = 2.2) +
  scale_colour_manual(values = OI, guide = "none") +
  scale_fill_manual(values = OI, guide = "none") +
  scale_y_discrete(labels = SHORT, expand = expansion(add = c(.7, 1.0))) +
  # WIDENED to 28-41. It used to start at 30, which squashed the very gap the panel
  # exists to show: these optima sit WELL below 37, not just shy of it.
  scale_x_continuous(breaks = seq(28, 40, 2)) +
  coord_cartesian(xlim = c(28, 41)) +
  annotate("text", x = (T_BODY + T_FEVER) / 2, y = length(ord1) + .85, size = 2.1,
           colour = RED, label = "febrile") +
  labs(x = "Temperature (°C)", y = NULL,
       title = "Apparent CUE optima lie below 37 °C",
       subtitle = "\u25cf CUE optimum (95% CrI)   \u25a1 growth optimum   \u00b7 isolates") +
  th + theme(axis.text.y = element_text(size = 7.5, face = "italic"),
             axis.line.y = element_blank(), axis.ticks.y = element_blank(),
             panel.grid.major.y = element_line(linewidth = .2, colour = "grey94"))

message("\nWriting figures:")

# TWO REBALANCED MAIN FIGURES (not a 3-way split of one continuous result).
# Figure 1 = the rate decoupling (MECHANISM): a growth turns over, b respiration
# does not, c the activation-energy gap. Figure 2 (assembled at the end, once the
# fever panel f2b exists) = CONSEQUENCES: a apparent-CUE curves, b CUE optima,
# c the 37->40 C change. The full relative-R/G cost curves (f2def) are a monotonic
# transform of the CUE curves [R/G = (1-CUE)/CUE] and move to the Supplement.
FIG1 <- (f1a | f1b) / f1c +
  patchwork::plot_layout(heights = c(1, 0.95), guides = "collect") +
  patchwork::plot_annotation(tag_levels = list(c("a", "b", "c"))) &
  theme(plot.tag = element_text(size = 10, face = "bold"),
        legend.position = "bottom", legend.box = "horizontal")
save_fig(FIG1, "FIG1_decoupling", 150)


# =============================================================================
# FIGURE 2 - THE BILL   (what the host actually costs)
# =============================================================================
#   a  the performance-cost plane: growth kept vs carbon burned, 37 -> 40 C
#   b  the bill, with credible intervals on BOTH quantities
#   c  which clade differences actually resolve  (8 of 10 - see below)
#
# NOTE. The "own cheapest temperature" that the carbon tax is measured against is
# exactly T_opt(CUE) from Figure 1e. Fig 1 establishes the quantities; Fig 2 spends
# them. Say that once in the caption and the two figures lock together without
# repeating a single panel.
#
# BOTH AXES OF PANEL a ARE RATIOS:
#   growth retained = G(T)/G(Topt)              -> lnB0 cancels
#   carbon tax      = [R/G](T) / min[R/G]       -> alpha and lnB0 cancel
# so neither depends on cell volume or inoculum density.
# =============================================================================

# ---- 2a: WHAT THE CARBON TAX IS ---------------------------------------------
# The tax was used as an axis in every other panel and DEFINED IN NONE of them -
# the reader had to reverse-engineer it from the caption. This panel is the
# definition, drawn:
#
#     carbon tax(T) = [R(T)/G(T)]  /  min_T [R/G]
#
# i.e. how many times more carbon this organism burns per unit of growth than it
# does at its own cheapest temperature. 1x = as cheap as it ever gets, and the
# temperature where that happens IS T_opt(CUE) from Figure 1e - so this panel is
# also the bridge between the two figures.
taxcv <- lapply(GRPS, function(g) {
  p <- dplyr::filter(P, Group == g)
  i <- if (nrow(p) > 1200) sample.int(nrow(p), 1200) else seq_len(nrow(p))
  M <- vapply(TG, function(t) tax(t, p$E[i], p$ER[i], p$Eh[i], p$Th[i]),
              numeric(length(i)))
  tibble::tibble(Group = g, T_C = TG,
    med = apply(M, 2, stats::median),
    lo  = apply(M, 2, stats::quantile, .025, names = FALSE),
    hi  = apply(M, 2, stats::quantile, .975, names = FALSE))
}) %>% dplyr::bind_rows() %>% dplyr::mutate(Group = factor(Group, levels = GRPS)) %>%
  dplyr::left_join(TMAX_G, by = "Group") %>%
  dplyr::filter(T_C <= Tmax_obs) %>% dplyr::select(-Tmax_obs)  # clip at last observed T

tax_min <- taxcv %>% dplyr::group_by(Group) %>%
  dplyr::slice_min(med, n = 1, with_ties = FALSE) %>% dplyr::ungroup()

f2def <- ggplot(taxcv, aes(T_C, med, colour = Group, fill = Group)) +
  annotate("rect", xmin = T_BODY, xmax = T_FEVER, ymin = -Inf, ymax = Inf,
           fill = RED, alpha = .07) +
  geom_vline(xintercept = c(T_BODY, T_FEVER), colour = RED,
             linetype = "22", linewidth = .32) +
  geom_hline(yintercept = 1, colour = "grey55", linewidth = .4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = .11) +
  geom_line(linewidth = .8) +
  geom_point(data = tax_min, aes(fill = Group), shape = 21, size = 2.2,
             colour = "white", stroke = .6) +
  scale_colour_manual(values = OI, labels = FULL) +
  scale_fill_manual(values = OI, guide = "none") +
  scale_x_continuous(breaks = seq(24, 44, 4)) +
  scale_y_log10(breaks = c(1, 1.5, 2, 3, 5, 8),
                labels = c("1x", "1.5x", "2x", "3x", "5x", "8x")) +
  coord_cartesian(xlim = c(T_MIN, T_MAX), ylim = c(.95, 8)) +
  annotate("text", x = T_BODY - 0.5, y = 7.2, hjust = 1, size = 2.1, colour = RED,
           label = "body 37 \u00b0C") +
  annotate("text", x = T_FEVER + 0.5, y = 7.2, hjust = 0, size = 2.1, colour = RED,
           label = "fever 40 \u00b0C") +
  labs(x = "Temperature (\u00b0C)",
       y = "Relative R/G,  (R/G)(T) / min(R/G)",
       title = "Temperature dependence of the respiration-to-growth ratio",
       subtitle = paste0("(R/G)(T) \u00f7 min(R/G): respiration per unit growth, relative to each ",
                         "taxon's own\nminimum (\u25cf = that minimum = T_opt(CUE), Fig 1e; ",
                         "1\u00d7 = each taxon's cheapest temperature).")) +
  guides(colour = guide_legend(ncol = 1)) +
  th + theme(legend.position = c(.28, .72))

# 2a. HOT BRANCH ONLY, and stopped at 42 C. Running the paths out to 44 made all
# five sprawl through the same corner and cross each other - you could not tell
# whose was whose. The story is 37 -> 40; 2 C of headroom is enough to show the
# direction of travel without the tangle.
TRAJ_END <- 42
traj <- lapply(GRPS, function(g) {
  p <- dplyr::filter(P, Group == g)
  i <- if (nrow(p) > 1200) sample.int(nrow(p), 1200) else seq_len(nrow(p))
  tp <- med(topt(p$E, p$Eh, p$Th)); if (!is.finite(tp)) tp <- T_MIN
  tg <- TG[TG >= tp & TG <= TRAJ_END]
  if (length(tg) < 2) tg <- TG[TG >= TRAJ_END - 2 & TG <= TRAJ_END]
  tibble::tibble(Group = g, T_C = tg,
    growth = vapply(tg, function(t) med(gpct(t, p$E[i], p$Eh[i], p$Th[i])), numeric(1)),
    tax    = vapply(tg, function(t) med(tax(t, p$E[i], p$ER[i], p$Eh[i], p$Th[i])), numeric(1)))
}) %>% dplyr::bind_rows() %>% dplyr::mutate(Group = factor(Group, levels = GRPS))

m37 <- S %>% dplyr::transmute(Group, x = growth37, y = tax37)
m40 <- S %>% dplyr::transmute(Group, x = growth40, y = tax40,
                              lo = tax40_lo, hi = tax40_hi,
                              glo = g40_lo,  ghi = g40_hi,
                              lab = SHORT[as.character(Group)])
mv  <- dplyr::left_join(m37, m40 %>% dplyr::select(Group, x1 = x, y1 = y), by = "Group")

# The panel uses clip = "off" so bold edge labels are not cut - but that also lets
# the faint trajectory lines spray past the top-right corner. Cap the trajectory
# data to the panel's y-limit so they stop AT the frame instead of over it.
Y_TOP <- max(m40$hi) * 1.10
traj  <- dplyr::filter(traj, tax <= Y_TOP)

.repel2 <- requireNamespace("ggrepel", quietly = TRUE)

f2a <- ggplot() +
  # the good corner. Label moved DOWN into the empty band just above the 1x line,
  # away from the cluster of open circles that was sitting on top of "no penalty".
  annotate("rect", xmin = 88, xmax = 102, ymin = .965, ymax = 1.12,
           fill = "#2E8B57", alpha = .07) +
  annotate("text", x = 95, y = 1.005, size = 2.0, colour = "#2E8B57",
           fontface = "italic", label = "full growth, no penalty") +
  geom_hline(yintercept = 1, colour = "grey80", linewidth = .3) +
  geom_path(data = traj, aes(growth, tax, colour = Group),
            linewidth = .45, alpha = .35) +
  geom_linerange(data = m40, aes(x = x, ymin = lo, ymax = hi, colour = Group),
                 linewidth = .5, alpha = .85) +
  geom_linerange(data = m40, aes(y = y, xmin = glo, xmax = ghi, colour = Group),
                 linewidth = .5, alpha = .85) +
  geom_point(data = m37, aes(x, y, colour = Group), shape = 21, fill = "white",
             size = 2.4, stroke = .8) +
  geom_segment(data = mv, aes(x = x, y = y, xend = x1, yend = y1, colour = Group),
               linewidth = .6, arrow = arrow(length = unit(4.5, "pt"), type = "closed")) +
  geom_point(data = m40, aes(x, y, fill = Group), shape = 21, size = 3.2,
             colour = "white", stroke = .65) +
  { if (.repel2)
      ggrepel::geom_text_repel(data = m40, aes(x, y, label = lab, colour = Group),
        size = 2.2, fontface = "bold", seed = 2, box.padding = .5,
        point.padding = .4, min.segment.length = .3, segment.size = .25,
        segment.alpha = .5, max.overlaps = Inf, show.legend = FALSE)
    else
      geom_text(data = m40, aes(x, y, label = lab, colour = Group),
        size = 2.2, fontface = "bold", hjust = 1.15, vjust = -1.0,
        show.legend = FALSE) } +
  scale_colour_manual(values = OI, guide = "none") +
  scale_fill_manual(values = OI, guide = "none") +
  scale_x_reverse(breaks = seq(100, 40, -10), labels = function(x) paste0(x, "%"),
                  expand = expansion(mult = c(.10, .13))) +   # room for edge labels
  scale_y_log10(breaks = c(1, 1.5, 2, 3, 4), labels = c("1×", "1.5×", "2×", "3×", "4×")) +
  # clip = "off": a bold repelled label (e.g. "C. auris I" at the high-growth edge)
  # was being cut to "C. auris" by the panel boundary. Let it spill into the margin.
  coord_cartesian(ylim = c(.97, Y_TOP), clip = "off") +
  labs(x = "Relative growth rate  (% of taxon max)  →  worse",
       y = "Relative respiratory cost (×)  →  worse",
       title = "What fever does: growth kept vs carbon burned",
       subtitle = paste0("○ 37 °C body   ● 40 °C fever (95% CrI)   arrow = the move fever forces\n",
                         "Cost is vs each taxon's own optimum (T_opt CUE, Fig 1e); both axes ratios")) + th

# ---- 2b: the bill, with intervals on BOTH quantities -------------------------
# The first version drew bare bars (no uncertainty) and encoded the two quantities
# inconsistently - growth loss as bar length, carbon cost as text. Everything else
# in the figure carries a credible interval; this should too. Two matched forests.
# Both halves are now 37 -> 40 RATIOS, and they are the numerator and denominator
# of the same fraction:
#       cost increase = (respiration at 40 / at 37)  /  (growth at 40 / at 37)
# The old left-hand panel was "growth lost, in percentage points". That is correct
# but it is a trip hazard - percentage points vs percent is exactly the kind of
# thing a referee queries - and it also broke the symmetry with the right-hand
# panel. Ratios on both sides: no subtraction, no units, nothing to explain.
# Short strip titles - the long versions were truncated inside the facets
# ("hrough the fever (growth", "n per unit growth ([R/G] a"). The formulae are
# already spelled out in the panel subtitle, so the strips only need to name the
# quantity.
LAB_KEPT <- "Growth-rate retention, 37\u219240 \u00b0C"
LAB_COST <- "R/G change, 37\u219240 \u00b0C"
bill <- dplyr::bind_rows(
  S %>% dplyr::transmute(Group, panel = LAB_KEPT,
                         v = growth_kept, lo = gk_lo, hi = gk_hi,
                         lab = sprintf("%.2f", growth_kept)),
  S %>% dplyr::transmute(Group, panel = LAB_COST,
                         v = fever_cost, lo = fc_lo, hi = fc_hi,
                         lab = sprintf("%.2f×", fever_cost))) %>%
  dplyr::mutate(panel = factor(panel, levels = c(LAB_KEPT, LAB_COST)),
                # dagger: Duo's 40 C values extrapolate beyond its last observed
                # temperature (38 C) - flagged on the value label + in the subtitle
                lab = ifelse(as.character(Group) == "Duo", "\u2020", lab))

# Isolate-level fever cost (tax40/tax37 per isolate) - faint points behind the
# group estimates on the cost panel (the "expose heterogeneity" layer).
# Which isolates actually had detectable growth at 40 C?  wells_by_isolate.csv is
# written by the same filter used for Fig 3 (fit_valid & has_curvature, after the
# 04 exclusions).  Isolates with 0/5 wells at 40 C still receive a model-implied
# ratio, but that ratio is identified by the TPC beyond observed positive growth,
# so it is drawn HOLLOW and excluded from the growth-positive subset marker.
grew40 <- tryCatch({
  readr::read_csv(file.path(tables_dir, "wells_by_isolate.csv"),
                  show_col_types = FALSE) %>%
    dplyr::rename(Isolate = 1) %>%
    dplyr::transmute(Isolate, grew = .data[["40"]] > 0)
}, error = function(e) NULL)

bill_iso <- tryCatch({
  ti <- readr::read_csv(file.path(tables_dir, "carbon_tax_isolate.csv"),
                        show_col_types = FALSE) %>%
    dplyr::filter(Group %in% GRPS)
  out <- ti %>% tidyr::pivot_wider(id_cols = c(Isolate, Group), names_from = T_C,
                                   values_from = tax) %>%
    dplyr::transmute(Isolate, Group,
                     panel = factor(LAB_COST, levels = c(LAB_KEPT, LAB_COST)),
                     v = `40` / `37`) %>%
    dplyr::filter(is.finite(v))
  if (!is.null(grew40)) out <- dplyr::left_join(out, grew40, by = "Isolate")
  if (!"grew" %in% names(out)) out$grew <- TRUE
  out %>% dplyr::mutate(grew = ifelse(is.na(grew), TRUE, grew))
}, error = function(e) NULL)

# Descriptive summary over ONLY the growth-positive isolates of each taxon. For
# C. auris clades all 3/3 isolates grew, so this sits on the all-isolate estimate;
# for Hae and para it is a single isolate, so it is drawn as a diamond WITHOUT an
# interval - a subset marker, not a replacement group estimate.
bill_pos <- if (!is.null(bill_iso)) {
  bill_iso %>% dplyr::filter(grew %in% TRUE) %>%
    dplyr::group_by(Group, panel) %>%
    dplyr::summarise(v = median(v), n = dplyr::n(), .groups = "drop")
} else NULL

ordb <- rev(as.character(GRPS))   # fixed taxonomic order, matches Fig 1c
bill$Group <- factor(bill$Group, levels = ordb)

# The reference line differs by panel: 1.0 = "no change" in BOTH, but it means
# "kept everything" on the left and "cost nothing extra" on the right.
ref1 <- tibble::tibble(panel = factor(c(LAB_KEPT, LAB_COST),
                                      levels = c(LAB_KEPT, LAB_COST)), x = 1)

# Duo's 40 C quantities (both columns) require extrapolation past its last
# temperature with quantifiable growth (38 C), so NO quantitative estimate is
# drawn for it here: its row carries an explicit annotation instead, and the
# model-conditional numbers are reported in the text/supplement only.
bill_est <- bill     %>% dplyr::filter(as.character(Group) != "Duo")
iso_est  <- if (!is.null(bill_iso)) bill_iso %>%
              dplyr::filter(as.character(Group) != "Duo") else NULL
pos_est  <- if (!is.null(bill_pos)) bill_pos %>%
              dplyr::filter(as.character(Group) != "Duo") else NULL
duo_note <- tibble::tibble(
  panel = factor(c(LAB_KEPT, LAB_COST), levels = c(LAB_KEPT, LAB_COST)),
  Group = factor("Duo", levels = ordb),
  x     = c(0.67, 1.52),
  lab   = c("no detectable growth at 40 \u00b0C",
            "(R/G) not data-supported at 40 \u00b0C"))

f2b <- ggplot(bill_est, aes(v, Group, colour = Group)) +
  geom_vline(data = ref1, aes(xintercept = x), colour = "grey55",
             linewidth = .4, inherit.aes = FALSE) +
  { if (!is.null(iso_est) && nrow(iso_est) > 0)
      geom_point(data = iso_est %>% dplyr::mutate(Group = factor(Group, levels = ordb)),
                 aes(shape = grew), size = 1.15, alpha = .55,
                 position = position_nudge(y = -.52))
    else NULL } +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 1.0, lineend = "round") +
  geom_point(aes(fill = Group), shape = 21, size = 2.8, colour = "white", stroke = .6) +
  { if (!is.null(pos_est) && nrow(pos_est) > 0)
      geom_point(data = pos_est %>% dplyr::mutate(Group = factor(Group, levels = ordb)),
                 shape = 23, size = 2.3, colour = "white", stroke = .5,
                 aes(fill = Group), position = position_nudge(y = -.26))
    else NULL } +
  { if (!is.null(pos_est) && nrow(pos_est) > 0)
      geom_text(data = pos_est %>% dplyr::filter(n < 3) %>%
                  dplyr::mutate(Group = factor(Group, levels = ordb)),
                aes(label = sprintf("%.2f\u00d7  ", v)), size = 1.85, hjust = 1,
                position = position_nudge(y = -.26), show.legend = FALSE)
    else NULL } +
  geom_text(aes(label = lab), vjust = -1.3, size = 2.2, fontface = "bold") +
  geom_text(data = duo_note, aes(x, Group, label = lab), inherit.aes = FALSE,
            size = 1.9, colour = "grey30", fontface = "italic", hjust = .5) +
  scale_x_continuous(trans = "log2", breaks = c(.25, .5, 1, 2, 4, 8),
                     labels = c("0.25\u00d7", "0.5\u00d7", "1\u00d7", "2\u00d7", "4\u00d7", "8\u00d7")) +
  facet_wrap(~ panel, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = OI, guide = "none") +
  scale_fill_manual(values = OI, guide = "none") +
  # SHORT labels (no region), not FULL. The wide italic "C. auris IV (South
  # America)" labels here were the widest thing in the left column, and patchwork
  # aligns panel a's left edge to this row - so c's long labels were pushing a's
  # whole plot area to the right, leaving the big empty gap beside a's y-title.
  # The regions are already in panel a's legend; c does not need to repeat them.
  scale_y_discrete(limits = ordb, labels = SHORT,
                   expand = expansion(add = c(.85, .95))) +
  labs(x = NULL, y = NULL,
       title = "Physiological change from 37 to 40 °C",
       subtitle = paste0("\u25cf hierarchical estimate   \u00b7 isolates   ",
                         "\u25c6 median among 40 \u00b0C growth-positive   \u25cb isolate without growth at 40 \u00b0C")) +
  th + theme(axis.text.y = element_text(size = 6.8, face = "italic"),
             axis.line.y = element_blank(), axis.ticks.y = element_blank(),
             panel.grid.major.y = element_line(linewidth = .2, colour = "grey94"),
             strip.background = element_blank(),
             strip.text = element_text(size = 6.6, face = "bold", colour = "grey20"),
             panel.spacing.x = unit(10, "pt"))

# ---- 2c: which clade differences resolve ------------------------------------
# Framing corrected. This used to read "Only N of 10 ... unresolved at n = 3",
# written when a normal-approximation sanity check wrongly said just 3 resolved.
# The REAL posterior resolves 8 of 10, because E/Eh/Th/E_R are strongly CORRELATED
# within a draw and that halves the interval on the derived quantity. It is a
# finding, not an apology - say so.
# COMPACT contrast labels ("IV vs para", not "C. auris IV vs C. parapsilosis").
# The long form was the widest text in the figure's left column, and patchwork
# reserves that width for EVERY panel above it - which is what pushed panel a's
# y-title far left and left the big empty gap. The subtitle names the species.
TINY <- c(Clade1 = "I", Clade2 = "II", Clade3 = "III", Clade4 = "IV",
          glab = "glab", para = "para", Hae = "hae", Duo = "duo")
CONp <- CON %>%
  dplyr::mutate(lab = paste0(TINY[a], " vs ", TINY[b])) %>%
  dplyr::arrange(ratio) %>%
  dplyr::mutate(lab = factor(lab, levels = lab))
.nres <- sum(CON$credible)

f2c <- ggplot(CONp, aes(ratio, lab)) +
  geom_vline(xintercept = 1, colour = "grey30", linewidth = .45) +
  geom_linerange(aes(xmin = lo, xmax = hi, colour = credible), linewidth = .85) +
  geom_point(aes(colour = credible), size = 2) +
  scale_colour_manual(values = c(`TRUE` = "grey15", `FALSE` = "grey72"), guide = "none") +
  scale_x_log10(breaks = c(.3, .5, .7, 1, 1.5, 2)) +
  labs(x = "Ratio of relative R/G at 40 \u00b0C   (< 1 = the first taxon pays less)", y = NULL,
       title = "Posterior pairwise contrasts in relative R/G at 40 \u00b0C",
       subtitle = paste0("I–IV = C. auris clades; hae = C. haemulonii, ",
                         "duo = C. duobushaemulonii.\n",
                         "Dark = 95% CrI excludes 1; light = not separable. Contrasts involving duo extrapolate beyond 38 \u00b0C.")) +
  th + theme(axis.text.y = element_text(size = 6.6),
             axis.line.y = element_blank(), axis.ticks.y = element_blank(),
             panel.grid.major.y = element_line(linewidth = .2, colour = "grey94"))

# ---- 2e: Topt vs fever cost, with a POSTERIOR correlation --------------------
# fever cost = [R/G](40) / [R/G](37): reference temperature cancels, fully invariant.
fever_cost_draw <- function(p) {
  exp((lnR(T_FEVER, p$ER) - lnG(T_FEVER, p$E, p$Eh, p$Th)) -
      (lnR(T_BODY,  p$ER) - lnG(T_BODY,  p$E, p$Eh, p$Th)))
}
brdg <- lapply(GRPS, function(g) {
  p  <- dplyr::filter(P, Group == g)
  tg <- topt(p$E, p$Eh, p$Th)
  fc <- fever_cost_draw(p)
  tibble::tibble(Group = g,
    Topt = med(tg),  Topt_lo = qq(tg, .025),  Topt_hi = qq(tg, .975),
    cost = med(fc),  cost_lo = qq(fc, .025),  cost_hi = qq(fc, .975))
}) %>% dplyr::bind_rows() %>% dplyr::mutate(Group = factor(Group, levels = GRPS))

# Posterior of Spearman's rho across the 5 taxa: at each draw, rank the taxa by
# Topt and by fever cost and correlate. This carries the parameter uncertainty
# into the correlation instead of pretending the 5 medians are fixed points.
.Tm <- vapply(GRPS, function(g) { p <- dplyr::filter(P, Group == g)
                                  topt(p$E, p$Eh, p$Th) }, numeric(nd))
.Cm <- vapply(GRPS, function(g) { p <- dplyr::filter(P, Group == g)
                                  fever_cost_draw(p) }, numeric(nd))
.spear <- function(a, b) { ra <- rank(a); rb <- rank(b)
  d <- ra - rb; 1 - 6 * sum(d^2) / (length(a) * (length(a)^2 - 1)) }
.rho <- vapply(seq_len(nd), function(j) {
  a <- .Tm[j, ]; b <- .Cm[j, ]
  if (any(!is.finite(a)) || any(!is.finite(b))) return(NA_real_)
  .spear(a, b) }, numeric(1))
.rho <- .rho[is.finite(.rho)]
rho_med <- stats::median(.rho); rho_lo <- qq(.rho, .025); rho_hi <- qq(.rho, .975)
rho_pneg <- mean(.rho < 0)
message(sprintf("\nFig 2e  Spearman(Topt, fever cost): rho = %.2f [%.2f, %.2f], P(rho<0) = %.3f",
                rho_med, rho_lo, rho_hi, rho_pneg))

.rlab <- sprintf("posterior Spearman ρ = %.2f\n[%.2f, %.2f]   P(ρ<0) = %s",
                 rho_med, rho_lo, rho_hi,
                 if (rho_pneg > .999) ">0.999" else sprintf("%.3f", rho_pneg))

f2e <- ggplot(brdg, aes(Topt, cost, colour = Group)) +
  annotate("rect", xmin = T_BODY - .08, xmax = T_BODY + .08, ymin = -Inf, ymax = Inf,
           fill = RED, alpha = .16) +
  annotate("text", x = T_BODY, y = max(brdg$cost_hi), size = 2.0, colour = RED,
           hjust = 1.1, vjust = 1, label = "body\n37 °C", lineheight = .9) +
  geom_linerange(aes(xmin = Topt_lo, xmax = Topt_hi), linewidth = .45, alpha = .55) +
  geom_linerange(aes(ymin = cost_lo, ymax = cost_hi), linewidth = .45, alpha = .55) +
  geom_point(aes(fill = Group), shape = 21, size = 3.4, colour = "white", stroke = .7) +
  { if (requireNamespace("ggrepel", quietly = TRUE))
      ggrepel::geom_text_repel(aes(label = SHORT[as.character(Group)]), size = 2.2,
        fontface = "bold", seed = 3, box.padding = .5, min.segment.length = .3,
        segment.size = .25, segment.alpha = .5, max.overlaps = Inf, show.legend = FALSE)
    else geom_text(aes(label = SHORT[as.character(Group)]), size = 2.2,
        fontface = "bold", vjust = -1.3, show.legend = FALSE) } +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.08, vjust = -0.7, size = 2.15,
           colour = "grey25", label = .rlab, lineheight = .95) +
  scale_colour_manual(values = OI, guide = "none") +
  scale_fill_manual(values = OI, guide = "none") +
  scale_y_log10(breaks = c(1.2, 1.5, 2, 2.5), labels = c("1.2×", "1.5×", "2×", "2.5×")) +
  labs(x = "Fitted growth T_opt under assay conditions (°C)",
       y = "Change in R/G, 40 vs 37 °C (×)",
       title = "Higher fitted growth optima, smaller modeled fever costs",
       subtitle = paste0("Exploratory comparative pattern (consistent with thermal ",
                         "restriction): both axes derive\nfrom the same fitted curves; ",
                         "four of seven taxa are C. auris clades.")) +
  th + theme(plot.margin = margin(4, 10, 3, 4))



# FOUR panels, not five. The performance-cost PLANE (f2a) told the same story as
# the bill (f2b) - "who pays what at fever" - twice. Five panels of different
# shapes crammed together looked overstuffed, so the plane goes to Supplementary
# and Fig 2 keeps the four that each say something distinct:
#   a  what "cost" IS            (cost curves + legend + formula)
#   c  the bill                  (growth kept, extra cost)
#   d  the Casadevall correlation (the thesis - a PROPER panel, not a dual axis)
#   e  which clade differences resolve
# f2a (the plane) is still built above; it is simply not assembled here. Add it
# back to a row if you ever want it in the main figure.
# Pair by SHAPE: the two square plots (cost curves, correlation) share the top
# row; the two wide forests (bill, contrasts) stack full-width below. Order of the
# `+` operands sets the tag letters, so: a cost curves, b correlation, c bill,
# d contrasts.
# patchwork::free() releases panel a from LEFT-edge alignment with the wide
# forests below it. Without it, the forests' y-axis labels reserve a wide left
# strip and panel a's y-title is flung to the far left of that strip, leaving a
# big empty gap. free() lets a keep its own natural margin. (patchwork >= 1.2;
# if unavailable, the terse contrast/bill labels above keep the gap small anyway.)
# RESTRUCTURED (reviewer pass): the main figure keeps only the two panels that
# carry the fever story - the R/G cost curves (a) and the 37->40 change (b).
# The Topt-vs-fever-cost scatter (exploratory, axes share a data source) and the
# exhaustive pairwise contrasts move to Supplementary as standalone figures.
# ============================================================================
# FIGURE 2 - CONSEQUENCES (assembled here because it needs f2b, defined above):
#   a  apparent-CUE curves           (f1d)
#   b  CUE optima vs body temp        (f1e)
#   c  37->40 C change: retention + dR/G  (f2b)
# The full relative-R/G cost curves (f2def) are a monotonic transform of the CUE
# curves and move to the Supplement (FIG_SUPP_relative_RG).
# ============================================================================
FIG2 <- f1d + f1e + f2b +
  patchwork::plot_layout(design = "AB\nCC", widths = c(0.58, 0.42),
                         heights = c(1, 0.95)) +
  patchwork::plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
if (exists("fig_keep_add")) fig_keep_add("FIG2_consequences")
save_fig(FIG2, "FIG2_consequences", 200)

# ---- Supplementary figures ---------------------------------------------------
# Relative-R/G cost curves (monotonic transform of Fig 2a); Topt-vs-fever-cost
# scatter and pairwise contrasts (exploratory, retired from main text).
if (exists("fig_keep_add")) fig_keep_add(c("FIG_SUPP_relative_RG",
                                           "FIG_SUPP_topt_vs_fever_cost",
                                           "FIG_SUPP_pairwise_contrasts"))
save_fig(f2def, "FIG_SUPP_relative_RG", 110, w_mm = 150)
save_fig(f2e, "FIG_SUPP_topt_vs_fever_cost", 95, w_mm = 120)
save_fig(f2c, "FIG_SUPP_pairwise_contrasts", 110, w_mm = 150)

# ---- Leave-one-group-out robustness of the Topt vs fever-cost association ----
# Recomputes the posterior Spearman correlation (same .Tm/.Cm draw matrices as
# above) dropping one group at a time, plus a version collapsing the four
# C. auris clades to their per-draw mean (phylogenetic-clustering check).
if (exists("fig_keep_add")) fig_keep_add("FIG_SUPP_loo_correlation")
.loo_rows <- lapply(GRPS, function(g) {
  keep <- setdiff(GRPS, g)
  r <- vapply(seq_len(nd), function(j) {
    a <- .Tm[j, keep]; b <- .Cm[j, keep]
    if (any(!is.finite(a)) || any(!is.finite(b))) return(NA_real_)
    .spear(a, b) }, numeric(1))
  r <- r[is.finite(r)]
  tibble::tibble(what = paste0("without ", GROUP_LABEL_1L[g]),
                 med = stats::median(r), lo = qq(r, .025), hi = qq(r, .975),
                 pneg = mean(r < 0))
})
.r_coll <- vapply(seq_len(nd), function(j) {
  au <- c("Clade1", "Clade2", "Clade3", "Clade4")
  a <- c(mean(.Tm[j, au]), .Tm[j, c("para", "Hae", "Duo")])
  b <- c(mean(.Cm[j, au]), .Cm[j, c("para", "Hae", "Duo")])
  if (any(!is.finite(a)) || any(!is.finite(b))) return(NA_real_)
  .spear(a, b) }, numeric(1))
.r_coll <- .r_coll[is.finite(.r_coll)]
loo_df <- dplyr::bind_rows(
  tibble::tibble(what = "all seven taxa", med = rho_med, lo = rho_lo, hi = rho_hi,
                 pneg = rho_pneg),
  dplyr::bind_rows(.loo_rows),
  tibble::tibble(what = "C. auris clades collapsed (n = 4 taxa)",
                 med = stats::median(.r_coll), lo = qq(.r_coll, .025),
                 hi = qq(.r_coll, .975), pneg = mean(.r_coll < 0))) %>%
  dplyr::mutate(what = factor(what, levels = rev(what)))
f_loo <- ggplot(loo_df, aes(med, what)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = .4) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = .9, colour = "grey30") +
  geom_point(size = 2.2, colour = "grey10") +
  geom_text(aes(label = sprintf("P(\u03c1<0) = %.2f", pneg)), vjust = -1.1,
            size = 2.0, colour = "grey35") +
  coord_cartesian(xlim = c(-1.05, 1.05)) +
  labs(x = "Posterior Spearman \u03c1 (growth T_opt vs 37\u219240 \u00b0C change in R/G)",
       y = NULL,
       title = "Leave-one-out robustness of the comparative association",
       subtitle = "Point = posterior median, line = 95% CrI; exploratory (taxa are phylogenetically structured).") +
  th + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
             panel.grid.major.y = element_line(linewidth = .2, colour = "grey94"))
save_fig(f_loo, "FIG_SUPP_loo_correlation", 90, w_mm = 140)



message("\n", strrep("=", 66))
message("Values: tables/fig_values.csv   Contrasts: tables/fig_contrasts.csv")
message("Two figures. Panels 1a/1b are ABSOLUTE rates (depend on cell volume +")
message("inoculum). Everything you CLAIM - Fig 1c/1d and ALL of Fig 2 - is a ratio")
message("or a slope, invariant to both. Say so in the legend.")
message(strrep("=", 66), "\n")
