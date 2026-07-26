# =============================================================================
# 10_carbon_tax.R  --  WHAT DOES IT COST TO LIVE IN A HUMAN?
# =============================================================================
#
# 13_thermal_economics.R found WHERE each taxon's carbon-economy optimum sits
# (Topt(CUE), all below 37 C). That is only half the question. An optimum at
# 34.4 C is irrelevant if the penalty at 37 C is 2%; it matters enormously if
# the penalty is 3x. This script measures the PENALTY.
#
# THE QUANTITY
# ------------
# Carbon cost per unit growth is R/G. Define the TAX as that cost at temperature
# T, relative to the cost at the cheapest temperature the taxon has:
#
#       Tax(T) = [ R(T)/G(T) ]  /  min_T [ R(T)/G(T) ]        (dimensionless, >= 1)
#
# "At 40 C this organism burns Tax(40) times more carbon per unit of growth than
#  it does at its own best temperature."
#
# WHY IT IS ASSUMPTION-FREE (the same trick as 13, and it goes further)
# --------------------------------------------------------------------
# Unknown scalings G -> aG (cell volume, carbon density) and R -> cR (inoculum):
#
#       Tax(T) = [cR(T) / aG(T)] / [cR(T0) / aG(T0)] = [R(T)/G(T)] / [R(T0)/G(T0)]
#
# a and c cancel. But note that in log space
#
#       ln(R/G) = (alpha - lnB0) + (E_R - E)*c*(1/Tref - 1/T) + log(1 + exp(u))
#
# the INTERCEPTS (alpha - lnB0) are constant in T, so they cancel in the ratio too.
# Tax(T) therefore depends ONLY on E, E_R, Eh and Th - four pure SHAPE parameters.
# Not on cell volume, not on carbon density, not on inoculum, not even on the
# fitted intercepts. There is nothing left in it to get wrong.
#
# INTERNAL CHECK: min_T [R/G] occurs where dlnR/dT = dlnG/dT, which is exactly
# Topt(CUE) from 13. The script asserts this - if it ever fails, the two scripts
# have drifted apart and BOTH are wrong.
#
# WHY IT BEATS Topt(CUE) AS A STATISTIC
# -------------------------------------
# Topt(CUE) captures only WHERE the optimum is. Tax captures where it is AND how
# steeply performance falls away from it. Two taxa can share an optimum and pay
# wildly different penalties at 37 C if one curve is flat and the other is sharp.
# In this dataset that distinction is the whole result: C. auris Clade I has an
# optimum no higher than C. parapsilosis, but a much FLATTER curve - so it pays
# far less to sit at human body temperature.
#
# CAVEATS (inherited from 13, unchanged)
# --------------------------------------
#   * Two-stage posterior: growth and respiration were fitted separately, so
#     pairing their draws assumes independence they do not strictly have. This
#     INFLATES the intervals (conservative), it does not shrink them.
#   * n = 3 isolates per clade. Clade-level contrasts have limited power.
#   * Tax at 44 C is a steep extrapolation of the SHAPE (though still inside the
#     measured temperature range) and glabrata's Eh is poorly identified, so its
#     44 C interval will be very wide. Do not lead with the 44 C column.
#
# INPUTS   models/bayes_growth_ss.rds, models/bayes_resp_arr.rds
#          models/bayes_data_growth.rds
# OUTPUTS  tables/carbon_tax_group.csv          (tax at 37 / 40 / 42 / 44 C)
#          tables/carbon_tax_isolate.csv
#          tables/carbon_tax_curves.csv
#          tables/carbon_tax_contrasts_fever.csv   (auris vs comparators at 40 C)
#          figures/FIG_carbon_tax.png / .pdf
#
# RUN AFTER 08_bayesian_models.R (13 is not required, but they should agree).
# =============================================================================

# C1: --file= FIRST. Under `Rscript scripts/<this>.R` neither sys.frame(1)$ofile
# nor rstudioapi resolves, so this fell back to getwd() and then died on
# "cannot open file .../config.R". Sourcing from run_all.R was unaffected,
# which is why the bug stayed hidden.
.this_dir <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) return(dirname(normalizePath(
    # R replaces every space in --file= with "~+~". A project path with a
    # space in it (this one has two) therefore comes back mangled, and
    # source() then fails on a path that does not exist.
    gsub("~+~", " ", sub("^--file=", "", fa[1]), fixed = TRUE),
    mustWork = FALSE)))
  d <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (length(d) == 0 || is.na(d) || !nzchar(d)) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
        nzchar(rstudioapi::getActiveDocumentContext()$path)) {
      d <- dirname(rstudioapi::getActiveDocumentContext()$path)
    } else d <- getwd()
  }
  d
})
source(file.path(.this_dir, "config.R"))

need <- c("posterior", "ggplot2", "dplyr", "tidyr", "tibble", "readr", "patchwork")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) {
  stop("Install first:\n  install.packages(c(",
       paste(sprintf('"%s"', miss), collapse = ", "), "))")
}
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble); library(patchwork)
})

CINV    <- 11604.51812
T_BODY  <- 37
T_FEVER <- 40
T_MIN_C <- 22
T_MAX_C <- 44
T_REPORT <- c(37, 40, 42, 44)

message("\n=== 10_carbon_tax.R ============================================")
message("What does it cost to live in a human?")
message("Tax(T) = [R/G](T) / min[R/G]. Depends only on E, E_R, Eh, Th.\n")

req <- function(p) { if (!file.exists(p)) stop("Not found: ", p,
                     "\nRun 08_bayesian_models.R first."); p }
fit_g      <- readRDS(req(file.path(models_dir, "bayes_growth_ss.rds")))
fit_r      <- readRDS(req(file.path(models_dir, "bayes_resp_arr.rds")))
growth_dat <- readRDS(req(file.path(models_dir, "bayes_data_growth.rds")))

dg <- posterior::as_draws_df(fit_g)
dr <- posterior::as_draws_df(fit_r)
n_draws <- min(nrow(dg), nrow(dr))
dg <- dg[seq_len(n_draws), ]; dr <- dr[seq_len(n_draws), ]
message(sprintf("Posterior draws: %d", n_draws))

GRPS <- STRAIN_GROUPS[STRAIN_GROUPS %in% unique(as.character(growth_dat$Group))]
iso_map <- growth_dat %>% dplyr::distinct(Isolate, Group) %>%
  dplyr::mutate(Isolate = as.character(Isolate), Group = as.character(Group))

# never call this `get` - see the comment in 13
draw_col <- function(d, nm) {
  if (!nm %in% names(d)) return(rep(NA_real_, nrow(d)))
  as.numeric(d[[nm]])
}

# =============================================================================
# 1) THE TAX, DRAW BY DRAW
# =============================================================================
# ln(R/G) up to an additive constant that cancels in the ratio.
ln_cost <- function(TC, E, ER, Eh, Th) {
  TK <- TC + 273.15
  u  <- Eh * CINV * (1 / Th - 1 / TK)
  # log1p(exp(u)) computed stably: for large u it is just u.
  sp <- ifelse(u > 30, u, log1p(exp(pmin(u, 30))))
  (ER - E) * CINV * (1 / 293.15 - 1 / TK) + sp
}

# Topt(CUE): where R/G is cheapest. Closed form (identical to 13).
topt_cue <- function(E, ER, Eh, Th) {
  s <- (E - ER) / Eh
  out <- rep(NA_real_, length(s))
  ok <- is.finite(s) & s > 0 & s < 1
  out[ok] <- 1 / (1 / Th[ok] - stats::qlogis(s[ok]) / (Eh[ok] * CINV)) - 273.15
  out
}

tax_draws <- function(E, ER, Eh, Th, TC) {
  t0 <- topt_cue(E, ER, Eh, Th)
  # If there is no interior optimum, the cheapest temperature in the measured
  # range is the lower bound - use it rather than returning NA and hiding the curve.
  t0[!is.finite(t0)] <- T_MIN_C
  t0 <- pmin(pmax(t0, T_MIN_C), T_MAX_C)
  exp(ln_cost(TC, E, ER, Eh, Th) - ln_cost(t0, E, ER, Eh, Th))
}

grp_par <- lapply(GRPS, function(g) {
  tibble::tibble(
    Group = g,
    E  = draw_col(dg, paste0("b_E_Group",  g)),
    Eh = draw_col(dg, paste0("b_Eh_Group", g)),
    Th = draw_col(dg, paste0("b_Th_Group", g)),
    ER = draw_col(dr, paste0("b_E_Group",  g)))
}) %>% dplyr::bind_rows()

# ---- INTERNAL CHECK: does min(R/G) really sit at Topt(CUE)? ------------------
# If 13 and 14 have drifted apart this catches it immediately.
.chk <- dplyr::filter(grp_par, Group == GRPS[1])[1, ]
.t0  <- topt_cue(.chk$E, .chk$ER, .chk$Eh, .chk$Th)
if (is.finite(.t0)) {
  .grid <- seq(T_MIN_C, T_MAX_C, by = 0.01)
  .num  <- .grid[which.min(ln_cost(.grid, .chk$E, .chk$ER, .chk$Eh, .chk$Th))]
  if (abs(.num - .t0) > 0.05) {
    stop(sprintf(paste0("CONSISTENCY CHECK FAILED: closed-form Topt(CUE) = %.3f but the\n",
                        "numerical minimum of R/G is at %.3f. 13 and 14 disagree - fix before\n",
                        "trusting anything downstream."), .t0, .num))
  }
  message(sprintf("Check passed: closed-form Topt(CUE) = %.2f C = numerical argmin of R/G (%.2f C).",
                  .t0, .num))
}

q   <- function(x, p) stats::quantile(x[is.finite(x)], p, names = FALSE)
med <- function(x) stats::median(x[is.finite(x)])

# ---- tax at the reported temperatures ---------------------------------------
grp_tax <- lapply(GRPS, function(g) {
  p <- dplyr::filter(grp_par, Group == g)
  lapply(T_REPORT, function(tc) {
    tx <- tax_draws(p$E, p$ER, p$Eh, p$Th, tc)
    tibble::tibble(Group = g, T_C = tc,
                   tax = med(tx), lo = q(tx, .025), hi = q(tx, .975))
  }) %>% dplyr::bind_rows()
}) %>% dplyr::bind_rows()

readr::write_csv(grp_tax, file.path(tables_dir, "carbon_tax_group.csv"))

message("\n--- CARBON TAX  (fold-increase in respiration per unit growth,")
message("    relative to each taxon's own cheapest temperature) ------------")
wide <- grp_tax %>%
  dplyr::mutate(cell = sprintf("%.2f [%.2f-%.2f]", tax, lo, hi)) %>%
  dplyr::select(Group, T_C, cell) %>%
  tidyr::pivot_wider(names_from = T_C, values_from = cell, names_prefix = "T")
print(as.data.frame(wide), row.names = FALSE)

# ---- full curves -------------------------------------------------------------
TGRID <- seq(T_MIN_C, T_MAX_C, by = 0.25)
tax_curves <- lapply(GRPS, function(g) {
  p <- dplyr::filter(grp_par, Group == g)
  idx <- if (nrow(p) > 1500) sample.int(nrow(p), 1500) else seq_len(nrow(p))
  M <- vapply(TGRID, function(tc)
    tax_draws(p$E[idx], p$ER[idx], p$Eh[idx], p$Th[idx], tc), numeric(length(idx)))
  tibble::tibble(Group = g, T_C = TGRID,
                 tax = apply(M, 2, stats::median),
                 lo  = apply(M, 2, stats::quantile, .025, names = FALSE),
                 hi  = apply(M, 2, stats::quantile, .975, names = FALSE))
}) %>% dplyr::bind_rows()
readr::write_csv(tax_curves, file.path(tables_dir, "carbon_tax_curves.csv"))

# ---- per isolate -------------------------------------------------------------
iso_tax <- lapply(seq_len(nrow(iso_map)), function(i) {
  iso <- iso_map$Isolate[i]; g <- iso_map$Group[i]
  re <- function(d, par) {
    nm <- sprintf("r_Isolate__%s[%s,Intercept]", par, iso)
    if (nm %in% names(d)) as.numeric(d[[nm]]) else rep(0, nrow(d))
  }
  E  <- draw_col(dg, paste0("b_E_Group",  g)) + re(dg, "E")
  Th <- draw_col(dg, paste0("b_Th_Group", g)) + re(dg, "Th")
  Eh <- draw_col(dg, paste0("b_Eh_Group", g))
  ER <- draw_col(dr, paste0("b_E_Group",  g)) + re(dr, "E")
  lapply(c(T_BODY, T_FEVER), function(tc) {
    tx <- tax_draws(E, ER, Eh, Th, tc)
    tibble::tibble(Isolate = iso, Group = g, T_C = tc,
                   tax = med(tx), lo = q(tx, .025), hi = q(tx, .975))
  }) %>% dplyr::bind_rows()
}) %>% dplyr::bind_rows() %>%
  dplyr::arrange(match(Group, STRAIN_GROUPS), Isolate, T_C)
readr::write_csv(iso_tax, file.path(tables_dir, "carbon_tax_isolate.csv"))

# =============================================================================
# 2) THE CONTRAST THAT MATTERS: C. auris vs the comparators, AT FEVER
# =============================================================================
# This is the claim the paper wants to make, so it gets tested directly rather
# than eyeballed off a figure. Ratio of taxes, draw by draw.
AURIS <- intersect(c("Clade1", "Clade2", "Clade3", "Clade4"), GRPS)
COMP  <- intersect(c("glab", "para"), GRPS)

tax_at <- function(g, tc) {
  p <- dplyr::filter(grp_par, Group == g)
  tax_draws(p$E, p$ER, p$Eh, p$Th, tc)
}

fever <- expand.grid(auris = AURIS, comp = COMP, stringsAsFactors = FALSE) %>%
  dplyr::rowwise() %>%
  dplyr::mutate({
    r <- tax_at(auris, T_FEVER) / tax_at(comp, T_FEVER)
    r <- r[is.finite(r)]
    tibble::tibble(
      ratio    = stats::median(r),
      lo       = stats::quantile(r, .025, names = FALSE),
      hi       = stats::quantile(r, .975, names = FALSE),
      P_cheaper = mean(r < 1))          # P(auris pays LESS than the comparator)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(credible = hi < 1) %>%
  dplyr::arrange(ratio)

readr::write_csv(fever, file.path(tables_dir, "carbon_tax_contrasts_fever.csv"))

message("\n--- At 40 C (fever): does each C. auris clade pay LESS than the")
message("    comparator species?  ratio < 1 means cheaper. -----------------")
for (i in seq_len(nrow(fever))) {
  f <- fever[i, ]
  message(sprintf("  %-8s vs %-5s  ratio %.2f [%.2f-%.2f]   P(cheaper) = %.3f  %s",
                  f$auris, f$comp, f$ratio, f$lo, f$hi, f$P_cheaper,
                  ifelse(f$credible, "<-- credible", "")))
}

# =============================================================================
# 3) FIGURE
# =============================================================================
OI <- c(Clade1 = "#0072B2", Clade2 = "#56B4E9", Clade3 = "#009E73",
        Clade4 = "#E69F00", glab = "#CC79A7", para = "#D55E00")
LAB <- c(Clade1 = "C. auris I", Clade2 = "C. auris II", Clade3 = "C. auris III",
         Clade4 = "C. auris IV", glab = "C. glabrata", para = "C. parapsilosis")

theme_pub <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 7, colour = "grey20"),
        axis.title = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3, colour = "grey20"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey20"),
        legend.text = element_text(size = 7), legend.title = element_blank(),
        legend.key.height = unit(8, "pt"), legend.key.width = unit(12, "pt"),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key = element_rect(fill = alpha("white", 0.85), colour = NA),
        plot.title = element_text(size = 8, face = "bold", hjust = 0),
        plot.margin = margin(4, 6, 4, 4))

tax_curves$Group <- factor(tax_curves$Group, levels = GRPS)
grp_tax$Group    <- factor(grp_tax$Group,    levels = GRPS)
iso_tax$Group    <- factor(iso_tax$Group,    levels = GRPS)

# ---- Panel A: the tax curves. log10 y - a "3x" penalty must LOOK 3x. --------
pA <- ggplot(tax_curves, aes(T_C, tax, colour = Group, fill = Group)) +
  annotate("rect", xmin = T_BODY, xmax = T_FEVER, ymin = -Inf, ymax = Inf,
           fill = "#B22222", alpha = 0.07) +
  geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.35) +
  geom_vline(xintercept = c(T_BODY, T_FEVER), colour = "#B22222",
             linetype = "22", linewidth = 0.32) +
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.11) +
  geom_line(linewidth = 0.65) +
  scale_colour_manual(values = OI, labels = LAB) +
  scale_fill_manual(values = OI, labels = LAB) +
  scale_x_continuous(breaks = seq(24, 44, 4)) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50),
                labels = c("1x", "2x", "5x", "10x", "20x", "50x")) +
  coord_cartesian(xlim = c(T_MIN_C, T_MAX_C), ylim = c(0.95, 50)) +
  annotate("text", x = T_BODY - 0.45, y = 33, angle = 90, hjust = 1, size = 2.2,
           colour = "#B22222", label = "body 37°C") +
  annotate("text", x = T_FEVER + 0.9, y = 33, angle = 90, hjust = 1, size = 2.2,
           colour = "#B22222", label = "fever 40°C") +
  annotate("text", x = T_MIN_C + 0.5, y = 1.03, hjust = 0, vjust = 0, size = 2.3,
           fontface = "italic", colour = "grey30",
           label = "1x = the taxon's own cheapest temperature") +
  labs(x = "Temperature (°C)",
       y = "Carbon tax   (respiration per unit growth,\nrelative to own optimum)",
       title = "The cost of being warm") +
  guides(fill = "none", colour = guide_legend(ncol = 1)) +
  theme_pub +
  theme(legend.position = c(0.20, 0.80))

# ---- Panel B: tax at 37 and 40 C, with isolates ----------------------------
ordg <- grp_tax %>% dplyr::filter(T_C == T_FEVER) %>%
  dplyr::arrange(tax) %>% dplyr::pull(Group) %>% as.character()

bd <- grp_tax %>% dplyr::filter(T_C %in% c(T_BODY, T_FEVER)) %>%
  dplyr::mutate(Group = factor(Group, levels = ordg),
                what  = factor(ifelse(T_C == T_BODY, "37 °C (body)", "40 °C (fever)"),
                               levels = c("37 °C (body)", "40 °C (fever)")))
bi <- iso_tax %>% dplyr::filter(T_C %in% c(T_BODY, T_FEVER)) %>%
  dplyr::mutate(Group = factor(Group, levels = ordg),
                what  = factor(ifelse(T_C == T_BODY, "37 °C (body)", "40 °C (fever)"),
                               levels = c("37 °C (body)", "40 °C (fever)")))

pB <- ggplot(bd, aes(x = tax, y = Group, colour = Group)) +
  geom_vline(xintercept = 1, colour = "grey40", linewidth = 0.35) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.9, lineend = "round") +
  geom_point(data = bi, shape = 16, size = 0.9, alpha = 0.5,
             position = position_nudge(y = 0.26)) +
  geom_point(aes(fill = Group), shape = 21, size = 2.6, colour = "white", stroke = 0.55) +
  geom_text(aes(label = sprintf("%.2f", tax)), vjust = 2.2, size = 2.0,
            show.legend = FALSE) +
  facet_wrap(~ what, nrow = 1) +
  scale_colour_manual(values = OI, guide = "none") +
  scale_fill_manual(values = OI, guide = "none") +
  scale_y_discrete(labels = LAB, expand = expansion(add = c(0.7, 0.7))) +
  scale_x_log10(breaks = c(1, 1.5, 2, 3, 5), labels = c("1x", "1.5x", "2x", "3x", "5x")) +
  labs(x = "Carbon tax (log scale)", y = NULL,
       title = "C. auris pays less to sit in a fever than either comparator") +
  theme_pub +
  theme(axis.text.y = element_text(size = 7.5, face = "italic"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey93"),
        strip.background = element_blank(),
        strip.text = element_text(size = 7.5, face = "bold", colour = "grey20"))

fig <- pA / pB +
  patchwork::plot_layout(heights = c(1, 0.95)) +
  patchwork::plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

W <- 183 / 25.4; H <- 150 / 25.4
ggsave(file.path(figures_dir, "FIG_carbon_tax.png"), fig,
       width = W, height = H, dpi = 600, bg = "white")
tryCatch({
  dev_pdf <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file.path(figures_dir, "FIG_carbon_tax.pdf"), fig,
         width = W, height = H, device = dev_pdf, bg = "white")
}, error = function(e) message("PDF not written: ", conditionMessage(e)))

message("\nWritten:")
message("  tables/carbon_tax_group.csv")
message("  tables/carbon_tax_isolate.csv")
message("  tables/carbon_tax_curves.csv")
message("  tables/carbon_tax_contrasts_fever.csv")
message("  figures/FIG_carbon_tax.png / .pdf")
message("=================================================================\n")
