# =============================================================================
# 17_fig2.R  --  FIGURE 2: the carbon cost of fever (and its four supplementary panels)
# =============================================================================
#   Rscript scripts/17_fig2.R        or  source() it in RStudio
# All data, draws and panel constructors live in fig_common.R (shared with 16_fig1.R).
# =============================================================================
.this_dir <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) return(dirname(normalizePath(
    gsub("~+~", " ", sub("^--file=", "", fa[1]), fixed = TRUE), mustWork = FALSE)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
})
source(file.path(.this_dir, "fig_common.R"))
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

