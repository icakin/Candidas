# =============================================================================
# 19_fig4.R  --  FIGURE 4: the etcGEM counterfactual falsification
# =============================================================================
#   Rscript scripts/19_fig4.R        or  source() it in RStudio
#
# Four panels, drawn from tables the etcGEM layer (gem/) writes; nothing is computed
# here that is not in those tables.
#   A  measured against etcGEM-predicted thermal performance, one facet per species.
#      C. auris is the calibration target (in-sample); the three relatives are predicted
#      from their genomes with every parameter frozen (out-of-sample).
#   B  the four proteomes' predicted enzyme Tm distributions on the axis the model
#      demands: the separation it needs is wider than the whole axis, so the four
#      species collapse onto one curve. Inset: the same curves magnified.
#   C  required interspecies separation (model inversion) against the sequence-predicted
#      paired difference and against a MEASURED congeneric benchmark (Walunjkar et al.
#      2025, Mol Biol Evol 42:msaf137: S. cerevisiae vs S. uvarum, 1.6 C over 827 pairs).
#   D  the dynamic greedy search: shifting the strongest bottleneck enzymes one by one by
#      -8 C never drives a relative's 40 C growth below detection.
#
# Inputs (regenerate with `python3 gem/26_fig4_tables.py`, which flattens the etcGEM JSON):
#   gem/tables/measured_tpc_honest.csv        17_build_measured_tpc.py   (panel A, measured)
#   gem/tables/etcgem_tpc_pred_fine.csv       18_build_etcgem_tpc.py     (panel A, predicted)
#   gem/tables/fig4/tm_predicted.csv          predicted Tm per enzyme     (panel B)
#   gem/tables/fig4/paired.csv                paired ortholog differences (panels B, C)
#   gem/tables/fig4/requirements.csv          model-required separations  (panels B, C)
#   gem/tables/fig4/dyn_sparse.csv            greedy-search trajectories  (panel D)
# Every printed number is pinned to its producing script in gem/FIG4_LOCKED.md.
#
# Needs only: dplyr, readr (base graphics). Same conventions as 18_fig3.R.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr) })

root <- if (file.exists("results/tables/fit_metrics.csv")) "." else
        if (file.exists("../results/tables/fit_metrics.csv")) ".." else
        path.expand("~/Desktop/Projects/Candidas")
.cfg <- file.path(root, "scripts/config.R")
if (file.exists(.cfg)) try(source(.cfg), silent = TRUE)
if (exists("fig_keep_add")) fig_keep_add("FIG4_etcgem_counterfactual")

GT  <- file.path(root, "gem", "tables")
FIG <- file.path(root, "results", "figures", "manuscript")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
rd <- function(...) read_csv(file.path(GT, ...), show_col_types = FALSE, progress = FALSE)

MEAS <- rd("measured_tpc_honest.csv")
PRED <- rd("etcgem_tpc_pred_fine.csv")
TMD  <- rd("fig4", "tm_predicted.csv")
PAIR <- rd("fig4", "paired.csv")
REQ  <- rd("fig4", "requirements.csv")
DYN  <- rd("fig4", "dyn_sparse.csv")

SP4  <- c("auris", "haemulonii", "duobushaemulonii", "parapsilosis")
RELS <- SP4[-1]
LAB  <- c(auris = "C. auris", haemulonii = "C. haemulonii",
          duobushaemulonii = "C. duobushaemulonii", parapsilosis = "C. parapsilosis")
# one species palette, identical in A, B and D; C uses no species colours
SPCOL <- c(auris = "#1d3f73", haemulonii = "#4d9dc0",
           duobushaemulonii = "#e08214", parapsilosis = "#7a5195")
C_PRED <- "#b0b0b0"; C_MEAS <- "#4f4f4f"; C_REQ <- "#a5453b"
RED <- "#a5453b"; INK <- "#1a1a1a"; SEC <- "#454545"; GREY <- "#7a7a7a"
DEG <- "°"; THR <- 0.05                      # detection threshold, h^-1

# ---- the numbers the panels quote -------------------------------------------
req_at <- function(k, thr = THR) REQ$median_required[REQ$param == k & REQ$threshold == thr]
reqd <- c(Tm = req_at("Tm"), Topt = req_at("Topt"))
if (any(is.na(reqd))) stop("no finite required separation at the ", THR, " h^-1 threshold; ",
                           "the figure cannot be drawn")
# the requirement is a model inversion: no sampling error, but it moves with the detection
# threshold, and at the strictest threshold no Tm shift reaches failure at all (unbounded)
req_range <- lapply(c(Tm = "Tm", Topt = "Topt"), function(k) {
  v <- REQ$median_required[REQ$param == k]
  list(lo = min(v, na.rm = TRUE), hi = max(v, na.rm = TRUE), unbounded = any(is.na(v)))
})
best <- function(k) { p <- PAIR[PAIR$param == k, ]; p[which.max(p$mean), ] }
bTm <- best("Tm"); bTo <- best("Topt")
dpair   <- bTm$mean                                  # largest paired dTm, auris's favour
foldgap <- reqd["Tm"] / dpair
MEAS_TM <- 1.6                                       # Walunjkar et al. 2025, 827 ortholog pairs

# ---- helpers ----------------------------------------------------------------
cx <- function(pt) pt / 12                           # par(ps = 12): cex in points
lw <- function(pt) pt / 0.75
clean <- function(sides = c(1, 2)) box(bty = "n")
tm_lab   <- expression(italic(T)[m]); topt_lab <- expression(italic(T)[opt])
# panel heading, placed in the outer margin of the current panel's figure region
heading <- function(txt, dy = 0.02) {
  u <- par("usr"); f <- par("fig"); pl <- par("plt")
  # x at the left edge of the panel's plot region, y just above it, in device coords
  par(xpd = NA)
  text(grconvertX(f[1] + pl[1] * (f[2] - f[1]), "ndc", "user"),
       grconvertY(f[3] + pl[4] * (f[4] - f[3]) + dy, "ndc", "user"),
       txt, adj = c(0, 0), cex = cx(9.6), font = 2, col = INK)
  par(xpd = FALSE)
}

draw_fig4 <- function() {
  par(ps = 12, family = "sans", xaxs = "i", yaxs = "i", lend = "butt", ljoin = "round")

  # ======================= A: measured vs predicted TPCs ======================
  # 2 x 2 facets inside the top-left quadrant
  ax <- c(0.085, 0.475); ay <- c(0.555, 0.875)
  fw <- (ax[2] - ax[1]) / 2; fh <- (ay[2] - ay[1]) / 2
  for (i in seq_along(SP4)) {
    sp <- SP4[i]; col <- (i - 1) %% 2; row <- (i - 1) %/% 2      # row 0 = top
    par(fig = c(ax[1] + col * fw + 0.012, ax[1] + (col + 1) * fw - 0.006,
                ay[2] - (row + 1) * fh + 0.028, ay[2] - row * fh - 0.012),
        mar = c(0, 0, 0, 0), new = (i > 1))
    plot.new(); plot.window(xlim = c(21, 45), ylim = c(-0.045, 0.92))
    abline(h = THR, lty = 3, lwd = lw(0.7), col = "#9a9a9a")
    pp <- PRED[PRED$species == sp, ]; pp <- pp[order(pp$T), ]
    lines(pp$T, pp$pred_mu, col = GREY, lwd = lw(1.3), lty = 2)
    mm <- MEAS[MEAS$species == sp & !is.na(MEAS$mu_honest), ]; mm <- mm[order(mm$T), ]
    lines(mm$T, mm$mu_honest, col = SPCOL[sp], lwd = lw(2.0))
    # every assayed temperature gets a point, zeros included: observed failures, not
    # missing data. The red cross means EVERY well was dead at that temperature.
    points(mm$T, mm$mu_honest, pch = 16, cex = 0.55, col = SPCOL[sp])
    dead <- mm[mm$n_alive == 0, ]
    if (nrow(dead)) points(dead$T, rep(0, nrow(dead)), pch = 4, cex = 0.9, lwd = lw(1.4), col = RED)
    axis(1, at = c(22, 30, 38, 44), labels = if (row == 1) c(22, 30, 38, 44) else FALSE,
         cex.axis = cx(7), padj = -1.6, tck = -0.03, lwd = lw(0.8), col = INK)
    axis(2, at = c(0, 0.4, 0.8), labels = if (col == 0) c("0", "0.4", "0.8") else FALSE,
         cex.axis = cx(7), las = 1, hadj = 0.6, tck = -0.03, lwd = lw(0.8), col = INK)
    box(bty = "l", lwd = lw(0.8), col = INK)
    par(xpd = NA)
    text(21, 0.92 + 0.02 * 0.965, bquote(italic(.(LAB[sp]))~.(if (sp == "auris") "  (calibration target)" else "")),
         adj = c(0, 0), cex = cx(7.6), col = INK)
    if (i == 1) {
      text(22.0, 0.885, "measured", adj = c(0, 1), cex = cx(7), col = SPCOL["auris"])
      text(22.0, 0.745, "etcGEM",   adj = c(0, 1), cex = cx(7), col = GREY)
    }
    if (i == 2) text(44.5, 0.885, "× = every well dead", adj = c(1, 1), cex = cx(7), col = RED)
    if (i == 3) {
      text(21 - 3.3, 0.92 + 0.10, expression("growth rate (h"^-1*")"), adj = c(0.5, 0), srt = 90, cex = cx(7.6))
      text(45 + 1.5, -0.045 - 0.34, paste0("temperature (", DEG, "C)"), adj = c(0.5, 1), cex = cx(7.6))
    }
    par(xpd = FALSE)
  }
  par(fig = c(ax, ay), mar = c(0, 0, 0, 0), new = TRUE); plot.new(); plot.window(0:1, 0:1)
  par(xpd = NA)
  text(0, 1.105, "A  The model fits its calibration target but misses the relatives’ collapse",
       adj = c(0, 0), cex = cx(9.6), font = 2, col = INK)
  par(xpd = FALSE)

  # ======================= B: predicted Tm distributions ======================
  par(fig = c(0.555, 0.975, 0.555, 0.875), mar = c(2.6, 2.4, 0, 0), mgp = c(1.5, 0.45, 0), new = TRUE)
  vals <- split(TMD$pred_tm, TMD$species)[SP4]
  meds <- sapply(vals, median); med_lo <- min(meds); med_hi <- max(meds)
  kde  <- lapply(vals, density, adjust = 1, n = 1024)
  X0 <- 44.0; X1 <- 44.0 + reqd["Tm"] + 10.0
  pk <- max(sapply(kde, function(k) max(k$y)))
  plot.new(); plot.window(xlim = c(X0, X1), ylim = c(0, pk * 1.50))
  for (sp in SP4) lines(kde[[sp]]$x, kde[[sp]]$y, col = SPCOL[sp], lwd = lw(1.6))
  yb <- pk * 1.30
  arrows(med_lo, yb, med_lo + reqd["Tm"], yb, code = 3, length = 0.06, angle = 25, lwd = lw(1.6), col = RED)
  text(med_lo + reqd["Tm"] / 2, yb + pk * 0.05,
       sprintf("separation the model requires: %.0f %sC", reqd["Tm"], DEG),
       adj = c(0.5, 0), cex = cx(7.8), col = RED, font = 2)
  text(X0 + (X1 - X0) * 0.435, pk * 0.235, adj = c(0, 0.5), cex = cx(7.2), col = INK,
       bquote(atop("mean paired "*Delta*italic(T)[m]*" ("*italic("C. auris")*" − "*italic(.(LAB[bTm$relative]))*")",
                   "= "*.(sprintf("%.2f", dpair))*" "*.(DEG)*"C  (95% CI "*.(sprintf("%.2f", bTm$lo))*"–"*
                   .(sprintf("%.2f", bTm$hi))*"; n = "*.(bTm$n)*")")))
  axis(1, cex.axis = cx(8), padj = -0.8, tck = -0.02, lwd = lw(0.8), col = INK)
  axis(2, cex.axis = cx(8), las = 1, hadj = 0.8, tck = -0.02, lwd = lw(0.8), col = INK)
  box(bty = "l", lwd = lw(0.8), col = INK)
  mtext(bquote("predicted enzyme "*italic(T)[m]*" ("*.(DEG)*"C)"), side = 1, line = 1.5, cex = cx(9))
  mtext("density", side = 2, line = 1.55, cex = cx(9))
  par(xpd = NA)
  text(X0, pk * 1.50 * 1.10, adj = c(0, 0), cex = cx(8.8), font = 2, col = INK,
       bquote(bold("B  Paired "*Delta*italic(T)[m]*" is directionally consistent but ~"*
                   .(sprintf("%.0f", foldgap))*"× below the fitted requirement")))
  par(xpd = FALSE)
  # inset: the same curves magnified. They stay superimposed.
  fB <- par("fig")
  par(fig = c(fB[1] + 0.415 * (fB[2] - fB[1]) + 0.02, fB[1] + 0.98 * (fB[2] - fB[1]),
              fB[3] + 0.30 * (fB[4] - fB[3]) + 0.03, fB[3] + 0.80 * (fB[4] - fB[3]) + 0.01),
      mar = c(1.2, 0.2, 0.9, 0.2), new = TRUE)
  xi <- seq(47, 63, length.out = 500)
  yi <- max(sapply(kde, function(k) max(approx(k$x, k$y, xi, rule = 2)$y)))
  plot.new(); plot.window(xlim = c(47, 63), ylim = c(0, yi * 1.05))
  rect(med_lo, 0, med_hi, yi * 1.05, col = adjustcolor("#333333", 0.18), border = NA)
  for (sp in SP4) lines(xi, approx(kde[[sp]]$x, kde[[sp]]$y, xi, rule = 2)$y, col = SPCOL[sp], lwd = lw(1.5))
  axis(1, at = seq(48, 62, 2), cex.axis = cx(7.2), padj = -1.4, tck = -0.03, lwd = lw(0.6), col = INK)
  axis(1, at = c(47, 63), labels = FALSE, tck = 0, lwd = lw(0.6), col = INK)
  mtext(sprintf("zoomed view — medians span %.2f %sC", med_hi - med_lo, DEG), side = 3, line = 0.15,
        cex = cx(6.8), col = SEC)
  legend("topright", inset = c(-0.02, -0.02), bty = "n", cex = cx(6.0), lwd = lw(1.5), seg.len = 1.0,
         y.intersp = 0.85, x.intersp = 0.5, col = SPCOL[SP4],
         legend = sapply(SP4, function(s) sprintf("%s (%.2f)", LAB[s], meds[s])), text.font = 3)

  # ======================= C: required vs predicted vs measured ===============
  par(fig = c(0.085, 0.475, 0.085, 0.405), mar = c(2.6, 5.2, 0, 0.4), mgp = c(1.5, 0.45, 0), new = TRUE)
  C_XLO <- 0.12
  pred <- c(Tm = dpair, Topt = bTo$mean)
  ys <- c(Tm = 1, Topt = 0)
  plot.new(); plot.window(xlim = log10(c(C_XLO, 55)), ylim = c(-0.62, 1.45))
  L10 <- log10
  for (k in c("Tm", "Topt")) {
    y <- ys[k]; b <- if (k == "Tm") bTm else bTo
    segments(L10(pred[k]), y, L10(reqd[k]), y, col = "#cfcfcf", lwd = lw(1.4), lty = 3)
    # sampling CI on the paired prediction (open to the left when it includes zero)
    if (b$lo <= 0) {
      segments(L10(C_XLO * 1.55), y, L10(b$hi), y, col = "#6f6f6f", lwd = lw(1.5))
      arrows(L10(C_XLO * 1.60), y, L10(C_XLO * 1.06), y, length = 0.05, angle = 25, col = "#6f6f6f", lwd = lw(1.3))
      text(L10(C_XLO * 1.60), y - 0.22, "interval includes zero", adj = c(0, 1), cex = cx(6.4),
           col = "#6f6f6f", font = 3)
      segments(L10(b$hi), y - 0.075, L10(b$hi), y + 0.075, col = "#6f6f6f", lwd = lw(1.5))
    } else {
      segments(L10(b$lo), y, L10(b$hi), y, col = "#6f6f6f", lwd = lw(1.5))
      segments(L10(c(b$lo, b$hi)), y - 0.075, L10(c(b$lo, b$hi)), y + 0.075, col = "#6f6f6f", lwd = lw(1.5))
    }
    # the requirement's range over detection thresholds (open when unbounded)
    rr <- req_range[[k]]; cr <- adjustcolor(C_REQ, 0.6)
    segments(L10(rr$lo), y, L10(rr$hi), y, col = cr, lwd = lw(1.5))
    segments(L10(rr$lo), y - 0.075, L10(rr$lo), y + 0.075, col = cr, lwd = lw(1.5))
    if (rr$unbounded) {
      arrows(L10(rr$hi), y, L10(rr$hi * 1.42), y, length = 0.05, angle = 25, col = cr, lwd = lw(1.5))
      par(xpd = NA)
      text(L10(rr$hi * 1.5), y - 0.155, sprintf("unbounded at %.2f h⁻¹", min(REQ$threshold)),
           adj = c(0.5, 1), cex = cx(7.2), col = C_REQ)
      par(xpd = FALSE)
    } else segments(L10(rr$hi), y - 0.075, L10(rr$hi), y + 0.075, col = cr, lwd = lw(1.5))
    points(L10(pred[k]), y, pch = 21, cex = 1.75, bg = C_PRED, col = "#222222", lwd = lw(0.8))
    points(L10(reqd[k]), y, pch = 21, cex = 1.75, bg = C_REQ,  col = "#222222", lwd = lw(0.8))
    text(L10(pred[k]), y + 0.13, sprintf("%.2f %sC", pred[k], DEG), adj = c(0.5, 0), cex = cx(8), col = SEC)
    text(L10(reqd[k]), y + 0.13, sprintf("%.0f %sC", reqd[k], DEG), adj = c(0.5, 0), cex = cx(8.5), col = RED, font = 2)
    if (k == "Tm") {
      points(L10(MEAS_TM), y, pch = 23, cex = 1.8, bg = C_MEAS, col = "#222222", lwd = lw(0.8))
      text(L10(MEAS_TM), y + 0.14, sprintf("%.1f %sC", MEAS_TM, DEG), adj = c(0.5, 0), cex = cx(8), col = C_MEAS, font = 2)
      text(L10(sqrt(MEAS_TM * reqd[k])), y - 0.21, sprintf("≈%.0f×", reqd[k] / MEAS_TM), adj = c(0.5, 1),
           cex = cx(10), font = 2)
      text(L10(MEAS_TM), y - 0.40, "S. cerevisiae / S. uvarum,\n8 °C apart in growth limit",
           adj = c(0.5, 1), cex = cx(6.4), col = SEC, font = 3)
    } else {
      text(L10(sqrt(pred[k] * reqd[k])), y - 0.20, sprintf("≈%.0f×", reqd[k] / pred[k]), adj = c(0.5, 1),
           cex = cx(10), font = 2)
    }
  }
  at <- c(0.2, 0.5, 1, 2, 5, 10, 20, 50)
  axis(1, at = L10(at), labels = at, cex.axis = cx(8), padj = -0.8, tck = -0.02, lwd = lw(0.8), col = INK)
  axis(2, at = c(0, 1), labels = FALSE, tck = 0, lwd = 0)
  par(xpd = NA)
  text(L10(C_XLO) - 0.06, 0, expression(enzyme~italic(T)[opt]), adj = c(1, 0.5), cex = cx(9.5))
  text(L10(C_XLO) - 0.06, 1, expression(enzyme~italic(T)[m]),   adj = c(1, 0.5), cex = cx(9.5))
  par(xpd = FALSE)
  axis(1, at = L10(c(C_XLO, 55)), labels = FALSE, tck = 0, lwd = lw(0.8), col = INK)
  mtext(paste0("interspecies separation (", DEG, "C, log)"), side = 1, line = 1.5, cex = cx(9))
  par(xpd = NA)
  .lg <- c("sequence-predicted (paired orthologs)", "model-required", "measured congeneric benchmark")
  legend(L10(C_XLO) - 0.05, 1.45 + 0.20, xjust = 0, yjust = 0, horiz = TRUE, bty = "n", cex = cx(7.4),
         pch = c(21, 21, 23), pt.bg = c(C_PRED, C_REQ, C_MEAS), col = "#222222", pt.cex = 1.3,
         x.intersp = 0.6, text.width = strwidth(.lg, cex = cx(7.4)) * 1.08, legend = .lg)
  text(L10(C_XLO) - 0.05, 1.45 + 0.52, "C  Required shifts exceed predictions and a measured benchmark",
       adj = c(0, 0), cex = cx(9.6), font = 2, col = INK)
  par(xpd = FALSE)

  # ======================= D: dynamic greedy search ===========================
  par(fig = c(0.555, 0.975, 0.085, 0.405), mar = c(2.6, 2.8, 0, 0.4), mgp = c(1.5, 0.45, 0), new = TRUE)
  nmax <- max(DYN$n)
  plot.new(); plot.window(xlim = c(0, nmax + 1), ylim = c(0, 0.78))
  for (rel in RELS) for (k in c("Tm", "Topt")) {
    tr <- DYN[DYN$relative == rel & DYN$param == k, ]; tr <- tr[order(tr$n), ]
    lines(tr$n, tr$mu40, col = SPCOL[rel], lwd = lw(1.9), lty = if (k == "Tm") 1 else 2)
    ev <- tr[seq(1, nrow(tr), by = 10), ]
    points(ev$n, ev$mu40, pch = if (k == "Tm") 16 else 17, cex = 0.7, col = SPCOL[rel])
  }
  abline(h = THR, lty = 2, lwd = lw(1.0), col = "#6f6f6f")
  text(nmax * 0.97, 0.075, "detection (must reach)", adj = c(1, 0), cex = cx(7.2), col = SEC)
  text((nmax + 1) / 2, 0.78 * 0.965, "no pair or triple in the beam search crossed detection",
       adj = c(0.5, 1), cex = cx(7.6), col = INK)
  axis(1, cex.axis = cx(8), padj = -0.8, tck = -0.02, lwd = lw(0.8), col = INK)
  axis(2, cex.axis = cx(8), las = 1, hadj = 0.8, tck = -0.02, lwd = lw(0.8), col = INK)
  box(bty = "l", lwd = lw(0.8), col = INK)
  mtext(paste0("cumulative no. of bottleneck enzymes shifted −8 ", DEG, "C"), side = 1, line = 1.5, cex = cx(9))
  mtext(bquote("predicted "*mu*" at 40 "*.(DEG)*"C (h"^-1*")"), side = 2, line = 1.6, cex = cx(9))
  par(xpd = NA)
  legend(-0.3, 0.78 + 0.02, xjust = 0, yjust = 0, horiz = TRUE, bty = "n", cex = cx(7.0), seg.len = 1.7,
         x.intersp = 0.5, text.width = c(0.16, 0.20, 0.16, 0.10, 0.10) * (nmax + 1),
         col = c(SPCOL[RELS], "#555555", "#555555"), lwd = lw(1.7), lty = c(1, 1, 1, 1, 2),
         legend = c(LAB[RELS], expression(italic(T)[m]*" shift"), expression(italic(T)[opt]*" shift")),
         text.font = c(3, 3, 3, 1, 1))
  text(-0.3, 0.78 + 0.145, "D  Targeted perturbations do not drive relatives below detection",
       adj = c(0, 0), cex = cx(9.6), font = 2, col = INK)
  par(xpd = FALSE)

  # ======================= title ==============================================
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE); plot.new(); plot.window(0:1, 0:1)
  text(0.5, 0.963, paste("Within the etcGEM, sequence-predicted enzyme thermal properties do not",
                         "reproduce the observed thermal divergence"),
       adj = c(0.5, 0.5), cex = cx(11.6), font = 2, col = INK)
}

# ---- write ------------------------------------------------------------------
ok_to_write <- function(f) !exists(".fig_keep_ok") || isTRUE(.fig_keep_ok(f))
png_path <- file.path(FIG, "FIG4_etcgem_counterfactual.png")
pdf_path <- file.path(FIG, "FIG4_etcgem_counterfactual.pdf")
W <- 11.5; H <- 8.6                                   # matches the Python original's canvas
.dev_type <- if (capabilities("aqua")) "quartz" else if (capabilities("cairo")) "cairo" else "Xlib"
if (ok_to_write(png_path)) {
  png(png_path, width = W, height = H, units = "in", res = 300, type = .dev_type, bg = "white")
  draw_fig4(); dev.off()
}
if (ok_to_write(pdf_path)) {
  if (capabilities("aqua")) quartz(file = pdf_path, type = "pdf", width = W, height = H) else
  if (capabilities("cairo")) cairo_pdf(pdf_path, width = W, height = H) else
    pdf(pdf_path, width = W, height = H)
  draw_fig4(); dev.off()
}
cat("wrote", png_path, "and", pdf_path, "\n")
cat(sprintf("  panel B: largest paired dTm %.3f C (%s, n = %d, p = %.1e) vs required %.2f C -> %.0fx\n",
            dpair, bTm$relative, bTm$n, bTm$p, reqd["Tm"], foldgap))
cat(sprintf("  panel C: Topt paired %.3f C vs required %.2f C -> %.0fx; measured benchmark %.1f C -> %.0fx\n",
            bTo$mean, reqd["Topt"], reqd["Topt"] / bTo$mean, MEAS_TM, reqd["Tm"] / MEAS_TM))
cat(sprintf("  panel D: greedy search to n = %d enzymes; min mu40 reached %.3f h^-1 (detection %.2f)\n",
            max(DYN$n), min(DYN$mu40), THR))
