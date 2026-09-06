# =============================================================================
# 17_fig3.R - FIGURE 3: high-temperature growth vs phylogenomic proximity
#
# Canonical generator. Reproduces FIG3_v8 exactly and replaces the earlier
# port of phylo/make_fig3.py (the febrile-persistence figure), which is retired.
#
# Everything in the figure is computed here from the pipeline tables. Nothing is
# read from a hand-curated exclusion list.
#
# DETECTION CRITERION. Growth is inferred from the dissolved-oxygen trace, not
# measured independently. A well is scored positive when BOTH hold:
#   (1) curvature   rt = r * fit-window length  >=  the cut, and
#   (2) signal      the well consumed >= MIN_O2_DRAWDOWN mg/L of O2 over the trace.
# Test (2) is the check declared in 07_oxygen_fits.R but never applied there; it
# is applied in 04_trim_selector.R only, as the interactive `no_signal` flag.
# Without it the fitter returns confident curvature on flat, dead wells, which is
# what produced the spurious non-monotonic rt at 42-44 C in earlier drafts.
# A well failing either test is a NEGATIVE, not a missing observation, so the
# denominators stay at every sampled isolate.
#
# THE CUT is the midpoint of the largest gap in the pooled isolate distribution.
# It is derived from these data, but it never sees the species labels, so it
# cannot have been chosen to favour the contrast. Panel C therefore also reports
# the counts at the pipeline's own pre-existing CURV_MIN_RT = 0.5.
#
# Needs only: dplyr, readr  (base graphics; no ggplot2 / patchwork / ape).
# Run:  Rscript scripts/17_fig3.R      or  source() it in RStudio.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr) })

root <- if (file.exists("results/tables/fit_metrics.csv")) "." else
        if (file.exists("../results/tables/fit_metrics.csv")) ".." else
        path.expand("~/Desktop/Projects/Candidas")

# The pipeline gates figure writes through a FIG_KEEP whitelist installed by
# config.R (it monkey-patches ggsave/pdf). Source it if present so the gate
# exists, then declare this figure. Harmless when run stand-alone.
.cfg <- file.path(root, "scripts/config.R")
if (file.exists(.cfg)) try(source(.cfg), silent = TRUE)
if (exists("fig_keep_add")) fig_keep_add("FIG3_discordance")

TAB <- file.path(root, "results", "tables")
FIG <- file.path(root, "results", "figures", "manuscript")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

MIN_O2_DRAWDOWN <- 2.0     # mg/L; matches 04_trim_selector.R
CURV_MIN_RT     <- 0.5     # the pipeline's pre-existing curvature criterion
TFOCUS          <- 40      # focal temperature for the contrast
TMAX            <- 44      # assay ceiling
R2_MIN          <- 0.90
N_MIN           <- 50

SPP  <- c("auris", "duobushaemulonii", "haemulonii", "parapsilosis")
TIP  <- c(auris = "auris_cladeI", duobushaemulonii = "duobushaemulonii",
          haemulonii = "haemulonii", parapsilosis = "parapsilosis")
LAB  <- c(auris = "C. auris", duobushaemulonii = "C. duobushaemulonii",
          haemulonii = "C. haemulonii", parapsilosis = "C. parapsilosis")

# ---- 1. per-well curvature, with the pipeline's technical QC -----------------
L   <- read_csv(file.path(TAB, "fit_coefficients_long.csv"), show_col_types = FALSE)
MET <- read_csv(file.path(TAB, "fit_metrics.csv"),           show_col_types = FALSE)
DER <- read_csv(file.path(TAB, "derived_N0_R_results_with_carbon.csv"), show_col_types = FALSE)

lut <- DER %>% distinct(OTU, otu_name)

species_of <- function(nm) dplyr::case_when(
  grepl("^Clade", nm) ~ "auris",
  grepl("^Hae",   nm) ~ "haemulonii",
  grepl("^Duo",   nm) ~ "duobushaemulonii",
  grepl("^para",  nm) ~ "parapsilosis",
  TRUE                ~ NA_character_)

W <- L %>%
  filter(parameter == "r") %>%
  mutate(Replicate = toupper(Replicate)) %>%
  select(-any_of(c("fit_valid", "keep", "has_curvature", "delta_aicc_curv"))) %>%
  left_join(lut, by = "OTU") %>%
  left_join(MET %>% mutate(Replicate = toupper(Replicate)) %>%
              select(T, OTU, Replicate, r2, n, r_at_bound, K_at_bound, fit_valid),
            by = c("T", "OTU", "Replicate")) %>%
  filter(fit_valid, !r_at_bound, !K_at_bound, r2 >= R2_MIN, n >= N_MIN) %>%
  mutate(rt = Estimate * T_end_min, sp = species_of(otu_name)) %>%
  filter(!is.na(sp))

# ---- 2. metabolic signal: O2 consumed from the trace peak onwards ------------
# Measured on the TRACE, not inside the fit window: "dead culture" is a property
# of the curve, not of whichever window the trimming step happened to pick.
O2 <- read_csv(file.path(TAB, "Oxygen_Data_Filtered.csv"), show_col_types = FALSE)
DD <- O2 %>%
  filter(is.finite(Oxygen), is.finite(Time)) %>%
  mutate(Replicate = toupper(Replicate)) %>%
  arrange(T, OTU, Replicate, Time) %>%
  group_by(T, OTU, Replicate) %>%
  summarise(drawdown = { pk <- which.max(Oxygen); Oxygen[pk] - min(Oxygen[pk:n()]) },
            .groups = "drop")

W <- W %>%
  left_join(DD, by = c("T", "OTU", "Replicate")) %>%
  mutate(live = !is.na(drawdown) & drawdown >= MIN_O2_DRAWDOWN,
         rt_live = ifelse(live, rt, 0))          # no signal -> no growth

# ---- 3. isolate level: best well per isolate per temperature -----------------
ISO <- W %>% group_by(sp, otu_name, T) %>%
  summarise(rt_live = max(rt_live), .groups = "drop")

F40 <- ISO %>% filter(T == TFOCUS)
A <- F40$rt_live[F40$sp == "auris"]; B <- F40$rt_live[F40$sp != "auris"]
n1 <- length(A); n2 <- length(B)

PMW  <- suppressWarnings(wilcox.test(A, B, exact = FALSE))$p.value
CLES <- mean(outer(A, B, function(a, b) (a > b) + 0.5 * (a == b)))

v   <- sort(F40$rt_live)
gi  <- which.max(diff(v))
THR <- (v[gi] + v[gi + 1]) / 2                   # label-blind, largest-gap cut

k1 <- sum(A >= THR); k2 <- sum(B >= THR)
PF <- fisher.test(matrix(c(k1, n1 - k1, k2, n2 - k2), 2, byrow = TRUE))$p.value
kS <- sum(B >= CURV_MIN_RT)
PFS <- fisher.test(matrix(c(n1, 0, kS, n2 - kS), 2, byrow = TRUE))$p.value

wilson <- function(k, n, z = qnorm(0.975)) {
  p <- k / n; d <- 1 + z^2 / n
  c((p + z^2 / (2 * n)) / d - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d,
    (p + z^2 / (2 * n)) / d + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d)
}
w1 <- wilson(k1, n1); w2 <- wilson(k2, n2)       # Newcombe hybrid score
RD   <- k1 / n1 - k2 / n2
RDLO <- RD - sqrt((k1 / n1 - w1[1])^2 + (w2[2] - k2 / n2)^2)
RDHI <- RD + sqrt((w1[2] - k1 / n1)^2 + (k2 / n2 - w2[1])^2)

# per-species summaries
cnt <- sapply(SPP, function(s) sum(F40$rt_live[F40$sp == s] >= THR))
tot <- sapply(SPP, function(s) sum(F40$sp == s))
fr  <- cnt / tot
# Per-ISOLATE highest temperature with growth, kept as the individual values.
# Panel B previously showed only their median, which reads as a species limit: it
# put 38 C against C. haemulonii while one C. haemulonii isolate in fact grew at
# 44 C, and one C. parapsilosis isolate at 42 C. Those isolates are now drawn.
iso_upper <- lapply(SPP, function(s) {
  isos <- unique(ISO$otu_name[ISO$sp == s])
  sapply(isos, function(o) {
    d <- ISO[ISO$otu_name == o & ISO$rt_live >= CURV_MIN_RT, ]
    if (nrow(d)) max(d$T) else 34 })
})
names(iso_upper) <- SPP
upper <- sapply(iso_upper, median)

# Wilson intervals on each species' fraction. With n = 2, 3, 3 in the relatives
# these are extremely wide, and drawing the fractions as bare points asserted a
# precision the data does not carry. The inference elsewhere in the figure (the
# risk difference, and the threshold-free rank test) already handles the small n
# correctly; it was only the graphic that overstated it.
ci <- lapply(SPP, function(s) wilson(cnt[[s]], tot[[s]]))
names(ci) <- SPP

# Isolate-level curvature at 40 C, per species: the PRIMARY, threshold-free
# result named in the caption, which the figure did not previously show.
CURV <- lapply(SPP, function(s) sort(F40$rt_live[F40$sp == s]))
names(CURV) <- SPP
CMAX <- max(unlist(CURV)) * 1.10

cat(sprintf("THR=%.3f  MW P=%.2e  CLES=%.2f  %d/%d vs %d/%d  Fisher P=%.4f\n",
            THR, PMW, CLES, k1, n1, k2, n2, PF))
cat(sprintf("RD=%.3f [%.3f, %.3f] | at rt > %.1f: %d/%d vs %d/%d, P=%.4f\n",
            RD, RDLO, RDHI, CURV_MIN_RT, n1, n1, kS, n2, PFS))
write_csv(F40, file.path(TAB, "fig3_isolate_40C.csv"))

# ---- 4. tree: minimal newick reader, then the induced 4-tip topology ---------
read_newick <- function(path) {
  s <- gsub("[[:space:]]", "", paste(readLines(path, warn = FALSE), collapse = ""))
  s <- sub(";$", "", s); i <- 1
  tok <- function() substr(s, i, i)
  parse_len <- function() {
    if (tok() != ":") return(0)
    i <<- i + 1; j <- i
    while (i <= nchar(s) && grepl("[0-9eE.+-]", tok())) i <<- i + 1
    as.numeric(substr(s, j, i - 1))
  }
  parse_name <- function() {
    j <- i
    while (i <= nchar(s) && !grepl("[,:();]", tok())) i <<- i + 1
    substr(s, j, i - 1)
  }
  parse_node <- function() {
    kids <- list()
    if (tok() == "(") {
      i <<- i + 1
      repeat {
        kids <- c(kids, list(parse_node()))
        if (tok() == ",") { i <<- i + 1; next }
        break
      }
      i <<- i + 1                      # ")"
    }
    nm <- parse_name(); bl <- parse_len()
    list(name = nm, brlen = bl, kids = kids)
  }
  parse_node()
}

# depth = root-to-node distance; invariant under pruning, so the induced tree
# keeps the true branch lengths without any rescaling.
annotate <- function(nd, depth = 0) {
  nd$depth <- depth + nd$brlen
  nd$kids <- lapply(nd$kids, annotate, depth = nd$depth)
  nd$tips <- if (!length(nd$kids)) nd$name else unlist(lapply(nd$kids, `[[`, "tips"))
  nd
}
# drop unsampled tips and collapse the resulting single-child nodes
induce <- function(nd, keep) {
  if (!length(nd$kids)) return(if (nd$name %in% keep) nd else NULL)
  k <- Filter(Negate(is.null), lapply(nd$kids, induce, keep = keep))
  if (!length(k)) return(NULL)
  if (length(k) == 1) return(k[[1]])
  nd$kids <- k; nd$tips <- unlist(lapply(k, `[[`, "tips")); nd
}

TR <- induce(annotate(read_newick(file.path(root, "phylo", "pg_rooted.nwk"))), TIP)
TR$depth <- 0                                     # draw from the induced root
rowy <- setNames(seq_along(SPP) - 1, TIP[SPP])
yof  <- function(nd) mean(rowy[nd$tips])
segs <- list()
walk <- function(nd) {
  for (ch in nd$kids) {
    segs[[length(segs) + 1]] <<- c(nd$depth, yof(nd), nd$depth, yof(ch))   # vertical
    segs[[length(segs) + 1]] <<- c(nd$depth, yof(ch), ch$depth, yof(ch))   # horizontal
    walk(ch)
  }
}
walk(TR)
SEG  <- do.call(rbind, segs)
# tip x positions, collected in one walk of the induced tree
tipx <- c()
collect <- function(nd) if (!length(nd$kids)) tipx[[nd$name]] <<- nd$depth else
  invisible(lapply(nd$kids, collect))
collect(TR)
XMAX <- max(unlist(tipx))

# ---- 5. the figure ----------------------------------------------------------
INK <- "#1a1a1a"; MUTED <- "#4d4d4d"; FAINT <- "#dedede"; GUIDE <- "#c4c4c4"
TEAL <- "#2f9aa6"; NAVY <- "#1d3f73"; RED <- "#a5453b"
RAMP <- colorRampPalette(c("#cfe9e6", "#7fc4c9", "#3f8fb5", "#1d3f73"))(101)
ramp <- function(p) RAMP[1 + round(100 * max(0, min(1, 0.10 + 0.90 * p)))]
DEG <- "\u00b0"; DAG <- "\u2020"

draw_fig3 <- function() {
  par(ps = 12, family = "sans", xaxs = "i", yaxs = "i")
  cx <- function(pt) pt / 12
  lw <- function(pt) pt / 0.75
  YL   <- c(3.75, -1.70)
  YTOP <- 3.12
  YN   <- -0.34          # axis numbers
  YT   <- -1.22          # axis title (reversed ylim: more negative is higher)
  YS   <- -0.84          # subtitle, between title and numbers
  PY   <- c(0.315, 0.845)

  hdr <- function(xmid, title, sub) {
    text(xmid, YT, title, cex = cx(6.4), col = MUTED, adj = c(0.5, 0))
    text(xmid, YS, sub,   cex = cx(5.4), col = "#767676", adj = c(0.5, 0))
  }

  # ---- A: phylogram
  par(fig = c(0.016, 0.104, PY[1], PY[2]), mar = c(0, 0, 0, 0), new = FALSE)
  plot.new(); plot.window(xlim = c(-0.02, XMAX * 1.02), ylim = YL)
  for (s in seq_len(nrow(SEG)))
    segments(SEG[s, 1], SEG[s, 2], SEG[s, 3], SEG[s, 4], col = "#6b6b6b", lwd = lw(1.1), lend = 1)
  segments(0, 3.48, 0.2, 3.48, col = INK, lwd = lw(1.6), lend = 1, xpd = NA)
  text(0, 3.72, expression(paste("0.2 subs site"^-1)), cex = cx(5.2), col = MUTED,
       adj = c(0, 1), xpd = NA)
  for (s in SPP)                                    # tip guides, so A and B read as one table
    segments(tipx[[TIP[s]]], rowy[TIP[s]], XMAX * 1.30, rowy[TIP[s]],
             col = GUIDE, lwd = lw(0.5), lty = "13", xpd = NA)

  # ---- B: fraction growing at 40 C, with interval
  par(fig = c(0.104, 0.520, PY[1], PY[2]), new = TRUE)
  plot.new(); plot.window(xlim = c(0, 10), ylim = YL)
  X0 <- 5.1; X1 <- 8.6; CN <- 9.6
  xf <- function(p) X0 + (X1 - X0) * p
  for (v in c(0, 0.5, 1)) {
    segments(xf(v), YN + 0.12, xf(v), YTOP, col = FAINT, lwd = lw(0.6))
    text(xf(v), YN, format(v), cex = cx(5.4), col = "#767676", adj = c(0.5, 0))
  }
  hdr(mean(c(X0, X1)), paste0("Isolates growing at ", TFOCUS, " ", DEG, "C"),
      "point, fraction; bar, 95% CI")
  text(CN, YN, "n/N", cex = cx(5.4), col = "#767676", adj = c(0.5, 0))
  for (s in SPP) {
    y <- rowy[TIP[s]]
    segments(0, y, 1.05, y, col = GUIDE, lwd = lw(0.5), lty = "13")
    segments(CN + 0.42, y, 10, y, col = GUIDE, lwd = lw(0.5), lty = "13")
    lb <- if (s == "parapsilosis") bquote(bolditalic(.(LAB[[s]]))~"(outgroup)")
          else bquote(bolditalic(.(LAB[[s]])))
    text(4.90, y, lb, cex = cx(6.2), col = INK, adj = c(1, 0.5))
    lo <- max(0, ci[[s]][1]); hi <- min(1, ci[[s]][2])
    segments(xf(lo), y, xf(hi), y, col = "#5f9aa3", lwd = lw(1.2), lend = 1)
    for (e in c(lo, hi)) segments(xf(e), y - 0.13, xf(e), y + 0.13, col = "#5f9aa3", lwd = lw(1.2))
    points(xf(fr[[s]]), y, pch = 21, cex = cx(10) * 1.05, bg = ramp(fr[[s]]),
           col = "white", lwd = lw(0.8))
    text(CN, y, sprintf("%d/%d", cnt[[s]], tot[[s]]), cex = cx(5.8), col = INK, adj = c(0.5, 0.5))
  }

  # ---- C: isolate curvature at 40 C -- the primary, threshold-free result
  par(fig = c(0.560, 0.988, PY[1], PY[2]), new = TRUE)
  plot.new(); plot.window(xlim = c(-1.30, CMAX), ylim = YL)
  for (v in 0:3) if (v <= CMAX) {
    segments(v, YN + 0.12, v, YTOP, col = FAINT, lwd = lw(0.6))
    text(v, YN, v, cex = cx(5.4), col = "#767676", adj = c(0.5, 0))
  }
  segments(THR, YN + 0.12, THR, YTOP, col = "#c39a95", lwd = lw(0.7), lty = "13")
  text(THR, YTOP + 0.34, bquote(italic(rt)*" = "*.(sprintf("%.2f", THR))*" reference"), cex = cx(5.0), col = "#9c6f69", adj = c(0.5, 0))
  hdr(CMAX / 2, expression("Isolate curvature at 40 "*degree*"C, "*italic(rt)),
      "every isolate; bar, median")
  SHORT <- c(auris = "auris", duobushaemulonii = "duobushaemulonii",
             haemulonii = "haemulonii", parapsilosis = "parapsilosis")
  for (s in SPP) {                                  # labels and guides, under the marks
    y <- rowy[TIP[s]]
    text(-1.26, y, bquote(italic(.(SHORT[[s]]))), cex = cx(5.2), col = MUTED, adj = c(0, 0.5))
    segments(-0.16, y, min(CURV[[s]]) - 0.055, y, col = GUIDE, lwd = lw(0.5), lty = "13")
  }
  for (s in SPP) {
    y <- rowy[TIP[s]]; u <- CURV[[s]]
    j <- if (length(u) > 1)
      seq(-0.21, 0.21, length.out = length(u))[rank(u, ties.method = "first")] else 0
    points(u, y + j, pch = 21, cex = cx(10) * 0.76,
           bg = if (s == "auris") NAVY else "#ffffff",
           col = if (s == "auris") "white" else "#3f6670",
           lwd = if (s == "auris") lw(0.8) else lw(1.15))
    segments(median(u), y - 0.20, median(u), y + 0.20, col = "#3a3a3a", lwd = lw(1.0), lend = 1)
  }

  # ---- D: effect size (secondary, count-based)
  par(fig = c(0.104, 0.520, 0.070, 0.245), mar = c(2.1, 0, 0, 0), new = TRUE)
  plot.new(); plot.window(xlim = c(-0.10, 1.03), ylim = c(-1.5, 1.15))
  abline(v = 0, col = "#dcdcdc", lwd = lw(0.6), lty = "44")
  segments(RDLO, 0, RDHI, 0, col = "#4a4a4a", lwd = lw(1.1))
  for (xx in c(RDLO, RDHI)) segments(xx, -0.15, xx, 0.15, col = "#4a4a4a", lwd = lw(1.1))
  points(RD, 0, pch = 21, bg = NAVY, col = "white", lwd = lw(0.7), cex = cx(8) * 1.25)
  axis(1, at = c(0, 0.5, 1.0), labels = c("0", "0.5", "1"), cex.axis = cx(5.4),
       col = "#bdbdbd", col.axis = MUTED, tck = -0.07, lwd = lw(0.7), mgp = c(0, 0.22, 0))
  mtext(bquote("Risk difference at "*.(TFOCUS)*" "*degree*"C ("*italic("C. auris")*
               " - pooled non-auris)"), side = 1, line = 1.05, cex = cx(5.8), col = INK)
  text(RD, 0.62, sprintf("+%.2f  (95%% CI +%.2f to +%.2f)", RD, RDLO, RDHI),
       cex = cx(5.8), col = INK, font = 2, adj = c(0.5, 0))

  # ---- notes: one callout, not caption prose
  par(fig = c(0.545, 0.992, 0.055, 0.250), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rect(0, 0.02, 1, 0.98, col = "#faf6f5", border = "#e0cdca", lwd = lw(0.6))
  text(0.035, 0.80, "Growth at 40 \u00b0C does not track relatedness",
       cex = cx(6.0), font = 2, col = RED, adj = c(0, 0.5))
  .dw <- strwidth(expression(italic(d)), cex = cx(5.6))
  text(0.035, 0.50, expression(italic(d)), cex = cx(5.6), col = INK, adj = c(0, 0.5))
  text(0.035 + .dw, 0.50,
       " = 0.33\u20130.35 subs/site to nearest relatives;  12/12 vs 3/8 isolates grew",
       cex = cx(5.6), col = INK, adj = c(0, 0.5))
  text(0.035, 0.20, bquote("Counts: Fisher "*italic(P)*" = 0.004;  curvature: Mann-Whitney "*
                           italic(P)*" = 0.001"),
       cex = cx(5.6), col = INK, adj = c(0, 0.5))

  # ---- panel labels and title
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.016, 0.882, "A  Phylogeny", cex = cx(7.2), font = 2, col = INK, adj = c(0, 0))
  text(0.190, 0.882, "B  Thermal phenotype aligned to phylogeny", cex = cx(7.2),
       font = 2, col = INK, adj = c(0, 0))
  text(0.560, 0.882, "C  Continuous curvature contrast", cex = cx(7.2), font = 2,
       col = INK, adj = c(0, 0))
  text(0.104, 0.268, sprintf("D  %d %sC contrast (count-based)", TFOCUS, DEG),
       cex = cx(7.2), font = 2, col = INK, adj = c(0, 0))
  text(0.5, 0.942, bquote(bold("High-temperature growth in ")*bolditalic("C. auris")*
                          bold(" is discordant with phylogenomic proximity")),
       cex = cx(8.2), col = INK, adj = c(0.5, 0))
}

# Devices. The PDF must carry UTF-8 (the assay-ceiling dagger and the degree
# sign): the default pdf() device is single-byte and silently substitutes a dot
# for them, so use cairo_pdf when it is available. cairo_pdf is not the device
# config.R gates, so ask the whitelist directly before writing.
ok_to_write <- function(f) !exists(".fig_keep_ok") || isTRUE(.fig_keep_ok(f))

png_path <- file.path(FIG, "FIG3_discordance.png")
pdf_path <- file.path(FIG, "FIG3_discordance.pdf")
# Designed at final print size (183 mm full page width), NOT at a large size to be
# reduced later: cx() emits literal points, so a figure drawn 13 in wide and shrunk to
# 7.2 in turns its 7 pt labels into 3.8 pt, below every journal's minimum.
W <- 7.2; H <- 4.35

# Device choice. On macOS R reports capabilities("cairo") TRUE but the cairo DLL
# needs XQuartz (libXrender); without it png(type="cairo") fails silently and the
# old file is left in place. Use the native quartz device there, cairo elsewhere.
.dev_type <- if (capabilities("aqua")) "quartz" else if (capabilities("cairo")) "cairo" else "Xlib"
if (ok_to_write(png_path)) {
  png(png_path, width = W, height = H, units = "in", res = 300, type = .dev_type, bg = "white")
  draw_fig3(); dev.off()
}
if (ok_to_write(pdf_path)) {
  if (capabilities("aqua")) quartz(file = pdf_path, type = "pdf", width = W, height = H) else
  if (capabilities("cairo")) cairo_pdf(pdf_path, width = W, height = H) else
    pdf(pdf_path, width = W, height = H)
  draw_fig3(); dev.off()
}
cat("wrote", png_path, "and", pdf_path, "\n")
