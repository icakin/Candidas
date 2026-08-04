# =============================================================================
# 16_supplementary_figures.R  --  etc-GEM SUPPLEMENTARY FIGURES (diagnostics + validation)
# =============================================================================
# The diagnostics behind 14_main_figures.R, in one supplementary:
#   a  predicted vs observed (calibrated) -- goodness of fit, 1:1 line
#   b  residuals vs temperature -- no systematic bias
#   c  identifiability: RMSE as each knob is swept round its optimum (a clear
#      valley = the data constrain that knob; the posterior stand-in)
#   d  differential-evolution vs Bayesian capacity -- the two calibration
#      methods agree, so the capacity result is not method-dependent
#
# Data: cauris_etcgem outputs/supp_data/. Style matches 14_main_figures.R.
#   Rscript 16_supplementary_figures.R   (needs ggplot2, patchwork)
# =============================================================================
# script directory: works under Rscript (--file=), source() and RStudio
.fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.d0 <- if (length(.fa)) dirname(normalizePath(sub("^--file=", "", .fa[1]), mustWork = FALSE)) else
         tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
if (length(.d0) == 0 || is.na(.d0) || !nzchar(.d0)) {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    .d0 <- dirname(rstudioapi::getActiveDocumentContext()$path) else .d0 <- getwd()
}
# project root: nearest ancestor of the script (or cwd) holding results/ and scripts/
.find_root <- function(p) {
  p <- normalizePath(p, mustWork = FALSE)
  for (i in 0:4) {
    if (dir.exists(file.path(p, "results")) && dir.exists(file.path(p, "scripts"))) return(p)
    p <- normalizePath(file.path(p, ".."), mustWork = FALSE)
  }
  NA_character_
}
.root <- .find_root(.d0)
if (is.na(.root)) .root <- .find_root(getwd())
if (is.na(.root)) stop("Cannot locate the project root (a folder containing results/ and scripts/).")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

OI <- c(I = "#0072B2", II = "#CC79A7", III = "#009E73", IV = "#E69F00")
SHORT <- c(I = "C. auris I", II = "C. auris II", III = "C. auris III", IV = "C. auris IV")
RED <- "#B22222"; GRPS <- c("I", "II", "III", "IV")
# species name italic, clade numeral roman, in every legend and axis
SHORT_P <- setNames(paste0("italic('C. auris')~", GRPS), GRPS)
lab_sp  <- function(x) { p <- SHORT_P[as.character(x)]; p[is.na(p)] <- paste0("'", x[is.na(p)], "'"); parse(text = p) }
T_BODY <- 37; T_FEVER <- 40; T_MIN <- 22; T_MAX <- 44
XLAB_T <- expression("Temperature ("*degree*"C)")
th <- theme_classic(base_size = 8) + theme(
  axis.text = element_text(size = 6.8, colour = "grey25"),
  axis.title = element_text(size = 7.5, colour = "grey10"),
  axis.line = element_line(linewidth = .3, colour = "grey30"),
  axis.ticks = element_line(linewidth = .3, colour = "grey30"),
  legend.text = element_text(size = 6.4), legend.title = element_blank(),
  legend.key.height = unit(7.5, "pt"), legend.key.width = unit(13, "pt"),
  legend.background = element_rect(fill = alpha("white", .85), colour = NA),
  plot.title = element_text(size = 8, face = "bold", hjust = 0, colour = "grey5"),
  plot.subtitle = element_text(size = 6.8, colour = "grey35", lineheight = 1.05),
  plot.tag = element_text(size = 10, face = "bold"),
  plot.title.position = "plot", plot.margin = margin(4, 6, 3, 4))
febrile <- function() list(
  annotate("rect", xmin = T_BODY, xmax = T_FEVER, ymin = -Inf, ymax = Inf, fill = RED, alpha = .07),
  geom_vline(xintercept = c(T_BODY, T_FEVER), colour = RED, linetype = "22", linewidth = .32))

cand <- c(file.path(.root, "cauris_etcgem", "strains", "eci_cauris", "outputs", "supp_data"),
          file.path(.root, "outputs", "supp_data"),
          file.path(.d0,  "supp_data"),
          file.path(.d0,  "..", "outputs", "supp_data"),
          file.path(.d0,  "..", "cauris_etcgem", "strains", "eci_cauris", "outputs", "supp_data"),
          "supp_data")
DD <- cand[which(vapply(cand, dir.exists, logical(1)))[1]]
if (is.na(DD)) stop("supp_data/ not found. Looked in:\n  ", paste(cand, collapse = "\n  "))
FD <- file.path(.root, "results", "figures", "manuscript"); dir.create(FD, recursive = TRUE, showWarnings = FALSE)
FD <- normalizePath(FD, mustWork = FALSE)
rd <- function(f) read.csv(file.path(DD, f))
save_fig <- function(f, name, h_mm, w_mm = 183) {
  ggsave(file.path(FD, paste0(name, ".png")), f, width = w_mm/25.4, height = h_mm/25.4,
         dpi = 600, bg = "white")
  tryCatch({ dv <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
    ggsave(file.path(FD, paste0(name, ".pdf")), f, width = w_mm/25.4, height = h_mm/25.4,
           device = dv, bg = "white") }, error = function(e) NULL)
  message("  ", name, ".png + .pdf")
}

# ---- a: predicted vs observed -----------------------------------------------
pvo <- rd("predicted_vs_observed.csv"); pvo$clade <- factor(pvo$clade, levels = GRPS)
r2 <- sapply(GRPS, function(cl) { d <- pvo[pvo$clade == cl, ]
  1 - sum((d$predicted - d$observed)^2) / sum((d$observed - mean(d$observed))^2) })
lim <- c(min(c(pvo$observed, pvo$predicted)) * .9, max(c(pvo$observed, pvo$predicted)) * 1.05)
fa <- ggplot(pvo, aes(observed, predicted, colour = clade)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "22", linewidth = .4) +
  geom_point(size = 1.4, alpha = .9) +
  annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = 2.2, parse = TRUE,
           fontface = "bold", colour = "grey15",
           label = sprintf("RMSE == %.3f~h^-1", sqrt(mean(pvo$residual^2)))) +
  scale_colour_manual(values = OI, labels = lab_sp) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = expression("Observed "*italic(mu)*" (h"^-1*")"), y = expression("Predicted (h"^-1*")"),
       title = "Predictions track measurements",
       subtitle = bquote("Points on 1:1 = perfect.  "*R^2*": "*.(paste(sprintf("%.2f", r2), collapse=" / "))),
       tag = "a") +
  guides(colour = guide_legend(ncol = 1)) + th + theme(legend.position = c(.82, .27))

# ---- b: residuals vs temperature --------------------------------------------
fb <- ggplot(pvo, aes(temp_C, residual, colour = clade)) +
  febrile() + geom_hline(yintercept = 0, colour = "grey55", linewidth = .4) +
  geom_line(linewidth = .5, alpha = .7) + geom_point(size = 1.2, alpha = .9) +
  scale_colour_manual(values = OI, guide = "none") +
  scale_x_continuous(breaks = seq(24, 44, 4)) +
  labs(x = XLAB_T, y = expression("Residual (h"^-1*")"), title = "Residuals unstructured across the working range",
       subtitle = paste0("No bias 24-40 \u00B0C; all four clades deviate together at the 22 and 44 \u00B0C limits.\n",
                         "Dashed lines and shading, 37-40 \u00B0C febrile window."), tag = "b") + th

# ---- c: identifiability (3 knobs) -------------------------------------------
prof <- rd("identifiability_profiles.csv"); opt <- rd("identifiability_optima.csv")
prof$clade <- factor(prof$clade, levels = GRPS); opt$clade <- factor(opt$clade, levels = GRPS)
knlab <- c(dTopt = "dTopt shift", dCp_scale = "dCp-scale", kcat_scale = "kcat-scale (capacity)")
prof$kn <- factor(knlab[prof$knob], levels = knlab); opt$kn <- factor(knlab[opt$knob], levels = knlab)
fc <- ggplot(prof, aes(value, rmse, colour = clade)) +
  geom_line(linewidth = .6) +
  geom_point(data = opt, aes(opt_value, opt_rmse, fill = clade), shape = 21, size = 1.9,
             colour = "white", stroke = .4) +
  facet_wrap(~ kn, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = OI, labels = lab_sp) + scale_fill_manual(values = OI, guide = "none") +
  coord_cartesian(ylim = c(0, quantile(prof$rmse, .95, na.rm = TRUE))) +
  labs(x = NULL, y = expression("RMSE (h"^-1*")"),
       title = "Each knob has a well-defined optimum (dot = fit)",
       subtitle = paste0("One-dimensional profiles, remaining knobs held at the fit; not a joint identifiability test.\n",
                         "In the posterior dCp-scale and kcat-scale are correlated (r = 0.43 to 0.70 by clade)."), tag = "c") +
  guides(colour = guide_legend(nrow = 1)) +
  th + theme(strip.background = element_blank(),
             strip.text = element_text(size = 6.8, face = "bold", colour = "grey20"),
             legend.position = "top", legend.justification = "right",
             panel.spacing.x = unit(8, "pt"))

# ---- d: DE vs Bayesian capacity agree ---------------------------------------
de <- rd("capacity_bootstrap.csv")[, c("clade", "capacity")]         # DE point
bs <- rd("bayes_summary.csv"); bs <- bs[bs$param == "kcat_scale", c("clade","median","lo","hi")]
cap <- merge(de, bs, by = "clade")
ordc <- cap$clade[order(cap$capacity)]; cap$clade <- factor(cap$clade, levels = ordc)
fd <- ggplot(cap) +
  geom_pointrange(aes(median, clade, xmin = lo, xmax = hi, colour = clade),
                  linewidth = .6, size = .5) +
  geom_point(aes(capacity, clade), shape = 4, size = 2.2, stroke = .9, colour = "grey25") +
  geom_text(aes(median, clade, label = sprintf("%.2f", median)), vjust = -1.05,
            size = 2.0, colour = "grey15") +
  scale_colour_manual(values = OI, guide = "none") + scale_y_discrete(labels = lab_sp) +
  labs(x = "Effective growth-capacity (kcat-scale)", y = NULL,
       title = "Two calibration methods agree on the ordering",
       subtitle = "Dot + bar = Bayesian (median, 95% CrI);  x = differential evolution (the two agree to within 0.02)", tag = "d") +
  th + theme(axis.text.y = element_text(size = 6.8),
             axis.line.y = element_blank(), axis.ticks.y = element_blank(),
             panel.grid.major.y = element_line(linewidth = .2, colour = "grey94"))

FIG <- (fa | fb) / fc / fd + patchwork::plot_layout(heights = c(1, .95, .5)) +
  patchwork::plot_annotation(
    title = "etc-GEM model diagnostics: fit, identifiability, and method robustness",
    theme = theme(plot.title = element_text(size = 9.5, face = "bold", colour = "grey5",
                                            margin = margin(b = 3)))) &
  theme(plot.tag = element_text(size = 10, face = "bold"))

message("Writing model supplementary:")
save_fig(FIG, "FIG_MODEL_SUPP", 190)
message("FIG_MODEL_SUPP -> ", FD)


# =============================================================================
# SUPPLEMENTARY: validation / elimination figure (merged in; formerly 16_supplementary_figures.R)
# Writes FIG_MODEL_SUPP_validation. Runs after FIG_MODEL_SUPP above; re-defined
# helpers below do not affect it. Reads etc-GEM supp_data + data/expression/.
# =============================================================================
# =============================================================================
# 16_supplementary_figures.R  --  SUPPLEMENTARY validation / elimination figure
# =============================================================================
# Supports Fig 3. Shows the capacity trait is (a) clade-level, (b) the only
# calibration knob that separates clades, (c) not in the >99%-identical coding
# sequence, and (d) not in bulk transcription (public Clade I vs II RNA-seq).
#   a  capacity by clade + variance components (ICC1)          [capacity_isolates.csv]
#   b  clade pairs resolved per knob (kcat_scale/dTopt/dCp)    [bayes_summary.csv]
#   c  cross-clade proteome identity (diamond RBH, whole prot.) [hard-coded result]
#   d  metabolic vs genome log2FC, Clade I vs II, GSE165762     [capacity_expression_CladeI_vs_II.csv]
# Style matches 20/21 (Okabe-Ito, theme_classic).  Rscript 16_supplementary_figures.R
# =============================================================================
# script directory: works under Rscript (--file=), source() and RStudio
.fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.d0 <- if (length(.fa)) dirname(normalizePath(sub("^--file=", "", .fa[1]), mustWork = FALSE)) else
         tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
if (length(.d0) == 0 || is.na(.d0) || !nzchar(.d0)) {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    .d0 <- dirname(rstudioapi::getActiveDocumentContext()$path) else .d0 <- getwd()
}
# project root: nearest ancestor of the script (or cwd) holding results/ and scripts/
.find_root <- function(p) {
  p <- normalizePath(p, mustWork = FALSE)
  for (i in 0:4) {
    if (dir.exists(file.path(p, "results")) && dir.exists(file.path(p, "scripts"))) return(p)
    p <- normalizePath(file.path(p, ".."), mustWork = FALSE)
  }
  NA_character_
}
.root <- .find_root(.d0)
if (is.na(.root)) .root <- .find_root(getwd())
if (is.na(.root)) stop("Cannot locate the project root (a folder containing results/ and scripts/).")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

OI <- c(I = "#0072B2", II = "#CC79A7", III = "#009E73", IV = "#E69F00")
SHORT <- c(I = "C. auris I", II = "C. auris II", III = "C. auris III", IV = "C. auris IV")
GRPS <- c("I", "II", "III", "IV"); RED <- "#B22222"
th <- theme_classic(base_size = 8) + theme(
  axis.text = element_text(size = 6.8, colour = "grey25"),
  axis.title = element_text(size = 7.5, colour = "grey10"),
  axis.line = element_line(linewidth = .3, colour = "grey30"),
  axis.ticks = element_line(linewidth = .3, colour = "grey30"),
  legend.position = "none",
  plot.title = element_text(size = 8, face = "bold", hjust = 0, colour = "grey5"),
  plot.subtitle = element_text(size = 6.8, colour = "grey35"),
  plot.title.position = "plot", plot.margin = margin(4, 6, 3, 4))

cand <- c(file.path(.root, "cauris_etcgem", "strains", "eci_cauris", "outputs", "supp_data"),
          file.path(.root, "outputs", "supp_data"),
          file.path(.d0,  "supp_data"),
          file.path(.d0,  "..", "outputs", "supp_data"),
          file.path(.d0,  "..", "cauris_etcgem", "strains", "eci_cauris", "outputs", "supp_data"),
          "supp_data")
DD <- cand[which(vapply(cand, dir.exists, logical(1)))[1]]
if (is.na(DD)) stop("supp_data/ not found. Looked in:\n  ", paste(cand, collapse = "\n  "))
FD <- file.path(.root, "results", "figures", "manuscript"); dir.create(FD, recursive = TRUE, showWarnings = FALSE)
FD <- normalizePath(FD, mustWork = FALSE)

# ---- a  capacity is a clade-level trait (variance components) ----------------
cap <- read.csv(file.path(DD, "capacity_isolates.csv"))
cap$clade <- factor(cap$clade, levels = GRPS)
av  <- summary(aov(capacity ~ clade, cap))[[1]]
msb <- av["clade", "Mean Sq"]; msw <- av["Residuals", "Mean Sq"]; ni <- 3
icc <- (msb - msw) / (msb + (ni - 1) * msw)
Fv  <- av["clade", "F value"]; pv <- av["clade", "Pr(>F)"]
ord <- names(sort(tapply(cap$capacity, cap$clade, mean)))
cap$cl2 <- factor(cap$clade, levels = ord)
fa <- ggplot(cap, aes(cl2, capacity, colour = clade)) +
  stat_summary(fun = mean, geom = "crossbar", width = .5, linewidth = .5, aes(colour = clade)) +
  geom_jitter(width = .12, size = 1.5) +
  scale_colour_manual(values = OI) + scale_x_discrete(labels = paste("Clade", ord)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.5, size = 2.1, colour = "grey15",
           label = sprintf("between-clade variance = %.0f%% (ICC1)\nANOVA F(3,8)=%.1f, p=%.1e", 100*icc, Fv, pv)) +
  labs(x = NULL, y = "effective growth-capacity (kcat-scale)",
       title = "Effective capacity is a clade-level trait (n = 12)", subtitle = "isolate-level", tag = "a") + th

# ---- b  only capacity separates clades (non-overlapping 95% CrI) -------------
bs <- read.csv(file.path(DD, "bayes_summary.csv"))
sep_count <- function(param) {
  s <- bs[bs$param == param, ]; rownames(s) <- s$clade; n <- 0L
  for (p in combn(GRPS, 2, simplify = FALSE))
    if (s[p[1], "hi"] < s[p[2], "lo"] || s[p[2], "hi"] < s[p[1], "lo"]) n <- n + 1L
  n
}
sepd <- data.frame(knob = c("kcat_scale", "dTopt", "dCp_scale"),
                   n = vapply(c("kcat_scale","dTopt","dCp_scale"), sep_count, integer(1)))
sepd$knob <- factor(sepd$knob, levels = sepd$knob)
fb <- ggplot(sepd, aes(knob, n, fill = knob == "kcat_scale")) +
  geom_col(width = .7) + geom_text(aes(label = paste0(n, "/6")), vjust = -.5, size = 2.6, fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "#0072B2", "FALSE" = "grey75")) +
  coord_cartesian(ylim = c(0, 6.3)) +
  labs(x = NULL, y = "clade pairs resolved (of 6)",
       title = "Only capacity separates the clades", subtitle = "non-overlapping 95% CrI", tag = "b") + th

# ---- c  coding sequence >99% identical (diamond RBH; hard-coded result) ------
idd <- data.frame(clade = c("II", "III", "IV"), pid = c(99.8, 99.4, 100.0))
idd$clade <- factor(idd$clade, levels = c("II","III","IV"))
fc <- ggplot(idd, aes(clade, pid, fill = clade)) +
  geom_col(width = .7) + geom_hline(yintercept = 99, linetype = "22", colour = "grey55", linewidth = .3) +
  geom_text(aes(label = sprintf("%.1f%%", pid)), vjust = -.5, size = 2.6, fontface = "bold") +
  scale_fill_manual(values = OI[c("II","III","IV")]) +
  scale_x_discrete(labels = paste("Clade", c("II","III","IV"), "vs I")) +
  coord_cartesian(ylim = c(95, 100.4)) +
  labs(x = NULL, y = "median protein identity (%)",
       title = "Coding sequence is >99% identical", subtitle = "diamond RBH, whole proteome (~5,400 prot.)", tag = "c") + th

# ---- d  not bulk transcription either (GSE165762, Clade I vs II) -------------
nf <- file.path(.root, "data", "expression", "capacity_expression_CladeI_vs_II.csv")
if (!file.exists(nf)) nf <- file.path(.d0, "..", "data", "expression", "capacity_expression_CladeI_vs_II.csv")
if (!file.exists(nf)) nf <- file.path(.d0, "capacity_expression_CladeI_vs_II.csv")
if (file.exists(nf)) {
  de <- read.csv(nf); de$grp <- ifelse(de$is_metab, "metabolic", "genome background")
  m <- de$logFC[de$is_metab]; o <- de$logFC[!de$is_metab]
  pw <- suppressWarnings(wilcox.test(m, o)$p.value)
  fd <- ggplot(de, aes(logFC, after_stat(density), colour = grp, fill = grp)) +
    geom_density(alpha = .25, linewidth = .5) +
    geom_vline(xintercept = 0, linetype = "22", colour = "grey55", linewidth = .3) +
    scale_colour_manual(values = c("genome background" = "grey55", "metabolic" = RED)) +
    scale_fill_manual(values = c("genome background" = "grey75", "metabolic" = RED)) +
    coord_cartesian(xlim = c(-4, 4)) +
    annotate("text", x = Inf, y = Inf, hjust = 1.04, vjust = 1.5, size = 2.0, colour = "grey15",
             label = sprintf("median metab %+.3f\nmedian bkgd  %+.3f\nWilcoxon p = %.2f (n.s.)",
                             median(m), median(o), pw)) +
    labs(x = "log2 fold change (Clade I vs II)", y = "density",
         title = "Not bulk transcription either", subtitle = "metabolic enzymes track the genome", tag = "d") +
    th + theme(legend.position = c(.02, .98), legend.justification = c(0,1),
               legend.title = element_blank(), legend.text = element_text(size = 6),
               legend.key.height = unit(7, "pt"), legend.background = element_rect(fill = NA, colour = NA))
} else {
  fd <- ggplot() + annotate("text", 0, 0, size = 2.4, colour = "grey40",
        label = "run 13_capacity_expression.R first\n(writes capacity_expression_CladeI_vs_II.csv)") +
        theme_void() + labs(tag = "d")
  message("NOTE: capacity_expression_CladeI_vs_II.csv not found next to this script; panel d is a placeholder.")
}

FIG <- (fa | fb) / (fc | fd) +
  patchwork::plot_annotation(
    title = "The effective-capacity axis is clade-structured, not explained by coding variation or bulk transcription",
    theme = theme(plot.title = element_text(size = 10, face = "bold", colour = "grey5", margin = margin(b = 3)))) &
  theme(plot.tag = element_text(size = 11, face = "bold"))

ggsave(file.path(FD, "FIG_MODEL_SUPP_validation.png"), FIG, width = 183/25.4, height = 150/25.4, dpi = 600, bg = "white")
tryCatch(ggsave(file.path(FD, "FIG_MODEL_SUPP_validation.pdf"), FIG, width = 183/25.4, height = 150/25.4,
                device = if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf, bg = "white"),
         error = function(e) NULL)
message("Wrote FIG_MODEL_SUPP_validation.(png|pdf) -> ", FD)
