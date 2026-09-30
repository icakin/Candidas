# =============================================================================
# 14_fig_model.R  (ARCHIVED)  --  the etc-GEM model figure of the CUT single-proteome
# etc-GEM. Split out of scripts/14_main_figures.R on 2026-09-06 when the repository
# was reorganised. It read archive/old_etcgem/supp_data/ and wrote FIG_MODEL.png and
# FIG_MODEL_SUPP_consistency.png (now in archive/old_etcgem/figures/). It depended on
# helpers defined in the top half of 14_main_figures.R (theme `th`, palette, `.d0`)
# and is kept for the record, not for running.
# =============================================================================
# =============================================================================
# FIGURE 3 - etc-GEM model figure   (merged in; formerly 14_main_figures.R)
# Different data source: reads the etc-GEM calibration outputs from
#   cauris_etcgem/strains/eci_cauris/outputs/supp_data/ (produced by the Python
#   pipeline). Figures 1-2 above are already written by this point; if supp_data
#   is absent the Fig 3 block below stops with a clear message and 1-2 still stand.
# Runs after Fig 1/2, so its re-defined colour keys (I-IV) do not affect them.
# =============================================================================
# =============================================================================
# 14_main_figures.R  --  THE etc-GEM MODEL FIGURE (paper Fig 3)
# =============================================================================
# One publication figure, everything the sequence-grounded etc-GEM shows:
#   a  BEFORE vs AFTER calibration, per clade -- genome-only (grey dashed) vs the
#      Bayesian posterior fit (colour + 95% band) over the measured points.
#   b  SEQUENCE CAN'T EXPLAIN IT -- the genome gives ONE curve for all clades
#      (enzymes ~99.5-100% identical), yet the measured clades differ ~2x.
#   c  WHAT DIFFERS = an effective growth-capacity axis -- posterior of the kcat-scale knob per clade
#      (an effective scalar, not identified enzyme abundance).
#   d  EFFECTIVE CAPACITY COVARIES WITH ECONOMICS -- capacity (fit to GROWTH only) covaries with
#      predicts the INDEPENDENTLY MEASURED fever tax (Fig 2), n = 12 isolates.
#   e  ...and the independently measured respiration. (The old circular
#      capacity-vs-peak-growth check now lives in the supplement, see below.)
#
# Bayesian fits come from joint_calibrate_bayes.py (emcee); the a-priori curve is
# the uncalibrated genome prediction. Data: cauris_etcgem outputs/supp_data/ and
# project tables/ (carbon_tax_isolate.csv, derived_N0_R_results_with_carbon.csv).
# Style matches 14_main_figures.R (Okabe-Ito, theme_classic).
#   Rscript 14_main_figures.R   (needs ggplot2, patchwork)
# =============================================================================
# Packaging: resolve the script directory from --file= FIRST. Under
# `Rscript scripts/<this>.R` neither rstudioapi nor sys.frame(1)$ofile resolves,
# so this fell back to getwd() and then died on "cannot open file .../config.R".
# Sourcing it from run_all.R was unaffected, which is why the bug stayed hidden.
.d0 <- local({
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
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

# ---- style (identical to 14_main_figures.R) --------------------------------
OI <- c(I = "#0072B2", II = "#CC79A7", III = "#009E73", IV = "#E69F00")
SHORT <- c(I = "C. auris I", II = "C. auris II", III = "C. auris III", IV = "C. auris IV")
RED <- "#B22222"; GRPS <- c("I", "II", "III", "IV")
T_BODY <- 37; T_FEVER <- 40; T_MIN <- 22; T_MAX <- 44
XLAB_T <- expression("Temperature ("*degree*"C)")
YLAB_MU <- expression("Specific growth rate, "*italic(mu)*" (h"^-1*")")
th <- theme_classic(base_size = 8) + theme(
  axis.text = element_text(size = 6.8, colour = "grey25"),
  axis.title = element_text(size = 7.5, colour = "grey10"),
  axis.line = element_line(linewidth = .3, colour = "grey30"),
  axis.ticks = element_line(linewidth = .3, colour = "grey30"),
  legend.text = element_text(size = 6.4), legend.title = element_blank(),
  legend.key.height = unit(7.5, "pt"), legend.key.width = unit(13, "pt"),
  legend.background = element_rect(fill = alpha("white", .85), colour = NA),
  legend.key = element_rect(fill = alpha("white", .85), colour = NA),
  plot.title = element_text(size = 8, face = "bold", hjust = 0, colour = "grey5"),
  plot.subtitle = element_text(size = 6.8, colour = "grey35", lineheight = 1.05),
  plot.tag = element_text(size = 10, face = "bold"),
  plot.title.position = "plot", plot.margin = margin(4, 6, 3, 4))
febrile <- function() list(
  annotate("rect", xmin = T_BODY, xmax = T_FEVER, ymin = -Inf, ymax = Inf,
           fill = RED, alpha = .07),
  geom_vline(xintercept = c(T_BODY, T_FEVER), colour = RED, linetype = "22", linewidth = .32))

# ---- locate data + output ---------------------------------------------------
cand <- c(file.path(.d0, "supp_data"), file.path(.d0, "..", "outputs", "supp_data"),
          file.path(.d0, "..", "cauris_etcgem", "strains", "eci_cauris", "outputs", "supp_data"),
          "supp_data")
DD <- cand[which(vapply(cand, dir.exists, logical(1)))[1]]
if (is.na(DD)) stop("supp_data/ not found.")
# project tables/ (isolate-level economics live here, not in supp_data)
tcand <- c(file.path(.d0, "..", "results", "tables"), file.path(.d0, "results", "tables"), file.path("results", "tables"))
TD <- tcand[which(vapply(tcand, dir.exists, logical(1)))[1]]
if (is.na(TD)) stop("tables/ not found (need carbon_tax_isolate.csv, derived_N0_R_results_with_carbon.csv).")
FD <- file.path(.d0, "..", "results", "figures", "manuscript"); dir.create(FD, recursive = TRUE, showWarnings = FALSE)
FD <- normalizePath(FD, mustWork = FALSE)
rd  <- function(f) read.csv(file.path(DD, f))
rdt <- function(f) read.csv(file.path(TD, f))
save_fig <- function(f, name, h_mm, w_mm = 183) {
  ggsave(file.path(FD, paste0(name, ".png")), f, width = w_mm/25.4, height = h_mm/25.4,
         dpi = 600, bg = "white")
  tryCatch({ dv <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
    ggsave(file.path(FD, paste0(name, ".pdf")), f, width = w_mm/25.4, height = h_mm/25.4,
           device = dv, bg = "white") }, error = function(e) NULL)
  message("  ", name, ".png + .pdf")
}
fac <- function(x) factor(SHORT[as.character(x)], levels = SHORT[GRPS])

cur  <- rd("calibration_curves.csv")            # clade,temp_C,apriori,calibrated
meas <- rd("calibration_measured.csv")          # clade,temp_C,observed
pred <- rd("bayes_predictive.csv")              # clade,temp_C,med,lo,hi
samp <- rd("bayes_posterior_samples.csv")       # clade,dTopt,dCp_scale,kcat_scale
summ <- rd("bayes_summary.csv")                 # clade,param,median,lo,hi
cur$clade  <- factor(cur$clade,  levels = GRPS)
meas$clade <- factor(meas$clade, levels = GRPS)
pred$clade <- factor(pred$clade, levels = GRPS)
samp$clade <- factor(samp$clade, levels = GRPS)
summ$clade <- factor(summ$clade, levels = GRPS)
cur$fc <- fac(cur$clade); meas$fc <- fac(meas$clade); pred$fc <- fac(pred$clade)

# =============================================================================
# a  BEFORE (a-priori) vs AFTER (Bayesian) -- faceted by clade
# =============================================================================
ymax <- max(c(pred$hi, meas$observed, cur$apriori)) * 1.14
fa <- ggplot() +
  febrile() +
  geom_line(data = cur, aes(temp_C, apriori), colour = "grey55",
            linetype = "22", linewidth = .6) +
  geom_ribbon(data = pred, aes(temp_C, ymin = lo, ymax = hi, fill = clade), alpha = .22) +
  geom_line(data = pred, aes(temp_C, med, colour = clade), linewidth = .8) +
  geom_line(data = cur, aes(temp_C, calibrated), colour = "white", linetype = "22", linewidth = 1.0) +
  geom_point(data = meas, aes(temp_C, observed), shape = 16, size = .9,
             colour = "grey20", alpha = .85) +
  facet_wrap(~ fc, nrow = 1) +
  scale_colour_manual(values = OI, guide = "none") +
  scale_fill_manual(values = OI, guide = "none") +
  scale_x_continuous(breaks = seq(24, 44, 8)) +
  coord_cartesian(xlim = c(T_MIN, T_MAX), ylim = c(0, ymax)) +
  labs(x = XLAB_T, y = YLAB_MU,
       title = "Before vs after calibration",
       subtitle = paste0("Grey dashed = a-priori (genome only).  White dashed = calibrated model.  Colour + band = Bayesian data fit (95% CrI). ",
                         "Points = measured.")) +
  th + theme(strip.background = element_blank(),
             strip.text = element_text(size = 7, face = "bold.italic", colour = "grey20"),
             panel.spacing.x = unit(7, "pt"))

# =============================================================================
# b  SEQUENCE CAN'T EXPLAIN IT -- one genome curve vs the measured spread
#    (near-identity independently confirmed: diamond RBH whole-proteome median
#     99.4-100% identity across clades I/II/III/IV.)
# =============================================================================
apri <- cur[cur$clade == "I", c("temp_C", "apriori")]   # identical across clades
fb <- ggplot() +
  febrile() +
  geom_line(data = apri, aes(temp_C, apriori), colour = "grey15", linewidth = 1.0) +
  geom_point(data = meas, aes(temp_C, observed, colour = clade), size = 1.2, alpha = .9) +
  annotate("text", x = T_MIN + 0.4, y = 0.30, hjust = 0, vjust = 1, size = 2.0,
           colour = "grey30", lineheight = .95,
           label = "black = the one\ngenome-derived curve") +
  scale_colour_manual(values = OI, guide = "none") +
  scale_x_continuous(breaks = seq(24, 44, 4)) +
  coord_cartesian(xlim = c(T_MIN, T_MAX), ylim = c(0, NA)) +
  labs(x = XLAB_T, y = YLAB_MU, title = "Sequence can't explain the clades",
       subtitle = "One genome curve; clades differ ~2-fold", tag = "b") +
  th

# =============================================================================
# c  CAPACITY POSTERIOR per clade (Bayesian) -- and it is the ONLY calibration
#    knob that separates clades (5/6 pairs; dTopt 1/6, dCp_scale 0/6).
# =============================================================================
kc <- summ[summ$param == "kcat_scale", ]
ordc <- as.character(kc$clade[order(kc$median)])
samp$cl2 <- factor(samp$clade, levels = ordc); kc$cl2 <- factor(kc$clade, levels = ordc)
fc <- ggplot() +
  geom_violin(data = samp, aes(kcat_scale, cl2, fill = clade), colour = NA,
              alpha = .35, orientation = "y", width = .9) +
  geom_pointrange(data = kc, aes(median, cl2, xmin = lo, xmax = hi, colour = clade),
                  linewidth = .6, size = .45) +
  geom_text(data = kc, aes(median, cl2, label = sprintf("%.2f", median)),
            vjust = -1.05, size = 2.1, colour = "grey15") +
  scale_fill_manual(values = OI, guide = "none") +
  scale_colour_manual(values = OI, guide = "none") +
  scale_y_discrete(labels = SHORT) +
  labs(x = "Effective growth-capacity (kcat-scale, posterior)", y = NULL,
       title = "What differs: effective capacity",
       subtitle = "posterior per clade; the only knob that separates clades", tag = "c") +
  th + theme(axis.text.y = element_text(size = 6.8, face = "italic"),
             axis.line.y = element_blank(), axis.ticks.y = element_blank(),
             panel.grid.major.y = element_line(linewidth = .2, colour = "grey94"))

# =============================================================================
# d + e  CAPACITY IS CONSEQUENTIAL OUT-OF-SAMPLE (isolate level, n = 12)
#   kcat_scale is the fitted HEIGHT knob, so capacity-vs-peak-growth is a
#   self-consistency check (r=0.96), NOT evidence -> moved to supplement.
#   The non-circular test: capacity (fit to GROWTH only) predicts the
#   INDEPENDENTLY MEASURED fever tax (d) and respiration (e).
# =============================================================================
capi <- read.csv(file.path(DD, "capacity_isolates.csv"))  # clade,isolate,capacity,peak
tax  <- rdt("carbon_tax_isolate.csv")                     # Isolate,Group,T_C,tax
der  <- rdt("derived_N0_R_results_with_carbon.csv")       # already carries otu_name

db  <- der[der$keep %in% TRUE & der$T >= 32 & der$T <= 36, ]      # growth-optimal band
resp <- aggregate(respiration_fgC_h ~ otu_name, db, mean); names(resp)[2] <- "resp_opt"
fev  <- tax[tax$T_C == T_FEVER, c("Isolate", "tax")]; names(fev) <- c("otu_name", "fever_tax40")

iso <- merge(capi, fev, by.x = "isolate", by.y = "otu_name")
iso <- merge(iso, resp, by.x = "isolate", by.y = "otu_name")
iso$clade <- factor(iso$clade, levels = GRPS)

corci <- function(x, y, B = 5000) {                              # r + bootstrap 95% CI + p
  r <- cor(x, y); p <- suppressWarnings(cor.test(x, y)$p.value); n <- length(x)
  bs <- replicate(B, { i <- sample.int(n, n, TRUE)
                       if (sd(x[i]) > 0 && sd(y[i]) > 0) cor(x[i], y[i]) else NA_real_ })
  ci <- quantile(bs, c(.025, .975), na.rm = TRUE)
  sprintf("r = %+.2f [%+.2f, %+.2f]\np = %s", r, ci[1], ci[2],
          if (p < 1e-3) formatC(p, format = "e", digits = 0) else sprintf("%.3f", p))
}
base_d <- function() list(
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey55",
              linewidth = .5, linetype = "22"),
  geom_point(aes(fill = clade), shape = 21, size = 2.3, colour = "white", stroke = .6),
  scale_colour_manual(values = OI, guide = "none"),
  scale_fill_manual(values = OI, guide = "none"),
  scale_x_continuous(expand = expansion(mult = c(.10, .10))))

fd <- ggplot(iso, aes(capacity, fever_tax40, colour = clade)) + base_d() +
  annotate("text", x = Inf, y = Inf, hjust = 1.04, vjust = 1.35, size = 1.95,
           colour = "grey15", lineheight = .95, label = corci(iso$capacity, iso$fever_tax40)) +
  labs(x = "Effective growth-capacity (kcat-scale)", y = "Fever tax at 40 °C (×)",
       title = "Effective capacity covaries with economics",
       subtitle = "cross-phenotype: higher capacity, cheaper fever (n = 12)", tag = "d") + th

fe <- ggplot(iso, aes(capacity, resp_opt, colour = clade)) + base_d() +
  annotate("text", x = Inf, y = Inf, hjust = 1.04, vjust = 1.35, size = 1.95,
           colour = "grey15", lineheight = .95, label = corci(iso$capacity, iso$resp_opt)) +
  labs(x = "Effective growth-capacity (kcat-scale)",
       y = expression("Respiration @32-36 "*degree*"C (fgC cell"^-1*" h"^-1*")"),
       title = "...and measured respiration",
       subtitle = "cross-phenotype: high capacity, less carbon burned", tag = "e") + th

# =============================================================================
# SUPPLEMENT: the demoted self-consistency check (capacity vs peak growth)
# =============================================================================
fs <- ggplot(iso, aes(capacity, peak, colour = clade)) + base_d() +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.5, size = 2.2,
           colour = "grey15", lineheight = .95, label = corci(iso$capacity, iso$peak)) +
  labs(x = "Effective growth-capacity (kcat-scale)",
       y = expression("Peak growth, "*italic(mu)[max]*" (h"^-1*")"),
       title = "Consistency check (NOT independent evidence)",
       subtitle = "capacity is the fitted height knob; this is expected, not a result") + th

# =============================================================================
# assemble: a full-width facet row on top; b | c | (d over e) below
# =============================================================================
fa <- fa + labs(tag = "a")
# three rows: (a) full-width small-multiples, (b|c) middle, (d|e) bottom.
# side-by-side d|e gives each panel real width, so no clipped subtitles / overlap.
FIG <- fa / (fb | fc) / (fd | fe) + patchwork::plot_layout(heights = c(1, 1, 1)) +
  patchwork::plot_annotation(
    title = "Near-identical sequences yield distinct growth curves captured by effective capacity",
    theme = theme(plot.title = element_text(size = 10.5, face = "bold", colour = "grey5",
                                            margin = margin(b = 3)))) &
  theme(plot.tag = element_text(size = 11, face = "bold"))

message("Writing model figure:")
save_fig(FIG, "FIG_MODEL", 215)
save_fig(fs,  "FIG_MODEL_SUPP_consistency", 70, w_mm = 90)
message("FIG_MODEL (+ supp) -> ", FD)
