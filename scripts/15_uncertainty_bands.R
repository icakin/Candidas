# =============================================================================
# 12_uncertainty_bands.R   (ADD-ON, non-destructive)
# =============================================================================
# Respiration and CUE versus temperature, per clade, with an HONEST uncertainty
# band that combines TWO sources:
#   1. replicate scatter        - the vial-to-vial spread (biological + measurement)
#   2. temperature equilibration - the onset-temperature effect from script 11,
#                                  small at mid T, large at the extremes
# combined in quadrature:  sd_total = sqrt(sd_replicate^2 + sd_temperature^2).
#
# The point estimate is the UNCHANGED pipeline value. Only the band is new.
# Growth is not shown here: it never uses N0, so it carries no temperature term.
# The carbon-side range (RQ, cell-carbon) is a SEPARATE systematic multiplier and
# is reported as a stated assumption, not folded into this per-point band.
#
# INPUT   results/tables/derived_N0_R_results_with_carbon.csv
#         results/tables/temperature_equilibration_percurve.csv   (from script 11)
# OUTPUT  results/figures/fig_respiration_uncertainty.png
#         results/figures/fig_cue_uncertainty.png
#         results/tables/respiration_cue_uncertainty.csv
# RUN     Rscript scripts/12_uncertainty_bands.R   (run 11 first)
# =============================================================================

.this_dir <- {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
}
source(file.path(.this_dir, "config.R"))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

pcurve <- file.path(tables_dir, "temperature_equilibration_percurve.csv")
if (!file.exists(pcurve)) stop("Run 11_temperature_equilibration_sensitivity.R first; missing ", pcurve)
pv <- readr::read_csv(pcurve, show_col_types = FALSE)

# ---- combine the two uncertainty sources per clade x temperature ------------
# work on log respiration and on CUE directly.
band <- pv %>%
  group_by(clade, T) %>%
  summarise(
    n            = n(),
    resp         = median(respiration_fgC_h),
    # replicate scatter of log respiration -> SE of the mean
    sd_rep_logR  = if (n() > 1) sd(log(respiration_fgC_h)) / sqrt(n()) else NA_real_,
    # temperature-equilibration component, as a fractional (log) shift
    sd_tmp_logR  = abs(median(resp_swing)) / 100,
    cue          = median(cue_pipe),
    sd_rep_cue   = if (n() > 1) sd(cue_pipe) / sqrt(n()) else NA_real_,
    sd_tmp_cue   = abs(median(cue_swing)) / 100 * median(cue_pipe),
    .groups = "drop") %>%
  mutate(
    sd_rep_logR = ifelse(is.na(sd_rep_logR), 0, sd_rep_logR),
    sd_rep_cue  = ifelse(is.na(sd_rep_cue),  0, sd_rep_cue),
    # quadrature combination
    sd_tot_logR = sqrt(sd_rep_logR^2 + sd_tmp_logR^2),
    resp_lo = resp * exp(-1.96 * sd_tot_logR),
    resp_hi = resp * exp( 1.96 * sd_tot_logR),
    sd_tot_cue = sqrt(sd_rep_cue^2 + sd_tmp_cue^2),
    cue_lo = pmax(0, cue - 1.96 * sd_tot_cue),
    cue_hi = pmin(1, cue + 1.96 * sd_tot_cue))
readr::write_csv(band, file.path(tables_dir, "respiration_cue_uncertainty.csv"))

# ---- figures: point estimate + combined band, faceted by clade --------------
resp_p <- ggplot(band, aes(T, resp)) +
  geom_ribbon(aes(ymin = resp_lo, ymax = resp_hi), fill = "firebrick", alpha = .20) +
  geom_line(colour = "firebrick", linewidth = .7) +
  geom_point(colour = "firebrick", size = 1.2) +
  facet_wrap(~ clade, scales = "free_y") +
  labs(title = "Respiration per cell vs temperature, with combined uncertainty",
       subtitle = "Band = replicate scatter + temperature-equilibration (quadrature). Point = pipeline value.",
       x = "temperature (C)", y = "respiration (fg C cell-1 h-1)") +
  theme_bw(base_size = 11)
ggsave(file.path(figures_dir, "fig_respiration_uncertainty.png"), resp_p,
       width = 10, height = 6.5, dpi = 150)

cue_p <- ggplot(band, aes(T, cue)) +
  geom_ribbon(aes(ymin = cue_lo, ymax = cue_hi), fill = "steelblue4", alpha = .20) +
  geom_line(colour = "steelblue4", linewidth = .7) +
  geom_point(colour = "steelblue4", size = 1.2) +
  facet_wrap(~ clade, scales = "free_y") +
  labs(title = "CUE vs temperature, with combined uncertainty",
       subtitle = "Band = replicate scatter + temperature-equilibration. Carbon assumptions (RQ, cell C) are a separate stated range.",
       x = "temperature (C)", y = "CUE") +
  theme_bw(base_size = 11)
ggsave(file.path(figures_dir, "fig_cue_uncertainty.png"), cue_p,
       width = 10, height = 6.5, dpi = 150)

cat("\nwrote fig_respiration_uncertainty.png, fig_cue_uncertainty.png, and the table.\n")
cat("median band width (95%) as +/- % of value:\n")
print(as.data.frame(band %>% group_by(T) %>%
  summarise(resp_pm = round(100 * (exp(1.96 * median(sd_tot_logR)) - 1), 1),
            .groups = "drop")), row.names = FALSE)
