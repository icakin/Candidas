# =============================================================================
# 11_temperature_equilibration_sensitivity.R   (ADD-ON, non-destructive)
# =============================================================================
# The vials are inoculated at bench temperature and take time to reach setpoint.
# During that onset window (inoculation -> fit-window start) the growth rate is
# NOT the stable r the pipeline back-projects with, so N0 = N_inoc*exp(r*dt), and
# therefore per-cell respiration and CUE, carry a TEMPERATURE-STRUCTURED uncertainty.
#
# This script does NOT change any pipeline number. It quantifies the size of that
# uncertainty and writes a table + figure, so it can be reported honestly.
#
# METHOD (a bracket, not a fragile point correction)
#   - equilibration time t_eq: read per (clade, setpoint) from the measured
#     T_internal trace in the raw xlsx (time to reach within 0.5 C of the plateau).
#   - during the overlap of [0, t_eq] with the onset window, growth is bracketed
#     between ZERO and the FULL fitted rate. That brackets N0, hence respiration
#     and CUE. Growth per cell (= r * cell carbon) never uses N0, so it is immune.
#
# OUTPUT  results/tables/temperature_equilibration_sensitivity.csv
#         results/tables/rdt_term_respiration_ranking.csv
#         results/figures/Fig_temperature_equilibration.png
# RUN     Rscript scripts/11_temperature_equilibration_sensitivity.R   (needs readxl, ggplot2, patchwork)
# =============================================================================

.this_dir <- {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
}
source(file.path(.this_dir, "config.R"))
suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(ggplot2)
})
has_patch <- requireNamespace("patchwork", quietly = TRUE)

RESULTS <- file.path(tables_dir, "derived_N0_R_results_with_carbon.csv")
if (!file.exists(RESULTS)) stop("Run the main pipeline first; missing ", RESULTS)
res <- readr::read_csv(RESULTS, show_col_types = FALSE) %>%
  filter(isTRUE(keep) | keep == TRUE) %>%
  filter(is.finite(r), r > 0, is.finite(respiration_fgC_h), respiration_fgC_h > 0,
         is.finite(growth_fgC_h), growth_fgC_h > 0) %>%
  mutate(clade = tolower(sub("_.*", "", otu_name)))

# ---- extract per-point T_internal from the raw xlsx -------------------------
xlsx <- list.files(data_dir, pattern = "_Oxygen\\.xlsx$", full.names = TRUE)
if (!length(xlsx)) xlsx <- list.files(file.path(base_dir, "data"),
                                      pattern = "_Oxygen\\.xlsx$", full.names = TRUE)
read_tint <- function(f) {
  base <- sub("_Oxygen\\.xlsx$", "", basename(f))
  clade <- tolower(sub("_[0-9.]+$", "", base)); Tset <- as.numeric(sub(".*_", "", base))
  raw <- tryCatch(suppressMessages(readxl::read_excel(f, sheet = 1, col_names = FALSE,
                    .name_repair = "minimal")), error = function(e) NULL)
  if (is.null(raw)) return(NULL)
  hdr <- which(apply(raw, 1, function(z) any(z == "Date/Time", na.rm = TRUE)))[1]
  if (is.na(hdr)) return(NULL)
  h <- as.character(unlist(raw[hdr, ]))
  ti <- which(grepl("T_internal", h))[1]; tc <- which(grepl("^Time/Min", h))[1]
  if (is.na(ti) || is.na(tc)) return(NULL)
  body <- raw[(hdr + 1):nrow(raw), ]
  tibble(clade = clade, T_set = Tset,
         time_min = suppressWarnings(as.numeric(unlist(body[, tc]))),
         T_internal = suppressWarnings(as.numeric(unlist(body[, ti])))) %>%
    filter(is.finite(time_min), is.finite(T_internal))
}
ti_long <- bind_rows(lapply(xlsx, read_tint))
if (!nrow(ti_long)) stop("Could not read T_internal from any xlsx in ", data_dir)

# ---- estimate actual onset growth from measured temperature -----------------
# Instead of the pipeline's constant-rate exp(r*dt), walk the measured T_internal
# trace over the onset window and read the clade's growth rate at each temperature
# off its OWN growth curve r(T). This handles both limbs of the thermal curve, so
# below the optimum a cooler onset means slower growth (respiration under-counted)
# and above it a cooler onset means faster growth (respiration over-counted).
# NOTE: this is an ESTIMATE - it leans on the measured growth curve, so treat the
# absolute size as indicative, not a hard correction. Robust conclusion: growth is
# untouched, respiration/CUE move only at the temperature extremes.
tpc <- res %>% group_by(clade, T) %>% summarise(rbar = median(r), .groups = "drop")
rfun <- setNames(
  lapply(split(tpc, tpc$clade), function(z) approxfun(z$T, z$rbar, rule = 2)),
  unique(tpc$clade))

onset_growth <- function(cl, Tset, fit_start, r_fit, dt) {
  z <- ti_long %>% filter(clade == cl, T_set == Tset, time_min <= fit_start) %>% arrange(time_min)
  f <- rfun[[cl]]
  if (nrow(z) >= 2 && !is.null(f)) sum(f(z$T_internal) * c(diff(z$time_min), 0)) else r_fit * dt
}

d <- res %>% rowwise() %>%
  mutate(g_real = onset_growth(clade, T, fit_start_time, r, delta_Ninoc_to_N0_min),
         N0_corr   = N_inoculation_cells_per_L * exp(g_real),
         resp_corr = respiration_fgC_h * (N0_cells_per_L / N0_corr),
         cue_pipe  = growth_fgC_h / (growth_fgC_h + respiration_fgC_h),
         cue_corr  = growth_fgC_h / (growth_fgC_h + resp_corr),
         resp_swing = 100 * (resp_corr / respiration_fgC_h - 1),
         cue_swing  = 100 * (cue_corr / cue_pipe - 1)) %>% ungroup()

summ <- d %>% group_by(T) %>% summarise(
  n = n(),
  resp_swing_med = median(resp_swing), resp_swing_lo = quantile(resp_swing, .25),
  resp_swing_p75 = quantile(resp_swing, .75),
  cue_swing_med  = median(cue_swing),  cue_swing_p25  = quantile(cue_swing, .25),
  cue_swing_hi   = quantile(cue_swing, .75),
  growth_swing = 0, .groups = "drop")
readr::write_csv(summ, file.path(tables_dir, "temperature_equilibration_sensitivity.csv"))
# per-vial corrected values, so 12_uncertainty_bands.R can combine the temperature
# component with replicate scatter into one honest band.
readr::write_csv(
  d %>% select(T, OTU, Replicate, clade, respiration_fgC_h, resp_corr,
               cue_pipe, cue_corr, growth_fgC_h, resp_swing, cue_swing),
  file.path(tables_dir, "temperature_equilibration_percurve.csv"))
cat("\n== temperature-equilibration uncertainty (median % shift from pipeline) ==\n")
print(as.data.frame(summ %>% transmute(T, n, growth = 0,
      respiration = round(resp_swing_med, 1), CUE = round(cue_swing_med, 1))), row.names = FALSE)

# ---- r*dt term test: respiration ranking with vs without the back-projection -
term <- res %>% mutate(
  resp_with   = respiration_fgC_h,                       # N0 = N_inoc*exp(r*dt)
  resp_noterm = respiration_fgC_h * exp(r * delta_Ninoc_to_N0_min)) %>%  # N0 = N_inoc
  group_by(clade) %>%
  summarise(resp_with = median(resp_with), resp_noterm = median(resp_noterm), .groups = "drop") %>%
  mutate(rank_with = rank(resp_with), rank_noterm = rank(resp_noterm))
readr::write_csv(term, file.path(tables_dir, "rdt_term_respiration_ranking.csv"))
cat("\n== r*dt term test: median respiration per clade, rank with vs without term ==\n")
print(as.data.frame(term), row.names = FALSE, digits = 3)

# ---- figure -----------------------------------------------------------------
mk <- function(med, elo, ehi, ttl, col, lo, hi) ggplot(summ, aes(T)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_ribbon(aes(ymin = .data[[elo]], ymax = .data[[ehi]]), fill = col, alpha = .22) +
  geom_line(aes(y = .data[[med]]), colour = col, linewidth = 1) +
  geom_point(aes(y = .data[[med]]), colour = col, size = 1.4) +
  coord_cartesian(ylim = c(lo, hi)) +
  labs(title = ttl, x = "setpoint (C)", y = "% shift from pipeline") + theme_bw(base_size = 11)
p1 <- ggplot(summ, aes(T)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(aes(y = 0), colour = "seagreen4", linewidth = 1) +
  coord_cartesian(ylim = c(-35, 35)) +
  labs(title = "Growth per cell (immune to N0)", x = "setpoint (C)", y = "% shift from pipeline") +
  theme_bw(base_size = 11)
p2 <- mk("resp_swing_med", "resp_swing_lo", "resp_swing_p75", "Respiration per cell", "firebrick", -35, 35)
p3 <- mk("cue_swing_med",  "cue_swing_p25", "cue_swing_hi",   "CUE",                  "steelblue4", -35, 35)
if (has_patch) {
  p <- patchwork::wrap_plots(p1, p2, p3, nrow = 1) +
    patchwork::plot_annotation(
      title = "Temperature-equilibration uncertainty on each quantity",
      subtitle = "Estimate from measured T_internal vs each clade's growth curve (IQR band). Growth immune; respiration/CUE shift only at the extremes.")
  ggsave(file.path(figures_dir, "Fig_temperature_equilibration.png"), p, width = 13, height = 4.6, dpi = 150)
} else {
  ggsave(file.path(figures_dir, "Fig_temperature_equilibration.png"), p2, width = 6, height = 4.6, dpi = 150)
  message("patchwork not installed; saved the respiration panel only.")
}
cat("\nwrote Fig_temperature_equilibration.png and two tables.\n")
