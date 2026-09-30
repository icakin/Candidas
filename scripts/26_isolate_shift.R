# =============================================================================
# 26_isolate_shift.R - is thermotolerance a translation of the curve or a reshaping?
#
# Isolate-level posteriors (group effect + isolate random effect, same draws for
# growth and respiration). For every isolate and every draw:
#   T_opt     growth optimum                       = topt(E, Eh, Th)
#   T_carbon  carbon-economy optimum (min of R/G)  = topt(E - E_R, Eh, Th)
#   gap       T_opt - T_carbon
#   T_2x      cost-doubling temperature on the high limb (as in 23_thermal_headroom.R)
# Across the 20 isolates, per draw:
#   slope of T_carbon on T_opt (1 = rigid translation), slope of gap on T_opt (0 = rigid),
#   between-isolate SD of gap, E, E_R, Eh and Th, so the reader can see WHICH parameter
#   carries the between-isolate variation in T_opt.
# Base R for the summaries (no summarise() name reuse, see 23 for why).
#
# Run:  Rscript scripts/26_isolate_shift.R
# OUT   results/tables/isolate_shift_draws.csv     (2000 draws x 20 isolates)
#       results/tables/isolate_shift_isolates.csv  (per-isolate medians + 95% CrI)
#       results/tables/isolate_shift_slopes.csv    (per-draw across-isolate statistics)
# =============================================================================
suppressPackageStartupMessages({ library(posterior) })
root <- if (file.exists("results/tables/fit_metrics.csv")) "." else
        if (file.exists("../results/tables/fit_metrics.csv")) ".." else path.expand("~/Desktop/Projects/Candidas")
source(file.path(root, "scripts/config.R"))
CINV <- 11604.51812; TREF <- 293.15; T_HI <- 70; ND <- 2000
models_dir <- file.path(root, "results/rds")
fit_g <- readRDS(file.path(models_dir, "bayes_growth_ss.rds")); fit_r <- readRDS(file.path(models_dir, "bayes_resp_arr.rds"))
gdat  <- readRDS(file.path(models_dir, "bayes_data_growth.rds"))
dg <- as_draws_df(fit_g); dr <- as_draws_df(fit_r)
nd <- min(nrow(dg), nrow(dr)); set.seed(1); keep <- sort(sample(seq_len(nd), min(ND, nd)))
dg <- dg[keep, ]; dr <- dr[keep, ]; nd <- length(keep)
iso <- unique(data.frame(Isolate = as.character(gdat$Isolate), Group = as.character(gdat$Group)))
