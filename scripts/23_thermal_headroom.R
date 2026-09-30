# =============================================================================
# 23_thermal_headroom.R - thermal headroom inside the host
#
# For EVERY posterior draw, independently:
#   T_min      = temperature minimising R/G for that draw  (= Sharpe-Schoolfield
#                optimum with E replaced by E - E_R; the carbon-economy optimum)
#   T_Mxcost   = the crossing ON THE HIGH-TEMPERATURE LIMB (T > T_min) at which R/G
#                reaches M times that draw's OWN minimum, for M = 1.5, 2 and 3.
#                Scale-free: R/G = (c/q) * K/r and the constant cancels in the ratio,
#                so no starting cell number, carbon quota, O2-to-C factor or
#                respiratory quotient enters. M = 2 is the main-text threshold; the
#                other two are the supplementary sensitivity test.
#                Draws with no crossing below T_HI are returned NA and are EXCLUDED
#                from the summaries rather than imputed or censored at the boundary;
#                n_ok reports how many remain (this bites only C. duobushaemulonii,
#                whose curve is steep enough that some draws never reach 3x).
#   headroom   = T_2xcost - 37 and T_2xcost - 40
#   T50        = the high-temperature crossing where growth falls to half its own
#                peak. Reported but flagged: for the C. auris clades this lands at
#                44-46 C, at or beyond the 44 C assay ceiling, so it is extrapolated.
#
# Maths identical to fig_common.R (CINV, TREF, lnG, ln_cost, topt).
#
# Run:  Rscript scripts/23_thermal_headroom.R
# OUT   results/tables/thermal_headroom_draws.csv   (per-draw, for the figure)
#       results/tables/thermal_headroom_summary.csv (medians + 95% CrI)
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(posterior); library(tibble) })

root <- if (file.exists("results/tables/fit_metrics.csv")) "." else
        if (file.exists("../results/tables/fit_metrics.csv")) ".." else
        path.expand("~/Desktop/Projects/Candidas")
source(file.path(root, "scripts/config.R"))

CINV <- 11604.51812; TREF <- 293.15
T_BODY <- 37; T_FEVER <- 40; T_LO <- 15; T_HI <- 70

fit_g <- readRDS(file.path(models_dir, "bayes_growth_ss.rds"))
fit_r <- readRDS(file.path(models_dir, "bayes_resp_arr.rds"))
gdat  <- readRDS(file.path(models_dir, "bayes_data_growth.rds"))
dg <- as_draws_df(fit_g); dr <- as_draws_df(fit_r)
nd <- min(nrow(dg), nrow(dr)); dg <- dg[seq_len(nd), ]; dr <- dr[seq_len(nd), ]
GRPS <- STRAIN_GROUPS[STRAIN_GROUPS %in% unique(as.character(gdat$Group))]
dcol <- function(d, nm) if (nm %in% names(d)) as.numeric(d[[nm]]) else rep(NA_real_, nrow(d))
message(sprintf("draws: %d | groups: %s", nd, paste(GRPS, collapse = ", ")))

lnG <- function(TC, E, Eh, Th) {
  TK <- TC + 273.15; u <- Eh * CINV * (1/Th - 1/TK)
  E * CINV * (1/TREF - 1/TK) - log1p(exp(pmin(u, 700)))
}
ln_cost <- function(TC, E, ER, Eh, Th) {          # log(R/G) up to a constant
  TK <- TC + 273.15; u <- Eh * CINV * (1/Th - 1/TK)
  (ER - E) * CINV * (1/TREF - 1/TK) + ifelse(u > 30, u, log1p(exp(pmin(u, 30))))
}
topt <- function(Es, Eh, Th) {
  s <- Es/Eh; if (!is.finite(s) || s <= 0 || s >= 1) return(NA_real_)
  1/(1/Th - stats::qlogis(s)/(Eh * CINV)) - 273.15
}
cross <- function(f, lo, hi) {                    # upward crossing of zero, or NA
  if (!is.finite(f(lo)) || !is.finite(f(hi))) return(NA_real_)
  if (f(lo) > 0 || f(hi) < 0) return(NA_real_)
  tryCatch(stats::uniroot(f, c(lo, hi), tol = 1e-4)$root, error = function(e) NA_real_)
}

out <- lapply(GRPS, function(g) {
  E  <- dcol(dg, paste0("b_E_Group",  g)); Eh <- dcol(dg, paste0("b_Eh_Group", g))
  Th <- dcol(dg, paste0("b_Th_Group", g)); ER <- dcol(dr, paste0("b_E_Group",  g))
  MULT <- c(1.5, 2, 3)
  res <- vapply(seq_len(nd), function(i) {
    Tmin <- topt(E[i] - ER[i], Eh[i], Th[i])
    Topt <- topt(E[i],         Eh[i], Th[i])
    if (!is.finite(Tmin) || !is.finite(Topt)) return(rep(NA_real_, 3 + length(MULT)))
    c0 <- ln_cost(Tmin, E[i], ER[i], Eh[i], Th[i])
    tM <- vapply(MULT, function(m)
      cross(function(TC) ln_cost(TC, E[i], ER[i], Eh[i], Th[i]) - c0 - log(m), Tmin, T_HI),
      numeric(1))
    g0 <- lnG(Topt, E[i], Eh[i], Th[i])
    t50 <- cross(function(TC) (g0 + log(0.5)) - lnG(TC, E[i], Eh[i], Th[i]), Topt, T_HI)
    c(Tmin, Topt, tM, t50)
  }, numeric(3 + length(MULT)))
  tibble(Group = g, .draw = seq_len(nd),
         T_min_cost   = res[1, ], T_opt_growth = res[2, ],
         T_1.5xcost   = res[3, ], T_2xcost = res[4, ], T_3xcost = res[5, ],
         T50          = res[6, ],
         headroom_body  = res[4, ] - T_BODY,
         headroom_fever = res[4, ] - T_FEVER)
}) %>% bind_rows()

write.csv(out, file.path(tables_dir, "thermal_headroom_draws.csv"), row.names = FALSE)

q <- function(x, p) stats::quantile(x[is.finite(x)], p, names = FALSE)
# NB: every summary term reads the raw draw column, never a term created earlier in the
# same summarise() call - dplyr makes later expressions see the new value, which silently
# collapses every interval onto its median.
summ <- out %>% group_by(Group) %>% summarise(
  n_ok             = sum(is.finite(T_2xcost)),
  T_2xcost_med     = median(T_2xcost, na.rm = TRUE),
  t2_lo            = q(T_2xcost, .025),
  t2_hi            = q(T_2xcost, .975),
  headroom_body_med = median(headroom_body, na.rm = TRUE),
  hb_lo            = q(headroom_body, .025),
  hb_hi            = q(headroom_body, .975),
  P_body_positive  = mean(headroom_body  > 0, na.rm = TRUE),
  headroom_fever_med = median(headroom_fever, na.rm = TRUE),
  hf_lo            = q(headroom_fever, .025),
  hf_hi            = q(headroom_fever, .975),
  P_fever_positive = mean(headroom_fever > 0, na.rm = TRUE),
  T50_med          = median(T50, na.rm = TRUE),
  t50_lo           = q(T50, .025),
  t50_hi           = q(T50, .975),
  .groups = "drop") %>% arrange(desc(T_2xcost_med))
write.csv(summ, file.path(tables_dir, "thermal_headroom_summary.csv"), row.names = FALSE)
print(as.data.frame(summ), row.names = FALSE, digits = 3)
message("\nwritten: thermal_headroom_draws.csv, thermal_headroom_summary.csv")
