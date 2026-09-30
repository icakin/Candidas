# =============================================================================
# 65_isolate_crossings.R - per-ISOLATE cost-multiplier crossings
#
# Same maths as 23_thermal_headroom.R, but each isolate's own curve from the SAME
# hierarchical fit: group fixed effect + that isolate's random effects on E, Th (growth)
# and E (respiration). Eh has no isolate term (group level only, as fitted).
# For every posterior draw: T_min (carbon-economy optimum), and the high-limb crossing
# where R/G reaches M x that draw's own minimum, M = 1.5, 2, 3. Scale-free as before.
# Used to compare the carbon-economy collapse with the observed growth limit per isolate
# (results/tables/lik_transition.csv) in scripts/66_safety_margin_lik.py.
#
# Run:  Rscript scripts/65_isolate_crossings.R          (needs the brms fits; ~minutes)
# OUT   results/tables/isolate_crossings_draws.csv, isolate_crossings_summary.csv
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(posterior); library(tibble) })
root <- if (file.exists("results/tables/fit_metrics.csv")) "." else
        if (file.exists("../results/tables/fit_metrics.csv")) ".." else
        path.expand("~/Desktop/Projects/Candidas")
source(file.path(root, "scripts/config.R"))
CINV <- 11604.51812; TREF <- 293.15; T_HI <- 70
fit_g <- readRDS(file.path(models_dir, "bayes_growth_ss.rds"))
fit_r <- readRDS(file.path(models_dir, "bayes_resp_arr.rds"))
gdat  <- readRDS(file.path(models_dir, "bayes_data_growth.rds"))
dg <- as_draws_df(fit_g); dr <- as_draws_df(fit_r)
nd <- min(nrow(dg), nrow(dr)); dg <- dg[seq_len(nd), ]; dr <- dr[seq_len(nd), ]
iso_map <- gdat %>% distinct(Isolate, Group) %>%
  mutate(Isolate = as.character(Isolate), Group = as.character(Group))
dcol <- function(d, nm) if (nm %in% names(d)) as.numeric(d[[nm]]) else rep(0, nrow(d))
re   <- function(d, par, iso) dcol(d, sprintf("r_Isolate__%s[%s,Intercept]", par, iso))
message(sprintf("draws: %d | isolates: %d", nd, nrow(iso_map)))
ln_cost <- function(TC, E, ER, Eh, Th) {
  TK <- TC + 273.15; u <- Eh * CINV * (1/Th - 1/TK)
  (ER - E) * CINV * (1/TREF - 1/TK) + ifelse(u > 30, u, log1p(exp(pmin(u, 30))))
}
topt <- function(Es, Eh, Th) {
  s <- Es/Eh; if (!is.finite(s) || s <= 0 || s >= 1) return(NA_real_)
  1/(1/Th - stats::qlogis(s)/(Eh * CINV)) - 273.15
}
cross <- function(f, lo, hi) {
  if (!is.finite(f(lo)) || !is.finite(f(hi))) return(NA_real_)
  if (f(lo) > 0 || f(hi) < 0) return(NA_real_)
  tryCatch(stats::uniroot(f, c(lo, hi), tol = 1e-4)$root, error = function(e) NA_real_)
}
MULT <- c(1.5, 2, 3)
out <- lapply(seq_len(nrow(iso_map)), function(k) {
  iso <- iso_map$Isolate[k]; g <- iso_map$Group[k]
  E  <- dcol(dg, paste0("b_E_Group",  g)) + re(dg, "E",  iso)
  Th <- dcol(dg, paste0("b_Th_Group", g)) + re(dg, "Th", iso)
  Eh <- dcol(dg, paste0("b_Eh_Group", g))
  ER <- dcol(dr, paste0("b_E_Group",  g)) + re(dr, "E",  iso)
  res <- vapply(seq_len(nd), function(i) {
    Tmin <- topt(E[i] - ER[i], Eh[i], Th[i]); Topt <- topt(E[i], Eh[i], Th[i])
    if (!is.finite(Tmin) || !is.finite(Topt)) return(rep(NA_real_, 2 + length(MULT)))
    c0 <- ln_cost(Tmin, E[i], ER[i], Eh[i], Th[i])
    tM <- vapply(MULT, function(m)
      cross(function(TC) ln_cost(TC, E[i], ER[i], Eh[i], Th[i]) - c0 - log(m), Tmin, T_HI), numeric(1))
    c(Tmin, Topt, tM)
  }, numeric(2 + length(MULT)))
  tibble(Isolate = iso, Group = g, .draw = seq_len(nd),
         T_min_cost = res[1, ], T_opt_growth = res[2, ],
         T_1.5xcost = res[3, ], T_2xcost = res[4, ], T_3xcost = res[5, ])
}) %>% bind_rows()
write.csv(out, file.path(tables_dir, "isolate_crossings_draws.csv"), row.names = FALSE)
q <- function(x, p) stats::quantile(x[is.finite(x)], p, names = FALSE)
summ <- out %>% group_by(Isolate, Group) %>% summarise(
  n_ok = sum(is.finite(T_3xcost)),
  T_opt_med = median(T_opt_growth, na.rm = TRUE),
  T_1.5x_med = median(T_1.5xcost, na.rm = TRUE), T_1.5x_lo = q(T_1.5xcost, .025), T_1.5x_hi = q(T_1.5xcost, .975),
  T_2x_med = median(T_2xcost, na.rm = TRUE), T_2x_lo = q(T_2xcost, .025), T_2x_hi = q(T_2xcost, .975),
  T_3x_med = median(T_3xcost, na.rm = TRUE), T_3x_lo = q(T_3xcost, .025), T_3x_hi = q(T_3xcost, .975),
  .groups = "drop")
write.csv(summ, file.path(tables_dir, "isolate_crossings_summary.csv"), row.names = FALSE)
print(as.data.frame(summ), row.names = FALSE, digits = 3)
message("\nwritten: isolate_crossings_draws.csv, isolate_crossings_summary.csv")
