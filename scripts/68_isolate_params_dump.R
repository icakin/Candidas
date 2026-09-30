# =============================================================================
# 68_isolate_params_dump.R - per-isolate and per-group posterior draws of every
# curve parameter, so the downstream Python analyses (69-73) can run without R.
#   growth:      lnB0, E, Eh, Th   (group effect + isolate random effect where fitted)
#   respiration: alpha, ER
# All 12,000 draws are written per isolate and per group (used by the Python figures) plus
# the full posterior medians. QUOTA (fg C per cell, per group) is written so that
# growth in fg C cell^-1 h^-1 can be converted to r (h^-1) exactly as in fig_common.R.
# Run:  Rscript scripts/68_isolate_params_dump.R
# OUT   results/tables/isolate_params_draws.csv, isolate_params_summary.csv,
#       results/tables/group_quota.csv
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(posterior); library(tibble) })
root <- if (file.exists("results/tables/fit_metrics.csv")) "." else
        if (file.exists("../results/tables/fit_metrics.csv")) ".." else
        path.expand("~/Desktop/Projects/Candidas")
source(file.path(root, "scripts/config.R"))
fit_g <- readRDS(file.path(models_dir, "bayes_growth_ss.rds"))
fit_r <- readRDS(file.path(models_dir, "bayes_resp_arr.rds"))
gdat  <- readRDS(file.path(models_dir, "bayes_data_growth.rds"))
dg <- as_draws_df(fit_g); dr <- as_draws_df(fit_r)
nd <- min(nrow(dg), nrow(dr)); dg <- dg[seq_len(nd), ]; dr <- dr[seq_len(nd), ]
iso_map <- gdat %>% distinct(Isolate, Group) %>%
  mutate(Isolate = as.character(Isolate), Group = as.character(Group))
dcol <- function(d, nm) if (nm %in% names(d)) as.numeric(d[[nm]]) else rep(0, nrow(d))
re   <- function(d, par, iso) dcol(d, sprintf("r_Isolate__%s[%s,Intercept]", par, iso))
keep <- seq_len(nd)   # all draws (was every 4th) so Python figures match the printed numbers
P <- lapply(seq_len(nrow(iso_map)), function(k) {
  iso <- iso_map$Isolate[k]; g <- iso_map$Group[k]
  tibble(Isolate = iso, Group = g, .draw = keep,
    lnB0  = (dcol(dg, paste0("b_lnB0_Group",  g)) + re(dg, "lnB0",  iso))[keep],
    E     = (dcol(dg, paste0("b_E_Group",     g)) + re(dg, "E",     iso))[keep],
    Eh    = (dcol(dg, paste0("b_Eh_Group",    g)))[keep],
    Th    = (dcol(dg, paste0("b_Th_Group",    g)) + re(dg, "Th",    iso))[keep],
    alpha = (dcol(dr, paste0("b_alpha_Group", g)) + re(dr, "alpha", iso))[keep],
    ER    = (dcol(dr, paste0("b_E_Group",     g)) + re(dr, "E",     iso))[keep])
}) %>% bind_rows()
G <- lapply(unique(iso_map$Group), function(g) tibble(Isolate = paste0("GROUP_", g), Group = g, .draw = keep,
    lnB0 = dcol(dg, paste0("b_lnB0_Group", g))[keep], E = dcol(dg, paste0("b_E_Group", g))[keep],
    Eh = dcol(dg, paste0("b_Eh_Group", g))[keep], Th = dcol(dg, paste0("b_Th_Group", g))[keep],
    alpha = dcol(dr, paste0("b_alpha_Group", g))[keep], ER = dcol(dr, paste0("b_E_Group", g))[keep])) %>% bind_rows()
out <- bind_rows(P, G)
write.csv(out, file.path(tables_dir, "isolate_params_draws.csv"), row.names = FALSE)
summ <- out %>% group_by(Isolate, Group) %>% summarise(across(c(lnB0, E, Eh, Th, alpha, ER),
  list(med = ~median(.x), lo = ~quantile(.x, .025), hi = ~quantile(.x, .975))), .groups = "drop")
write.csv(summ, file.path(tables_dir, "isolate_params_summary.csv"), row.names = FALSE)
.der <- readr::read_csv(derived_csv, show_col_types = FALSE)
.q <- .der %>% filter(is.finite(growth_fgC_h), is.finite(growth_C_per_C_h), growth_fgC_h > 0, growth_C_per_C_h > 0) %>%
  mutate(OTU = as.integer(OTU), q = growth_fgC_h / growth_C_per_C_h) %>%
  left_join(readr::read_csv(file.path(tables_dir, "otu_names.csv"), show_col_types = FALSE) %>%
              transmute(OTU = as.integer(OTU), Group = as.character(group)), by = "OTU") %>%
  group_by(Group) %>% summarise(quota_fgC = median(q), .groups = "drop")
write.csv(.q, file.path(tables_dir, "group_quota.csv"), row.names = FALSE)
print(as.data.frame(summ %>% select(Isolate, E_med, Eh_med, Th_med, ER_med)), row.names = FALSE, digits = 4)
message(sprintf("\nwritten: isolate_params_draws.csv (%d rows), isolate_params_summary.csv, group_quota.csv", nrow(out)))
