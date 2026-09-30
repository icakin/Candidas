# ============================================================================
# REFIT UNDER AUTOMATIC REPROCESSING  —  run in YOUR R (brms 2.23.0, rstan 2.32.7)
#   Rscript reprocess/40_refit_brms.R
# Step 1 reproduces the LEGACY posterior from the legacy data (sanity check).
# Step 2 fits the SAME model to the automatic supported-growth wells.
# Nothing in the legacy pipeline is overwritten; all output goes to reprocess/.
# ============================================================================
suppressPackageStartupMessages({library(dplyr); library(readr); library(brms); library(posterior)})
# --- locate the project root robustly (Rscript-safe; no sys.frame) ----------
.args <- commandArgs(trailingOnly = FALSE)
.fa   <- grep("^--file=", .args, value = TRUE)
.sdir <- if (length(.fa)) dirname(normalizePath(sub("^--file=", "", .fa[1]))) else getwd()
root  <- normalizePath(file.path(.sdir, ".."), mustWork = FALSE)
if (!dir.exists(file.path(root, "reprocess"))) root <- normalizePath(getwd())
if (!dir.exists(file.path(root, "reprocess")))
  stop("Cannot find the project root. Run from the Candidas folder: Rscript reprocess/40_refit_brms.R")
RP <- file.path(root, "reprocess"); TB <- file.path(root, "results", "tables")
stopifnot(file.exists(file.path(RP, "dataset_growth.csv")),
          file.exists(file.path(root, "results", "rds", "bayes_data_growth.rds")))
message("project root: ", root)
CINV <- 11604.51812; TREF <- 293.15
ITER <- 4000; WARM <- 1000; CH <- 4; SEED <- 20260921

# Formula literals must match the LEGACY formula exactly (brms does not
# substitute R variables into non-linear formulas). Verbatim from bayes_growth_ss.rds:
bf_ss <- brms::bf(
  y | se(se_y, sigma = TRUE) ~ lnB0 - (E * 11604.51812) * ((1/TK) - (1/293.15)) - log(1 + exp((Eh * 11604.51812) * ((1/Th) - (1/TK)))),
  lnB0 ~ 0 + Group + (1 | Isolate),
  E    ~ 0 + Group + (1 | Isolate),
  Eh   ~ 0 + Group,
  Th   ~ 0 + Group + (1 | Isolate),
  nl = TRUE)

fit_one <- function(dat, tag) {
  dat <- dat %>% dplyr::mutate(Group = factor(Group), Isolate = factor(Isolate)) %>%
                 dplyr::filter(is.finite(y), is.finite(se_y), se_y > 0, is.finite(TK))
  Th_mu <- max(dat$TK) + 2
  pri <- c(brms::prior_string(sprintf("normal(%.6f, 2)", mean(dat$y)), nlpar="lnB0"),
           brms::prior_string("normal(0.65, 0.5)", nlpar="E",  lb="0"),
           brms::prior_string("normal(3, 2)",      nlpar="Eh", lb="0"),
           brms::prior_string(sprintf("normal(%.6f, 5)", Th_mu), nlpar="Th", lb="273.15"),
           brms::prior_string("student_t(3, 0, 1)",   class="sd", nlpar="lnB0"),
           brms::prior_string("student_t(3, 0, 0.5)", class="sd", nlpar="E"),
           brms::prior_string("student_t(3, 0, 2)",   class="sd", nlpar="Th"),
           brms::prior_string("exponential(1)", class="sigma"))
  f <- brms::brm(bf_ss, data=dat, prior=pri, iter=ITER, warmup=WARM, chains=CH,
                 seed=SEED, cores=CH, control=list(adapt_delta=0.95, max_treedepth=12),
                 refresh=200)
  saveRDS(f, file.path(RP, sprintf("refit_growth_%s.rds", tag)))
  message(sprintf("[%s] n=%d  saved", tag, nrow(dat))); f
}

# ---- 1) LEGACY reproduction ------------------------------------------------
legacy <- readRDS(file.path(root,"results","rds","bayes_data_growth.rds"))
message(sprintf("legacy rows: %d", nrow(legacy)))
f_leg <- fit_one(legacy %>% dplyr::select(y, se_y, TK, Group, Isolate), "legacyrepro")

# ---- 2) AUTOMATIC supported-growth ----------------------------------------
CELL <- legacy %>% dplyr::distinct(otu_name, cell_carbon_fg)
auto <- readr::read_csv(file.path(RP,"dataset_growth.csv"), show_col_types=FALSE) %>%
  dplyr::left_join(CELL, by="otu_name") %>%
  dplyr::filter(is.finite(r_auto_h), r_auto_h > 0, is.finite(r_se_h), r_se_h > 0) %>%
  dplyr::mutate(
    TK      = T + 273.15,
    growth  = r_auto_h * cell_carbon_fg,     # fg C / cell / h, same construction as legacy
    y       = log(growth),
    se_y    = r_se_h / r_auto_h,             # delta method: sd(log r) = se(r)/r
    Group   = group, Isolate = otu_name) %>%
  dplyr::select(y, se_y, TK, Group, Isolate)
message(sprintf("automatic rows: %d", nrow(auto)))
f_auto <- fit_one(auto, "automatic")

# ---- 3) posterior DIFFERENCES, not just overlapping intervals ---------------
grp <- sort(unique(as.character(legacy$Group)))
dl <- posterior::as_draws_df(f_leg); da <- posterior::as_draws_df(f_auto)
n <- min(nrow(dl), nrow(da))
topt <- function(E,Eh,Th){ s<-E/Eh; ifelse(s>0 & s<1, 1/(1/Th - qlogis(s)/(Eh*CINV)) - 273.15, NA) }
res <- do.call(rbind, lapply(grp, function(g){
  gE  <- function(d,p) as.numeric(d[[paste0("b_",p,"_Group",g)]])[1:n]
  To_l <- topt(gE(dl,"E"), gE(dl,"Eh"), gE(dl,"Th"))
  To_a <- topt(gE(da,"E"), gE(da,"Eh"), gE(da,"Th"))
  d    <- To_a - To_l
  data.frame(Group=g,
    Topt_legacy=median(To_l,na.rm=TRUE), Topt_auto=median(To_a,na.rm=TRUE),
    diff_med=median(d,na.rm=TRUE),
    diff_lo=quantile(d,.025,na.rm=TRUE), diff_hi=quantile(d,.975,na.rm=TRUE),
    P_diff_gt0=mean(d>0,na.rm=TRUE),
    E_legacy=median(gE(dl,"E")), E_auto=median(gE(da,"E")),
    dE_med=median(gE(da,"E")-gE(dl,"E")),
    Th_legacy=median(gE(dl,"Th"))-273.15, Th_auto=median(gE(da,"Th"))-273.15,
    dTh_med=median(gE(da,"Th")-gE(dl,"Th")))
}))
readr::write_csv(res, file.path(RP,"refit_growth_comparison.csv"))
print(res, digits=4)
message("\nWritten: reprocess/refit_growth_{legacyrepro,automatic}.rds and refit_growth_comparison.csv")
message("NOTE: this is the GROWTH model only. The respiration model and Figures 2-3 follow once this is checked.")
