# =============================================================================
# 22_clinical_alignment.R - does the measured thermal carbon economy line up with
# independent clinical standing?
#
# TWO LEVELS, TESTED SEPARATELY. They are not pooled and the second is not extra n.
#
#   SPECIES LEVEL. The clinical variable is the WHO fungal priority pathogens list
#   (2022), taken as published and fixed BEFORE any physiology was examined:
#       C. auris             critical        rank 3
#       C. parapsilosis      high            rank 2
#       C. haemulonii        not listed      rank 1
#       C. duobushaemulonii  not listed      rank 1
#   The list is built from case fatality, incidence, sequelae and antifungal
#   resistance. It contains nothing about growth or respiration, which is what
#   makes any alignment with our measurements non-trivial.
#
#   CLADE LEVEL. Within C. auris, clade II is reported as overwhelmingly otic
#   rather than invasive: 57 of 61 clade II isolates from ear infections
#   (Chow et al. 2019, J Clin Microbiol, doi 10.1128/JCM.00007-19), a pattern the
#   authors call uncharacteristic of the other clades and caution may partly
#   reflect sampling of non-sterile sites.
#
# WHAT IS COMPUTED. Everything is a posterior quantity carried draw by draw from
# the same fits used in Figs 1 and 2, so parameter uncertainty propagates instead
# of five medians being treated as fixed points. With four species and ties in the
# ranking, a rank correlation is weak, so the direct ordering probabilities are
# the primary output and the Spearman posterior is reported alongside it.
#
# THIS IS AN OBSERVATION, NOT A TEST. Four species, shared ancestry, no correction
# for the many ways clinical importance and thermal physiology could covary.
# Wording in the manuscript must say so.
#
# Run:  Rscript scripts/22_clinical_alignment.R
# OUT   results/tables/clinical_alignment_species.csv
#       results/tables/clinical_alignment_clades.csv
#       results/tables/clinical_alignment_summary.txt
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tibble); library(posterior) })

root <- if (file.exists("results/tables/fit_metrics.csv")) "." else
        if (file.exists("../results/tables/fit_metrics.csv")) ".." else
        path.expand("~/Desktop/Projects/Candidas")
source(file.path(root, "scripts/config.R"))

CINV <- 11604.51812; TREF <- 293.15; T_BODY <- 37; T_FEVER <- 40; T_MIN <- 22; T_MAX <- 44

fit_g <- readRDS(file.path(models_dir, "bayes_growth_ss.rds"))
fit_r <- readRDS(file.path(models_dir, "bayes_resp_arr.rds"))
gdat  <- readRDS(file.path(models_dir, "bayes_data_growth.rds"))
dg <- as_draws_df(fit_g); dr <- as_draws_df(fit_r)
nd <- min(nrow(dg), nrow(dr)); dg <- dg[seq_len(nd), ]; dr <- dr[seq_len(nd), ]
GRPS <- STRAIN_GROUPS[STRAIN_GROUPS %in% unique(as.character(gdat$Group))]
dcol <- function(d, nm) if (nm %in% names(d)) as.numeric(d[[nm]]) else rep(NA_real_, nrow(d))

# --- identical maths to fig_common.R -----------------------------------------
lnG <- function(TC, E, Eh, Th) { TK <- TC + 273.15; u <- Eh * CINV * (1/Th - 1/TK)
  E * CINV * (1/TREF - 1/TK) - log1p(exp(pmin(u, 700))) }
lnR <- function(TC, ER) ER * CINV * (1/TREF - 1/(TC + 273.15))
topt <- function(Eslope, Eh, Th) { s <- Eslope/Eh; o <- rep(NA_real_, length(s))
  k <- is.finite(s) & s > 0 & s < 1
  o[k] <- 1/(1/Th[k] - stats::qlogis(s[k])/(Eh[k]*CINV)) - 273.15; o }
fever_cost <- function(E, ER, Eh, Th)
  exp((lnR(T_FEVER, ER) - lnG(T_FEVER, E, Eh, Th)) -
      (lnR(T_BODY,  ER) - lnG(T_BODY,  E, Eh, Th)))
growth_kept <- function(E, Eh, Th) exp(lnG(T_FEVER, E, Eh, Th) - lnG(T_BODY, E, Eh, Th))

M <- function(g) {
  E <- dcol(dg, paste0("b_E_Group", g)); Eh <- dcol(dg, paste0("b_Eh_Group", g))
  Th <- dcol(dg, paste0("b_Th_Group", g)); ER <- dcol(dr, paste0("b_E_Group", g))
  list(cost = fever_cost(E, ER, Eh, Th), kept = growth_kept(E, Eh, Th),
       Tcue = topt(E - ER, Eh, Th), Topt = topt(E, Eh, Th))
}
D <- setNames(lapply(GRPS, M), GRPS)

AURIS <- intersect(c("Clade1","Clade2","Clade3","Clade4"), GRPS)
pool  <- function(field) Reduce(`+`, lapply(AURIS, function(g) D[[g]][[field]])) / length(AURIS)
SPP <- list(auris = list(cost = pool("cost"), kept = pool("kept"),
                         Tcue = pool("Tcue"), Topt = pool("Topt")),
            para = D[["para"]], Hae = D[["Hae"]], Duo = D[["Duo"]])
WHO <- c(auris = 3, para = 2, Hae = 1, Duo = 1)   # critical / high / unlisted
nm  <- names(WHO)

# --- species level ------------------------------------------------------------
cm <- sapply(nm, function(s) SPP[[s]]$cost)       # draws x species
km <- sapply(nm, function(s) SPP[[s]]$kept)
ok <- apply(is.finite(cm) & is.finite(km), 1, all)
cm <- cm[ok, , drop = FALSE]; km <- km[ok, , drop = FALSE]
rho_c <- apply(cm, 1, function(v) suppressWarnings(cor(v, WHO, method = "spearman")))
rho_k <- apply(km, 1, function(v) suppressWarnings(cor(v, WHO, method = "spearman")))

P_auris_cheapest <- mean(apply(cm, 1, function(v) which.min(v) == 1))
P_auris_keeps    <- mean(apply(km, 1, function(v) which.max(v) == 1))
P_order          <- mean(cm[,"auris"] < cm[,"para"] & cm[,"para"] < cm[,"Hae"] &
                         cm[,"Hae"]   < cm[,"Duo"])
P_listed_lt_un   <- mean(pmax(cm[,"auris"], cm[,"para"]) < pmin(cm[,"Hae"], cm[,"Duo"]))

spec <- tibble(species = nm, who_rank = as.integer(WHO),
  who_group = c("critical","high","not listed","not listed"),
  fever_cost = apply(cm, 2, median),
  fc_lo = apply(cm, 2, quantile, .025), fc_hi = apply(cm, 2, quantile, .975),
  growth_kept = apply(km, 2, median),
  gk_lo = apply(km, 2, quantile, .025), gk_hi = apply(km, 2, quantile, .975))
write.csv(spec, file.path(tables_dir, "clinical_alignment_species.csv"), row.names = FALSE)

# --- clade level --------------------------------------------------------------
cc <- sapply(AURIS, function(g) D[[g]]$cost); kk <- sapply(AURIS, function(g) D[[g]]$kept)
tt <- sapply(AURIS, function(g) D[[g]]$Topt)
okc <- apply(is.finite(cc) & is.finite(kk) & is.finite(tt), 1, all)
cc <- cc[okc, , drop = FALSE]; kk <- kk[okc, , drop = FALSE]; tt <- tt[okc, , drop = FALSE]
i2 <- match("Clade2", AURIS)
P2_costliest <- mean(apply(cc, 1, which.max) == i2)
P2_least_kept<- mean(apply(kk, 1, which.min) == i2)
P2_coolest   <- mean(apply(tt, 1, which.min) == i2)
clade <- tibble(clade = AURIS,
  fever_cost = apply(cc, 2, median), fc_lo = apply(cc, 2, quantile, .025),
  fc_hi = apply(cc, 2, quantile, .975),
  growth_kept = apply(kk, 2, median), Topt = apply(tt, 2, median))
write.csv(clade, file.path(tables_dir, "clinical_alignment_clades.csv"), row.names = FALSE)

# --- report -------------------------------------------------------------------
out <- c(
"CLINICAL ALIGNMENT  (observation, not a test: 4 species, shared ancestry)",
"",
"SPECIES LEVEL  vs WHO fungal priority pathogens list 2022",
capture.output(print(as.data.frame(spec), row.names = FALSE, digits = 3)),
"",
sprintf("  posterior Spearman(fever cost, WHO rank)   = %.2f [%.2f, %.2f], P(rho<0) = %.3f",
        median(rho_c), quantile(rho_c,.025), quantile(rho_c,.975), mean(rho_c < 0)),
sprintf("  posterior Spearman(growth kept, WHO rank)  = %.2f [%.2f, %.2f], P(rho>0) = %.3f",
        median(rho_k), quantile(rho_k,.025), quantile(rho_k,.975), mean(rho_k > 0)),
sprintf("  P(C. auris has the lowest fever cost of the four)      = %.3f", P_auris_cheapest),
sprintf("  P(C. auris retains the most growth at 40 C)            = %.3f", P_auris_keeps),
sprintf("  P(cost ordering auris < para < Hae < Duo exactly)      = %.3f", P_order),
sprintf("  P(both WHO-listed species cheaper than both unlisted)  = %.3f", P_listed_lt_un),
"",
"CLADE LEVEL  vs the otic-not-invasive phenotype of clade II (Chow et al. 2019)",
capture.output(print(as.data.frame(clade), row.names = FALSE, digits = 3)),
"",
sprintf("  P(clade II has the highest fever cost of the four)     = %.3f", P2_costliest),
sprintf("  P(clade II retains the least growth at 40 C)           = %.3f", P2_least_kept),
sprintf("  P(clade II has the coolest growth optimum)             = %.3f", P2_coolest),
"",
"The WHO list is built from case fatality, incidence, sequelae and resistance.",
"It uses no thermal or metabolic information, so an alignment is not circular.",
"It is still four species with shared ancestry: report as a pattern to be tested",
"with more taxa, never as evidence of cause.")
writeLines(out, file.path(tables_dir, "clinical_alignment_summary.txt"))
cat(paste(out, collapse = "\n"), "\n")
