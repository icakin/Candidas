# =============================================================================
# C11 PART A - establish what K is, and derive the balance quantity
# =============================================================================
# Everything here is a CHECK against the committed pipeline output. Nothing is
# assumed: each algebraic step is verified numerically on all usable series.
# =============================================================================

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
       value = TRUE)[1])), "00_common.R"))

d <- c11_load()

hdr <- function(x) message("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74))
res <- list()
note <- function(key, value, detail = "") {
  res[[length(res) + 1]] <<- tibble::tibble(check = key,
                                            value = as.character(value),
                                            detail = detail)
  message(sprintf("  %-52s %s", key, value))
}

hdr("PART A.1  What is K?")

# --- the fit target ----------------------------------------------------------
src <- readLines(file.path(SCRIPTS, "07_oxygen_fits.R"), warn = FALSE)
fts <- grep("^FIT_TO_SPLINE\\s*<-", src, value = TRUE)
note("07_oxygen_fits.R fit-target switch", trimws(fts[1]),
     sprintf("line %d", grep("^FIT_TO_SPLINE\\s*<-", src)[1]))
stopifnot(grepl("FALSE", fts[1]))

# The model itself, from config.R:
#     resp_model(r, K, t, O2_0) = O2_0 + (K/r) * (1 - exp(r*t))
# so   dO2/dt = -K * exp(r*t)   and   dO2/dt |_{t=0} = -K.
# Oxygen is mg/L (config.R's data dictionary) and Time is minutes, so K is a
# VOLUMETRIC rate in mg O2 per litre per minute. No normalisation is applied
# anywhere, so no solubility (O2ref) term enters -- the difference from the
# sister repository's D8, where the fit is to a normalised trace.
note("model", "O2(t) = O2_0 + (K/r)(1 - exp(r t))", "config.R resp_model()")
note("units of K", "mg O2 / L / min  (volumetric)", "-dO2/dt at t = 0")
note("units of r", "1 / min", "")
note("units of K / r", "mg O2 / L", "the drawdown amplitude of the fitted curve")
note("solubility / O2ref term required?", "NO",
     "fit is to raw Oxygen (FIT_TO_SPLINE = FALSE), so K is already volumetric")

hdr("PART A.2  The algebra, verified on the committed table")

# 07 forms respiration by integrating the model over the fit window and dividing
# by the biomass integral:
#     C_tot           = (K/r) (e^{r T_end} - 1)                 [mg O2 / L]
#     biomass_integral= N0 (e^{r T_end} - 1) / r                [cell min / L]
#     R_O2_mg_cell_min= C_tot / biomass_integral = K / N0
# The (e^{r T_end} - 1) factors and the r cancel EXACTLY, so the window length
# does not enter. Check it.
e_RK <- with(d, max(abs(R_O2_mg_cell_min - K / N0_cells_per_L) /
                      abs(R_O2_mg_cell_min), na.rm = TRUE))
note("max rel. error   R_O2_mg_cell_min = K / N0", format(e_RK, digits = 3),
     sprintf("%d series", nrow(d)))

e_N0 <- with(d, max(abs(N0_cells_per_L -
                          N_inoculation_cells_per_L * exp(r * delta_Ninoc_to_N0_min)) /
                      abs(N0_cells_per_L), na.rm = TRUE))
note("max rel. error   N0 = N_inoc exp(r delta)", format(e_N0, digits = 3), "")

# Growth:      G = r * q * 60                    q = cell_carbon_fg  [fg C/cell]
# Respiration: R = c * K / N0                    c = MG_TO_FG * O2_TO_C_MASS * RQ * MIN_TO_H
# so           CUE = G/(G+R) = 1 / (1 + (c/(60 q)) * K/(r N0))
cue_pred <- with(d, 1 / (1 + (C11_C / (60 * cell_carbon_fg)) *
                           (K / (r * N0_cells_per_L))))
e_CUE <- max(abs(cue_pred - d$CUE) / abs(d$CUE), na.rm = TRUE)
note("max rel. error   CUE = 1/(1 + (c/60q) K/(r N0))", format(e_CUE, digits = 3),
     sprintf("c = %.6e", C11_C))

stopifnot(e_RK < 1e-12, e_N0 < 1e-12, e_CUE < 1e-12)

hdr("PART A.3  The balance quantity")

# CUE = 1 / (1 + a * B). CUE is maximised exactly where B is minimised, for any
# POSITIVE a that does not depend on temperature. Substituting N0:
#
#     CUE = 1 / (1 + [ c / (60 q N_inoc) ] * K / (r e^{r delta}) )
#             \_________ a _____________/   \______ B ________/
#
# q (cell_carbon_fg) is per ISOLATE and constant across temperature; N_inoc is
# likewise per isolate and constant across temperature; c is a pure constant.
# So a is a temperature-independent positive scalar within a taxon and DROPS OUT
# of the argmin entirely. What does NOT drop out is delta: it is per series and
# varies with temperature.
#
#   B_bp(T)   = K / (r * exp(r * delta))     the published pipeline's balance
#   B_free(T) = K / r                        with the back-projection removed
#                                            (delta = 0, i.e. N0 = N_inoc)
#
# B_free is the fully scale-free quantity: no N0, no carbon quota, no O2-to-C
# factor, no RQ, no inoculation density.
d <- d %>% mutate(B_bp = K / (r * exp(r * delta_Ninoc_to_N0_min)),
                  B_free = K / r)

# a is temperature-independent within an isolate: check it directly.
a_var <- d %>%
  mutate(a = C11_C / (60 * cell_carbon_fg * N_inoculation_cells_per_L)) %>%
  group_by(OTU) %>%
  summarise(n_T = n_distinct(T), n_a = n_distinct(signif(a, 12)), .groups = "drop")
note("isolates whose 'a' varies with temperature",
     sum(a_var$n_a > 1), sprintf("of %d isolates; a = c/(60 q N_inoc)", nrow(a_var)))
stopifnot(all(a_var$n_a == 1))

# and therefore argmax(CUE) == argmin(B_bp), series by series -> check on the
# taxon x temperature aggregate that the pipeline's own CUE ranks inversely to B_bp
chk <- d %>% group_by(group, T) %>%
  summarise(cue = exp(mean(log(CUE))), b = exp(mean(log(B_bp))), .groups = "drop") %>%
  group_by(group) %>%
  summarise(argmax_cue = T[which.max(cue)], argmin_B = T[which.min(b)], .groups = "drop")
note("taxa where argmax(pipeline CUE) != argmin(B_bp)",
     sum(chk$argmax_cue != chk$argmin_B), paste(C11_TAXA, collapse = ", "))

out <- dplyr::bind_rows(res)
readr::write_csv(out, file.path(C11_TAB, "partA_checks.csv"))
readr::write_csv(chk, file.path(C11_TAB, "partA_argmin_equivalence.csv"))
readr::write_csv(d %>% select(T, OTU, Replicate, group, otu_label, r, K,
                              delta_Ninoc_to_N0_min, cell_carbon_fg,
                              N_inoculation_cells_per_L, N0_cells_per_L,
                              CUE, growth_fgC_h, respiration_fgC_h, B_bp, B_free),
                file.path(C11_TAB, "partA_series.csv"))

hdr("PART A done")
message("  wrote ", file.path(C11_TAB, "partA_checks.csv"))
message("  wrote ", file.path(C11_TAB, "partA_argmin_equivalence.csv"))
message("  wrote ", file.path(C11_TAB, "partA_series.csv"), "  (", nrow(d), " series)")
