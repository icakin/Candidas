# =============================================================================
# C11 PART E - figures
# =============================================================================
# Writes into reports/C11_scale_free/figures/. config.R's figure whitelist is
# declared through fig_keep_add() in 00_common.R, which is the entry point added
# in C10 precisely so a script can write a figure the base list does not know.
# =============================================================================

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
       value = TRUE)[1])), "00_common.R"))

E    <- readr::read_csv(file.path(C11_TAB, "partB_activation_energies.csv"), show_col_types = FALSE)
opt  <- readr::read_csv(file.path(C11_TAB, "partB_modelfree_optimum.csv"),   show_col_types = FALSE)
cur  <- readr::read_csv(file.path(C11_TAB, "partB_balance_curves.csv"),      show_col_types = FALSE)
gap  <- readr::read_csv(file.path(C11_TAB, "partC_gap_decomposition.csv"),   show_col_types = FALSE)

lab  <- setNames(vapply(C11_TAXA, c11_label, character(1)), C11_TAXA)
ordf <- function(g) factor(lab[as.character(g)], levels = unname(lab))

th <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93", colour = NA),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(colour = "grey35", size = 8.6))
RED <- "#B22222"

# ---- Fig 1: activation energies --------------------------------------------
f1d <- dplyr::bind_rows(
  E %>% transmute(taxon = ordf(group), quantity = "E_r  (growth, Sharpe-Schoolfield)",
                  est = E_r, lo = E_r_lo, hi = E_r_hi),
  E %>% transmute(taxon = ordf(group), quantity = "E_K  (respiration, Arrhenius, scale-free)",
                  est = E_K, lo = E_K_lo, hi = E_K_hi))
p1 <- ggplot(f1d, aes(est, taxon, colour = quantity)) +
  geom_vline(xintercept = 0.65, linetype = "22", colour = "grey55") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .18,
                 position = position_dodge(width = .5), linewidth = .5) +
  geom_point(size = 2.2, position = position_dodge(width = .5)) +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  labs(x = "activation energy  E  (eV)", y = NULL, colour = NULL,
       title = "Activation energies of growth and of the volumetric O2 uptake rate",
       subtitle = paste("E_K is fitted to K itself - no N0, no carbon quota, no O2-to-C factor, no RQ.",
                        "\nDashed line: the 0.65 eV metabolic-theory value. Bars: bootstrap over isolates.")) +
  th + theme(legend.position = "bottom")
ggsave(file.path(C11_FIG, "c11_fig1_activation_energies.png"), p1,
       width = 8.2, height = 3.8, dpi = 200, bg = "white")

# ---- Fig 2: the balance quantity against temperature ------------------------
f2 <- cur %>% mutate(taxon = ordf(group))
mins <- opt %>% transmute(taxon = ordf(group), T_opt, B = NA_real_)
mins$B <- vapply(seq_len(nrow(mins)), function(i) {
  s <- dplyr::filter(f2, taxon == mins$taxon[i]); exp(min(s$logB))
}, numeric(1))
p2 <- ggplot(f2, aes(T, B)) +
  annotate("rect", xmin = 37, xmax = 40, ymin = -Inf, ymax = Inf, fill = RED, alpha = .08) +
  geom_vline(xintercept = 37, colour = RED, linetype = "22", linewidth = .4) +
  geom_line(colour = "grey30", linewidth = .5) +
  geom_point(size = 1.2, colour = "grey30") +
  geom_point(data = mins, aes(T_opt, B), colour = "#0072B2", size = 2.8) +
  geom_vline(data = mins, aes(xintercept = T_opt), colour = "#0072B2",
             linetype = "solid", linewidth = .35, alpha = .6) +
  facet_wrap(~taxon, nrow = 1, scales = "free_y") +
  scale_y_log10() +
  labs(x = expression("Temperature ("*degree*"C)"),
       y = expression("balance  K/r  (mg "*O[2]*" "*L^-1*")"),
       title = "The balance quantity K/r, whose minimum is the CUE optimum",
       subtitle = paste("CUE is maximised exactly where K/r is minimised; every conversion constant sits in a factor that cancels.",
                        "\nBlue: the model-free minimum. Red band: 37-40 C, body to fever.")) + th
ggsave(file.path(C11_FIG, "c11_fig2_balance_vs_T.png"), p2,
       width = 10.5, height = 3.4, dpi = 200, bg = "white")

# ---- Fig 3: model-free vs fitted -------------------------------------------
f3 <- gap %>%
  transmute(taxon = ordf(group),
            `model-free, scale-free (K/r)` = T_free,
            `model-free, published N0`     = T_bp,
            `model-free, equilibration-corrected` = T_eq,
            `published (fitted Bayesian CUE)` = T_pub) %>%
  tidyr::pivot_longer(-taxon, names_to = "estimate", values_to = "T_opt") %>%
  mutate(estimate = factor(estimate, levels = c(
    "model-free, scale-free (K/r)", "model-free, published N0",
    "model-free, equilibration-corrected", "published (fitted Bayesian CUE)")))
lims <- gap %>% transmute(taxon = ordf(group), lo = T_free_lo, hi = T_free_hi,
                          estimate = factor("model-free, scale-free (K/r)",
                                            levels = levels(f3$estimate)), T_opt = T_free)
p3 <- ggplot(f3, aes(T_opt, taxon, colour = estimate, shape = estimate)) +
  annotate("rect", xmin = 37, xmax = Inf, ymin = -Inf, ymax = Inf, fill = RED, alpha = .07) +
  geom_vline(xintercept = 37, colour = RED, linetype = "22") +
  geom_errorbarh(data = lims, aes(xmin = lo, xmax = hi), height = .14, linewidth = .45) +
  geom_point(size = 2.4, position = position_dodge(width = .45)) +
  scale_colour_manual(values = c("#0072B2", "#009E73", "#E69F00", "grey25")) +
  scale_shape_manual(values = c(16, 17, 15, 4)) +
  labs(x = expression("T"[opt]*"(CUE)  ("*degree*"C)"), y = NULL,
       colour = NULL, shape = NULL,
       title = "Model-free versus fitted CUE optimum",
       subtitle = paste("Bars: bootstrap 95% interval on the scale-free estimate.",
                        "Shaded: at or above body temperature.")) +
  th + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
ggsave(file.path(C11_FIG, "c11_fig3_modelfree_vs_fitted.png"), p3,
       width = 8.2, height = 4.2, dpi = 200, bg = "white")

# ---- Fig 4: gap decomposition ----------------------------------------------
f4 <- gap %>%
  transmute(taxon = ordf(group),
            `N0 back-projection` = gap_N0,
            `functional form`    = gap_form,
            `solubility / units` = gap_solubility) %>%
  tidyr::pivot_longer(-taxon, names_to = "source", values_to = "degrees") %>%
  mutate(source = factor(source, levels = c("N0 back-projection", "functional form",
                                            "solubility / units")))
tot <- gap %>% transmute(taxon = ordf(group), degrees = gap_total)
p4 <- ggplot(f4, aes(degrees, taxon, fill = source)) +
  geom_vline(xintercept = 0, colour = "grey45") +
  geom_col(width = .6) +
  geom_point(data = tot, aes(degrees, taxon), inherit.aes = FALSE,
             shape = 18, size = 3, colour = "black") +
  scale_fill_manual(values = c("#0072B2", "#E69F00", "grey75")) +
  labs(x = expression("contribution to  T"[opt]^{"model-free"}*" - T"[opt]^{"published"}*"  ("*degree*"C)"),
       y = NULL, fill = NULL,
       title = "What separates the model-free optimum from the published one",
       subtitle = paste("Negative = the model-free optimum is COOLER, so the margin to 37 C widens.",
                        "\nDiamond: the total. The solubility term is exactly zero (K is volumetric).")) +
  th + theme(legend.position = "bottom")
ggsave(file.path(C11_FIG, "c11_fig4_gap_decomposition.png"), p4,
       width = 8.2, height = 3.8, dpi = 200, bg = "white")

message("\nwrote 4 figures to ", C11_FIG)
for (f in list.files(C11_FIG, full.names = TRUE)) message("  ", basename(f))
