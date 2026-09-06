# =============================================================================
# 16_fig1.R  --  FIGURE 1: growth and respiration decouple with temperature
# =============================================================================
#   Rscript scripts/16_fig1.R        or  source() it in RStudio
# All data, draws and panel constructors live in fig_common.R (shared with 17_fig2.R).
# =============================================================================
.this_dir <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) return(dirname(normalizePath(
    gsub("~+~", " ", sub("^--file=", "", fa[1]), fixed = TRUE), mustWork = FALSE)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
})
source(file.path(.this_dir, "fig_common.R"))
FIG1 <- (f1a | f1b) / f1c +
  patchwork::plot_layout(heights = c(1, 0.95), guides = "collect") +
  patchwork::plot_annotation(tag_levels = list(c("a", "b", "c"))) &
  theme(plot.tag = element_text(size = 10, face = "bold"),
        legend.position = "bottom", legend.box = "horizontal")
save_fig(FIG1, "FIG1_decoupling", 150)


