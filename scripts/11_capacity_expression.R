#!/usr/bin/env Rscript
# =====================================================================
# Does low ENZYME CAPACITY show up as low METABOLIC-GENE EXPRESSION?
# Test on the largest capacity gap in the model: Clade II (kcat_scale 0.41)
# vs Clade I (0.61), using the only public C. auris clade-resolved RNA-seq
# at a controlled condition: Bing et al. 2021 (GEO GSE165762, 37C YPD).
#
# Prediction from the etc-GEM: Clade II runs ~40% less effective enzyme
# capacity than Clade I. If capacity is an EXPRESSION/allocation trait,
# the model's 928 central-metabolic enzyme genes should be coordinately
# LOWER in Clade II than Clade I (positive log2FC for Clade I vs II),
# and shifted BELOW the genome background (that specificity is the real
# test -- Clade II grows slower so many genes move; capacity predicts the
# metabolic enzymes move MORE than the genome-wide shift).
#
# Samples (GSE165762): Clade I = 470140, 470147, 470154 (3 isolates x 3 reps);
#                      Clade II = CBS10913 (CBS109, 3 reps). Gene IDs = CJI97_*.
# Your model's enzyme_gene_map.csv `gene_id` column is ALSO CJI97_* -> direct join.
#
# RUN (from anywhere):
#   Rscript scripts/11_capacity_expression.R
# Needs ONLY edgeR (Bioconductor). GEOquery is NOT required anymore.
#
# C1 EXECUTION FIXES (no analysis change):
#   * the hard-coded /Users/ilgazcakin/... fallback is gone; every path is now
#     derived from base_dir (config.R), so this runs on any machine;
#   * the mid-script setwd() is gone. It made every relative path below depend
#     on where R happened to be started, and left the caller's working directory
#     silently changed when run_all.R source()d this file. All paths are explicit;
#   * outputs honour CANDIDAS_EXPRESSION_OUT so a non-destructive re-run does not
#     write into data/.
# The GEO count matrix is COMMITTED at data/expression/GSE165762/ - this script
# does not need to download anything. The download is kept only as a fallback.
# =====================================================================

.this_dir <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) return(dirname(normalizePath(
    # R replaces every space in --file= with "~+~". A project path with a
    # space in it (this one has two) therefore comes back mangled, and
    # source() then fails on a path that does not exist.
    gsub("~+~", " ", sub("^--file=", "", fa[1]), fixed = TRUE),
    mustWork = FALSE)))
  d <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (length(d) == 0 || is.na(d) || !nzchar(d)) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
        nzchar(rstudioapi::getActiveDocumentContext()$path)) {
      d <- dirname(rstudioapi::getActiveDocumentContext()$path)
    } else d <- getwd()
  }
  d
})
source(file.path(.this_dir, "config.R"))   # base_dir, EXPRESSION_DIR, EXPRESSION_OUT_DIR

suppressMessages({
  if (!requireNamespace("edgeR", quietly=TRUE)) stop("install edgeR: BiocManager::install('edgeR')")
  library(edgeR)
})

# --- inputs: all under data/expression/, resolved from base_dir --------------
GENELIST <- file.path(EXPRESSION_DIR, "model_metabolic_genes_CJI97.txt")
if (!file.exists(GENELIST))
  stop("EXPRESSION_GENELIST_MISSING: ", GENELIST, "\n",
       "  This file is committed to the repository; a missing copy means an\n",
       "  incomplete checkout. Re-clone or restore data/expression/.", call. = FALSE)
met_genes <- readLines(GENELIST)

OUT_DIR <- EXPRESSION_OUT_DIR
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Gene list : %s\nOutput dir: %s\nModel metabolic enzyme genes: %d\n",
            GENELIST, OUT_DIR, length(met_genes)))

# ---- 1. the GEO processed count matrix (COMMITTED - no download needed) -----
cf <- file.path(EXPRESSION_DIR, "GSE165762", "GSE165762_Raw_counts.txt.gz")
if (file.exists(cf)) {
  cat("Using committed GEO count matrix (no download):\n  ", cf, "\n")
} else {
  dir.create(dirname(cf), showWarnings=FALSE, recursive=TRUE)
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165762/suppl/GSE165762_Raw_counts.txt.gz"
  cat("Committed copy absent. Downloading", url, "...\n"); options(timeout=1200)
  download.file(url, cf, mode="wb")
}
raw <- read.table(gzfile(cf), header=TRUE, sep="\t", row.names=1,
                  check.names=FALSE, stringsAsFactors=FALSE)
cat(sprintf("Raw count table: %d rows x %d columns\n", nrow(raw), ncol(raw)))
cat("Columns as deposited:\n"); print(colnames(raw))
raw <- raw[grepl("^CJI97_", rownames(raw)), , drop=FALSE]     # keep gene rows only

# ---- 2. assign each column to a clade -------------------------------
clade_of <- function(x){
  if (grepl("470140|470147|470154", x)) return("CladeI")
  if (grepl("CBS109|CBS10913",     x)) return("CladeII")
  NA_character_
}
grp <- vapply(colnames(raw), clade_of, character(1))
if (!all(c("CladeI","CladeII") %in% grp) && ncol(raw) == 12) {
  warning("Column names lack strain tokens; falling back to GEO deposit order ",
          "(cols 1-9 = Clade I 470140/147/154 x3reps; cols 10-12 = Clade II CBS10913). ",
          "VERIFY the printed columns match this order before trusting the result.")
  grp <- c(rep("CladeI",9), rep("CladeII",3))
}
cat("\nColumn -> clade:\n"); print(data.frame(column=colnames(raw), clade=grp), row.names=FALSE)
raw <- raw[, !is.na(grp), drop=FALSE]; grp <- grp[!is.na(grp)]
stopifnot(all(c("CladeI","CladeII") %in% grp))

# ---- 3. count matrix -----------------------------------------------
M <- as.matrix(raw); mode(M) <- "numeric"
cat(sprintf("Count matrix: %d genes x %d samples (CladeI=%d, CladeII=%d)\n",
            nrow(M), ncol(M), sum(grp=="CladeI"), sum(grp=="CladeII")))

# ---- 4. edgeR: Clade I vs Clade II per-gene log2FC ------------------
d <- DGEList(counts=M, group=grp)
keep <- rowSums(cpm(d) > 1) >= 3; d <- d[keep,,keep.lib.sizes=FALSE]
d <- calcNormFactors(d); d <- estimateDisp(d)
et <- exactTest(d, pair=c("CladeII","CladeI"))   # +logFC = higher in Clade I = lower in low-capacity Clade II
tt <- topTags(et, n=Inf)$table; tt$gene <- rownames(tt)
tt$is_metab <- tt$gene %in% met_genes

# ---- 5. THE TEST ---------------------------------------------------
mfc <- tt$logFC[tt$is_metab]; ofc <- tt$logFC[!tt$is_metab]
w <- wilcox.test(mfc, ofc)
cat("\n================ RESULT ================\n")
cat(sprintf("metabolic enzyme genes expressed & tested: %d\n", length(mfc)))
cat(sprintf("median log2FC  metabolic  (CladeI vs II): %+.3f\n", median(mfc)))
cat(sprintf("median log2FC  genome background        : %+.3f\n", median(ofc)))
cat(sprintf("metabolic-vs-background (specificity): Wilcoxon W=%.0f, p=%.2e\n", w$statistic, w$p.value))
cat(sprintf("fraction of metabolic genes higher in Clade I: %.1f%%\n", 100*mean(mfc>0)))
cat("\nREAD IT LIKE THIS:\n")
cat("  Capacity=expression is SUPPORTED if metabolic median > background median\n")
cat("  AND the Wilcoxon specificity p < 0.05 (metabolic enzymes down in Clade II\n")
cat("  by MORE than the genome-wide slow-growth shift). If metabolic == background,\n")
cat("  that's just the growth-rate effect, not a capacity-allocation signal.\n")

.out_csv <- file.path(OUT_DIR, "capacity_expression_CladeI_vs_II.csv")
.out_pdf <- file.path(OUT_DIR, "capacity_expression_barcode.pdf")
write.csv(tt[,c("gene","logFC","logCPM","PValue","FDR","is_metab")],
          .out_csv, row.names=FALSE)
pdf(.out_pdf, width=7, height=4)
o <- order(tt$logFC)
plot(tt$logFC[o], type="n", xlab="genes ranked by log2FC (CladeI vs II)",
     ylab="log2 fold change", main="Metabolic enzymes (red) vs genome")
abline(h=0, col="grey70")
points(which(!tt$is_metab[o]), tt$logFC[o][!tt$is_metab[o]], pch=".", col="grey75")
points(which(tt$is_metab[o]),  tt$logFC[o][tt$is_metab[o]],  pch="|", col="red", cex=0.6)
dev.off()
cat("\nWrote:\n  ", .out_csv, "\n  ", .out_pdf, "\n")
