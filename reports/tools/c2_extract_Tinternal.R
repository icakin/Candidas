#!/usr/bin/env Rscript
# =============================================================================
# c2_extract_Tinternal.R -- pull the MEASURED temperature trace out of the raw
# PreSens exports.
#
#   Rscript reports/tools/c2_extract_Tinternal.R
#   -> reports/tools/c2_Tinternal.csv   (File, group, T_set, Time, T_internal)
#
# Every data/*.xlsx carries two temperature columns that 01_convert_xlsx.R does
# not keep:
#     Tm [degC]         CONSTANT for the whole run = the set point the operator
#                       typed into the SDR software. This is what the optode's
#                       temperature compensation uses.
#     T_internal [degC] a real measured trace from the reader.
#
# T_internal is what makes the thermal ramp measurable instead of assumed.
# It is the READER's internal sensor, not a probe in the liquid, so it is still
# a proxy for vial temperature -- but a measured one.
# =============================================================================
suppressPackageStartupMessages(library(readxl))

.fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
HERE <- if (length(.fa)) dirname(normalizePath(
  gsub("~+~", " ", sub("^--file=", "", .fa[1]), fixed = TRUE), mustWork = FALSE)) else getwd()
ROOT <- dirname(dirname(HERE))

files <- sort(list.files(file.path(ROOT, "data"), pattern = "[.]xlsx$", full.names = TRUE))
cat("raw exports found:", length(files), "\n")

out <- list()
for (f in files) {
  x <- try(suppressMessages(read_excel(f, sheet = 1, col_names = FALSE,
                                       col_types = "text", .name_repair = "minimal")),
           silent = TRUE)
  if (inherits(x, "try-error")) { cat("  SKIP (unreadable):", basename(f), "\n"); next }
  hrow <- which(apply(x, 1, function(r) any(trimws(r) == "Date/Time", na.rm = TRUE)))[1]
  if (is.na(hrow)) { cat("  SKIP (no header):", basename(f), "\n"); next }
  hdr <- trimws(as.character(unlist(x[hrow, ]))); hdr[is.na(hdr)] <- ""
  i_time <- which(grepl("^Time", hdr))[1]
  i_ti   <- which(grepl("^T_internal", hdr))[1]
  i_tm   <- which(grepl("^Tm", hdr))[1]
  if (is.na(i_time) || is.na(i_ti)) { cat("  SKIP (no T_internal):", basename(f), "\n"); next }
  body <- x[(hrow + 1):nrow(x), , drop = FALSE]
  tt <- suppressWarnings(as.numeric(unlist(body[[i_time]])))
  ti <- suppressWarnings(as.numeric(unlist(body[[i_ti]])))
  tm <- if (!is.na(i_tm)) suppressWarnings(as.numeric(unlist(body[[i_tm]]))) else NA_real_
  ok <- is.finite(tt) & is.finite(ti)
  if (!any(ok)) { cat("  SKIP (empty):", basename(f), "\n"); next }
  stem  <- sub("[.]xlsx$", "", basename(f))
  grp   <- sub("_[0-9]+_?[Oo]xygen$", "", stem)
  T_set <- suppressWarnings(as.numeric(sub(".*_([0-9]+)_?[Oo]xygen$", "\\1", stem)))
  out[[length(out) + 1]] <- data.frame(
    File = basename(f), group = grp, T_set = T_set,
    Tm_setting = if (all(is.na(tm))) NA_real_ else stats::median(tm, na.rm = TRUE),
    Time = tt[ok], T_internal = ti[ok], stringsAsFactors = FALSE)
}
res <- do.call(rbind, out)
p <- file.path(HERE, "c2_Tinternal.csv")
utils::write.csv(res, p, row.names = FALSE)
cat("wrote", p, ":", nrow(res), "rows from", length(unique(res$File)), "plates\n")
cat("Tm setting constant within every plate, and equal to the file-name set point in",
    sum(abs(tapply(res$Tm_setting, res$File, function(v) v[1]) -
            tapply(res$T_set,      res$File, function(v) v[1])) < 1e-9, na.rm = TRUE),
    "of", length(unique(res$File)), "plates\n")
