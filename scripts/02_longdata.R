# =============================================================================
# 02_longdata.R - Convert raw oxygen CSVs from wide to long format
# =============================================================================
# Reads:  data/*_Oxygen.csv         (e.g. Clade1_22_Oxygen.csv, glab_30_Oxygen.csv)
#         Each CSV (one plate = one group at one temperature) has columns:
#             Time, T, OTU<c>_R1, ..., OTU<c>_R5   for the plate's 3 isolates
#         i.e. one column per (isolate, replicate) combination: 3 isolates x 5
#         replicates = 15 oxygen series per file. OTU codes are unique isolate
#         codes (1..18) assigned by 01_convert_xlsx.R.
#         Temperature is carried in the T column (the converter fills it from the
#         measured Tm) and, as a fallback, in the number before "_Oxygen" in the
#         file name.
#
# Writes:
#   tables/Oxygen_All_Long.csv
#
# Output columns:
#   File, Time, T, OTU, Replicate, Oxygen
#
# The wide isolate/replicate columns are split into a numeric OTU (isolate code
# 1..18) and a character Replicate ("R1".."R5"). Downstream scripts key every
# series on (T, OTU, Replicate).
# =============================================================================


# =============================================================================
# 0) Source shared config
# =============================================================================

# Packaging: resolve the script directory from --file= FIRST. Under
# `Rscript scripts/<this>.R` neither rstudioapi nor sys.frame(1)$ofile resolves,
# so this fell back to getwd() and then died on "cannot open file .../config.R".
# Sourcing it from run_all.R was unaffected, which is why the bug stayed hidden.
.this_dir <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) return(dirname(normalizePath(
    # R replaces every space in --file= with "~+~", so a project path that
    # contains a space comes back mangled and source() then fails on a path
    # that does not exist. Un-mangle it before normalising.
    gsub("~+~", " ", sub("^--file=", "", fa[1]), fixed = TRUE),
    mustWork = FALSE)))
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
})

source(file.path(.this_dir, "config.R"))


# =============================================================================
# 1) Load required packages
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)


# =============================================================================
# 2) Check directories
# =============================================================================

message("Working directory: ", getwd())
message("Reading data from: ", data_dir)
message("Writing tables to: ", tables_dir)

if (!dir.exists(data_dir)) {
  stop("Data directory does not exist: ", data_dir)
}

if (!dir.exists(tables_dir)) {
  dir.create(tables_dir, recursive = TRUE)
  message("Created tables directory: ", tables_dir)
}


# =============================================================================
# 3) Find oxygen CSV files
# =============================================================================
# Match any "*_Oxygen.csv" file (e.g. Clade1_22_Oxygen.csv), tolerating a stray
# space before the .csv extension. These are produced by 01_convert_xlsx.R.
# =============================================================================

all_files <- list.files(
  path = data_dir,
  recursive = TRUE,
  full.names = TRUE
)

files <- all_files[
  stringr::str_detect(
    basename(all_files),
    stringr::regex(
      "_Oxygen\\s*\\.csv\\s*$",
      ignore_case = TRUE
    )
  )
]

if (length(files) == 0) {
  stop(
    paste(
      "No matching *_Oxygen.csv files found in", data_dir,
      "\nExpected names like Clade1_22_Oxygen.csv (run 01_convert_xlsx.R first).",
      "\nCheck file names with: list.files(data_dir, recursive = TRUE)"
    )
  )
}

message("Found files:")
print(files)


# =============================================================================
# 4) Function to read one oxygen file
# =============================================================================

read_one <- function(path) {

  df <- readr::read_csv(path, show_col_types = FALSE)

  # Clean column names (strip a UTF-8 BOM and surrounding whitespace).
  names(df) <- stringr::str_replace_all(names(df), "^\\xef\\xbb\\xbf", "")
  names(df) <- stringr::str_trim(names(df))

  # Extract temperature from filename, e.g. "Clade1_22_Oxygen.csv" -> 22
  # (the number just before "_Oxygen"). Fallback only: the converter already
  # fills the T column, so this is used only when T is missing.
  temp_from_name <- stringr::str_match(
    basename(path),
    stringr::regex("(\\d+)\\s*_Oxygen", ignore_case = TRUE)
  )[, 2]

  temp_from_name <- suppressWarnings(
    as.numeric(stringr::str_trim(temp_from_name))
  )

  # Add T if missing or fully NA
  if (!"T" %in% names(df)) {
    df <- df %>%
      mutate(T = temp_from_name)
  } else if (all(is.na(df$T))) {
    df <- df %>%
      mutate(T = temp_from_name)
  }

  # Add source filename
  df <- df %>%
    mutate(File = basename(path), .before = 1)

  return(df)
}


# =============================================================================
# 5) Read and bind all raw files
# =============================================================================

raw <- purrr::map_dfr(files, read_one)

message("\nRaw data preview:")
print(head(raw))

message("\nRaw columns:")
print(names(raw))


# =============================================================================
# 6) Define isolate/replicate column pattern
# =============================================================================
# Expected oxygen columns: "OTU<c>_R<r>", e.g. OTU1_R1, OTU18_R5 (c = isolate
# code 1..18, r = replicate 1..5).
# =============================================================================

series_pattern <- "^OTU(\\d+)_R(\\d+)$"

series_cols <- names(raw)[
  stringr::str_detect(names(raw), series_pattern)
]

if (length(series_cols) == 0) {
  stop(
    paste(
      "No isolate/replicate columns matched the expected pattern:",
      series_pattern,
      "\nExpected column names like OTU1_R1, OTU2_R3, ..., OTU18_R5.",
      "\nAvailable columns are:",
      paste(names(raw), collapse = ", ")
    )
  )
}

message("\nMatched otu/replicate columns:")
print(series_cols)


# =============================================================================
# 7) Convert from wide to long format and split OTU / Replicate
# =============================================================================
# pivot_longer collapses every OTUc_Rr column into a single "series" name,
# then we split that name into a numeric OTU and a character Replicate.
# =============================================================================

long_data <- raw %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(series_cols),
    names_to = "series",
    values_to = "Oxygen"
  ) %>%
  tidyr::drop_na(Oxygen) %>%
  dplyr::mutate(
    OTU     = as.integer(stringr::str_match(series, series_pattern)[, 2]),
    Replicate = paste0("R", stringr::str_match(series, series_pattern)[, 3])
  ) %>%
  dplyr::select(-series)


# =============================================================================
# 7b) Apply data-scope filters from config.R
# =============================================================================
# EXCLUDE_TEMPS drops the listed temperatures (degC); EXCLUDE_OTUS drops the
# listed otus. NULL / empty for either means "drop none".
# =============================================================================

if (exists("EXCLUDE_TEMPS") && !is.null(EXCLUDE_TEMPS) && length(EXCLUDE_TEMPS) > 0) {
  before <- nrow(long_data)
  long_data <- long_data %>% filter(!T %in% as.numeric(EXCLUDE_TEMPS))
  message(sprintf(
    "Filter EXCLUDE_TEMPS = {%s}: %d / %d rows retained",
    paste(EXCLUDE_TEMPS, collapse = ", "), nrow(long_data), before
  ))
  if (nrow(long_data) == 0) {
    stop("No rows left after EXCLUDE_TEMPS filter. Check temperatures in config.R.")
  }
}

if (exists("EXCLUDE_OTUS") && !is.null(EXCLUDE_OTUS) && length(EXCLUDE_OTUS) > 0) {
  before <- nrow(long_data)
  long_data <- long_data %>% filter(!OTU %in% as.integer(EXCLUDE_OTUS))
  message(sprintf(
    "Filter EXCLUDE_OTUS = {%s}: %d / %d rows retained",
    paste(EXCLUDE_OTUS, collapse = ", "), nrow(long_data), before
  ))
  if (nrow(long_data) == 0) {
    stop("No rows left after EXCLUDE_OTUS filter. Check otus in config.R.")
  }
}


# =============================================================================
# 8) Final tidy output
# =============================================================================

long_data <- long_data %>%
  arrange(File, T, OTU, Replicate, Time) %>%
  dplyr::select(dplyr::any_of(c(
    "File",
    "Time",
    "T",
    "OTU",
    "Replicate",
    "Oxygen"
  )))


# =============================================================================
# 9) Save long-format oxygen table
# =============================================================================

out_table <- file.path(tables_dir, "Oxygen_All_Long.csv")

readr::write_csv(long_data, out_table)

message("\nOutput written: ", out_table)


# =============================================================================
# 10) Final checks
# =============================================================================

message("\nLong data preview:")
print(head(long_data))

message("\nLong data structure:")
dplyr::glimpse(long_data)

message("\nCounts by File, T, OTU, Replicate:")
print(
  long_data %>%
    count(File, T, OTU, Replicate),
  n = Inf
)

message("\nDone: 02_longdata.R completed successfully.")
