# =============================================================================
# 01_convert_xlsx.R - Interactive PreSens (24-well) xlsx -> tidy CSV converter
# =============================================================================
# Turns the raw PreSens oxygen exports in data/*.xlsx into the wide CSVs that
# 02_longdata.R expects (columns: Time, T, <OTU>_R1, <OTU>_R2, ...).
#
# PLATE LAYOUT for THIS project (Candida temperature x isolate):
#   Each Excel file is ONE 24-well plate at ONE temperature for ONE clade/species
#   GROUP (the group is the file-name prefix, e.g. "Clade1", "glab", "para").
#       - rows A, B, C  = 3 different ISOLATES of that group
#       - columns 1..5  = 5 replicate wells per isolate
#       - row D and column 6 are empty ("No Sensor") and are ignored.
#   So each file becomes 3 isolates x 5 replicates = 15 oxygen series.
#
# In the app you:
#   - tick the files to convert (all recognised groups pre-ticked);
#   - type an ISOLATE NAME for each row A/B/C of every group (free text - e.g.
#     the strain/isolate ID). Every isolate gets a unique code (OTU 1..18) from
#     its (group, row); your typed names are saved to tables/otu_names.csv and
#     shown in every result table and plot;
#   - set the TEMPERATURE (degC) per file (auto-filled from the plate's measured
#     Tm column; edit if needed).
#
# It writes one "<originalname>_Oxygen.csv" per file into data/.
#
# RUN:
#   RStudio: open this file -> click "Run App".
#   Terminal: Rscript scripts/01_convert_xlsx.R --app      <- the --app flag is REQUIRED
#             (without it the file just defines the app and returns, so that
#              run_all.R can source it without blocking on a Shiny server)
#
# You only need to run this once, before 02_longdata.R. After the CSVs exist,
# run the normal pipeline (run_all.R -> 02 -> 03 -> 06).
# =============================================================================

# ---- Locate the project and load shared registry ----------------------------
# --file= is checked FIRST (C1): under `Rscript scripts/01_convert_xlsx.R`
# neither sys.frame(1)$ofile nor rstudioapi is available, so this used to fall
# through to getwd() and then fail on "cannot open file .../config.R".
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
    } else {
      d <- getwd()
    }
  }
  d
})

# config.R provides STRAIN_GROUPS, GROUP_SPECIES, PLATE_ROWS, N_REPLICATES,
# canonical_group() and isolate_code(), so the code assignment is defined in one
# place and matches the rest of the pipeline.
source(file.path(.this_dir, "config.R"))

# ---- Packages ---------------------------------------------------------------
need <- c("shiny", "readxl", "readr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) {
  stop("Install missing package(s) first:\n  install.packages(c(",
       paste(sprintf('"%s"', miss), collapse = ", "), "))")
}
library(shiny)

# =============================================================================
# Core conversion logic (usable without Shiny)
# =============================================================================

# The group for a file is the file-name prefix before the trailing
# "_<temperature>_Oxygen" (case-insensitive), mapped to its canonical name.
group_from_name <- function(fname) {
  stem <- sub("\\.xlsx?$", "", basename(fname), ignore.case = TRUE)
  prefix <- sub("_[0-9]+\\s*_?[Oo]xygen$", "", stem)      # strip _<temp>_Oxygen
  prefix <- sub("_[0-9]+$", "", prefix)                    # strip any trailing _<temp>
  canonical_group(prefix)
}

# Temperature guessed from the file name: the number just before "_Oxygen".
temp_from_name <- function(fname) {
  m <- regmatches(fname, regexpr("(\\d+)(?=\\s*_[Oo]xygen)", fname, perl = TRUE))
  if (length(m) == 1 && nzchar(m)) as.numeric(m) else NA_real_
}

# Read one PreSens sheet: returns Time vector, a data.frame of well traces (one
# column per well, named e.g. "A1"), the well labels, and the plate's measured
# Tm (median of the Tm column).
read_presens_sheet <- function(path) {
  raw <- readxl::read_excel(
    path, sheet = 1, col_names = FALSE,
    col_types = "text", .name_repair = "minimal"
  )
  raw <- as.data.frame(raw, stringsAsFactors = FALSE)

  is_hdr <- apply(raw, 1, function(r) any(trimws(r) == "Date/Time", na.rm = TRUE))
  hdr_row <- which(is_hdr)[1]
  if (is.na(hdr_row)) {
    stop("Could not find a 'Date/Time' header row in ", basename(path))
  }

  hdr <- trimws(as.character(unlist(raw[hdr_row, ])))
  hdr[is.na(hdr)] <- ""

  time_col  <- which(grepl("^Time", hdr))[1]         # "Time/Min."
  well_cols <- which(grepl("^[A-Da-d][1-6]$", hdr))  # A1 .. D6
  tm_col    <- which(grepl("^Tm", hdr))[1]           # "Tm [C]" measured temp
  if (is.na(time_col) || length(well_cols) == 0) {
    stop("Could not locate the Time and well columns in ", basename(path))
  }

  wells <- toupper(hdr[well_cols])

  dat <- raw[(hdr_row + 1):nrow(raw), , drop = FALSE]
  time_vals <- suppressWarnings(as.numeric(dat[[time_col]]))
  keep <- is.finite(time_vals)
  dat <- dat[keep, , drop = FALSE]
  time_vals <- time_vals[keep]

  # "No Sensor" (and any other non-numeric) cells become NA here.
  well_df <- as.data.frame(
    lapply(well_cols, function(ci) suppressWarnings(as.numeric(dat[[ci]]))),
    stringsAsFactors = FALSE
  )
  names(well_df) <- wells

  tm_val <- if (!is.na(tm_col)) {
    suppressWarnings(stats::median(as.numeric(dat[[tm_col]]), na.rm = TRUE))
  } else NA_real_

  list(time = time_vals, wells = well_df, well_labels = wells, tm = tm_val)
}

# Assemble the tidy wide frame for one plate (one group).
#   group       : canonical group name (e.g. "Clade1").
#   n_rep       : number of replicate columns to read (default N_REPLICATES = 5).
#   temperature : numeric degC written into the T column.
# Rows A/B/C are the 3 isolates; each becomes OTU<code>_R<rep> columns, where
# code = isolate_code(group, row) is unique across the whole dataset.
build_output <- function(sheet, group, n_rep = N_REPLICATES, temperature) {
  wl   <- sheet$well_labels
  rowL <- toupper(substr(wl, 1, 1))
  colN <- as.integer(substr(wl, 2, nchar(wl)))

  get_well <- function(rl, cn) {
    idx <- which(rowL == rl & colN == cn)
    if (length(idx) == 0) return(NULL)
    sheet$wells[[idx]]
  }

  out <- data.frame(Time = sheet$time, T = temperature, check.names = FALSE)

  for (rep in seq_len(n_rep)) {
    for (rl in PLATE_ROWS) {                 # A, B, C = 3 isolates
      code <- isolate_code(group, rl)
      nm   <- paste0("OTU", code, "_R", rep)
      col  <- get_well(rl, rep)              # column `rep` = replicate R<rep>
      out[[nm]] <- if (is.null(col)) NA_real_ else col
    }
  }
  out
}

# Convert one file end-to-end and (optionally) write it next to the source.
convert_one <- function(path, group, n_rep = N_REPLICATES, temperature,
                        write = TRUE) {
  sheet <- read_presens_sheet(path)
  out <- build_output(sheet, group, n_rep, temperature)
  if (write) {
    out_path <- sub("\\.xlsx?$", ".csv", path, ignore.case = TRUE)
    readr::write_csv(out, out_path)
    attr(out, "out_path") <- out_path
  }
  out
}

# =============================================================================
# Shiny app
# =============================================================================

xlsx_files <- list.files(data_dir, pattern = "\\.xlsx?$", ignore.case = TRUE,
                         full.names = FALSE)
file_groups <- vapply(xlsx_files, group_from_name, character(1))
known_files <- xlsx_files[!is.na(file_groups)]   # files whose prefix is a known group
preselect   <- known_files

ui <- fluidPage(
  titlePanel("PreSens 24-well xlsx -> CSV converter (Candida: temperature x isolate)"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      checkboxGroupInput(
        "files", "Files to convert:",
        choices = xlsx_files, selected = preselect
      ),
      helpText("Each file is one plate: rows A/B/C = 3 isolates, columns 1-5 = ",
               "replicates. Files whose name does not start with a known group ",
               "are left unticked."),
      numericInput("nrep", "Number of replicates (columns):",
                   value = N_REPLICATES, min = 1, max = 6, step = 1),
      tags$hr(),
      strong("Isolate name for each row, per group:"),
      helpText("Type the isolate ID for rows A, B and C of each group. Each ",
               "isolate gets a unique code (OTU 1-18); your names are saved to ",
               "tables/otu_names.csv and shown in every table and plot."),
      uiOutput("strain_inputs"),
      tags$hr(),
      strong("Temperature (degC) per file:"),
      helpText("Auto-filled from the plate's measured Tm column; edit if needed."),
      uiOutput("temp_inputs"),
      tags$hr(),
      actionButton("convert", "Convert to CSV", class = "btn-primary"),
      tags$br(), tags$br(),
      verbatimTextOutput("status")
    ),
    mainPanel(
      width = 8,
      h4("Isolate code map"),
      tableOutput("name_map"),
      h4("Preview (first selected file)"),
      tableOutput("preview"),
      helpText("CSV columns are OTU<code>_R1 ... OTU<code>_R5 (internal codes), ",
               "grouped by replicate. Codes map to your names via ",
               "tables/otu_names.csv.")
    )
  )
)

server <- function(input, output, session) {

  # Groups present among the currently-selected files, in registry order.
  groups_present <- reactive({
    req(input$files)
    g <- unique(stats::na.omit(vapply(input$files, group_from_name, character(1))))
    STRAIN_GROUPS[STRAIN_GROUPS %in% g]
  })

  # One text box per (group, row A/B/C), defaulting to "<group>_iso<1..3>".
  output$strain_inputs <- renderUI({
    gp <- groups_present()
    if (length(gp) == 0) return(helpText("Select at least one recognised file."))
    do.call(tagList, lapply(gp, function(g) {
      tagList(
        tags$b(sprintf("%s  (%s)", g, unname(GROUP_SPECIES[g]))),
        lapply(seq_along(PLATE_ROWS), function(i) {
          rl <- PLATE_ROWS[i]
          textInput(
            paste0("strain_", g, "_", rl),
            label = sprintf("row %s -> OTU%d:", rl, isolate_code(g, rl)),
            value = sprintf("%s_iso%d", g, i)
          )
        }),
        tags$hr()
      )
    }))
  })

  # One numeric box per selected file, defaulting to the measured Tm (falling
  # back to the temperature parsed from the file name).
  output$temp_inputs <- renderUI({
    req(input$files)
    lapply(input$files, function(f) {
      guess <- tryCatch(read_presens_sheet(file.path(data_dir, f))$tm,
                        error = function(e) NA_real_)
      if (is.na(guess)) guess <- temp_from_name(f)
      numericInput(paste0("temp_", make.names(f)), label = f, value = guess)
    })
  })

  # The typed display name for (group, row); blank falls back to the default.
  strain_display <- reactive({
    gp <- groups_present()
    rows <- list()
    for (g in gp) for (i in seq_along(PLATE_ROWS)) {
      rl <- PLATE_ROWS[i]
      v  <- input[[paste0("strain_", g, "_", rl)]]
      nm <- if (is.null(v) || !nzchar(trimws(v))) sprintf("%s_iso%d", g, i) else trimws(v)
      rows[[length(rows) + 1]] <- data.frame(
        OTU      = isolate_code(g, rl),
        otu_name = nm,
        group    = g,
        species  = unname(GROUP_SPECIES[g]),
        stringsAsFactors = FALSE
      )
    }
    if (length(rows) == 0) return(NULL)
    do.call(rbind, rows)
  })

  name_map_df <- reactive({
    d <- strain_display()
    req(!is.null(d))
    d[order(d$OTU), , drop = FALSE]
  })

  output$name_map <- renderTable({
    d <- name_map_df()
    req(!is.null(d))
    data.frame(
      `CSV code` = paste0("OTU", d$OTU),
      Group      = d$group,
      `Your name`= d$otu_name,
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })

  output$preview <- renderTable({
    req(input$files)
    f <- input$files[1]
    g <- group_from_name(f)
    validate(need(!is.na(g), "First selected file has an unrecognised group prefix."))
    sheet <- read_presens_sheet(file.path(data_dir, f))
    tval  <- input[[paste0("temp_", make.names(f))]]
    head(build_output(sheet, g, input$nrep, tval), 6)
  }, digits = 2)

  observeEvent(input$convert, {
    req(input$files)

    names_msg <- tryCatch({
      readr::write_csv(name_map_df(), file.path(tables_dir, "otu_names.csv"))
      paste0("Wrote name map: ", file.path(tables_dir, "otu_names.csv"))
    }, error = function(e) paste0("WARNING: could not write otu_names.csv: ",
                                  conditionMessage(e)))

    msgs <- character(0)
    for (f in input$files) {
      g <- group_from_name(f)
      if (is.na(g)) {
        msgs <- c(msgs, sprintf("SKIP %s : unrecognised group prefix", f))
        next
      }
      path <- file.path(data_dir, f)
      tval <- input[[paste0("temp_", make.names(f))]]
      res <- tryCatch({
        out <- convert_one(path, g, input$nrep, tval, write = TRUE)
        sprintf("OK  %s [%s] -> %s  (%d rows, %d series, T=%s)",
                f, g, basename(attr(out, "out_path")),
                nrow(out), ncol(out) - 2, tval)
      }, error = function(e) sprintf("FAIL %s : %s", f, conditionMessage(e)))
      msgs <- c(msgs, res)
    }
    output$status <- renderText(paste(c(msgs, names_msg), collapse = "\n"))
  })
}

# ---- headless guard (C1) ----------------------------------------------------
# This file used to call runApp() whenever it was NOT interactive, which is
# exactly the case under Rscript - so `Rscript scripts/run_all.R` blocked
# forever on a Shiny server nobody could see. Now: sourcing this file always
# DEFINES ui/server and returns; the app launches only when a human asked for
# it. Nothing is written unless the app runs, so the committed converter output
# is never touched by an unattended run.
#   RStudio               -> "Run App", or just source() (interactive)
#   Terminal, on purpose  -> Rscript scripts/01_convert_xlsx.R --app
#                            (or CANDIDAS_RUN_APP=1 Rscript scripts/01_convert_xlsx.R)
#   run_all.R / any batch -> defines objects, launches nothing
.candidas_run_app <- function() {
  if (isTRUE(getOption("candidas.headless")))            return(FALSE)
  if (identical(Sys.getenv("CANDIDAS_HEADLESS"), "1"))   return(FALSE)
  if (identical(Sys.getenv("CANDIDAS_RUN_APP"), "1"))    return(TRUE)
  if ("--app" %in% commandArgs(trailingOnly = TRUE))     return(TRUE)
  interactive()
}

if (.candidas_run_app()) {
  if (interactive()) {
    shinyApp(ui, server)
  } else {
    message("Launching converter app at http://127.0.0.1:7788 ...")
    runApp(shinyApp(ui, server), host = "127.0.0.1", port = 7788,
           launch.browser = TRUE)
  }
} else {
  message("01_convert_xlsx.R: headless - ui/server defined, app NOT launched, ",
          "nothing written.\n  To use it: Rscript scripts/01_convert_xlsx.R --app")
}
