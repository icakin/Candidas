# =============================================================================
# 05_cell_sizes.R - Enter the average cell VOLUME of each isolate
# =============================================================================
# Candida isolates differ in size, and cell size sets the CARBON PER CELL,
# which scales absolute growth (fg C / h) and CUE (it does NOT change
# respiration_fgC_h or the temperature-response shapes / activation energies).
#
# You enter the average cell VOLUME (in cubic micrometres) for each isolate, plus
# a carbon density. The app computes carbon per cell = volume x density and saves:
#
#     tables/otu_cell_sizes.csv   (columns: OTU, cell_volume_um3,
#                                  carbon_density_fg_per_um3, cell_carbon_fg)
#
# Isolate DISPLAY NAMES are set in 01_convert_xlsx.R (Step 1), not here; this app
# shows them next to each isolate code if tables/otu_names.csv exists.
#
# 07_oxygen_fits.R reads that file and uses each isolate's own cell carbon. Any
# isolate not in the file falls back to the single global cell size in config.R,
# so the pipeline still runs even if you skip this step.
#
# RUN (like the other apps):
#   RStudio: open this file -> "Run App".
#   Terminal: Rscript scripts/05_cell_sizes.R
#
# You can run this any time before 07_oxygen_fits.R (it does not depend on the
# other steps). Re-run 06 afterwards to apply new sizes.
# =============================================================================

# ---- Locate the project ----------------------------------------------------
.script_dir <- local({
  d <- tryCatch({
    if (requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable() &&
        nzchar(rstudioapi::getActiveDocumentContext()$path)) {
      dirname(rstudioapi::getActiveDocumentContext()$path)
    } else NA_character_
  }, error = function(e) NA_character_)
  if (length(d) == 0 || is.na(d) || !nzchar(d)) {
    d <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  }
  if (length(d) == 0 || is.na(d) || !nzchar(d)) d <- getwd()
  normalizePath(d, mustWork = FALSE)
})

base_dir   <- dirname(.script_dir)
tables_dir <- file.path(base_dir, "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(tables_dir, "otu_cell_sizes.csv")

# ---- Packages --------------------------------------------------------------
need <- c("shiny", "readr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) {
  stop("Install missing package(s) first:\n  install.packages(c(",
       paste(sprintf('"%s"', miss), collapse = ", "), "))")
}
library(shiny)

# ---- Which isolates? Read them from the long data if available, else default -
otu_ids_default <- 1:18
long_csv <- file.path(tables_dir, "Oxygen_All_Long.csv")
if (file.exists(long_csv)) {
  .lg <- tryCatch(readr::read_csv(long_csv, show_col_types = FALSE),
                  error = function(e) NULL)
  if (!is.null(.lg) && "OTU" %in% names(.lg)) {
    ids <- suppressWarnings(sort(unique(as.integer(.lg$OTU))))
    ids <- ids[is.finite(ids)]
    if (length(ids)) otu_ids_default <- ids
  }
}

# ---- Isolate display names (from 01_convert_xlsx.R), for nicer labels --------
otu_name_of <- function(o) paste0("OTU", o)
names_csv <- file.path(tables_dir, "otu_names.csv")
if (file.exists(names_csv)) {
  .nm <- tryCatch(readr::read_csv(names_csv, show_col_types = FALSE),
                  error = function(e) NULL)
  if (!is.null(.nm) && all(c("OTU", "otu_name") %in% names(.nm))) {
    .map <- stats::setNames(as.character(.nm$otu_name), as.character(as.integer(.nm$OTU)))
    otu_name_of <- function(o) {
      nm <- unname(.map[as.character(o)])
      ifelse(is.na(nm) | !nzchar(nm), paste0("OTU", o), nm)
    }
  }
}

# ---- Preload existing values -----------------------------------------------
DENSITY_DEFAULT <- 100   # fg C / um^3 (Candida/yeast default from config.R)
VOLUME_DEFAULT  <- 1     # um^3 (placeholder; enter your measured volume)
prev <- NULL
if (file.exists(out_csv)) {
  prev <- tryCatch(readr::read_csv(out_csv, show_col_types = FALSE),
                   error = function(e) NULL)
}
get_prev <- function(otu, col, default) {
  if (is.null(prev) || !all(c("OTU", col) %in% names(prev))) return(default)
  hit <- prev[[col]][as.integer(prev$OTU) == otu]
  if (length(hit) && is.finite(suppressWarnings(as.numeric(hit[1])))) as.numeric(hit[1]) else default
}

# ---- UI --------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Average cell volume per isolate"),
  sidebarLayout(
    sidebarPanel(
      width = 5,
      helpText("Enter the average cell VOLUME (cubic micrometres, um^3) for ",
               "each strain. Carbon density converts volume to carbon per cell."),
      numericInput("density",
                   "Carbon density (fg C / um^3), applied to all OTUs:",
                   value = get_prev(otu_ids_default[1], "carbon_density_fg_per_um3",
                                    DENSITY_DEFAULT),
                   min = 1, step = 10),
      tags$hr(),
      strong("Average cell volume (um^3):"),
      uiOutput("vol_inputs"),
      tags$hr(),
      actionButton("save", "Save cell sizes", class = "btn-primary"),
      tags$br(), tags$br(),
      verbatimTextOutput("status")
    ),
    mainPanel(
      width = 7,
      h4("Preview"),
      tableOutput("preview"),
      helpText("cell_carbon_fg = cell_volume_um3 x carbon density. This is what ",
               "07_oxygen_fits.R uses to turn growth rate into fg C / h.")
    )
  )
)

# ---- Server ----------------------------------------------------------------
server <- function(input, output, session) {

  otu_ids <- otu_ids_default

  output$vol_inputs <- renderUI({
    lapply(otu_ids, function(o) {
      fluidRow(
        column(4, tags$b(sprintf("OTU%d — %s", o, otu_name_of(o)))),
        column(8, numericInput(paste0("v_", o), "volume (um^3)",
                               value = get_prev(o, "cell_volume_um3", VOLUME_DEFAULT),
                               min = 0.001, step = 0.1))
      )
    })
  })

  size_table <- reactive({
    dens <- if (is.null(input$density) || !is.finite(input$density)) DENSITY_DEFAULT else input$density
    rows <- lapply(otu_ids, function(o) {
      v <- input[[paste0("v_", o)]]
      if (is.null(v) || !is.finite(v) || v <= 0) return(NULL)
      data.frame(OTU = o,
                 cell_volume_um3 = round(v, 4),
                 carbon_density_fg_per_um3 = dens,
                 cell_carbon_fg = round(v * dens, 4),
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  })

  output$preview <- renderTable({ size_table() }, digits = 3)

  observeEvent(input$save, {
    df <- size_table()
    if (is.null(df) || nrow(df) == 0) {
      output$status <- renderText("Nothing to save - enter at least one volume.")
      return()
    }
    ok <- tryCatch({ readr::write_csv(df, out_csv); TRUE },
                   error = function(e) { output$status <- renderText(paste("FAILED:", conditionMessage(e))); FALSE })
    if (isTRUE(ok)) {
      output$status <- renderText(paste0(
        "Saved ", nrow(df), " OTU cell size(s) to:\n", out_csv,
        "\nNow re-run 07_oxygen_fits.R to apply them."))
    }
  })
}

if (interactive()) {
  shinyApp(ui, server)
} else {
  message("Launching cell-size app at http://127.0.0.1:7799 ...")
  runApp(shinyApp(ui, server), host = "127.0.0.1", port = 7799, launch.browser = TRUE)
}
