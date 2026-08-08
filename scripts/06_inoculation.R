# =============================================================================
# 06_inoculation.R - Inoculation density per clade / species group
# =============================================================================
# WHY THIS EXISTS
#   Everything was inoculated at the same OD. But OD is a BIOMASS measure (light
#   scattering), not a cell count - so at the SAME OD you get MORE small cells
#   than large ones. C. glabrata cells (~1-4 um) are much smaller than
#   C. parapsilosis (~2.5-4 x 5-8 um), so an identical OD means glabrata went in
#   at a considerably higher CELL density.
#
#   A single global N_inoculation_cells_per_L (as in config.R) is therefore wrong
#   by construction. This app lets you set it per group.
#
# WHY IT MATTERS
#   N0 = N_inoc x exp(r x delta), and respiration = K / N0. So respiration scales
#   as 1/N_inoc: get the inoculum wrong for a species and its per-cell respiration
#   - and hence its CUE - is wrong by the same factor. Growth is NOT affected
#   (it depends on cell carbon, not cell number).
#
#   Note the useful consequence: if you set N_inoc consistently from OD (i.e.
#   N_inoc proportional to 1/cell_volume), then cell volume CANCELS out of CUE -
#   growth scales with volume and respiration scales with volume too. Your CUE
#   then rests on the inoculum BIOMASS (which the OD gives you) rather than on a
#   guessed cell size.
#
# TWO WAYS TO ENTER IT
#   1. FROM OD (recommended): give the inoculation OD, and for each group the
#      conversion factor "cells per mL at OD 1.0". Smaller cells => bigger factor.
#      N_inoc [cells/L] = OD x factor x 1000
#   2. DIRECTLY: type the cells/L for each group (e.g. from a haemocytometer count
#      - by far the best option if you have it).
#
# Writes: tables/otu_inoc.csv  (one row per isolate; 07_oxygen_fits.R reads it and
#         uses each isolate's own N_inoc instead of the global constant).
#
# RUN:  RStudio: open this file -> "Run App".   Then re-run 07_oxygen_fits.R.
#       Terminal: Rscript scripts/06_inoculation.R --app   <- --app is REQUIRED
#       (run_all.R source()s this file headless; it defines the app and returns)
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

need <- c("shiny", "readr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) {
  stop("Install first:\n  install.packages(c(",
       paste(sprintf('"%s"', miss), collapse = ", "), "))")
}
library(shiny)

out_csv <- file.path(tables_dir, "otu_inoc.csv")

# ---- isolates + groups -------------------------------------------------------
nm_csv <- file.path(tables_dir, "otu_names.csv")
if (!file.exists(nm_csv)) stop("Not found: ", nm_csv, "\nRun 01_convert_xlsx.R first.")
NM <- readr::read_csv(nm_csv, show_col_types = FALSE)
if (!"group" %in% names(NM)) stop("otu_names.csv has no `group` column - re-run 01.")
NM <- data.frame(OTU = as.integer(NM$OTU), isolate = as.character(NM$otu_name),
                 group = as.character(NM$group), stringsAsFactors = FALSE)
GROUPS <- unique(NM$group)

# Cell sizes (if 05_cell_sizes.R has been run) - used only to SUGGEST a factor.
VOL <- stats::setNames(rep(NA_real_, length(GROUPS)), GROUPS)
sz_csv <- file.path(tables_dir, "otu_cell_sizes.csv")
if (file.exists(sz_csv)) {
  .sz <- tryCatch(readr::read_csv(sz_csv, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(.sz) && "cell_volume_um3" %in% names(.sz)) {
    for (g in GROUPS) {
      o <- NM$OTU[NM$group == g]
      v <- suppressWarnings(as.numeric(.sz$cell_volume_um3[as.integer(.sz$OTU) %in% o]))
      v <- v[is.finite(v)]
      if (length(v)) VOL[[g]] <- stats::median(v)
    }
  }
}

# Reference: the textbook "OD 1.0 = 1-3e7 cells/mL" is calibrated on S. cerevisiae
# / C. albicans, which are ~50-65 um3 cells. Cells per OD scale as 1/volume, so a
# smaller cell gives a LARGER factor.
REF_VOL    <- 55      # um3, the cell the textbook factor refers to
REF_FACTOR <- 2e7     # cells per mL at OD 1.0 for a REF_VOL cell

suggest_factor <- function(g) {
  v <- VOL[[g]]
  if (!is.finite(v) || v <= 0) return(REF_FACTOR)
  signif(REF_FACTOR * (REF_VOL / v), 3)
}

# ---- previous values ---------------------------------------------------------
prev <- if (file.exists(out_csv)) {
  tryCatch(readr::read_csv(out_csv, show_col_types = FALSE), error = function(e) NULL)
} else NULL
prev_val <- function(g, col, default) {
  if (is.null(prev) || !all(c("group", col) %in% names(prev))) return(default)
  v <- suppressWarnings(as.numeric(prev[[col]][prev$group == g]))
  v <- v[is.finite(v)]
  if (length(v)) v[1] else default
}
OD_DEFAULT <- prev_val(GROUPS[1], "od", 0.0005)


ui <- fluidPage(
  titlePanel("Inoculation density per clade / species"),
  sidebarLayout(
    sidebarPanel(
      width = 5,
      radioButtons("mode", "How do you want to set it?",
                   choices = c("From OD (same OD, different cell sizes)" = "od",
                               "Directly, in cells per litre" = "direct"),
                   selected = "od"),
      conditionalPanel(
        "input.mode == 'od'",
        numericInput("od", "Inoculation OD (the same for every group):",
                     value = OD_DEFAULT, min = 0, step = 0.0001),
        helpText(tags$b("Cells per mL at OD 1.0"), " for each group. OD measures ",
                 "BIOMASS, so smaller cells give a LARGER value. Suggestions below ",
                 "are scaled from a ", REF_VOL, " um3 reference cell at ",
                 format(REF_FACTOR, scientific = TRUE), " cells/mL/OD, using the ",
                 "volumes you set in 05_cell_sizes.R."),
        uiOutput("factor_inputs")
      ),
      conditionalPanel(
        "input.mode == 'direct'",
        helpText("Type the inoculation density directly (cells per LITRE). ",
                 "A haemocytometer count is by far the best source."),
        uiOutput("direct_inputs")
      ),
      tags$hr(),
      actionButton("save", "Save inoculation densities", class = "btn-primary"),
      tags$br(), tags$br(),
      verbatimTextOutput("status")
    ),
    mainPanel(
      width = 7,
      h4("What will be used"),
      tableOutput("preview"),
      helpText(tags$b("Respiration scales as 1 / N_inoc"), ", so these values set ",
               "the absolute level of respiration and CUE. Growth is NOT affected ",
               "(it depends on cell carbon, not cell number). Re-run ",
               "07_oxygen_fits.R after saving."),
      tags$hr(),
      h4("Current global value in config.R"),
      verbatimTextOutput("global")
    )
  )
)

server <- function(input, output, session) {

  output$global <- renderText(
    sprintf("N_inoculation_cells_per_L = %s  (applied to every isolate until you save here)",
            format(N_inoculation_cells_per_L, scientific = TRUE)))

  output$factor_inputs <- renderUI({
    lapply(GROUPS, function(g) {
      lab <- if (is.finite(VOL[[g]])) sprintf("%s  (cell volume %.1f um3):", g, VOL[[g]])
             else sprintf("%s:", g)
      numericInput(paste0("f_", g), lab,
                   value = prev_val(g, "cells_per_mL_per_OD", suggest_factor(g)),
                   min = 0, step = 1e6)
    })
  })

  output$direct_inputs <- renderUI({
    lapply(GROUPS, function(g) {
      numericInput(paste0("d_", g), sprintf("%s  (cells per litre):", g),
                   value = prev_val(g, "N_inoc_cells_per_L", N_inoculation_cells_per_L),
                   min = 0, step = 1e6)
    })
  })

  vals <- reactive({
    rows <- lapply(GROUPS, function(g) {
      if (identical(input$mode, "od")) {
        f  <- input[[paste0("f_", g)]]
        od <- input$od
        if (is.null(f) || is.null(od) || !is.finite(f) || !is.finite(od)) return(NULL)
        n <- od * f * 1000                      # cells/mL -> cells/L
        data.frame(group = g, od = od, cells_per_mL_per_OD = f,
                   N_inoc_cells_per_L = n, stringsAsFactors = FALSE)
      } else {
        n <- input[[paste0("d_", g)]]
        if (is.null(n) || !is.finite(n)) return(NULL)
        data.frame(group = g, od = NA_real_, cells_per_mL_per_OD = NA_real_,
                   N_inoc_cells_per_L = n, stringsAsFactors = FALSE)
      }
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (!length(rows)) return(NULL)
    do.call(rbind, rows)
  })

  output$preview <- renderTable({
    v <- vals(); req(!is.null(v))
    ref <- N_inoculation_cells_per_L
    data.frame(
      Group = v$group,
      `Cell volume (um3)` = ifelse(is.finite(VOL[v$group]),
                                   round(VOL[v$group], 1), NA),
      `N_inoc (cells/L)`  = formatC(v$N_inoc_cells_per_L, format = "e", digits = 2),
      `vs config global`  = sprintf("%.2fx", v$N_inoc_cells_per_L / ref),
      `Respiration changes by` = sprintf("%.2fx", ref / v$N_inoc_cells_per_L),
      check.names = FALSE)
  })

  observeEvent(input$save, {
    v <- vals()
    if (is.null(v)) { output$status <- renderText("Nothing to save."); return() }
    out <- merge(NM, v, by = "group", all.x = TRUE)
    out <- out[order(out$OTU), c("OTU", "isolate", "group", "od",
                                 "cells_per_mL_per_OD", "N_inoc_cells_per_L")]
    ok <- tryCatch({ readr::write_csv(out, out_csv); TRUE },
                   error = function(e) { output$status <- renderText(
                     paste("FAILED:", conditionMessage(e))); FALSE })
    if (isTRUE(ok)) {
      output$status <- renderText(paste0(
        "Saved ", nrow(out), " isolates to:\n", out_csv,
        "\n\nNow re-run 07_oxygen_fits.R - it will use each isolate's own N_inoc.",
        "\nGrowth will not change; respiration and CUE will."))
    }
  })
}

# ---- headless guard --------------------------------------------------------
# THIS is the one that makes run_all.R hang: run_all.R source()s 06 in its
# loop and under Rscript the old `else` branch called runApp(), so the master
# runner blocked on a Shiny server forever and never reached 07 onwards.
# results/tables/otu_inoc.csv is a committed INPUT read by 07.
# Sourcing this file now always DEFINES ui/server and returns; the app launches
# only when a human asked for it, so nothing it owns is ever regenerated by an
# unattended run.
#   RStudio               -> "Run App", or just source() (interactive)
#   Terminal, on purpose  -> Rscript scripts/06_inoculation.R --app
#                            (or CANDIDAS_RUN_APP=1 Rscript scripts/06_inoculation.R)
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
    message("Launching inoculation app at http://127.0.0.1:7801 ...")
    runApp(shinyApp(ui, server), host = "127.0.0.1", port = 7801, launch.browser = TRUE)
  }
} else {
  message("06_inoculation.R: headless - ui/server defined, app NOT launched, ",
          "nothing written (results/tables/otu_inoc.csv is read as data by 07).\n",
          "  To use it: Rscript scripts/06_inoculation.R --app")
}
