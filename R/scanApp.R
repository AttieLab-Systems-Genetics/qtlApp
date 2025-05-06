#' Scan App Module
#'
#' @param id shiny identifier
#' @param main_par reactive list with selected_dataset, LOD_thr and which_trait
#' @param import reactive list with file_directory and markers
#'
#' @importFrom DT DTOutput renderDT
#' @importFrom shiny actionButton h4 moduleServer nearPoints NS plotOutput
#'             reactive reactiveValues renderPlot renderUI req setProgress shinyApp
#'             uiOutput withProgress
#' @importFrom bslib card card_header page_sidebar sidebar
#' @importFrom shinycssloaders withSpinner
#' @importFrom dplyr across mutate where
#' @importFrom stringr str_split
#'
#' @export
scanApp <- function() {
  ui <- bslib::page_sidebar(
    title = "Test Scan",
    sidebar = bslib::sidebar("side_panel",
      mainParInput("main_par"),    # "selected_dataset", "LOD_thr"
      mainParUI("main_par"),       # "which_trait", "selected_chr"
      peakInput("peak"),           # "which_peak"
      peakUI("peak"),              # MOVED: UI for the allele effects plot
      downloadInput("download"),   # downloadButton, filename
      downloadOutput("download")), # plot_table, inputs for Plots or Tables
    bslib::card(
      bslib::card_header("LOD profile"),
      scanOutput("scan_list")),
    bslib::card( # Card for the peaks table
      bslib::card_header("Peaks Table"),
      peakOutput("peak")       # ADDED: UI for the peaks datatable
    ),
    bslib::card(
      bslib::card_header("Clicked Peak on Scan"), 
      scanUI("scan_list"))
  )
  server <- function(input, output, session) {
      import <- importServer("import")
      main_par <- mainParServer("main_par", import)
      scan_list <- scanServer("scan_list", main_par, import)
      peak_list <- peakServer("peak", main_par, import)
      merged_list <- mergeServer("merged_list", scan_list, peak_list)
      # ** allele plot not working **
      downloadServer("download", merged_list)
  }
  shiny::shinyApp(ui = ui, server = server)
}
#' @rdname scanApp
#' @export
scanServer <- function(id, main_par, import) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Find selected_trait; special handling for `genes` or `isoforms`.
    selected_trait <- shiny::reactive({
      shiny::req(import(), main_par$selected_dataset())
      trait <- shiny::req(main_par$which_trait())
      # if gene or isoform then trait is `symbol_id``.
      trait_type <- get_trait_type(import(), main_par$selected_dataset())
      if(trait_type %in% c("genes", "isoforms")) {
        trait <- stringr::str_remove(main_par$which_trait(), "_.*")
      }
      trait
    })
    scans <- shiny::reactive({
      shiny::req(main_par$selected_dataset(), selected_trait())
      scans <- shiny::withProgress(
        message = paste("scan of", selected_trait(), "in progress"),
        value = 0, {
          shiny::setProgress(1)
          suppressMessages(
            trait_scan(import()$file_directory,
              main_par$selected_dataset(), selected_trait()))
        }
      )
    })
    scan_table <- shiny::reactive({
      shiny::req(scans(), main_par$which_trait(), main_par$LOD_thr())
      QTL_plot_visualizer(
        scans(), main_par$which_trait(), main_par$LOD_thr(), import()$markers)
    })
    scan_table_chr <- shiny::reactive({
      shiny::req(scan_table(), main_par$selected_chr())
      if (main_par$selected_chr() == "All") {
        scan_table()
      } else {
        dplyr::filter(scan_table(), .data$chr == main_par$selected_chr())
      }
    })
    scan_plot <- shiny::reactive({
      shiny::req(scan_table_chr(), main_par$LOD_thr(), main_par$selected_chr())
      ggplot_qtl_scan(scan_table_chr(), main_par$LOD_thr(), main_par$selected_chr())
    })
    output$scan_plot <- shiny::renderUI({
      shiny::req(scan_plot())
      shiny::plotOutput(ns("render_plot"), click = ns("plot_click")) |>
        shinycssloaders::withSpinner(type = 8, color = "#3498db")
    })
    output$render_plot <- shiny::renderPlot({
      shiny::req(scan_plot())
      scan_plot()
    })
    # See also plotly clicked_data in `scanlyApp()`.
    output$plot_click <-  DT::renderDT({
      shiny::req(scan_plot(), scan_table_chr(), main_par$selected_chr(), input$plot_click)
      xvar <- "position"
      if (main_par$selected_chr() == "All") xvar <- "BPcum"
      out <- shiny::nearPoints(scan_table_chr(), input$plot_click,
        xvar = xvar, yvar = "LOD",
        threshold = 10, maxpoints = 1, addDist = TRUE)
      dplyr::mutate(out, dplyr::across(dplyr::where(is.numeric), \(x) signif(x, 4)))
    })
    # The `file_name()` is used in `downloadServer()` for plot and table file names.
    file_name <- shiny::reactive({
      instanceID <- shiny::req(main_par$which_trait())
      if(shiny::req(main_par$selected_chr()) != "All") {
        instanceID <- paste0(instanceID, "_chr", main_par$selected_chr())
      }
      paste("scan", instanceID, sep = "_")
    })
    # Return `scan_list` = reactiveValues containing elements `filename`, `tables` and `plots`.
    # The tables and plots are reactiveValues with reactives `scan_table_chr` and `scan_plot`.
    # Access `file_name()` as `scan_list$filename()` and `scan_plot()` as `scan_list$plots$scan()`.
    # View plot names as `names(scan_list$plots)`.
    shiny::reactiveValues(
      filename = file_name,
      tables = shiny::reactiveValues(
        scan = scan_table_chr),
      plots  = shiny::reactiveValues(
        scan = scan_plot)
    )
  })
}
#' @rdname scanApp
#' @export
scanUI <- function(id) {
  ns <- shiny::NS(id)
  DT::DTOutput(ns("plot_click"))
}
#' @rdname scanApp
#' @export
scanOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("scan_plot"))
}