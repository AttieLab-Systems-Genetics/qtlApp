#' Scan App Module
#'
#' @param id shiny identifier
#' @param main_par reactive list with selected_dataset, LOD_thr and which_trait
#' @param import reactive list with file_directory and markers
#'
#' @importFrom DT DTOutput renderDT
#' @importFrom shiny actionButton h4 moduleServer nearPoints NS plotOutput
#'             reactive renderPlot renderUI req setProgress shinyApp
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
      mainParInput("main_par"), # "selected_dataset", "LOD_thr"
      mainParUI("main_par")    # "which_trait"
    ),
    scanOutput("scan_table")
  )
  server <- function(input, output, session) {
      import <- importServer("import")
      main_par <- mainParServer("main_par", import)
      scanServer("scan_table", main_par, import)
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
      shiny::req(import(), main_par$selected_dataset)
      trait <- shiny::req(main_par$which_trait)
      # if gene or isoform then trait is `symbol_id``.
      trait_type <- get_trait_type(import(), main_par$selected_dataset)
      if(trait_type %in% c("genes", "isoforms")) {
        trait <- stringr::str_remove(main_par$which_trait, "_.*")
      }
      trait
    })
    scans <- shiny::reactive({
      shiny::req(main_par$selected_dataset, selected_trait())
      scans <- shiny::withProgress(
        message = paste("scan of", selected_trait(), "in progress"),
        value = 0, {
          shiny::setProgress(1)
          suppressMessages(
            trait_scan(import()$file_directory,
              main_par$selected_dataset, selected_trait()))
        }
      )
    })
    scan_table <- shiny::reactive({
      shiny::req(scans(), main_par$which_trait, main_par$LOD_thr)
      QTL_plot_visualizer(
        scans(), main_par$which_trait, main_par$LOD_thr, import()$markers)
    })
    scan_plot <- shiny::reactive({
      shiny::req(scan_table())
      ggplot_qtl_scan(scan_table())
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
      shiny::req(scan_table(), input$plot_click)
      out <- shiny::nearPoints(scan_table(), input$plot_click,
        xvar = "BPcum", yvar = "LOD",
        threshold = 10, maxpoints = 1, addDist = TRUE)
      dplyr::mutate(out, dplyr::across(dplyr::where(is.numeric), \(x) signif(x, 4)))
    })
    # Return
    scan_table
  })
}
#' @rdname scanApp
#' @export
scanOutput <- function(id) {
    ns <- shiny::NS(id)
    bslib::card(
      bslib::card_header("LOD profile"),
      shiny::uiOutput(ns("scan_plot")),
      DT::DTOutput(ns("plot_click"))
    )
}
