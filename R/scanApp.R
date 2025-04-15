#' Scan App Module
#'
#' @param id shiny identifier
#' @param import reactive list with file_directory and markers
#' @param main_par reactive list with selected_dataset, LOD_thr and which_trait
#'
#' @importFrom DT DTOutput renderDT
#' @importFrom shiny actionButton moduleServer nearPoints NS plotOutput reactive renderPlot
#'             req setProgress shinyApp withProgress
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
      mainParUI("main_par"),    # "which_trait"
      scanUI("scan")            # "scan" actionButton
    ),
    scanOutput("scan")
  )
  server <- function(input, output, session) {
      import <- importServer("import")
      main_par <- mainParServer("main_par", import)
      scanServer("scan", main_par, import)
  }
  shiny::shinyApp(ui = ui, server = server)

}
#' @rdname scanApp
#' @export
scanServer <- function(id, main_par, import) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    chosen_trait <- shiny::reactive({
      shiny::req(main_par$which_trait)
      stringr::str_split(
        stringr::str_split(main_par$which_trait, pattern=" [(]")[[1]][2],
        pattern="[)]")[[1]][1]
    })
        
    # create the scans only when `scan` button clicked
    shiny::observeEvent(input$scan, {
      shiny::req(main_par$selected_dataset, chosen_trait())
      shiny::req(main_par$which_trait, main_par$LOD_thr)
      scans <- shiny::withProgress(
        message = paste("scan of", chosen_trait(), "in progress"),
        value = 0, {
          shiny::setProgress(1)
          trait_scan(import()$file_directory, main_par$selected_dataset, chosen_trait())
        })
      # ** check this code **
      scan_plot <- ggplot_qtl_scan(
        QTL_plot_visualizer(scans, main_par$which_trait, main_par$LOD_thr, import()$markers))
      output$scan_plot <- shiny::renderPlot({
        scan_plot[[1]]
      })
      output$scan_points <-  DT::renderDT({
        shiny::req(input$plot_click)
        out <- shiny::nearPoints(scan_plot[[2]], input$plot_click, xvar = "BPcum",
          yvar = "LOD", threshold = 10, maxpoints = 1, addDist = TRUE)
        dplyr::mutate(out, dplyr::across(dplyr::where(is.numeric), signif, 4))
      })
    })
  })
}
#' @rdname scanApp
#' @export
scanUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::actionButton(ns("scan"), "Show the LOD scan")
}
#' @rdname scanApp
#' @export
scanOutput <- function(id) {
    ns <- shiny::NS(id)
    bslib::card(
      bslib::card_header("LOD profile"),
      shiny::plotOutput(ns("scan_plot"), click = ns("plot_click")) |>
        shinycssloaders::withSpinner(color="#0dc5c1"),
      DT::DTOutput(ns("scan_points"))
    )
}
