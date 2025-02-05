#' Scan module
#'
#' @param id shiny identifier
#' @param file_directory data frame with file directory information
#' @param annotation_list list with annotation information
#' @param markers list object with marker information
#'
#' @importFrom DT DTOutput renderDT
#' @importFrom shiny actionButton moduleServer nearPoints NS plotOutput reactive renderPlot
#'             req setProgress shinyApp withProgress
#' @importFrom bslib card card_header page_sidebar sidebar
#' @importFrom shinycssloaders withSpinner
#' @importFrom dplyr across mutate where
#' @importFrom stringr str_split
#' @export
scanServer <- function(id, main_par, file_directory, annotation_list, markers) {
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
      shiny::req(main_par$which_trait, input$LOD_thr)
      scans <- shiny::withProgress(
        message = paste("scan of", chosen_trait(), "in progress"),
        value = 0, {
          shiny::setProgress(1)
          trait_scan(file_directory, main_par$selected_dataset, chosen_trait())
        })
      scan_plot <- QTL_plot_visualizer(scans, main_par$which_trait, input$LOD_thr, markers)
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
scanInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::sliderInput(ns("LOD_thr"),
    label = "LOD threshold for evaluation",
    min = 4,
    max = 20,
    value = 7.5,
    round = TRUE)
}
scanUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::actionButton(ns("scan"), "Show the LOD scan")
}
scanOutput <- function(id) {
    ns <- shiny::NS(id)
    bslib::card(
      bslib::card_header("LOD profile"),
      shiny::plotOutput(ns("scan_plot"), click = ns("plot_click")) |>
        shinycssloaders::withSpinner(color="#0dc5c1"),
      DT::DTOutput(ns("scan_points"))
    )
}
scanApp <- function() {
  source("qtlSetup.R")
  source("mainParServer.R")
  source("traitServer.R")
  ui <- bslib::page_sidebar(
    title = "Test Scan",
    sidebar = bslib::sidebar("side_panel",
      mainParInput("main_par"),
      mainParUI("main_par"),
      scanInput("scan"),
      scanUI("scan")
    ),
    scanOutput("scan")
  )
  server <- function(input, output, session) {
      main_par <- mainParServer("main_par", file_directory, annotation_list)
      scanServer("scan", main_par, file_directory, annotation_list, markers)
  }
  shiny::shinyApp(ui = ui, server = server)

}
