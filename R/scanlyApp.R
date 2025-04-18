#' Plotly Scan App Module
#'
#' @param id shiny identifier
#' @param import reactive list with file_directory and markers
#' @param main_par reactive list with selected_dataset, LOD_thr and which_trait
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
scanlyApp <- function() {
  ui <- bslib::page_sidebar(
    title = "Test Scan",
    sidebar = bslib::sidebar("side_panel",
      mainParInput("main_par"), # "selected_dataset", "LOD_thr"
      mainParUI("main_par")    # "which_trait"
    ),
    scanlyOutput("scanly")
  )
  server <- function(input, output, session) {
      import <- importServer("import")
      main_par <- mainParServer("main_par", import)
      scan_plot <- scanServer("scan", main_par, import)
      scanlyServer("scanly", main_par, scan_plot)
  }
  shiny::shinyApp(ui = ui, server = server)
}
#' @rdname scanApp
#' @export
scanlyServer <- function(id, main_par, scan_plot) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Plotly plot
    output$scan_plot <- shiny::renderUI({
      shiny::req(scan_plot())
      plotly::plotlyOutput(ns("render_plot"), click = ns("plot_click")) |>
        shinycssloaders::withSpinner(type = 8, color = "#3498db")
    })
    output$render_plot <- plotly::renderPlotly({
      shiny::req(scan_plot())
      scan_plot()
    })
    # Handle clicked points display
    observeEvent(event_data("plotly_click", source = "scan_plot"), {
      clicked_data(event_data("plotly_click", source = "scan_plot"))
    })
    # Create the interactive plotly plot
    output$scan_plot <- renderPlotly({
      req(scan_plot(), input$which_trait)
      ggplotly_qtl_scan()
    })
  })
}
#' @rdname scanApp
#' @export
scanlyOutput <- function(id) {
    ns <- shiny::NS(id)
    bslib::card(
      bslib::card_header("LOD profile"),
      shiny::uiOutput(ns("scan_plot")),
      shiny::h4("Clicked point info does not work at this time"),
      DT::DTOutput(ns("clicked_point_info"))
    )
}
