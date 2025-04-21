#' Plotly Scan App Module
#' 
#' This app has a point-click feature not yet implemented.
#' See `kalynn_R/latest_app_kalynn/app.R`.
#' Also not sure if `selected_chr` should be here or elsewhere.
#'
#' @param id shiny identifier
#' @param main_par reactive list with selected_dataset, LOD_thr and which_trait
#' @param scan_table reactive list with scan table
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
      mainParUI("main_par"),    # "which_trait"
      scanlyInput("scanly").    # "selected_chr"
    ),
    scanlyOutput("scanly")
  )
  server <- function(input, output, session) {
      import <- importServer("import")
      main_par <- mainParServer("main_par", import)
      scan_table <- scanServer("scan_table", main_par, import)
      peak_table <- peakServer("peak_table", main_par, import)
      scanlyServer("scanly", main_par, scan_table, peak_table)
  }
  shiny::shinyApp(ui = ui, server = server)
}
#' @rdname scanApp
#' @export
scanlyServer <- function(id, main_par, scan_table) {
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
    output$scan_plot <- renderPlotly({
      req(scan_table(), peak_table(), input$selected_chr)
      ggplotly_qtl_scan(scan_table(), peak_table(), input$selected_chr)
    })
    output$plot_click <- DT::renderDT({
      shiny::req(scan_table(), input$plot_click)
      out <- shiny::nearPoints(scan_table(), input$plot_click,
                               xvar = "BPcum", yvar = "LOD",
                               threshold = 10, maxpoints = 1, addDist = TRUE)
      dplyr::mutate(out,
        dplyr::across(dplyr::where(is.numeric), \(x) signif(x, 4)))
    })
  })
}
#' @rdname scanApp
#' @export
scanlyInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::selectInput(ns("selected_chr"), "Zoom to Chromosome:",
    choices = c("All",
                "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                "11", "12", "13", "14", "15", "16", "17", "18", "19",
                "X", "Y", "M"),
    selected = "All",
    width = "150px")
}
#' @rdname scanApp
#' @export
scanlyOutput <- function(id) {
    ns <- shiny::NS(id)
    bslib::card(
      bslib::card_header("LOD profile"),
      shiny::uiOutput(ns("scan_plot")),
      DT::DTOutput(ns("clicked_point_info"))
    )
}
