#' Plotly Scan App Module
#' 
#' This app has a point-click feature not yet implemented.
#' See `kalynn_R/latest_app_kalynn/app.R`.
#' Also not sure if `selected_chr` should be here or elsewhere.
#'
#' @param id shiny identifier
#' @param main_par reactive list with selected_dataset, LOD_thr and which_trait
#' @param scan_table reactive dataframe from scanServer
#' @param peak_table reactive dataframe from peakServer
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
      scanlyInput("scanly")     # "selected_chr"
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
scanlyServer <- function(id, main_par, scan_table, peak_table) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Plotly plot
    output$scan_plot <- shiny::renderUI({
      shiny::req(scan_plot())
      plotly::plotlyOutput(ns("render_plot")) |>
        shinycssloaders::withSpinner(type = 8, color = "#3498db")
    })
    scanly_plot <- shiny::reactive({
      req(scan_table(), peak_table(), input$selected_chr)
      # Check if the selected chromosome is "All"
      if (input$selected_chr == "All") {
        # Create a plot for all chromosomes
        ggplotly_qtl_scan(scan_table(), peak_table(), source = "scanly_plot")
      } else {
        # Create a plot for the selected chromosome
        ggplotly_qtl_scan(scan_table(), peak_table(), input$selected_chr, source = "scanly_plot")
      }
    })
    output$render_plot <- plotly::renderPlotly({
      req(scanly_plot())
    })
    # Handle clicked points display ** not quite right yet **
    click_data <- shiny::reactive({
      # Want to use the plotly click event data if present.
      # But when peak_table() changes, reset peak info
      shiny::req(peak_table(), main_par$LOD_thr)
      ordered_peaks <- highest_peaks(peak_table(), main_par$LOD_thr)
      # See Kalynn's app for how to get the peak info
      # and display it in the plotly click event.
      # ** look at output$clicked_point_info in app.R. **
      plotly::event_data("plotly_click", source = "scanly_plot")
    })
    shiny::observeEvent(click_data(), {
      if (!is.null(click_data)) {
        # Process the click data
          selected_x <- click_data()$x
          selected_y <- click_data()$y
          # Perform actions based on the clicked data
      }
    })
    output$click_info <- renderPrint({
      shiny::req(click_data())
      if (!is.null(click_data())) {
          paste("Clicked point: x =", click_data()$x, ", y =", click_data()$y)
      } else {
        "No point clicked yet."
      }
    })
    # Add observer for plotly double click event
    shiny::observeEvent(plotly::event_data("plotly_doubleclick", source = "scan_plot"), {
      if(input$selected_chr != "All") {
        updateSelectInput(session, "selected_chr", selected = "All")
      }
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
