#' Plotly Scan App Module
#' 
#'
#' @param id shiny identifier
#' @param main_par reactive list with selected_dataset, LOD_thr and which_trait
#' @param scan_list reactive list from scanServer
#' @param peak_list reactive list from peakServer
#'
#' @importFrom DT datatable DTOutput renderDT
#' @importFrom shiny actionButton h4 moduleServer nearPoints NS plotOutput
#'             reactive reactiveVal renderPlot renderUI req setProgress shinyApp
#'             uiOutput withProgress
#' @importFrom bslib card card_header page_sidebar sidebar
#' @importFrom shinycssloaders withSpinner
#' @importFrom dplyr across mutate where
#' @importFrom stringr str_split
#' @importFrom htmltools tags
#' @importFrom plotly event_data plotlyOutput renderPlotly
#'
#' @export
scanlyApp <- function() {
  ui <- bslib::page_sidebar(
    title = "Test Plotly Scan",
    sidebar = bslib::sidebar("side_panel",
      mainParInput("main_par"), # "selected_dataset", "LOD_thr"
      mainParUI("main_par")),   # "which_trait","selected_chr"
    bslib::card(
      bslib::card_header("LOD profile"),
      scanlyOutput("scanly")),
    bslib::card(
      bslib::card_header("Clicked Peak"),
      scanlyUI("scanly"))
  )
  server <- function(input, output, session) {
      import <- importServer("import")
      main_par <- mainParServer("main_par", import)
      scan_list <- scanServer("scan_list", main_par, import)
      peak_list <- peakServer("peak_list", main_par, import)
      scanlyServer("scanly", main_par, scan_list, peak_list)
  }
  shiny::shinyApp(ui = ui, server = server)
}
#' @rdname scanApp
#' @export
scanlyServer <- function(id, main_par, scan_list, peak_list) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Plotly plot
    output$scan_plot <- shiny::renderUI({
      plotly::plotlyOutput(ns("render_plot")) |>
        shinycssloaders::withSpinner(type = 8, color = "#3498db")
    })
    # Extract `scan_table` and `scan_plot` from `scan_list`.
    scan_table <- shiny::reactive({
      shiny::req(scan_list$tables$scan())
      scan_list$tables$scan()
    })
    scan_plot <- shiny::reactive({
      shiny::req(scan_list$plots$scan())
      scan_list$plots$scan()
    })
    peak_table <- shiny::reactive({
      shiny::req(peak_list$tables$peak())
      peak_list$tables$peak()
    })
    scanly_plot <- shiny::reactive({

      shiny::req(main_par$selected_chr())
      scan_object <- list(plot = shiny::req(scan_plot()), table = shiny::req(scan_table()))
      ggplotly_qtl_scan(scan_object, peak_table(), main_par$selected_chr(), "scanly_plot")
    })
    output$render_plot <- plotly::renderPlotly({
      shiny::req(scanly_plot())
      scanly_plot()
    })
    # Peak information based on plotly click.
    # `stable_peak()` stores the stable peak information for table display.
    # `stable_peak()` is not reactive but is updated by reactives.
    # `stable_peak()` changes when `max_peak()` or `plotly_click()` changes.
    # `max_peak()` has maximum peak when reactives change.
    # `plotly_click()` has click event data from plotly.
    # ** want stable_peak() to include existing peak table and clicked point **
    # ** want to show table at start and when selected_chr change **
    stable_peak <- shiny::reactiveVal(NULL)
    plotly_click <- shiny::reactive({
      plotly::event_data("plotly_click", source = "scanly_plot")
    })
    shiny::observeEvent(shiny::req(plotly_click()), {
      shiny::req(peak_table(), scan_table(), which_peak())
      out <- peak_info(peak_table(), scan_table(), which_peak(), plotly_click())
      stable_peak(out)
      out
    })
    which_peak <- shiny::reactive({
      shiny::req(peak_table(), main_par$LOD_thr)
      ordered_peaks <- highest_peaks(peak_table(), main_par$LOD_thr)
      ordered_peaks$marker[1]
    })
    max_peak <- shiny::reactive({
      shiny::req(peak_table(), scan_table(), which_peak(), main_par$selected_chr)
      out <- peak_info(peak_table(), scan_table(), which_peak())
      stable_peak(out)
      out
    })
    output$clicked_point_info <- DT::renderDT({
      shiny::req(stable_peak(), scan_plot(), scan_table(), which_peak(), plotly_click(),
        main_par$selected_chr)
      DT::datatable(
        stable_peak(),
        options = list(dom = 't', ordering = FALSE, pageLength = 1),
        rownames = FALSE,
        class = 'compact hover',
        caption = htmltools::tags$caption(
          style = 'caption-side: top; text-align: left; color: #2c3e50; font-weight: bold; font-size: 14px;',
          'Peak Information and Strain Effects'))
    })
    output$high_peak <- shiny::renderUI({
      list(
        shiny::renderPrint({
          if(shiny::isTruthy(plotly_click())) {
            paste("Clicked point: x =", plotly_click()$x, ", y =", plotly_click()$y)
          } else {
            "No point clicked yet."
          }
        }),
        DT::DTOutput(ns("clicked_point_info")))
    })
    shiny::observeEvent(
      plotly::event_data("plotly_doubleclick", source = "scanly_plot"), {
      if(main_par$selected_chr() != "All") {
        shiny::updateSelectInput(session, "selected_chr", selected = "All")
      }
    })
  })
}
#' @rdname scanApp
#' @export
scanlyUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("high_peak"))
}
#' @rdname scanApp
#' @export
scanlyOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("scan_plot"))
}
