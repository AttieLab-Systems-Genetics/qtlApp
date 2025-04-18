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
      req(plot_base(), input$which_trait)
      
      # Get the base plot and data
      plot_result <- plot_base()
      p <- plot_result$p
      plot_data <- plot_result$data
      
      # Show the highest LOD peak
      peaks_info <- peak_finder(file_directory, input$selected_dataset, input$which_trait)
      if (nrow(peaks_info) > 0) {
          peaks_info <- peaks_info %>%
              arrange(desc(lod)) %>%
              slice(1)
          
          peak_point <- plot_data %>%
              filter(markers == peaks_info$marker)
          
          if (nrow(peak_point) > 0) {
              # Add the red diamond at the peak
              if (input$selected_chr == "All") {
                  p <- p + geom_point(data = peak_point,
                                     aes(x = BPcum, y = LOD),
                                     color = "red",
                                     size = 3,
                                     shape = 20)
              } else {
                  p <- p + geom_point(data = peak_point,
                                     aes(x = position, y = LOD),
                                     color = "red",
                                     size = 3,
                                     shape = 20)
              }
          }
      }
      
      # Create formatted trait text for plot title using the official gene symbol
      trait_text <- paste0("<b style='font-size: 24px;'>", official_gene_symbol(), "</b>")
      
      # Create subtitle with peak information
      subtitle <- if(nrow(peaks_info) > 0) {
          chr_label <- if(peaks_info$chr %in% c(20,21,22)) {
              c("X","Y","M")[peaks_info$chr-19]
          } else {
              peaks_info$chr
          }
          paste0(
              "<span style='font-size: 16px;'>",
              "<b>Peak Marker:</b> ", peaks_info$marker,
              " (Chr", chr_label, ":", round(peaks_info$pos, 2), " Mb) | ",
              "<b>LOD:</b> ", round(peaks_info$lod, 2),
              "</span>"
          )
      } else {
          "<span style='font-size: 16px; color: #7f8c8d;'>No significant peaks</span>"
      }

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
