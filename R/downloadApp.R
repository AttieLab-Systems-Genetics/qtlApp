#' Download App Module
#'
#' @param id shiny identifier
#' @param import reactive list with file_directory and annotation_list
#' @param trait_type character string as reactive object
#'
#' @importFrom DT DTOutput renderDT
#' @importFrom shiny moduleServer NS reactive renderText req selectInput shinyApp
#'             textOutput
#' @importFrom bslib page
#' 
#' @export
downloadApp <- function() {
  ui <- bslib::page_sidebar(
    title = "Test Plotly Scan",
    sidebar = bslib::sidebar("side_panel",
      mainParInput("main_par"),   # "selected_dataset", "LOD_thr"
      mainParUI("main_par"),      # "which_trait"
      scanlyInput("scanly"),      # "selected_chr" ** this needs to move to mainParApp
      downloadInput("download")), # "download"
    bslib::card(
      bslib::card_header("LOD profile"),
      scanOutput("scan_table"))
  )
  server <- function(input, output, session) {
    import <- importServer("import")
    main_par <- mainParServer("main_par", import)
    scan_table <- scanServer("scan_table", main_par, import)
    peak_table <- peakServer("peak_table", main_par, import)
    downloadServer("download", main_par, peak_table, scan_table)
  }
  shiny::shinyApp(ui = ui, server = server)
}
#' @rdname downloadApp
#' @export
downloadServer <- function(id, main_par, peak_table, scan_table, plot_result) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    output$downloadInput <- shiny::renderUI({
      list(
        shiny::downloadButton("download_qtl_plot", "Download QTL Plot"),
        shiny::selectInput(ns("format"), "Format", choices = c("PNG", "PDF"))
      )
    })
    output$download_qtl_plot <- shiny::downloadHandler(
      filename = function() {
        shiny::req(input$format)
        # Include chromosome info in filename if specific chromosome is selected
        chr_suffix <- if(main_par$selected_chr != "All") paste0("_chr", main_par$selected_chr) else ""
        paste0("lod_plot_", main_par$which_trait, chr_suffix, "_", format(Sys.time(), "%Y%m%d"),
          ".", input$format)
      },
      content = function(file) {
        # Get the base plot - this already has the chromosome filtering applied
        # Get peak information for the current trait
        peaks_info <- dplyr::filter(peak_table(), .data$lod >= main_par$LOD_thr) |>
          dplyr::arrange(dplyr::desc(.data$lod))
        # Filter peaks to only show those in the selected chromosome if applicable
        if(main_par$selected_chr != "All") {
          chr_num <- chr_XYM(main_par$selected_chr)
          peaks_info <- dplyr::filter(peaks_info, .data$chr == chr_num)
        }
        # Get the highest peak (if any)
        if(nrow(peaks_info) > 0) {
          peaks_info <- peaks_info %>% slice(1)
          # Find the peak point in the filtered data
          peak_point <- dplyr::filter(scan_table, .data$markers == peaks_info$marker)
          if (nrow(peak_point) > 0) {
            # Add the red diamond at the peak
            if (main_par$selected_chr == "All") {
              xvar <- "BPcum"
            } else {
              xvar <- "position"
            }
            plot_result <- plot_result + 
              ggplot2::geom_point(data = peak_point,
                ggplot2::aes_string(x = xvar, y = "LOD"),
                color = "red", size = 3, shape = 18)
          }
        }
      }
    )
  })
}
#' @rdname downloadApp
#' @export
downloadInput <- function(id, main_par, import) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("downloadInput"))
}
