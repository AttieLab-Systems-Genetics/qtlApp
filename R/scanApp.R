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
      downloadInput("download"),   # downloadButton, filename
      downloadOutput("download")), # plot_table, inputs for Plots or Tables
    bslib::card(
      bslib::card_header("LOD profile"),
      scanOutput("scan_list")),
    bslib::card(
      bslib::card_header("Clicked Peak"),
      scanUI("scan_list"))
  )
  server <- function(input, output, session) {
      import <- importServer("import")
      main_par <- mainParServer("main_par", import)
      scan_list <- scanServer("scan_list", main_par, import)
      peak_list <- peakServer("peak_list", main_par, import)
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
    scan_table_chr <- shiny::reactive({
      # Evaluate reactive before req or use
      sel_chr <- main_par$selected_chr()
      shiny::req(scan_table(), sel_chr)
      if (sel_chr == "All") { 
        scan_table()
      } else {
        # Ensure chr column exists before filtering
        st <- scan_table()
        if("chr" %in% colnames(st)){
            dplyr::filter(st, .data$chr == sel_chr)
        } else {
            warning("scan_table_chr: 'chr' column not found in scan_table()")
            st # Return unfiltered table if 'chr' missing
        }
      }
    })
    scan_plot <- shiny::reactive({
      # Evaluate reactive before req or use
      sel_chr <- main_par$selected_chr()
      lod_thr <- main_par$LOD_thr()
      shiny::req(scan_table_chr(), lod_thr, sel_chr)
      ggplot_qtl_scan(scan_table_chr(), lod_thr, sel_chr)
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
      # Evaluate reactive before req or use
      sel_chr <- main_par$selected_chr()
      shiny::req(scan_plot(), scan_table_chr(), sel_chr, input$plot_click)
      xvar <- "position"
      if (sel_chr == "All") xvar <- "BPcum" 
      out <- shiny::nearPoints(scan_table_chr(), input$plot_click,
        xvar = xvar, yvar = "LOD",
        threshold = 10, maxpoints = 1, addDist = TRUE)
      dplyr::mutate(out, dplyr::across(dplyr::where(is.numeric), \(x) signif(x, 4)))
    })
    # The `file_name()` is used in `downloadServer()` for plot and table file names.
    file_name <- shiny::reactive({
      instanceID <- shiny::req(main_par$which_trait())
      # Evaluate selected_chr reactive before comparing
      if(shiny::req(main_par$selected_chr()) != "All") { 
        instanceID <- paste0(instanceID, "_chr", main_par$selected_chr())
      }
      paste("scan", instanceID, sep = "_")
    })
    # Return `scan_list` = reactiveValues containing elements `filename`, `tables` and `plots`.
    # The tables and plots are reactiveValues with reactives `scan_table_chr` and `scan_plot`.
    # Access `file_name()` as `scan_list$filename()` and `scan_plot()` as `scan_list$plots$scan()`.
    # View plot names as `names(scan_list$plots)`.
    
    # --- MODIFICATION: Return a REGULAR list containing reactives --- 
    return_value <- list(
      filename = file_name,           # The filename reactive
      scan_table = scan_table_chr,    # The scan table reactive
      scan_plot = scan_plot          # The scan plot reactive
    )
    message("--- scanServer returning object of class: list ---")
    return(return_value)
    # --------------------------------------------------------
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
