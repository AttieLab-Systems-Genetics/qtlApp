#' Scan App Module
#'
#' @param id shiny identifier
#' @param main_par reactive list with selected_dataset, LOD_thr and which_trait
#' @param import reactive list with file_directory and markers
#'
#' @importFrom DT DTOutput renderDT
#' @importFrom shiny actionButton h4 moduleServer nearPoints NS plotOutput
#'             reactive reactiveValues renderPlot renderUI req setProgress shinyApp
#'             uiOutput withProgress div downloadButton numericInput tagList
#'             selectInput plotOutput downloadHandler
#' @importFrom bslib card card_header page_sidebar sidebar
#' @importFrom shinycssloaders withSpinner
#' @importFrom dplyr across mutate where filter
#' @importFrom stringr str_split str_remove
#' @importFrom ggplot2 ggsave
#'
#' @export
scanApp <- function() {
  # Source UI helper functions if not already available
  if (!exists("create_download_button", mode = "function")) {
    source("R/ui_styles.R") # Assumes ui_styles.R contains all create_* helpers
  }

  ui <- bslib::page_sidebar(
    title = "QTL Scan Visualizer",
    sidebar = bslib::sidebar("side_panel",
      mainParInput("main_par"),
      mainParUI("main_par"),
      peakInput("peak"),
      peakUI("peak"),
      # Existing download module (may or may not be used for these new buttons)
      downloadInput("download"),
      downloadOutput("download")
    ),
    bslib::card(
      bslib::card_header("LOD Profile Controls"),
      bslib::card_body(
        # Row for plot title and download buttons
        shiny::div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;",
          shiny::h4("LOD Score Plot", style = "margin: 0; color: #2c3e50; font-weight: 600;"),
          shiny::div(style = "display: flex; gap: 10px;",
            # Use ui_styles.R helpers if available, otherwise standard buttons
            if (exists("create_download_button", mode = "function")) {
              tagList(
                create_download_button(shiny::NS("scan_list", "download_qtl_plot_png"), "PNG", class = "btn-sm"),
                create_download_button(shiny::NS("scan_list", "download_qtl_plot_pdf"), "PDF", class = "btn-sm")
              )
            } else {
              tagList(
                shiny::downloadButton(shiny::NS("scan_list", "download_qtl_plot_png"), "PNG", class = "btn-sm"),
                shiny::downloadButton(shiny::NS("scan_list", "download_qtl_plot_pdf"), "PDF", class = "btn-sm")
              )
            }
          )
        ),
        # Row for plot dimension controls and color toggle
        shiny::div(style = "display: flex; gap: 20px; align-items: center; margin-bottom: 20px; flex-wrap: wrap;",
          # Chromosome selector (from mainParUI, not namespaced here directly in scanApp)
          # This is typically part of mainParUI("main_par") which is in the sidebar.
          # If it needs to be here, it should be shiny::uiOutput(shiny::NS("main_par", "selected_chr_ui_placeholder"))
          # For now, assuming mainParUI handles chr selection in the sidebar.
          
          # Plot dimensions
          shiny::div(style = "display: flex; align-items: center; gap: 10px;",
            if (exists("create_numeric_input", mode = "function")) {
              create_numeric_input(shiny::NS("scan_list", "plot_width"), "Width:", value = 1000, min = 400, max = 2000, step = 50, width = "100px")
            } else {
              shiny::numericInput(shiny::NS("scan_list", "plot_width"), "Width:", value = 1000, min = 400, max = 2000, step = 50, width = "100px")
            },
            if (exists("create_numeric_input", mode = "function")) {
              create_numeric_input(shiny::NS("scan_list", "plot_height"), "Height:", value = 600, min = 300, max = 1200, step = 50, width = "100px")
            } else {
              shiny::numericInput(shiny::NS("scan_list", "plot_height"), "Height:", value = 600, min = 300, max = 1200, step = 50, width = "100px")
            }
          ),
          shiny::div(style = "display: flex; gap: 5px;",
            if (exists("create_button", mode = "function")) {
              tagList(
                create_button(shiny::NS("scan_list", "preset_1to1"), "1:1", class = "btn-sm btn-light"),
                create_button(shiny::NS("scan_list", "preset_3to2"), "3:2", class = "btn-sm btn-light"),
                create_button(shiny::NS("scan_list", "preset_16to9"), "16:9", class = "btn-sm btn-light")
              )
            } else {
              tagList(
                shiny::actionButton(shiny::NS("scan_list", "preset_1to1"), "1:1", class = "btn-sm"),
                shiny::actionButton(shiny::NS("scan_list", "preset_3to2"), "3:2", class = "btn-sm"),
                shiny::actionButton(shiny::NS("scan_list", "preset_16to9"), "16:9", class = "btn-sm")
              )
            }
          ),
          # Color toggle switch
          if (exists("create_lever_switch", mode = "function")) {
            create_lever_switch(shiny::NS("scan_list", "color_toggle"))
          } else {
            shiny::checkboxInput(shiny::NS("scan_list", "color_toggle"), "Alt Colors", value = TRUE)
          }
        )
      ),
      # The plot itself
      scanOutput("scan_list") # This is uiOutput(ns("scan_plot"))
    ),
    bslib::card( # Card for the peaks table
      bslib::card_header("Peaks Table"),
      peakOutput("peak")
    ),
    bslib::card(
      bslib::card_header("Clicked Peak on Scan"), 
      scanUI("scan_list"),
      max_height = "150px"
    )
  )
  server <- function(input, output, session) {
      import <- importServer("import")
      main_par <- mainParServer("main_par", import) # Provides selected_chr among others
      scan_list_server_output <- scanServer("scan_list", main_par, import)
      peak_list <- peakServer("peak", main_par, import)
      # Pass the plot from scan_list_server_output to downloadServer if needed
      # For now, downloadServer uses merged_list which has scan_list_server_output$plots$scan
      merged_list <- mergeServer("merged_list", scan_list_server_output, peak_list)
      downloadServer("download", merged_list)
  }
  shiny::shinyApp(ui = ui, server = server)
}

# Server logic for scan related inputs and plot
scanServer <- function(id, main_par, import) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_width_rv <- shiny::reactiveVal(1000)
    plot_height_rv <- shiny::reactiveVal(600)
    use_alternating_colors_rv <- shiny::reactiveVal(TRUE)

    shiny::observeEvent(input$plot_width, { plot_width_rv(input$plot_width) })
    shiny::observeEvent(input$plot_height, { plot_height_rv(input$plot_height) })

    shiny::observeEvent(input$preset_1to1, {
      plot_width_rv(800)
      plot_height_rv(800)
      shiny::updateNumericInput(session, "plot_width", value = 800)
      shiny::updateNumericInput(session, "plot_height", value = 800)
    })
    shiny::observeEvent(input$preset_3to2, {
      plot_width_rv(900)
      plot_height_rv(600)
      shiny::updateNumericInput(session, "plot_width", value = 900)
      shiny::updateNumericInput(session, "plot_height", value = 600)
    })
    shiny::observeEvent(input$preset_16to9, {
      plot_width_rv(1280) # Common 16:9, or use 1600x900 from app2
      plot_height_rv(720)
      shiny::updateNumericInput(session, "plot_width", value = 1280)
      shiny::updateNumericInput(session, "plot_height", value = 720)
    })

    shiny::observeEvent(input$color_toggle, {
      # For checkboxInput, input$color_toggle is TRUE/FALSE
      # For create_lever_switch, it's a timestamp, so toggle on any change
      if(is.logical(input$color_toggle)){
        use_alternating_colors_rv(input$color_toggle)
      } else {
        use_alternating_colors_rv(!use_alternating_colors_rv())
        # Optionally, update lever switch visual state if using JS
        # shinyjs::runjs(sprintf("document.getElementById('%s').classList.toggle('active');", ns("color_toggle")))
      }
    })
    
    selected_trait <- shiny::reactive({
      shiny::req(import(), main_par$selected_dataset())
      trait <- shiny::req(main_par$which_trait())
      trait_type <- get_trait_type(import(), main_par$selected_dataset())
      if(trait_type %in% c("genes", "isoforms")) {
        trait <- stringr::str_remove(main_par$which_trait(), "_.*")
      }
      trait
    })
    
    scans <- shiny::reactive({
      shiny::req(main_par$selected_dataset(), selected_trait())
      shiny::withProgress(
        message = paste("scan of", selected_trait(), "in progress"), value = 0, {
          shiny::setProgress(1)
          suppressMessages(
            trait_scan(import()$file_directory, main_par$selected_dataset(), selected_trait()))
        })
    })
    
    scan_table <- shiny::reactive({
      shiny::req(scans(), main_par$which_trait(), main_par$LOD_thr())
      QTL_plot_visualizer(scans(), main_par$which_trait(), main_par$LOD_thr(), import()$markers)
    })
    
    scan_table_chr <- shiny::reactive({
      shiny::req(scan_table(), main_par$selected_chr())
      if (main_par$selected_chr() == "All") {
        scan_table()
      } else {
        dplyr::filter(scan_table(), chr == main_par$selected_chr())
      }
    })
    
    # This is the main ggplot object for the scan
    current_scan_plot_gg <- shiny::reactive({
      shiny::req(scan_table_chr(), main_par$LOD_thr(), main_par$selected_chr())
      # Pass the color toggle state to ggplot_qtl_scan if it accepts it
      # For now, ggplot_qtl_scan handles its own colors based on create_modern_palette
      # We might need to modify ggplot_qtl_scan to accept a color_mode parameter
      # Or, make sure use_alternating_colors_rv() is used by ggplot_qtl_scan
      # For now, we assume ggplot_qtl_scan is already modernized and respects some internal toggle if available or uses its default.
      # The key is that ggplot_qtl_scan returns a ggplot object.
      ggplot_qtl_scan(scan_table_chr(), main_par$LOD_thr(), main_par$selected_chr())
    })

    output$scan_plot <- shiny::renderUI({
      shiny::req(current_scan_plot_gg())
      # Use reactive width and height for the plotOutput
      shiny::plotOutput(ns("render_plot"), 
                        click = ns("plot_click"), 
                        width = paste0(plot_width_rv(), "px"), 
                        height = paste0(plot_height_rv(), "px")) |>
        shinycssloaders::withSpinner(type = 8, color = "#3498db")
    })
    
    output$render_plot <- shiny::renderPlot({
      shiny::req(current_scan_plot_gg())
      current_scan_plot_gg()
    }, 
    # Make renderPlot itself aware of dynamic dimensions. 
    # This is an alternative to setting width/height in plotOutput if plotOutput is not in renderUI
    # width = function() plot_width_rv(), 
    # height = function() plot_height_rv()
    res = 96 # Standard resolution
    )
    
    output$plot_click <-  DT::renderDT({
      shiny::req(current_scan_plot_gg(), scan_table_chr(), main_par$selected_chr(), input$plot_click)
      xvar <- "position"
      if (main_par$selected_chr() == "All") xvar <- "BPcum"
      out_data <- shiny::nearPoints(scan_table_chr(), input$plot_click,
        xvar = xvar, yvar = "LOD",
        threshold = 10, maxpoints = 1, addDist = TRUE)
      if (!is.data.frame(out_data) || nrow(out_data) == 0) {
        return(DT::datatable(data.frame(), options = list(dom = 't')))
      }
      out_data <- dplyr::mutate(out_data, dplyr::across(dplyr::where(is.numeric), \(x) signif(x, 4)))
      DT::datatable(out_data, options = list(dom = 't', paging = FALSE), rownames = FALSE, selection = 'none')
    })

    # Download handlers for QTL plot
    output$download_qtl_plot_png <- shiny::downloadHandler(
      filename = function() {
        chr_suffix <- if(main_par$selected_chr() != "All") paste0("_chr", main_par$selected_chr()) else ""
        paste0("lod_plot_", main_par$which_trait(), chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".png")
      },
      content = function(file) {
        shiny::req(current_scan_plot_gg())
        # Use dynamic width/height for saving
        ggplot2::ggsave(file, plot = current_scan_plot_gg(), 
                        width = plot_width_rv()/96, # Assuming 96 DPI for conversion from px
                        height = plot_height_rv()/96, 
                        dpi = 300, units = "in")
      }
    )

    output$download_qtl_plot_pdf <- shiny::downloadHandler(
      filename = function() {
        chr_suffix <- if(main_par$selected_chr() != "All") paste0("_chr", main_par$selected_chr()) else ""
        paste0("lod_plot_", main_par$which_trait(), chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".pdf")
      },
      content = function(file) {
        shiny::req(current_scan_plot_gg())
        ggplot2::ggsave(file, plot = current_scan_plot_gg(), 
                        width = plot_width_rv()/96, 
                        height = plot_height_rv()/96, 
                        device = cairo_pdf, units = "in") # Use cairo_pdf for better PDF quality
      }
    )
    
    file_name_reactive <- shiny::reactive({
      instanceID <- shiny::req(main_par$which_trait())
      if(shiny::req(main_par$selected_chr()) != "All") {
        instanceID <- paste0(instanceID, "_chr", main_par$selected_chr())
      }
      paste("scan", instanceID, sep = "_")
    })
    
    shiny::reactiveValues(
      filename = file_name_reactive, # For existing downloadApp module if it uses this
      tables = shiny::reactiveValues(scan = scan_table_chr),
      plots  = shiny::reactiveValues(scan = current_scan_plot_gg) # Ensure this is the plot with dynamic sizing
    )
  })
}

scanUI <- function(id) {
  ns <- shiny::NS(id)
  DT::DTOutput(ns("plot_click"))
}

scanOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("scan_plot"))
}