#' Scan App Module
#'
#' @param id shiny identifier
#' @param main_par reactive list with selected_dataset, LOD_thr and which_trait
#' @param import reactive list with file_directory and markers
#'
#' @importFrom DT DTOutput renderDT datatable
#' @importFrom shiny actionButton h4 moduleServer NS plotOutput reactive reactiveVal
#'             reactiveValues renderPlot renderUI req setProgress shinyApp selectInput
#'             uiOutput withProgress div downloadButton numericInput tagList
#'             observeEvent updateNumericInput downloadHandler
#' @importFrom plotly plotlyOutput renderPlotly ggplotly event_data layout config
#' @importFrom bslib card card_header page_sidebar sidebar
#' @importFrom shinycssloaders withSpinner
#' @importFrom dplyr across mutate where filter select
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
              create_numeric_input(shiny::NS("scan_list", "plot_width"), "Width:", value = 1200, min = 400, max = 2000, step = 50, width = "100px")
            } else {
              shiny::numericInput(shiny::NS("scan_list", "plot_width"), "Width:", value = 1200, min = 400, max = 2000, step = 50, width = "100px")
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

    plot_width_rv <- shiny::reactiveVal(1200)
    plot_height_rv <- shiny::reactiveVal(600)
    use_alternating_colors_rv <- shiny::reactiveVal(TRUE)
    clicked_plotly_point_details_rv <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$plot_width, { plot_width_rv(input$plot_width) }, ignoreNULL = TRUE)
    shiny::observeEvent(input$plot_height, { plot_height_rv(input$plot_height) }, ignoreNULL = TRUE)

    shiny::observeEvent(input$preset_1to1, {
      shiny::updateNumericInput(session, "plot_width", value = 800)
      shiny::updateNumericInput(session, "plot_height", value = 800)
    })
    shiny::observeEvent(input$preset_3to2, {
      shiny::updateNumericInput(session, "plot_width", value = 900)
      shiny::updateNumericInput(session, "plot_height", value = 600)
    })
    shiny::observeEvent(input$preset_16to9, {
      shiny::updateNumericInput(session, "plot_width", value = 1280)
      shiny::updateNumericInput(session, "plot_height", value = 720)
    })

    shiny::observeEvent(input$color_toggle, {
      if(is.logical(input$color_toggle)){
        use_alternating_colors_rv(input$color_toggle)
      } else {
        use_alternating_colors_rv(!use_alternating_colors_rv())
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
      current_scan_table <- scan_table()
      # Ensure current_scan_table$chr is numeric for filtering if selected_chr is numeric
      # QTL_plot_visualizer should ensure $chr is numeric
      selected_chromosome <- main_par$selected_chr()
      if (selected_chromosome == "All") {
        current_scan_table
      } else {
        # Convert selected_chromosome to numeric to match chr column type (X=20, Y=21, M=22)
        # This conversion should ideally happen in mainParServer or selected_chr reactive
        sel_chr_num <- selected_chromosome 
        if(selected_chromosome == "X") sel_chr_num <- 20
        if(selected_chromosome == "Y") sel_chr_num <- 21
        if(selected_chromosome == "M") sel_chr_num <- 22
        sel_chr_num <- as.numeric(sel_chr_num)

        dplyr::filter(current_scan_table, chr == sel_chr_num)
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

    output$scan_plot_ui_render <- shiny::renderUI({
      shiny::req(current_scan_plot_gg())
      # Use reactive width and height for the plotOutput
      plotly::plotlyOutput(ns("render_plotly_plot"), 
                        width = paste0(plot_width_rv(), "px"), 
                        height = paste0(plot_height_rv(), "px")) |>
        shinycssloaders::withSpinner(type = 8, color = "#3498db")
    })
    
    output$render_plotly_plot <- plotly::renderPlotly({
      shiny::req(current_scan_plot_gg())
      plt <- plotly::ggplotly(current_scan_plot_gg(), 
                              source = ns("qtl_scan_plotly"), 
                              tooltip = c("x", "y", "chr")) # Basic tooltip, can be customized
      
      plt <- plt %>%
        plotly::layout(
          dragmode = "zoom", # Enables box zoom by default dragging
          hovermode = "closest",
          title = list( # Keep title minimal or remove if already handled by card header
            text = NULL 
          ),
          xaxis = list(title = current_scan_plot_gg()$labels$x), # Get x-axis label from ggplot
          yaxis = list(title = current_scan_plot_gg()$labels$y)  # Get y-axis label from ggplot
        ) %>%
        plotly::config(
          displaylogo = FALSE,
          # Plotly modebar by default includes zoom tools. Removing some less used ones.
          modeBarButtonsToRemove = c("select2d", "lasso2d", "hoverClosestCartesian", 
                                     "hoverCompareCartesian", "toggleSpikelines", "sendDataToCloud")
        )
      plt
    })

    shiny::observeEvent(plotly::event_data("plotly_click", source = ns("qtl_scan_plotly")), {
      ev_data <- plotly::event_data("plotly_click", source = ns("qtl_scan_plotly"))
      if (!is.null(ev_data) && length(ev_data) > 0) {
        # pointNumber is 0-indexed for the trace clicked
        # We need to ensure scan_table_chr() is available and gg_data is derived correctly
        # ggplotly can restructure data, so mapping back can be tricky if multiple traces exist
        # For a single trace plot, pointNumber directly maps to the row index (0-indexed)
        
        # Simplistic approach assuming pointNumber is reliable and plot data matches scan_table_chr()
        # This might need refinement if aes(color=chr) creates multiple traces in ggplotly
        clicked_idx <- ev_data$pointNumber[1] + 1 
        
        current_data <- scan_table_chr()
        if(clicked_idx <= nrow(current_data)){
          selected_point_df <- current_data[clicked_idx, , drop = FALSE]
          
          # Select relevant columns for display, similar to nearPoints output
          # Ensure 'markers', 'chr', 'position', 'LOD' exist
          cols_to_select <- c("markers", "chr", "position", "LOD")
          # Filter out missing columns before selecting
          cols_to_select <- cols_to_select[cols_to_select %in% names(selected_point_df)]
          
          if(length(cols_to_select) > 0){
            selected_point_df <- dplyr::select(selected_point_df, dplyr::all_of(cols_to_select))
            selected_point_df <- dplyr::mutate(selected_point_df, dplyr::across(dplyr::where(is.numeric), \(x) signif(x, 4)))
            clicked_plotly_point_details_rv(selected_point_df)
          } else {
            clicked_plotly_point_details_rv(data.frame(Info = "Selected point data columns not found."))
          }
        } else {
           clicked_plotly_point_details_rv(data.frame(Info = "Clicked point index out of bounds."))
        }
      } else {
        clicked_plotly_point_details_rv(NULL)
      }
    })
    
    output$plot_click_dt <- DT::renderDT({
      details <- clicked_plotly_point_details_rv()
      if (is.null(details) || nrow(details) == 0) {
        return(DT::datatable(data.frame(Info = "Click on the plot to see point details."), 
                             options = list(dom = 't'), rownames = FALSE))
      }
      DT::datatable(details, options = list(dom = 't', paging = FALSE), rownames = FALSE, selection = 'none')
    })

    # Download handlers for QTL plot
    output$download_qtl_plot_png <- shiny::downloadHandler(
      filename = function() {
        trait_name <- main_par$which_trait() %||% "plot"
        chr_suffix <- if(main_par$selected_chr() != "All") paste0("_chr", main_par$selected_chr()) else ""
        paste0("lod_plot_", trait_name, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".png")
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
        trait_name <- main_par$which_trait() %||% "plot"
        chr_suffix <- if(main_par$selected_chr() != "All") paste0("_chr", main_par$selected_chr()) else ""
        paste0("lod_plot_", trait_name, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".pdf")
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
      trait_name <- main_par$which_trait() %||% "scan"
      instanceID <- trait_name
      selected_chromosome <- main_par$selected_chr()
      if(!is.null(selected_chromosome) && selected_chromosome != "All") {
        instanceID <- paste0(instanceID, "_chr", selected_chromosome)
      }
      paste("scan", instanceID, sep = "_")
    })
    
    # `%||%` operator for cleaner NULL coalescing
    `%||%` <- function(a, b) if (!is.null(a)) a else b

    shiny::reactiveValues(
      filename = file_name_reactive, # For existing downloadApp module if it uses this
      tables = shiny::reactiveValues(scan = scan_table_chr),
      plots  = shiny::reactiveValues(scan = current_scan_plot_gg) # Ensure this is the plot with dynamic sizing
    )
  })
}

scanUI <- function(id) {
  ns <- shiny::NS(id)
  DT::DTOutput(ns("plot_click_dt"))
}

scanOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("scan_plot_ui_render"))
}