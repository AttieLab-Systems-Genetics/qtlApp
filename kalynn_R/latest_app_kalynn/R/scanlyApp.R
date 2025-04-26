#' Scanly Module Output UI
#'
#' Creates the plotly output container.
#'
#' @param id Module ID.
#'
#' @importFrom shiny NS uiOutput
#' @importFrom plotly plotlyOutput
#' @importFrom shinycssloaders withSpinner
#' @export
scanlyOutput <- function(id) {
  ns <- shiny::NS(id)
  # Wrap plotlyOutput in uiOutput to prevent errors during initialization
  shiny::uiOutput(ns("scan_plot_ui")) 
}

#' Scanly Module Info Table UI
#'
#' Creates the DT output container for clicked/selected point info.
#'
#' @param id Module ID.
#'
#' @importFrom shiny NS tagList h5
#' @importFrom DT DTOutput
#' @export
scanlyUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h5("Selected Point Information", 
              style = "margin: 0 0 10px 0; color: #2c3e50; font-size: 14px;"),
    DT::DTOutput(ns("clicked_point_info"))
  )
}


#' Scanly Module Server
#'
#' Handles the interactive QTL scan plot (Plotly) and clicked point info display.
#'
#' @param id Module ID.
#' @param main_par Reactive containing main parameters.
#' @param scan_object Reactive containing the ggplot object (`plot`) and 
#'        the underlying data (`table`) for the current scan.
#' @param peak_mod_reactives Reactive containing outputs from the peak module,
#'        specifically `peak_table` and `selected_peak`.
#' @param official_trait_symbol Reactive containing the official trait symbol string.
#'
#' @importFrom shiny moduleServer NS reactive reactiveVal observeEvent req isolate
#' @importFrom plotly renderPlotly event_data
#' @importFrom DT renderDT datatable
#' @importFrom htmltools tags
#' @importFrom shinycssloaders withSpinner
#' @export
scanlyServer <- function(id, main_par, scan_object, peak_mod_reactives, official_trait_symbol) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Render the Plotly plot UI
    output$scan_plot_ui <- shiny::renderUI({
        plotly::plotlyOutput(ns("scan_plot_render")) %>%
            shinycssloaders::withSpinner(type = 8, color = "#3498db")
    })

    # Reactive for the Plotly plot object
    scanly_plot_obj <- reactive({
      shiny::req(scan_object(), 
                 scan_object()$plot, 
                 scan_object()$table,
                 main_par()$selected_chr,
                 peak_mod_reactives()$peak_table, 
                 official_trait_symbol()) # Need official symbol for title
      
      message("scanlyServer: Generating plotly object.")
      
      # Use helper function to convert ggplot to plotly
      ggplotly_qtl_scan(
        scan_plot = scan_object()$plot, 
        plot_data = scan_object()$table,
        peak_table = peak_mod_reactives()$peak_table,
        official_trait_symbol = official_trait_symbol(),
        selected_chr = main_par()$selected_chr,
        # Pass plot dimensions if they are inputs elsewhere (e.g., in main app UI)
        # plot_width = input$plot_width, # Example if width/height inputs exist
        # plot_height = input$plot_height,
        source = ns("scan_plot_render") # Use namespaced source ID
      )
    })

    # Render the Plotly plot
    output$scan_plot_render <- plotly::renderPlotly({
      scanly_plot_obj() # Render the reactive plotly object
    })

    # --- Clicked/Selected Point Info Table --- 
    
    # Reactive value to store click event data
    plotly_click_data <- reactiveVal(NULL)
    
    # Observe plotly click events
    observeEvent(plotly::event_data('plotly_click', source = ns("scan_plot_render")), {
        message("scanlyServer: Plotly click event detected.")
        plotly_click_data(plotly::event_data('plotly_click', source = ns("scan_plot_render")))
    })
    
    # Clear click data when the trait changes to avoid showing old clicks
    observeEvent(main_par()$which_trait, {
        message("scanlyServer: Trait changed, clearing click data.")
        plotly_click_data(NULL)
    })

    # Reactive for the data frame to display in the DT table
    info_table_data <- reactive({
        # Ensure dependencies are met
        shiny::req(main_par()$selected_chr,
                 scan_object()$table, # Need plot data for click context
                 peak_mod_reactives()$peak_table, # Need peak table
                 peak_mod_reactives()$selected_peak) # Need selected peak
                 
        message("scanlyServer: Calculating info table data.")
        
        # Use the peak_info helper function
        peak_info(
            peak_table = peak_mod_reactives()$peak_table, 
            plot_data = scan_object()$table, 
            selected_peak_marker = peak_mod_reactives()$selected_peak, 
            click_event = plotly_click_data(), # Use the reactive value for click event
            selected_chr = main_par()$selected_chr
        )
    })

    # Render the DT table
    output$clicked_point_info <- DT::renderDT({
      df <- info_table_data()
      if (is.null(df)) {
          # Return an empty table or a message if no data
          return(DT::datatable(data.frame(Info = "Click plot or select peak"), 
                             options = list(dom = 't', ordering = FALSE), 
                             rownames = FALSE, 
                             caption = htmltools::tags$caption(
                                style = 'caption-side: top; text-align: left; color: #2c3e50; font-weight: bold; font-size: 14px;',
                                'Peak Information and Strain Effects'
                             )))
      }
      
      message("scanlyServer: Rendering info table.")
      DT::datatable(
        df,
        options = list(dom = 't', ordering = FALSE, pageLength = 1),
        rownames = FALSE,
        class = 'compact hover',
        caption = htmltools::tags$caption(
          style = 'caption-side: top; text-align: left; color: #2c3e50; font-weight: bold; font-size: 14px;',
          'Peak Information and Strain Effects'
        )
      )
    })
    
    # --- Double Click Logic --- 
    # Note: Resetting zoom via double-click directly in plotly is tricky.
    # A common workaround is to observe the double-click and update an input 
    # (like selected_chr) that controls the plot view. This requires the input
    # to be controllable by the server (e.g., using updateSelectInput in the main app
    # or passing a session object to the module - latter is complex).
    # For now, we just detect the double-click.
    observeEvent(plotly::event_data("plotly_doubleclick", source = ns("scan_plot_render")), {
        message("scanlyServer: Plotly double-click event detected.")
        # Ideally, this would trigger an update to main_par$selected_chr -> "All"
        # This requires communication back to the parent module or main app.
    })

    # Return value (optional, could return click data if needed elsewhere)
    return(reactive({ 
        list(click = plotly_click_data()) 
    }))
  })
} 