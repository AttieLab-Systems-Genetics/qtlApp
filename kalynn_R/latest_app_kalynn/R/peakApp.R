#' Peak Module Input UI
#'
#' Creates UI for selecting a peak.
#'
#' @param id Module ID.
#'
#' @importFrom shiny NS selectizeInput
#' @export
peakInput <- function(id) {
  ns <- shiny::NS(id)
  list(
    shiny::selectizeInput(ns("which_peak"), "Choose peak",
                 choices = NULL, # Server-side updated
                 multiple = FALSE,
                 options = list(
                     placeholder = 'Select a peak...',
                     onInitialize = I('function() { this.setValue(""); }')
                 ))
  )
}

#' Peak Module Plot UI
#'
#' Creates UI for displaying the allele effects plot.
#'
#' @param id Module ID.
#'
#' @importFrom shiny NS plotOutput uiOutput
#' @importFrom shinycssloaders withSpinner
#' @export
peakUI <- function(id) {
  ns <- shiny::NS(id)
  # Use uiOutput to render the plot dynamically based on data availability
  shiny::uiOutput(ns("allele_effects_ui")) 
}

#' Peak Module Output UI (Placeholder/Example)
#'
#' Example structure for a peak table output, if needed separately.
#' Currently not used in the main app structure based on app.R
#'
#' @param id Module ID.
#'
#' @importFrom DT DTOutput
#' @export
peakOutput <- function(id){
    ns <- shiny::NS(id)
    # Example: DT::DTOutput(ns("peak_table_output"))
    NULL # Returning NULL as it's not currently placed in app.R
}


#' Peak Module Server
#'
#' Handles logic for finding peaks, updating peak selection, and generating allele plots.
#'
#' @param id Module ID.
#' @param main_par Reactive containing main parameters (selected_dataset, which_trait, LOD_thr).
#' @param import_data Reactive containing imported data (file_directory, annotation_list).
#'
#' @importFrom shiny moduleServer NS observeEvent reactive req renderPlot renderUI updateSelectizeInput reactiveVal
#' @importFrom shinycssloaders withSpinner
#' @export
peakServer <- function(id, main_par, import_data) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive to get the actual selected trait name (handles genes/isoforms)
    selected_trait_name <- reactive({
      shiny::req(main_par()$which_trait, main_par()$selected_dataset)
      get_selected_trait(import_data()$file_directory, 
                         main_par()$which_trait, 
                         main_par()$selected_dataset)
    })

    # Reactive to find peaks for the selected trait
    peak_table_full <- reactive({
      shiny::req(main_par()$selected_dataset, selected_trait_name())
      message("peakServer: Fetching peaks for trait: ", selected_trait_name(), " in dataset: ", main_par()$selected_dataset)
      # Need peak_finder function (assuming it's loaded/sourced)
      peak_finder(import_data()$file_directory, 
                  main_par()$selected_dataset, 
                  selected_trait_name())
    })
    
    # Reactive for filtered and ordered peaks based on LOD threshold
    peaks_filtered_ordered <- reactive({ 
        shiny::req(peak_table_full(), main_par()$LOD_thr)
        highest_peaks(peak_table_full(), main_par()$LOD_thr)
    })

    # Update peak selection dropdown when filtered peaks change
    observeEvent(peaks_filtered_ordered(), {
        peaks_data <- peaks_filtered_ordered()
        
        if (!is.null(peaks_data) && nrow(peaks_data) > 0) {
            # Create named vector for dropdown choices
            marker_choices <- peaks_data$marker
            names(marker_choices) <- paste0(peaks_data$marker,
                                         " (Chr", chr_XYM(peaks_data$chr), ": ", round(peaks_data$pos, 2),
                                         " Mb, LOD: ", round(peaks_data$lod, 2), ")")
            
            # Keep current selection if it's still valid, otherwise select top peak
            current_selection <- shiny::isolate(input$which_peak)
            selected_val <- if (!is.null(current_selection) && current_selection %in% marker_choices) {
                 current_selection
             } else {
                 marker_choices[1] # Default to highest LOD peak
             }
            
            shiny::updateSelectizeInput(session, "which_peak",
                choices = marker_choices,
                selected = selected_val,
                server = FALSE # Usually few enough peaks to keep client-side
            )
             message("peakServer: Updated peak choices. Selected: ", selected_val)
        } else {
            # Clear choices and selection if no peaks meet threshold
            shiny::updateSelectizeInput(session, "which_peak",
                choices = character(0),
                selected = character(0)
            )
            message("peakServer: Cleared peak choices (no peaks >= LOD threshold).")
        }
    })

    # Reactive value to store the generated allele plot
    strain_plot_obj <- reactiveVal(NULL)

    # Render the allele effects plot UI (which contains the plot output)
    output$allele_effects_ui <- shiny::renderUI({
        shiny::plotOutput(ns("allele_effects_plot")) %>% 
            shinycssloaders::withSpinner(type = 8, color = "#3498db")
    })

    # Generate and render the allele effects plot
    output$allele_effects_plot <- shiny::renderPlot({
        # Requires selected peak and the full peak table
        shiny::req(input$which_peak, peak_table_full())
        message("peakServer: Rendering allele plot for peak: ", input$which_peak)
        
        # Pivot the data for the selected peak
        peak_long <- pivot_peaks(peak_table_full(), input$which_peak)
        
        # Generate plot using helper
        # Need ggplot_alleles function (assuming loaded/sourced)
        plot_alleles <- ggplot_alleles(
            peak_long = peak_long, 
            trait_name = selected_trait_name(),
            peak_marker = input$which_peak
            )
        
        # Store the plot for download handler (if download is handled here or needs access)
        strain_plot_obj(plot_alleles)
        
        # Return the plot to render
        plot_alleles
    })

    # Return reactive values needed by other modules (e.g., download)
    return(reactive({
      list(
        peak_table = peak_table_full(), # Return the full table for info display potentially
        selected_peak = input$which_peak, 
        allele_plot = strain_plot_obj() # Pass the generated plot object
      )
    }))
  })
} 