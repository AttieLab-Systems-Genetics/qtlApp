#' Peak App Module
#'
#' @param id shiny identifier
#' @param import reactive list with file_directory
#'
#' @importFrom DT datatable DTOutput renderDT
#' @importFrom shiny moduleServer NS observeEvent plotOutput reactive renderPlot
#'             renderText req selectizeInput setProgress shinyApp textOutput
#'             updateSelectizeInput withProgress
#' @importFrom bslib card card_header page_sidebar sidebar
#' @importFrom shinycssloaders withSpinner
#' @importFrom stringr str_split
#' @importFrom reshape2 melt
#' @importFrom ggplot2 aes element_blank element_line element_text geom_hline
#'             geom_point ggplot labs scale_color_manual theme theme_bw
#' 
#' @export
peakApp <- function() {
  ui <- bslib::page_sidebar(
    title = "Test Peak",
    sidebar = bslib::sidebar("side_panel",
      mainParInput("main_par"),    # "selected_dataset", "LOD_thr"
      mainParUI("main_par"),       # "which_trait"
      peakInput("peak"),           # "which_peak"
      downloadInput("download"),   # downloadButton, filename
      downloadOutput("download")), # plot_table, inputs for Plots or Tables
    bslib::card(
      bslib::card_header("Peaks"),
      peakOutput("peak")),
    bslib::card(
      bslib::card_header("Strain effects"),
      peakUI("peak"))
  )
  server <- function(input, output, session) {
      import <- importServer("import")
      main_par <- mainParServer("main_par", import)
      peak_list <- peakServer("peak", main_par, import)
      downloadServer("download", peak_list)
  }
  shiny::shinyApp(ui = ui, server = server)
}
#' @rdname peakApp
#' @export
peakServer <- function(id, main_par, import) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    shiny::observe({
      message("--- peakServer Debug ---")
      message(paste("Selected Dataset (from main_par):", shiny::isolate(main_par$selected_dataset())))
      message(paste("Which Trait (from main_par):", shiny::isolate(main_par$which_trait())))
      imp_data <- shiny::isolate(import())
      if(!is.null(imp_data) && !is.null(imp_data$file_directory)){
        message(paste("Imported file_directory rows:", nrow(imp_data$file_directory)))
      } else {
        message("Imported file_directory is NULL or not fully available yet.")
      }
      message("------------------------")
    })

    chosen_trait <- shiny::reactive({
      shiny::req(import(), main_par$which_trait(), main_par$selected_dataset()) 
      # Get selected trait and log it
      trait_val <- get_selected_trait(import(),
                         main_par$which_trait(), main_par$selected_dataset())
      message(paste("peakServer chosen_trait_reactive: Calculated chosen_trait as:", trait_val))
      trait_val
    })
    peak_table <- shiny::reactive({
      shiny::req(main_par$selected_dataset(), chosen_trait())
      
      # For logging
      current_dataset <- main_par$selected_dataset()
      current_chosen_trait <- chosen_trait()
      message(paste("peak_table_reactive: Attempting peak_finder for trait '", current_chosen_trait, "' in dataset '", current_dataset, "'"))
      
      # Call peak_finder (temporarily remove suppressMessages for debugging)
      found_peaks <- peak_finder(import()$file_directory, current_dataset) 
      message(paste("peak_table_reactive: peak_finder returned", nrow(found_peaks), "total peaks for dataset before trait subset."))
      
      # Subset for the chosen trait
      # Ensure 'trait' column exists in found_peaks before subsetting
      if("trait" %in% colnames(found_peaks)){
        result <- subset(found_peaks, trait == current_chosen_trait)
        message(paste("peak_table_reactive: After subsetting for trait '", current_chosen_trait, "', found", nrow(result), "peaks."))
      } else {
        message(paste("peak_table_reactive: 'trait' column not found in peaks returned by peak_finder. Cannot subset for trait '", current_chosen_trait, "'. Returning all found peaks for dataset."))
        result <- found_peaks # Or an empty data frame: result <- found_peaks[0,]
      }
      result
    })

    # Observer to update which_peak choices when peak_table changes
    shiny::observeEvent(peak_table(), {
      peaks <- peak_table()
      message("--- peakApp.R: Observer for which_peak (using updateSelectInput) ---")
      current_selected_peak <- shiny::isolate(input$which_peak) 

      if (is.null(peaks)) {
        message("peaks object is NULL. Clearing choices.")
        shiny::updateSelectInput(session, "which_peak", label = "Choose peak (no data)", 
                                    choices = character(0), 
                                    selected = character(0)) 
      } else {
        message(paste("peaks data frame has", nrow(peaks), "rows."))
        if (nrow(peaks) > 0) {
          message(paste("Columns in peaks data frame:", paste(colnames(peaks), collapse=", ")))
          if ("marker" %in% colnames(peaks) && length(peaks$marker) > 0) {
            peak_markers <- as.character(peaks$marker)
            choices <- stats::setNames(peak_markers, peak_markers) 
            
            message(paste("Generated choices for dropdown:", paste(names(choices), "=", choices, collapse="; ")))
            
            new_selected_value <- NULL
            if (length(choices) > 0) {
              if (!is.null(current_selected_peak) && current_selected_peak %in% choices) {
                new_selected_value <- current_selected_peak
              } else {
                new_selected_value <- choices[[1]] 
              }
            }
            message(paste("New selected value will be:", new_selected_value))

            shiny::updateSelectInput(session, "which_peak", 
                                      label = "Choose peak", 
                                      choices = choices, 
                                      selected = new_selected_value) 
            message("Updated 'which_peak' dropdown with actual peak markers.")
          } else {
            message("'marker' column NOT FOUND or empty in peaks data frame! Clearing choices.")
            shiny::updateSelectInput(session, "which_peak", label = "Choose peak (marker N/A)", 
                                        choices = character(0), 
                                        selected = character(0))
          }
        } else {
          message("No peaks found (nrow is 0). Updating dropdown to 'No peaks found'.")
          no_peaks_choice <- list("No peaks found" = "no_peaks_found_val")
          shiny::updateSelectInput(session, "which_peak", label = "Choose peak (none found)", 
                                      choices = no_peaks_choice, 
                                      selected = "no_peaks_found_val")
        }
      }
      message("--- End peakApp.R: Observer for which_peak ---")
    }, ignoreNULL = FALSE, priority = 1) 

    # Table ov peaks info ------------------------------------------------------
    output$peak_table <- DT::renderDT({DT::datatable(
      peak_table(),
      options = list(paging = TRUE,    ## paginate the output
                     pageLength = 5,   ## number of rows to output for each page
                     scrollX = TRUE,   ## enable scrolling on X axis
                     scrollY = TRUE,   ## enable scrolling on Y axis
                     autoWidth = TRUE, ## use smart column width handling
                     server = TRUE,    ## use client-side processing
                     dom = 'Bfrtip',
                     buttons = c('csv', 'excel'),
                     columnDefs = list(
                       list(targets = '_all', className = 'dt-center'),
                       list(targets = c(0, 8, 9), visible = FALSE))
      ),
      extensions = 'Buttons',
      selection = 'single',            ## enable selection of a single row
      filter = 'bottom',               ## include column filters at the bottom
      rownames = TRUE)                 ##  show row numbers/names
    })
 
    output$allele_effects_plot <- shiny::renderPlot({
      shiny::req(allele_plot())
      allele_plot()
    })
    allele_plot <- shiny::reactive({
      shiny::req(peak_table(), input$which_peak)
      peak <- pivot_peaks(peak_table(), input$which_peak)
      ggplot_alleles(peak)
    })
    file_name <- shiny::reactive({
      # Evaluate which_trait reactive
      trait_val <- shiny::req(main_par$which_trait())
      peak_val <- shiny::req(input$which_peak)
      instanceID <- paste(trait_val, peak_val, sep = "_") 
      
      # Debug before comparison
      sel_chr_val <- main_par$selected_chr()
      
      # Defensive check before comparison
      if (is.atomic(sel_chr_val) && is.character(sel_chr_val) && length(sel_chr_val) == 1) {
        # Evaluate selected_chr reactive and compare                        
        if(shiny::req(sel_chr_val) != "All") { 
          # Evaluate selected_chr reactive for paste
          instanceID <- paste0(instanceID, "_chr", sel_chr_val)
        }
      } else {
        warning("peakServer file_name: sel_chr_val is not a single character string!")
      }
      
      paste("peak", instanceID, sep = "_")
    })
    # Return `peak_list` = reactiveValues containing elements `filename`, `tables` and `plots`.
    shiny::reactiveValues(
      filename = file_name,
      tables = shiny::reactiveValues(
        peak = peak_table),
      plots  = shiny::reactiveValues(
        alleles = allele_plot)
    )
    # ------------------------------------------------------------------

  })
}
#' @rdname peakApp
#' @export
peakInput <- function(id) {
  ns <- shiny::NS(id)
  list(
    shiny::helpText("Choose a peak to see the strain effects.",
      "This only applies to the additive scans."),
    shiny::selectInput(ns("which_peak"),
      label = "Choose peak (selectInput Test)",
      choices = character(0),
      multiple = FALSE
    )
  )
}
#' @rdname peakApp
#' @export
peakUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::plotOutput(ns("allele_effects_plot"))
}
#' @rdname peakApp
#' @export
peakOutput <- function(id) {
  ns <- shiny::NS(id)
  DT::DTOutput(ns("peak_table"))
}
