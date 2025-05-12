#' Peak App Module
#'
#' @param id shiny identifier
#' @param import reactive list with file_directory
#'
#' @importFrom DT datatable DTOutput renderDT
#' @importFrom shiny moduleServer NS observeEvent plotOutput reactive renderPlot
#'             renderText req selectizeInput setProgress shinyApp textOutput
#'             updateSelectizeInput withProgress h4 div
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
peakServer <- function(id, main_par, import_data, peaks_cache) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    chosen_trait <- shiny::reactive({
      shiny::req(import_data(), main_par$which_trait(), main_par$selected_dataset()) 
      trait_val <- get_selected_trait(import_data(), main_par$which_trait(), main_par$selected_dataset())
      trait_val
    })
    
    peak_table <- shiny::reactive({
      shiny::req(main_par$selected_dataset(), chosen_trait()) 
      current_dataset <- main_par$selected_dataset()
      current_chosen_trait <- chosen_trait() 
      current_trait_type <- get_trait_type(import_data(), current_dataset)
      result <- peak_finder(import_data()$file_directory, 
                            current_dataset, 
                            selected_trait = current_chosen_trait, 
                            trait_type = current_trait_type,
                            cache_env = peaks_cache)
      result
    })

    shiny::observeEvent(peak_table(), {
      peaks <- peak_table()
      current_selected_peak <- shiny::isolate(input$which_peak) 
      if (is.null(peaks)) {
        shiny::updateSelectInput(session, "which_peak", label = "Choose peak (no data)", 
                                    choices = character(0), 
                                    selected = character(0)) 
      } else {
        if (nrow(peaks) > 0) {
          if ("marker" %in% colnames(peaks) && length(peaks$marker) > 0) {
            peak_markers <- as.character(peaks$marker)
            unique_peak_markers <- unique(peak_markers)
            choices <- stats::setNames(unique_peak_markers, unique_peak_markers) 
            new_selected_value <- NULL
            if (length(choices) > 0) {
              if (!is.null(current_selected_peak) && current_selected_peak %in% choices) {
                new_selected_value <- current_selected_peak
              } else {
                new_selected_value <- choices[[1]] 
              }
            }
            shiny::updateSelectInput(session, "which_peak", 
                                      label = "Choose peak", 
                                      choices = choices, 
                                      selected = new_selected_value) 
          } else {
            shiny::updateSelectInput(session, "which_peak", label = "Choose peak (marker N/A)", 
                                        choices = character(0), 
                                        selected = character(0))
          }
        } else {
          no_peaks_choice <- list("No peaks found" = "no_peaks_found_val")
          shiny::updateSelectInput(session, "which_peak", label = "Choose peak (none found)", 
                                      choices = no_peaks_choice, 
                                      selected = "no_peaks_found_val")
        }
      }
    }, ignoreNULL = FALSE, priority = 1)

    # Table of peaks info ------------------------------------------------------
    output$peak_table <- DT::renderDT({
      peaks_data <- peak_table() # Get the reactive data
      shiny::req(peaks_data) # Temporarily simplified req for debugging
      # shiny::req(peaks_data, nrow(peaks_data) > 0) # Original req
      
      # Define base columns to hide
      cols_to_hide_base <- c("ci_lo", "ci_hi")
      # Conditionally add gene_id if it exists
      cols_to_hide <- cols_to_hide_base
      if("gene_id" %in% colnames(peaks_data)){
          cols_to_hide <- c(cols_to_hide, "gene_id")
      }
      
      # Check which columns to hide actually exist in the current data
      cols_that_exist_and_should_be_hidden <- intersect(cols_to_hide, colnames(peaks_data))
      
      # Build columnDefs dynamically
      colDefs <- list(list(targets = '_all', className = 'dt-center'))
      if(length(cols_that_exist_and_should_be_hidden) > 0){
          colDefs <- c(colDefs, list(list(targets = cols_that_exist_and_should_be_hidden, visible = FALSE)))
      }
      
      DT::datatable(
          peaks_data, 
          options = list(paging = FALSE,
                         scrollX = TRUE,
                         scrollY = TRUE,
                         autoWidth = TRUE,
                         dom = 'Bt',
                         buttons = list(
                           list(extend = 'csv', className = 'btn-sm'),
                           list(extend = 'excel', className = 'btn-sm')
                         ),
                         columnDefs = colDefs # Use dynamically generated colDefs
          ),
          extensions = 'Buttons',
          selection = 'single',
          filter = 'none',
          rownames = TRUE
        )
    })
 
    output$allele_effects_plot <- shiny::renderPlot({
      shiny::req(allele_plot())
      allele_plot()
    })
    allele_plot <- shiny::reactive({
      # Use the already filtered peak_table()
      shiny::req(peak_table(), input$which_peak)
      # pivot_peaks expects the 'marker' column and A-H columns, which peak_finder should provide
      peak_long <- pivot_peaks(peak_table(), input$which_peak) 
      
      # Call ggplot_alleles without the extra arguments
      ggplot_alleles(peak_long)
    })
    file_name <- shiny::reactive({
      # Evaluate chosen_trait reactive (gene symbol)
      trait_val <- shiny::req(chosen_trait())
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
  # Source UI styling functions if not already loaded
  if (!exists("create_select_input", mode = "function")) {
    source("R/ui_styles.R")
  }
  
  ns <- shiny::NS(id)
  
  # Use modern styling if available, otherwise use standard controls
  if (exists("create_select_input", mode = "function")) {
    create_well_panel(
      h4("Strain Effects", style = "color: #2c3e50; margin-bottom: 15px;"),
      div(style = "color: #7f8c8d; margin-bottom: 15px;",
        "Select a peak to see strain effects (Only for additive scans)."
      ),
      div(style = "margin-bottom: 20px; position: relative; z-index: 2;",
        create_select_input(ns("which_peak"),
          label = NULL,
          choices = NULL,
          multiple = FALSE,
          options = list(
            placeholder = 'Select a peak...',
            onInitialize = I('function() { this.setValue(""); }')
          )
        )
      )
    )
  } else {
    list(
      shiny::helpText("Choose a peak to see the strain effects.",
        "This only applies to the additive scans."),
      shiny::selectInput(ns("which_peak"),
        label = "Choose peak",
        choices = character(0),
        multiple = FALSE
      )
    )
  }
}
#' @rdname peakApp
#' @export
peakUI <- function(id) {
  # Source UI styling functions if not already loaded
  if (!exists("create_plot_output", mode = "function")) {
    source("R/ui_styles.R")
  }
  
  ns <- shiny::NS(id)
  
  # Use modern styling if available, otherwise use standard output
  if (exists("create_plot_output", mode = "function")) {
    div(style = "margin-top: 20px; position: relative; z-index: 0;",
      shiny::plotOutput(ns("allele_effects_plot"), height = "400px") %>%
        shinycssloaders::withSpinner(type = 8, color = "#3498db", proxy.height = "400px")
    )
  } else {
    shiny::plotOutput(ns("allele_effects_plot"))
  }
}
#' @rdname peakApp
#' @export
peakOutput <- function(id) {
  ns <- shiny::NS(id)
  DT::DTOutput(ns("peak_table"))
}
