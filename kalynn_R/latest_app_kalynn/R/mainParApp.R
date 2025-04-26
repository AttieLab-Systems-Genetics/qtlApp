#' Main Parameter Module UI (Inputs)
#'
#' Creates UI elements for primary data selection inputs.
#'
#' @param id Module ID.
#' 
#' @importFrom shiny NS uiOutput sliderInput
#' @export
mainParInput <- function(id) {
  ns <- shiny::NS(id)
  list(
    # Dataset Selection (moved from app.R wellPanel)
    shiny::uiOutput(ns("selected_dataset_ui")),
    # LOD Threshold (moved from app.R LOD Plot tab)
    shiny::sliderInput(ns("LOD_thr"),
              label = shiny::h4("LOD Threshold", style = "color: #2c3e50; margin-bottom: 15px;"),
              min = 4, max = 120, value = 7, step = 0.5,
              ticks = TRUE)
  )
}

#' Main Parameter Module UI (Outputs/Dependent Inputs)
#'
#' Creates UI elements whose choices/values might depend on mainParInput.
#'
#' @param id Module ID.
#'
#' @importFrom shiny NS selectizeInput selectInput
#' @export
mainParUI <- function(id) {
  ns <- shiny::NS(id)
  list(
    # Trait Search (moved from app.R wellPanel)
    shiny::selectizeInput(ns("which_trait"),
                 label = shiny::h4("Search Trait", style = "color: #2c3e50; margin-bottom: 15px;"),
                 choices = NULL, # Server-side updated
                 selected = NULL,
                 multiple = FALSE,
                 options = list(
                     placeholder = 'Search gene symbol (e.g., Gnai3)',
                     maxOptions = 7,
                     # create = FALSE, # Allow creation? Example had it FALSE
                     maxItems = 1
                 )),
    # Chromosome Selector (moved from app.R LOD Plot tab)
    shiny::selectInput(ns("selected_chr"), "Zoom to Chromosome:",
                choices = c("All", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                          "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y", "M"),
                selected = "All",
                width = "150px"
    )
  )
}

#' Main Parameter Module Server
#'
#' Handles server logic for main parameter inputs, like updating trait choices.
#'
#' @param id Module ID.
#' @param import_data Reactive list containing `file_directory` and `annotation_list`.
#'
#' @importFrom shiny moduleServer NS observeEvent reactive req renderUI selectizeInput updateSelectizeInput isolate
#' @importFrom utils head
#' @export
mainParServer <- function(id, import_data) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Render the dataset selection UI dynamically
    output$selected_dataset_ui <- shiny::renderUI({
      shiny::req(import_data()$file_directory)
      choices <- unique(import_data()$file_directory$group)
      # Select first dataset by default if input is not already set
      selected_val <- shiny::isolate(input$selected_dataset) # Avoid triggering redraw if already set
      if(is.null(selected_val) || selected_val == "") {
          selected_val <- choices[1]
      }
      
      shiny::selectizeInput(ns("selected_dataset"),
          label = shiny::h4("Dataset Selection", style = "color: #2c3e50; margin-bottom: 15px;"),
          choices = choices,
          selected = selected_val, 
          multiple = FALSE,
          options = list(
              placeholder = 'Select a dataset...',
              onInitialize = I(sprintf('function() { this.setValue("%s"); }', selected_val)) # Ensure default selection
          ))
    })
    
    # Observe dataset selection to update trait choices
    shiny::observeEvent(input$selected_dataset, {
      shiny::req(import_data(), input$selected_dataset)
      
      # Get choices using the helper function
      trait_choices <- get_trait_choices(
          import_data()$file_directory, 
          import_data()$annotation_list, 
          input$selected_dataset
      )
      
      # Get initial subset for display if choices are numerous
      # This mirrors the behavior in the original app.R
      display_choices <- if (length(trait_choices) > 100) {
          utils::head(trait_choices, 100)
      } else {
          trait_choices
      }

      # Update the selectize input for traits
      shiny::updateSelectizeInput(session, "which_trait",
                                  choices = display_choices, 
                                  selected = "", # Clear selection when dataset changes
                                  server = TRUE) # Enable server-side searching if many options
                                  
       message("Updated trait choices for dataset: ", input$selected_dataset)
    }, ignoreNULL = TRUE, ignoreInit = TRUE) # ignoreNULL ensures it runs on first selection, ignoreInit prevents run on startup before UI exists

    # Return the reactive inputs for other modules to use
    return(reactive({
        list(
            selected_dataset = input$selected_dataset,
            which_trait = input$which_trait,
            LOD_thr = input$LOD_thr,
            selected_chr = input$selected_chr,
            # Include dependencies needed for downstream reactivity
            import_ready = !is.null(import_data())
        )
    }))
  })
} 