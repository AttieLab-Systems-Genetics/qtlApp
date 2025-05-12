#' Main Parameter App Module
#' 
#' @param id shiny identifier
#' @param import reactive list with file_directory and annotation_list
#'
#' @importFrom shiny moduleServer reactive renderPrint req selectizeInput
#'             shinyApp sliderInput updateSelectizeInput verbatimTextOutput
#' @importFrom bslib page
#' 
#' @export
mainParApp <- function(id) {
  ui <- bslib::page(
    title = "Test MainPar Module",
    mainParInput("main_par"), # "selected_dataset", "LOD_thr"
    mainParUI("main_par"),    # "which_trait", "which_chr"
    mainParOutput("main_par")
  )
  server <- function(input, output, session) {
    import <- importServer("import")
    mainParServer("main_par", import)
  }
  shiny::shinyApp(ui = ui, server = server)
}
#' @rdname mainParApp
#' @export
mainParServer <- function(id, import) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Populate Dataset Category choices
    output$dataset_category <- shiny::renderUI({
      shiny::req(import(), import()$file_directory)
      # Assuming file_directory now has a 'dataset_category' column
      if (!("dataset_category" %in% colnames(import()$file_directory))){
        warning("mainParServer: 'dataset_category' column missing in file_directory.")
        return(shiny::selectInput(ns("dataset_category"), label = "Dataset Category:", choices = "Error: category column missing"))
      }
      categories <- sort(unique(import()$file_directory$dataset_category))
      selected_category <- if (length(categories) > 0) categories[1] else NULL
      
      if (exists("create_select_input", mode = "function")) {
        create_select_input(ns("dataset_category"), label = "Dataset Category:", 
                            choices = categories, selected = selected_category, multiple = FALSE)
      } else {
        shiny::selectInput(ns("dataset_category"), label = "Dataset Category:", 
                           choices = categories, selected = selected_category, multiple = FALSE)
      }
    })

    # Populate Dataset choices based on selected category
    output$selected_dataset <- shiny::renderUI({
      shiny::req(import(), import()$file_directory, input$dataset_category)
      if (!("dataset_category" %in% colnames(import()$file_directory)) || !("group" %in% colnames(import()$file_directory))){
         warning("mainParServer: Required columns missing for dataset filtering.")
         return(shiny::selectInput(ns("selected_dataset"), label = "Choose a dataset:", choices = "Error: columns missing"))
      }
      filtered_files <- dplyr::filter(import()$file_directory, .data$dataset_category == input$dataset_category)
      choices <- sort(unique(filtered_files$group))
      selected <- if (length(choices) > 0) choices[1] else NULL
      
      if (exists("create_select_input", mode = "function")) {
        create_select_input(ns("selected_dataset"), label = "Choose a dataset:",
                            choices = choices, selected = selected, multiple = FALSE)
      } else {
        shiny::selectInput(ns("selected_dataset"), label = "Choose a dataset:",
                           choices = choices, selected = selected, multiple = FALSE)
      }
    })

    # Update trait choices (existing logic)
    shiny::observeEvent(shiny::req(input$selected_dataset), {
      shiny::req(import())
      current_ds <- input$selected_dataset
      message(paste("mainParServer: observeEvent for trait update. Dataset:", current_ds)) # DEBUG
      
      # Freeze input$which_trait while we update it to prevent premature reactions
      shiny::freezeReactiveValue(input, "which_trait")
      message("mainParServer: Froze input$which_trait.") # DEBUG

      # First, clear the existing selection
      shiny::updateSelectizeInput(session, "which_trait", choices = character(0), selected = character(0), 
                                  options = list(placeholder = 'Loading traits...')) 
      message("mainParServer: Cleared which_trait via updateSelectizeInput.") # DEBUG

      choices <- get_trait_choices(import(), current_ds)
      message(paste("mainParServer: Generated trait choices for", current_ds, ":", paste(head(choices), collapse=", "))) # DEBUG 
      
      new_selected_trait <- NULL
      if (!is.null(choices) && length(choices) > 0) {
        new_selected_trait <- choices[1] 
        message(paste("mainParServer: For dataset", current_ds, ", new_selected_trait will be:", new_selected_trait)) # DEBUG
      } else {
        message(paste("mainParServer: No trait choices found for dataset", current_ds))
      }
      
      # Update with new choices and select the first one
      # No need to freeze again if this is the last update to which_trait in this observer block
      shiny::updateSelectizeInput(session, "which_trait",
        choices = choices, 
        selected = new_selected_trait, 
        options = list(placeholder = 'Search or select trait...', maxItems = 1, maxOptions = 10), 
        server = TRUE 
      )
      message(paste("mainParServer: updateSelectizeInput for 'which_trait' called. Choices sent:", paste(head(choices),collapse=", "), "Attempted to select:", new_selected_trait)) # DEBUG
    })

    # Reactive for currently selected trait based on UI input, but validated against current dataset
    current_selected_trait <- shiny::reactive({
      shiny::req(input$selected_dataset) # Require selected_dataset first
      
      current_ds <- input$selected_dataset
      current_trait_input <- input$which_trait # Value from UI
      
      message(paste("mainParServer: current_selected_trait CALC START. Dataset:", current_ds, "UI which_trait:", current_trait_input)) # DEBUG

      valid_choices <- get_trait_choices(import(), current_ds)
      message(paste("mainParServer: current_selected_trait valid_choices for", current_ds, ":", paste(head(valid_choices), collapse=", "))) # DEBUG

      final_trait_to_return <- NULL

      if (!is.null(valid_choices) && length(valid_choices) > 0) {
        if (!is.null(current_trait_input) && current_trait_input %in% valid_choices) {
          final_trait_to_return <- current_trait_input
          message(paste("mainParServer: UI which_trait (", current_trait_input, ") is VALID. Returning it.")) # DEBUG
        } else {
          final_trait_to_return <- valid_choices[1] # Default to first valid choice
          message(paste("mainParServer: UI which_trait (", current_trait_input, ") is INVALID or NULL. Defaulting to first valid choice:", final_trait_to_return)) # DEBUG
        }
      } else {
        message(paste("mainParServer: No valid_choices for dataset", current_ds, ". Returning NULL from current_selected_trait.")) # DEBUG
      }
      
      # This req ensures we don't proceed with a NULL trait if no valid options exist.
      shiny::req(final_trait_to_return, cancelOutput = TRUE) # cancelOutput will prevent downstream errors if this is NULL
      message(paste("mainParServer: current_selected_trait CALC END. Returning:", final_trait_to_return)) # DEBUG
      return(final_trait_to_return)
    })

    # Show returned values.
    output$returns <- shiny::renderPrint({
      cat("dataset_category =", input$dataset_category,
          "\nselected_dataset =", input$selected_dataset,
          "\nwhich_trait =", input$which_trait,
          "\nselected_chr =", input$selected_chr,
          "\nLOD_thr =", input$LOD_thr)
    })
    
    # Return reactive expressions for inputs
    return(
      list(
        dataset_category = shiny::reactive(input$dataset_category),
        selected_dataset = shiny::reactive(input$selected_dataset),
        which_trait = current_selected_trait, # Use the new validated reactive
        selected_chr = shiny::reactive(input$selected_chr),
        LOD_thr = shiny::reactive(input$LOD_thr)
      )
    )
  })
}
#' @rdname mainParApp
#' @export
mainParInput <- function(id) {
  # Source UI styling functions if not already loaded
  if (!exists("create_slider_input", mode = "function")) {
    source("R/ui_styles.R")
  }
  
  ns <- shiny::NS(id)
  
  if (exists("create_select_input", mode = "function")) {
    list(
      shiny::uiOutput(ns("dataset_category")), # UI for new category dropdown
      shiny::uiOutput(ns("selected_dataset")), # Existing dataset dropdown
      create_slider_input(ns("LOD_thr"),
        label = "LOD threshold for evaluation",
        min = 4, max = 20, value = 7.5, step = 0.5)
    )
  } else {
    list(
      shiny::uiOutput(ns("dataset_category")),
      shiny::uiOutput(ns("selected_dataset")), 
      shiny::sliderInput(ns("LOD_thr"),
        label = "LOD threshold for evaluation",
        min = 4, max = 20, value = 7.5, step = 0.5)
    )
  }
}
#' @rdname mainParApp
#' @export
mainParUI <- function(id) {
  # Source UI styling functions if not already loaded
  if (!exists("create_select_input", mode = "function")) {
    source("R/ui_styles.R")
  }
  
  ns <- shiny::NS(id)
  
  # Use modern styling if available, otherwise use standard controls
  if (exists("create_select_input", mode = "function")) {
    list(
      create_select_input(ns("which_trait"),
        label = "Choose the trait",
        choices = NULL,
        selected = NULL,
        multiple = FALSE,
        options = list(
          placeholder = 'Search gene symbol...',
          maxItems = 1,
          maxOptions = 7
        )),
      create_select_input(ns("selected_chr"), 
        label = "Zoom to Chromosome:",
        choices = c("All",
                  "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                  "11", "12", "13", "14", "15", "16", "17", "18", "19",
                  "X", "Y", "M"),
        selected = "All",
        width = "150px")
    )
  } else {
  list(
    shiny::selectizeInput(ns("which_trait"),
      label = "Choose the trait",
      choices = NULL,
      multiple = FALSE,
        options = list(
          placeholder = 'Search gene symbol...',
          maxItems = 1,
          maxOptions = 7
        )),
      shiny::selectInput(ns("selected_chr"), 
        label = "Zoom to Chromosome:",
      choices = c("All",
                  "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                  "11", "12", "13", "14", "15", "16", "17", "18", "19",
                  "X", "Y", "M"),
      selected = "All",
        width = "150px")
    )
  }
}
#' @rdname mainParApp
#' @export
mainParOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::verbatimTextOutput(ns("returns"))
}
