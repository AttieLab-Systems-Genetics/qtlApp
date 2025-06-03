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
  # MAIN APP WRAPPER: Creates standalone Shiny app for testing the module
  ui <- bslib::page(
    title = "Test MainPar Module",
    mainParInput("main_par"), # Input controls: LOD threshold slider + dataset dropdowns
    mainParUI("main_par"), # Core UI: trait selector + chromosome dropdown
    mainParOutput("main_par") # Debug output: shows current parameter values
  )
  server <- function(input, output, session) {
    import <- importServer("import") # Get data import functionality
    mainParServer("main_par", import) # Connect main parameter server logic
  }
  shiny::shinyApp(ui = ui, server = server)
}

#' @rdname mainParApp
#' @export
mainParServer <- function(id, import) {
  # MAIN SERVER LOGIC: Handles all reactive behavior for parameter selection
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ====== DATASET CATEGORY DROPDOWN ======
    # Creates dropdown for selecting dataset categories
    output$dataset_category <- shiny::renderUI({
      shiny::req(import(), import()$file_directory) # Wait for import data

      # ERROR HANDLING: Check if required column exists
      if (!("dataset_category" %in% colnames(import()$file_directory))) {
        warning("mainParServer: 'dataset_category' column missing in file_directory.")
        return(shiny::selectInput(ns("dataset_category"),
          label = "Dataset Category:",
          choices = "Error: category column missing"
        ))
      }

      # GET CHOICES: Extract unique categories and sort them
      categories <- sort(unique(import()$file_directory$dataset_category))
      selected_category <- if (length(categories) > 0) categories[1] else NULL

      # RENDER INPUT: Use custom styling if available, otherwise default Shiny
      if (exists("create_select_input", mode = "function")) {
        create_select_input(ns("dataset_category"),
          label = "Dataset Category:",
          choices = categories, selected = selected_category, multiple = FALSE
        )
      } else {
        shiny::selectInput(ns("dataset_category"),
          label = "Dataset Category:",
          choices = categories, selected = selected_category, multiple = FALSE
        )
      }
    })

    # ====== SPECIFIC DATASET DROPDOWN ======
    # Creates dropdown for selecting specific datasets within chosen category
    output$selected_dataset <- shiny::renderUI({
      shiny::req(import(), import()$file_directory, input$dataset_category) # Wait for category selection

      # ERROR HANDLING: Check if required columns exist
      if (!("dataset_category" %in% colnames(import()$file_directory)) ||
        !("group" %in% colnames(import()$file_directory))) {
        warning("mainParServer: Required columns missing for dataset filtering.")
        return(shiny::selectInput(ns("selected_dataset"),
          label = "Choose a dataset:",
          choices = "Error: columns missing"
        ))
      }

      # FILTER DATA: Get datasets matching selected category
      filtered_files <- dplyr::filter(
        import()$file_directory,
        .data$dataset_category == input$dataset_category
      )
      choices <- sort(unique(filtered_files$group)) # Extract unique dataset groups
      selected <- if (length(choices) > 0) choices[1] else NULL

      # RENDER INPUT: Use custom styling if available
      if (exists("create_select_input", mode = "function")) {
        create_select_input(ns("selected_dataset"),
          label = "Choose a dataset:",
          choices = choices, selected = selected, multiple = FALSE
        )
      } else {
        shiny::selectInput(ns("selected_dataset"),
          label = "Choose a dataset:",
          choices = choices, selected = selected, multiple = FALSE
        )
      }
    })

    # ====== TRAIT SELECTION UPDATE LOGIC ======
    # Updates trait choices when dataset changes (server-side selectize for performance)
    shiny::observeEvent(shiny::req(input$selected_dataset), {
      shiny::req(import()) # Ensure import data is available
      current_ds <- input$selected_dataset

      # RESET TRAIT INPUT: Prevent conflicts during update
      shiny::freezeReactiveValue(input, "which_trait") # Temporarily freeze reactive
      shiny::updateSelectizeInput(session, "which_trait",
        choices = character(0),
        selected = character(0),
        options = list(placeholder = "Loading traits...")
      )

      # GET NEW TRAIT CHOICES: Based on selected dataset
      choices <- get_trait_choices(import(), current_ds) # External helper function
      new_selected_trait <- NULL
      if (!is.null(choices) && length(choices) > 0) {
        new_selected_trait <- choices[1] # Auto-select first trait
      }

      # UPDATE TRAIT DROPDOWN: With server-side processing for large lists
      shiny::updateSelectizeInput(session, "which_trait",
        choices = choices,
        selected = new_selected_trait,
        options = list(
          placeholder = "Search or select trait...",
          maxItems = 1, # Only allow single selection
          maxOptions = 10
        ), # Limit displayed options for performance
        server = TRUE # Enable server-side processing
      )
    })

    # ====== VALIDATED TRAIT REACTIVE ======
    # Ensures selected trait is always valid for current dataset
    current_selected_trait <- shiny::reactive({
      shiny::req(input$selected_dataset) # Wait for dataset selection
      current_ds <- input$selected_dataset
      current_trait_input <- input$which_trait

      # VALIDATE TRAIT: Check if current selection is valid for dataset
      valid_choices <- get_trait_choices(import(), current_ds)
      final_trait_to_return <- NULL

      if (!is.null(valid_choices) && length(valid_choices) > 0) {
        if (!is.null(current_trait_input) && current_trait_input %in% valid_choices) {
          final_trait_to_return <- current_trait_input # Keep current if valid
        } else {
          final_trait_to_return <- valid_choices[1] # Fallback to first valid option
        }
      }

      shiny::req(final_trait_to_return, cancelOutput = TRUE) # Cancel if no valid trait
      return(final_trait_to_return)
    })

    # ====== DEBUG OUTPUT ======
    # Shows current values of all parameters (for development/testing)
    output$returns <- shiny::renderPrint({
      cat(
        "dataset_category =", input$dataset_category,
        "\nselected_dataset =", input$selected_dataset,
        "\nwhich_trait =", input$which_trait,
        "\nselected_chr =", input$selected_chr,
        "\nLOD_thr =", input$LOD_thr
      )
    })

    # ====== RETURN VALUES ======
    # Makes all parameter values available to parent modules
    return(
      list(
        dataset_category = shiny::reactive(input$dataset_category),
        selected_dataset = shiny::reactive(input$selected_dataset),
        which_trait = current_selected_trait, # Use validated reactive (not raw input)
        selected_chr = shiny::reactive(input$selected_chr),
        LOD_thr = shiny::reactive(input$LOD_thr)
      )
    )
  })
}

#' @rdname mainParApp
#' @export
mainParInput <- function(id) {
  # INPUT UI: Creates parameter input controls (dataset selection + LOD threshold)

  # SOURCE STYLING: Load custom UI functions if available
  if (!exists("create_slider_input", mode = "function")) {
    source("R/ui_styles.R")
  }

  ns <- shiny::NS(id) # Create namespaced IDs for module

  if (exists("create_select_input", mode = "function")) {
    # CUSTOM STYLED INPUTS: Use enhanced UI functions
    list(
      shiny::uiOutput(ns("dataset_category")), # Dynamic dropdown (rendered in server)
      shiny::uiOutput(ns("selected_dataset")), # Dynamic dropdown (rendered in server)
      create_slider_input(ns("LOD_thr"), # Static slider with custom styling
        label = "LOD threshold for evaluation",
        min = 4, max = 20, value = 7.5, step = 0.5
      )
    )
  } else {
    # DEFAULT SHINY INPUTS: Fallback to standard Shiny widgets
    list(
      shiny::uiOutput(ns("dataset_category")),
      shiny::uiOutput(ns("selected_dataset")),
      shiny::sliderInput(ns("LOD_thr"), # Standard Shiny slider
        label = "LOD threshold for evaluation",
        min = 4, max = 20, value = 7.5, step = 0.5
      )
    )
  }
}

#' @rdname mainParApp
#' @export
mainParUI <- function(id) {
  # MAIN UI: Creates core selection controls (trait + chromosome)

  # SOURCE STYLING: Load custom UI functions if available
  if (!exists("create_select_input", mode = "function")) {
    source("R/ui_styles.R")
  }

  ns <- shiny::NS(id) # Create namespaced IDs

  if (exists("create_select_input", mode = "function")) {
    # CUSTOM STYLED INPUTS
    list(
      # TRAIT SELECTOR: Searchable dropdown with server-side processing
      create_select_input(ns("which_trait"),
        label = "Choose the trait",
        choices = NULL, # Populated dynamically in server
        selected = NULL,
        multiple = FALSE,
        options = list(
          placeholder = "Search gene symbol...",
          maxItems = 1, # Single selection only
          maxOptions = 7 # Limit displayed options
        )
      ),

      # CHROMOSOME SELECTOR: Static dropdown for genomic regions
      create_select_input(ns("selected_chr"),
        label = "Zoom to Chromosome:",
        choices = c(
          "All", # View all chromosomes
          "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
          "11", "12", "13", "14", "15", "16", "17", "18", "19",
          "X", "Y", "M"
        ), # Standard mouse/human chromosomes
        selected = "All",
        width = "150px"
      )
    )
  } else {
    # DEFAULT SHINY INPUTS: Fallback versions
    list(
      # TRAIT SELECTOR: Standard selectize input
      shiny::selectizeInput(ns("which_trait"),
        label = "Choose the trait",
        choices = NULL,
        multiple = FALSE,
        options = list(
          placeholder = "Search gene symbol...",
          maxItems = 1,
          maxOptions = 7
        )
      ),

      # CHROMOSOME SELECTOR: Standard select input
      shiny::selectInput(ns("selected_chr"),
        label = "Zoom to Chromosome:",
        choices = c(
          "All",
          "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
          "11", "12", "13", "14", "15", "16", "17", "18", "19",
          "X", "Y", "M"
        ),
        selected = "All",
        width = "150px"
      )
    )
  }
}

#' @rdname mainParApp
#' @export
mainParOutput <- function(id) {
  # OUTPUT UI: Creates debug display for current parameter values
  ns <- shiny::NS(id)
  shiny::verbatimTextOutput(ns("returns")) # Text output showing all current selections
}
