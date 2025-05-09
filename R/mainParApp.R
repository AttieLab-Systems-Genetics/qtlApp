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
    # Select `selected_dataset` = `group`.
    output$selected_dataset <- shiny::renderUI({
      shiny::req(import())
      choices <- unique(import()$file_directory$group)
      selected <- choices[1]
      shiny::selectizeInput(ns("selected_dataset"), label = "Choose a dataset",
        choices = choices, selected = selected, multiple = FALSE)
    })
    # Update trait choices.
    shiny::observeEvent(shiny::req(input$selected_dataset), {
     
      choices <- get_trait_choices(import(), input$selected_dataset)
      shiny::updateSelectizeInput(session, "which_trait",
        choices = choices, options = list(maxItems = 1, maxOptions = 5), server = TRUE)
    })
    # Show returned values.
    output$returns <- shiny::renderPrint({
      # Use input values directly for display if needed
      cat("selected_dataset =", input$selected_dataset,
          "\nwhich_trait =", input$which_trait,
          "\nselected_chr =", input$selected_chr,
          "\nLOD_thr =", input$LOD_thr)
    })
    
    # Return reactive expressions for inputs
    return(
      list(
        selected_dataset = shiny::reactive(input$selected_dataset),
        which_trait = shiny::reactive(input$which_trait),
        selected_chr = shiny::reactive(input$selected_chr),
        LOD_thr = shiny::reactive(input$LOD_thr)
        # Add other inputs if they need to be returned reactively
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
  
  # Use modern styling if available, otherwise use standard controls
  if (exists("create_slider_input", mode = "function")) {
    list(
      shiny::uiOutput(ns("selected_dataset")),
      create_slider_input(ns("LOD_thr"),
        label = "LOD threshold for evaluation",
        min = 4, max = 20, value = 7.5, step = 0.5)
    )
  } else {
  list(
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
