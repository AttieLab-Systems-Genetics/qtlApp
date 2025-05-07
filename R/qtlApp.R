#' QTL App
#' 
#' @param id shiny identifier
#' @param import reactive list with file_directory, annotation_list and markers
#' 
#' @importFrom shiny helpText moduleServer NS shinyApp tags
#' @importFrom bslib page_sidebar sidebar bs_theme
#' @importFrom shinyjs useShinyjs
#' 
#' @export
qtlApp <- function() {
  # Source UI styling functions if not already loaded
  if (!exists("create_fluid_page", mode = "function")) {
    source("R/ui_styles.R")
  }
  
  # Use modern UI if styling functions are available
  if (exists("create_fluid_page", mode = "function") && 
      exists("create_title_panel", mode = "function")) {
    
    ui <- create_fluid_page(
      shinyjs::useShinyjs(),
      
      # Add custom CSS styling
      tags$head(tags$style(custom_css)),
      
      # Modern title panel
      create_title_panel(
        "Pre-scanned QTL Visualizer for Diet DO Study",
        "Interactive visualization tool for QTL analysis"
      ),
      
      # Main content
      create_fluid_row(
        # Left sidebar
        create_column(3,
          create_well_panel(
            qtlInput("qtl")
          )
        ),
        
        # Main content area
        create_column(9,
          qtlOutput("qtl")
        )
      )
    )
    
  } else {
    # Fallback to default styling if modern UI functions aren't available
  ui <- bslib::page_sidebar(
    title = "Pre-scanned QTL visualizer, implemented for Diet DO study",
    shinyjs::useShinyjs(),
    sidebar = bslib::sidebar("side_panel",
      qtlInput("qtl")),
    qtlOutput("qtl")
  )
  }
  
  server <- function(input, output, session) {
    qtlServer("qtl")
  }
  
  shiny::shinyApp(ui = ui, server = server)
}
#' @rdname qtlApp
#' @export
qtlServer <- function(id) {
    shiny::moduleServer(id, function(input, output, session) {
      ns <- session$ns
      import <- importServer("import")
      main_par <- mainParServer("main_par", import)
      scan_list <- scanServer("scan_list", main_par, import)
      peak_table <- peakServer("peak_table", main_par, import)
      scanlyServer("scanly", main_par, scan_list, peak_table)
  })
}
#' @rdname qtlApp
#' @export
qtlInput <- function(id) {
  ns <- shiny::NS(id)
  list(
    shiny::helpText("Select your dataset, trait to show, and other options"),
    mainParInput(ns("main_par")), # "group", "LOD_thr"
    mainParUI(ns("main_par")),    # "which_trait", "selected_chr"
    peakInput(ns("peak_table")))  # "which_peak", "alleles" actionButton
}
#' @rdname qtlApp
#' @export
qtlOutput <- function(id) {
  ns <- shiny::NS(id)
  list(
    scanlyOutput(ns("scanly")),
    scanlyUI(ns("scanly")),
    peakOutput(ns("peak_table"))
  )
}
