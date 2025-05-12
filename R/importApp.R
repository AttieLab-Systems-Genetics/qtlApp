#' Import Import App Module
#'
#' @param id shiny identifier
#' @param file_directory data frame with file directory information
#' @param annotation_list list with annotation information
#' @param set_type character string as reactive object
#'
#' @importFrom DT DTOutput renderDT
#' @importFrom shiny moduleServer NS reactive req selectInput shinyApp
#' @importFrom bslib page_sidebar sidebar
#' @importFrom readxl read_excel
#' @importFrom tools file_ext
#' 
#' @export
importApp <- function() {
  ui <- bslib::page_sidebar(
    title = "Test Import Module",
    sidebar = bslib::sidebar("Choices",
      importUI("import")),
    importOutput("import")
  )
  server <- function(input, output, session) {
    importServer("import")
  }
  shiny::shinyApp(ui = ui, server = server)
}
#' @rdname importApp
#' @export
importServer <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    # Call the updated import_data function to load all necessary data
    loaded_data <- import_data()
    
    # Create caches - Call the (now empty) create_cache function 
    # just in case other code expects it to be called, though it does nothing.
    create_cache()
    
    # Return the loaded data list wrapped in a reactive
    shiny::reactive(loaded_data)
  })
}
#' @rdname importApp
#' @export
importUI <- function(id) {
  # Return empty list as UI is handled by mainApp
  list()
}
#' @rdname importApp
#' @export
importOutput <- function(id) {
  # Return empty list as output is handled by mainApp
  list()
}
