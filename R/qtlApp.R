#' QTL App
#' 
#' @param id shiny identifier
#' @param import reactive list with file_directory, annotation_list and markers
#' 
#' @importFrom shiny helpText moduleServer NS shinyApp
#' @importFrom bslib page_sidebar sidebar
#' @importFrom shinyjs useShinyjs
#' 
#' @export
qtlApp <- function() {
  ui <- bslib::page_sidebar(
    title = "Pre-scanned QTL visualizer, implemented for Diet DO study",
    shinyjs::useShinyjs(),
    sidebar = bslib::sidebar("side_panel",
      qtlInput("qtl")),
    qtlOutput("qtl")
  )
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
      peak_list <- peakServer("peak_list", main_par, import)
      scanlyServer("scanly", main_par, scan_list, peak_list)
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
    peakInput(ns("peak_list")))  # "which_peak", "alleles" actionButton
}
#' @rdname qtlApp
#' @export
qtlOutput <- function(id) {
  ns <- shiny::NS(id)
  list(
    scanlyOutput(ns("scanly")),
    scanlyUI(ns("scanly")),
    peakUI(ns("peak_list")),
    peakOutput(ns("peak_list"))
  )
}
