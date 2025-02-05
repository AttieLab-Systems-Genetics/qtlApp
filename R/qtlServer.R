#' QTL Server
#' 
#' @param id shiny identifier
#' @param file_directory data frame with file directory information
#' @param annotation_list list with annotation information
#' @param markers list object with marker information
#' 
#' @importFrom shiny helpText moduleServer NS shinyApp
#' @importFrom bslib page_sidebar sidebar
#' @importFrom shinyjs useShinyjs
#' @export
qtlServer <- function(id, file_directory, annotation_list, markers) {
    shiny::moduleServer(id, function(input, output, session) {
      ns <- session$ns
      main_par <- mainParServer("main_par", file_directory, annotation_list)
      scanServer("scan", main_par, file_directory, annotation_list, markers)
      peakServer("peak", main_par, file_directory)
  })
}
qtlInput <- function(id) {
  ns <- shiny::NS(id)
  list(
    shiny::helpText("Select your dataset, trait to show, and other options"),
    mainParInput(ns("main_par")), # "selected_dataset"
    scanInput(ns("scan")),        # "LOD_thr"
    mainParUI(ns("main_par")),    # "which_trait"
    scanUI(ns("scan")),           # "scan" actionButton
    peakInput(ns("peak")))        # "which_peak", "alleles" actionButton
}
#' @export
#' @rdname qtlServer
qtlOutput <- function(id) {
  ns <- shiny::NS(id)
  list(
    scanOutput(ns("scan")),
    peakOutput(ns("peak"))
  )
}
#' @export
#' @rdname qtlServer
qtlApp <- function() {
  source("qtlSetup.R")
  source("mainParServer.R")
  source("traitServer.R")
  source("scanServer.R")
  source("peakServer.R")
  ui <- bslib::page_sidebar(
    title = "Pre-scanned QTL visualizer, implemented for Diet DO study",
    shinyjs::useShinyjs(),
    sidebar = bslib::sidebar("side_panel",
      qtlInput("qtl")),
    qtlOutput("qtl")
  )
  server <- function(input, output, session) {
    qtlServer("qtl", file_directory, annotation_list, markers)
  }
  shiny::shinyApp(ui = ui, server = server)
}
