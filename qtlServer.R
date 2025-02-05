#' QTL Server
#' 
#' @param id shiny identifier
#' @param file_directory file directory list object
#' @param markers genome markers
#' @param annotation_list annotation list object
#' 
#' @importFrom shiny actionButton helpText moduleserver nearPoints NS observeEvent plotOutput
#'             reactive reactiveValues renderPlot renderPrint renderText req selectizeInput
#'             shinyApp sliderInput tagList updateSelectizeInput verbatimTextOutput
#' @importFrom stringr str_split
#' @importFrom DT datatable DTOutput renderDT
#' @importFrom reshape2 melt
#' @importFrom ggplot2 aes element_blank element_text element_line geom_hline geom_point ggplot
#'             labs scale_color_manual theme theme_bw
#' @importFrom bslib card card_header page_sidebar sidebar
#' @importFrom shinycssloaders withSpinner
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
  shiny::tagList(
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
