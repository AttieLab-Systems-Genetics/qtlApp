#' Download App
#'
#' @param id identifier for shiny reactive
#' @param prefix static prefix for filename
#' @param main_par input parameters from calling routine
#' @param download_list reactiveValues with  postfix,plot,table
#' @return nothing 
#'
#' @importFrom shiny column downloadButton downloadHandler fluidRow
#'             moduleServer NS radioButtons reactive renderUI req
#'             textAreaInput uiOutput
#' @importFrom utils write.csv    
#' @importFrom grDevices dev.off pdf
#' @export
downloadApp <- function(id) {
  ui <- bslib::page_sidebar(
    title = "Test Download",
    sidebar = bslib::sidebar("side_panel",
      mainParInput("main_par"), # "selected_dataset", "LOD_thr"
      mainParUI("main_par"),    # "which_trait", "selected_chr"
      downloadInput("download")), # plot_table
    bslib::card(
      bslib::card_header("LOD profile"),
      scanOutput("scan_list")),
    bslib::card(
      bslib::card_header("Download Info"),
      downloadOutput("download"))
  )
  server <- function(input, output, session) { 
    import <- importServer("import")
    main_par <- mainParServer("main_par", import)
    scan_list <- scanServer("scan_list", main_par, import)
    downloadServer("download", scan_list)
  }
  shiny::shinyApp(ui, server)
}
#' @rdname downloadApp
#' @export
downloadServer <- function(id, download_list) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    output$downloads <- shiny::renderUI({
      plot_table <- shiny::req(input$plot_table)
      shiny::downloadButton(ns(plot_table), plot_table)
    })
    output$filename <- renderUI({
      filename <- paste0(shiny::req(download_list$panel()), "_",
                         shiny::req(download_list$postfix()))
      shiny::textAreaInput(ns("filename"), "File Prefix:", filename)
    })
    output$Plots <- shiny::downloadHandler(
      filename = function() paste0(shiny::req(input$filename), ".pdf"),
      content = function(file) {
        # ** this is not working even thought plot seems fine **
        grDevices::pdf(file, width = 9,
                       height = shiny::req(download_list$height()))
        print(shiny::req(download_list$plot()))
        grDevices::dev.off()
      },
      contentType = "application/pdf")
    output$Tables <- shiny::downloadHandler(
      filename = function() paste0(shiny::req(input$filename), ".csv"),
      content = function(file) {
        table <- shiny::req(download_list$table())
        utils::write.csv(table, file, row.names = FALSE)
      })
  })
}
#' @rdname downloadApp
#' @export
downloadInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::radioButtons(ns("plot_table"), "", c("Plots","Tables"), "Plots",
    inline = TRUE)
}
#' @rdname downloadApp
#' @export
downloadOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h5("Download:"),
    shiny::fluidRow(
      shiny::column(3, shiny::uiOutput(ns("downloads"))),
      shiny::column(9, shiny::uiOutput(ns("filename")))))
}
