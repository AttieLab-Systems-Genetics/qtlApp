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
    ns <- session$ns
    #import <- read.csv(system.file("data/import.csv", package = "qtlApp"))
    import <- read.csv("data/import.csv")
    out <- list()
    for(i in 1:nrow(import)) {
      ext <- tools::file_ext(import$filename[i])
      switch(ext,
        csv =   out[[import$object[i]]] <- read.csv(import$filename[i]),
        rds =   out[[import$object[i]]] <- readRDS(import$filename[i]),
        xlsx =  out[[import$object[i]]] <- readxl::read_excel(import$filename[i])
      )
    }
    out$file_directory$group <- paste0(out$file_directory$diet, " ", 
      out$file_directory$trait_compartment, " ",
      out$file_directory$trait_type, ", ", 
      out$file_directory$scan_type)

    output$import_data <- shiny::renderUI({
      shiny::selectInput(ns("import_data"),
        label = "Choose a filename",
        choices = names(out))
      })
    annots <- names(out$annotation_list)
    output$annots <- shiny::renderUI({
      if(shiny::req(input$import_data) == 'annotation_list') {
        shiny::selectInput(ns("annotation_list"),
          label = "Choose a dataset",
          choices = annots)
        }
    })
    output$show_data <- DT::renderDT({
        object <- shiny::req(input$import_data)
        if(object == 'annotation_list') {
          df <- out[[object]][[shiny::req(input$annotation_list)]]
        } else {
          df <- out[[object]]
        }
        df
      }, 
      options = list(paging = TRUE,    ## paginate the output
        pageLength = 5,  ## number of rows to output for each page
        lengthMenu = c(5, 10, 25, 50, 100), ## number of rows to display
        autoWidth = TRUE,
        scrollX = TRUE,
        scrollY = "400px",
        searching = FALSE,
        ordering = FALSE,
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      )
    )
    # Return
    shiny::reactive(out)
  })
}
#' @rdname importApp
#' @export
importUI <- function(id) {
    ns <- shiny::NS(id)
    list(
        shiny::uiOutput(ns("import_data")),
        shiny::uiOutput(ns("annots")))
}
#' @rdname importApp
#' @export
importOutput <- function(id) {
    ns <- shiny::NS(id)
    DT::DTOutput(ns('show_data'))
}
