#' Main Parameter App Module
#' 
#' @param id shiny identifier
#' @param import reactive list with file_directory and annotation_list
#'
#' @importFrom shiny moduleServer reactive renderPrint req selectizeInput
#'             shinyApp updateSelectizeInput verbatimTextOutput
#' @importFrom bslib page
#' 
#' @export
mainParApp <- function(id) {
  ui <- bslib::page(
    title = "Test MainPar Module",
    mainParInput("main_par"),
    mainParUI("main_par"),
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

    # Select `group` = `selected_dataset`.
    output$group <- shiny::renderUI({
      shiny::req(import())
      choices <- import()$file_directory$group
      shiny::selectizeInput(ns("group"), label = "Choose a data group",
        choices = choices, selected = choices[1], multiple = FALSE)
    })
    trait_type <- shiny::reactive({
      get_trait_type(shiny::req(import()), input$group)
    })
    trait_list <- shiny::reactive({
      shiny::req(trait_type())
      get_trait_list(import(), trait_type())
    })
    # ** Want to join symbol with either gene.id or transcript.id **
    trait_id <- shiny::reactive({
      switch(shiny::req(trait_type()), # ** Expand for new data types ** #
        genes    = "gene.id",
        isoforms = "transcript.id",
        clinical = "data_name")
    })

    # Update trait choices.
    shiny::observeEvent(shiny::req(trait_type(), trait_list(), trait_id()), {
      choices <- paste0("(", trait_list()[[trait_id()]], ")")
      annotation_list <- shiny::req(import())$annotation_list
      if(trait_type() == "clinical") {
        choices <- paste(annotation_list$clinical$ID, choices)
      } else {
        choices <- paste(annotation_list$genes$symbol, choices)
      }
      shiny::updateSelectizeInput(session, "which_trait",
        choices = choices, options = list(maxItems = 1, maxOptions = 5), server = TRUE)
    })

    output$returns <- shiny::renderPrint({
      cat("group =", input$group,
          "\nwhich_trait =", input$which_trait)
    })
    # Return.
    input
  })
}
#' @rdname mainParApp
#' @export
mainParInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("group"))
}
#' @rdname mainParApp
#' @export
mainParUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::selectizeInput(ns("which_trait"),
    label = "Choose the trait",
    choices = NULL,
    multiple = FALSE,
    options = list(placeholder = 'Search...'))
}
#' @rdname mainParApp
#' @export
mainParOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::verbatimTextOutput(ns("returns"))
}
