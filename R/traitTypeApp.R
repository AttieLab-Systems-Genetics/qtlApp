#' Trait Type App Module
#'
#' @param id shiny identifier
#' @param import reactive list with file_directory and annotation_list
#' @param trait_type character string as reactive object
#'
#' @importFrom DT DTOutput renderDT
#' @importFrom shiny moduleServer NS reactive renderText req selectInput shinyApp
#'             textOutput
#' @importFrom bslib page
#' 
#' @export
traitTypeApp <- function() {
  ui <- bslib::page(
    title = "Test Trait Type Module",
    mainParInput("main_par"),
    shiny::textOutput("print")
  )
  server <- function(input, output, session) {
    import <- importServer("import")
    main_par <- mainParServer("main_par", import)
    trait_type <- traitTypeServer("trait_type", main_par, import)
    output$print <- shiny::renderText({
      shiny::req(trait_type())
    })
  }
  shiny::shinyApp(ui = ui, server = server)
}
#' @rdname traitTypeApp
#' @export
traitTypeServer <- function(id, main_par, import) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Return
    shiny::reactive({
      get_trait_type(shiny::req(import()), main_par$group)
    })
  })
}
