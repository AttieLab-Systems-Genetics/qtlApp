#' Trait App Module
#'
#' @param id shiny identifier
#' @param import reactive list with file_directory and annotation_list
#' @param trait_type character string as reactive object
#'
#' @importFrom DT DTOutput renderDT
#' @importFrom shiny moduleServer NS reactive req selectInput shinyApp
#' @importFrom bslib card card_header page_sidebar sidebar
#' 
#' @export
traitApp <- function() {
    ui <- bslib::page_sidebar(
        title = "Test Trait List",
        sidebar = bslib::sidebar("side_panel",
            shiny::uiOutput("group")),
        traitUI("trait"),
        traitOutput("trait")
    )
    server <- function(input, output, session) {
        import <- importServer("import")
        traitServer("trait", import, trait_type)

        # Choose a dataset
        output$group <- shiny::renderUI({
            shiny::req(import())
            file_directory <- import()$file_directory
            shiny::selectInput("group",
                label = "Choose a dataset to display",
                choices = unique(file_directory$group))
        })
        # Display data for all of the options
        trait_type <- shiny::reactive({
            get_trait_type(import(), input$group)
        })
    }
    shiny::shinyApp(ui = ui, server = server)
}
#' @rdname traitApp
#' @export
traitServer <- function(id, import, trait_type) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Display data for all of the options
        output$available_data <- DT::renderDT({
            shiny::req(import())$file_directory
        })
        # Display trait list
        output$trait_list <- DT::renderDT({
            shiny::req(trait_list())
        })
        # Get trait list
        trait_list <- shiny::reactive({
            shiny::req(trait_type())
            get_trait_list(import(), trait_type())
        })
        output$trait_show <- shiny::renderUI({
            shiny::req(trait_list(), trait_type())
            bslib::card(
                bslib::card_header(paste("Trait list", trait_type())),
                DT::DTOutput(ns('trait_list')))
        })

        # Return
        trait_list
    })
}
#' @rdname traitApp
#' @export
traitUI <- function(id) {
    ns <- shiny::NS(id)
    bslib::card(
        bslib::card_header("Datasets available"),
        DT::DTOutput(ns('available_data'))
    )
}
#' @rdname traitApp
#' @export
traitOutput <- function(id) {
    ns <- shiny::NS(id)
    shiny::uiOutput(ns("trait_show"))
}
