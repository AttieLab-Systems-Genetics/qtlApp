#' Trait module
#'
#' @param id shiny identifier
#' @param file_directory data frame with file directory information
#' @param annotation_list list with annotation information
#' @param set_type character string as reactive object
#'
#' @importFrom DT DTOutput renderDT
#' @importFrom shiny moduleServer NS reactive req selectInput shinyApp
#' @importFrom bslib card card_header page_sidebar sidebar
#' @export
traitServer <- function(id, file_directory, annotation_list, set_type) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Display data for all of the options--------------------------------------------------------
        output$available_data <- DT::renderDT({
            file_directory
        })
        # Display trait list
        output$trait_list <- DT::renderDT({
            shiny::req(trait_list())
        })
        
        trait_list <- shiny::reactive({
            tlist <- switch(shiny::req(set_type()),
                Genes    = data.frame(annotation_list$genes),
                Isoforms = data.frame(annotation_list$isoforms),
                Clinical = data.frame(annotation_list$clinical))
            assign("test_annot", tlist, envir = .GlobalEnv)
            tlist
        })
        # Return
        trait_list
    })
}
traitUI <- function(id) {
    ns <- shiny::NS(id)
    bslib::card(
        bslib::card_header("Datasets available"),
        DT::DTOutput(ns('available_data'))
    )
}
traitOutput <- function(id) {
    ns <- shiny::NS(id)
    bslib::card(
        bslib::card_header("Trait list"),
        DT::DTOutput(ns('trait_list'))
    )
}
#' @export
#' @rdname traitServer
traitApp <- function() {
    source("qtlSetup.R")
    ui <- bslib::page_sidebar(
        title = "Test Trait List",
        sidebar = bslib::sidebar("side_panel",
            shiny::selectInput("selected_dataset",
                label = "Choose a dataset to display",
                choices = unique(file_directory$group))),
        traitUI("trait"),
        traitOutput("trait")
    )
    server <- function(input, output, session) {
        set_type <- shiny::reactive({
            shiny::req(input$selected_dataset)
            type <- subset(file_directory, group==input$selected_dataset)$trait_type[1]
            assign("set_test", type, envir = .GlobalEnv)
            type
        })
        traitServer("trait", file_directory, annotation_list, set_type)
    }
    shiny::shinyApp(ui = ui, server = server)
}