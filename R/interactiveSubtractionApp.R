# R/interactiveSubtractionApp.R

#' Interactive Subtraction Analysis Module
#'
#' This module provides functionality to compute and visualize the difference
#' between interactive and additive QTL scans (interactive - additive).
#' Users can select either sex or diet as the interactive covariate.
#'
#' @param id Module ID
#' @export

library(shiny)
library(dplyr)
library(data.table)

#' Interactive Subtraction Input UI
#'
#' @param id Module ID
#' @export
interactiveSubtractionInput <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::h5("Interactive Analysis", style = "color: #2c3e50; margin-bottom: 10px; font-weight: bold;"),
        shiny::selectInput(
            ns("interaction_type"),
            label = "Interaction Type:",
            choices = c(
                "None (Additive Only)" = "none",
                "Sex Interaction" = "sex",
                "Diet Interaction" = "diet"
            ),
            selected = "none",
            width = "100%"
        ),
        shiny::conditionalPanel(
            condition = "input.interaction_type != 'none'",
            ns = ns,
            shiny::div(
                style = "margin-top: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
                shiny::p("Computing difference: Interactive - Additive",
                    style = "font-size: 12px; color: #6c757d; margin: 0;"
                ),
                shiny::checkboxInput(
                    ns("cache_results"),
                    label = "Cache computed results",
                    value = TRUE
                )
            )
        )
    )
}

#' Interactive Subtraction Server
#'
#' @param id Module ID
#' @param import_reactives Reactive list containing file_directory and annotation_list
#' @param main_par Reactive list with selected_dataset, which_trait, etc.
#' @export
interactiveSubtractionServer <- function(id, import_reactives, main_par) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Cache environment for storing computed differences
        subtraction_cache <- new.env(parent = emptyenv())

        # Return reactive values for use by other modules
        return(list(
            interaction_type = shiny::reactive(input$interaction_type),
            cache_env = subtraction_cache
        ))
    })
}
