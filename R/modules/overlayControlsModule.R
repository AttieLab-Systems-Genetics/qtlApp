#' Overlay Controls Module
#'
#' Renders overlay toggles for additive scans and exposes their states as
#' reactives to be consumed by the scan plot module.
#'
#' @param id Module ID
#' @param selected_dataset_group_reactive Reactive returning current dataset group name
#'
#' @return A list of three reactives: diet_toggle, sex_toggle, sex_diet_toggle
#' @importFrom shiny moduleServer NS uiOutput renderUI checkboxInput reactive
#' @importFrom htmltools div
#' @export
overlayControlsUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("overlay_toggles"))
}

#' @rdname overlayControlsUI
#' @export
overlayControlsServer <- function(id, selected_dataset_group_reactive) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$overlay_toggles <- shiny::renderUI({
      dataset_group <- selected_dataset_group_reactive()
      shiny::req(dataset_group)

      # Only show for HC_HF datasets that have interactive options
      if (!grepl("^HC_HF", dataset_group, ignore.case = TRUE)) {
        return(NULL)
      }

      toggles <- list()

      # Diet overlay availability
      if (any(grepl("Genes|Lipids|Clinical|Metabolites", dataset_group, ignore.case = TRUE))) {
        toggles <- c(toggles, list(
          shiny::checkboxInput(ns("overlay_diet"), "Overlay Diet", FALSE, width = "auto")
        ))
      }

      # Sex overlay availability
      if (any(grepl("Genes|Lipids|Clinical|Metabolites", dataset_group, ignore.case = TRUE))) {
        toggles <- c(toggles, list(
          shiny::checkboxInput(ns("overlay_sex"), "Overlay Sex", FALSE, width = "auto")
        ))
      }

      # Sex x Diet overlay availability
      if (any(grepl("Clinical|Lipid|Metabolites", dataset_group, ignore.case = TRUE))) {
        toggles <- c(toggles, list(
          shiny::checkboxInput(ns("overlay_sex_diet"), "Overlay Sex x Diet", FALSE, width = "auto")
        ))
      }

      if (length(toggles) == 0) return(NULL)

      htmltools::div(
        style = "display: flex; gap: 8px; align-items: center; margin-left: 15px; font-size: 12px; line-height: 1.1; white-space: nowrap; flex-wrap: nowrap;",
        toggles
      )
    })

    return(list(
      diet_toggle = shiny::reactive(isTRUE(input$overlay_diet)),
      sex_toggle = shiny::reactive(isTRUE(input$overlay_sex)),
      sex_diet_toggle = shiny::reactive(isTRUE(input$overlay_sex_diet))
    ))
  })
}


