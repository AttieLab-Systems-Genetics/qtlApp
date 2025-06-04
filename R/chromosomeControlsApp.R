#' Chromosome Controls Module
#'
#' @param id shiny identifier
#' @export
chromosomeControlsApp <- function() {
    ui <- bslib::page_sidebar(
        title = "Test Chromosome Controls",
        sidebar = bslib::sidebar(
            "side_panel",
            chromosomeControlsInput("chr_controls")
        ),
        chromosomeControlsOutput("chr_controls")
    )
    server <- function(input, output, session) {
        chromosomeControlsServer("chr_controls")
    }
    shiny::shinyApp(ui = ui, server = server)
}

#' @rdname chromosomeControlsApp
#' @export
chromosomeControlsServer <- function(id) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Observer for "Zoom to Chromosome" button
        shiny::observeEvent(input$zoom_to_chr, {
            selected_chr <- input$selected_chr
            if (!is.null(selected_chr) && selected_chr != "All") {
                message(paste("Chromosome Controls: Zooming to chromosome:", selected_chr))
                shiny::showNotification(
                    paste("Zoomed to chromosome", selected_chr),
                    type = "message",
                    duration = 2
                )
            } else {
                shiny::showNotification(
                    "Please select a specific chromosome to zoom to",
                    type = "warning",
                    duration = 3
                )
            }
        })

        # Observer for "Show All Chromosomes" button
        shiny::observeEvent(input$reset_chr_view, {
            message("Chromosome Controls: Resetting to show all chromosomes")
            shiny::updateSelectInput(session, "selected_chr", selected = "All")
            shiny::showNotification(
                "Showing all chromosomes",
                type = "message",
                duration = 2
            )
        })

        # Return the selected chromosome reactive
        return(list(
            selected_chr = shiny::reactive(input$selected_chr)
        ))
    })
}

#' @rdname chromosomeControlsApp
#' @export
chromosomeControlsInput <- function(id) {
    ns <- shiny::NS(id)

    div(
        style = "margin-bottom: 15px; background: #f8f9fa; padding: 10px; border-radius: 6px; border: 1px solid #bdc3c7;",
        div(
            style = "display: flex; align-items: center; gap: 15px; flex-wrap: wrap;",
            div(
                style = "flex: 1; min-width: 200px;",
                h6("🔍 Chromosome View", style = "color: #2c3e50; margin: 0 0 8px 0; font-weight: bold;"),
                shiny::selectInput(
                    ns("selected_chr"),
                    label = NULL,
                    choices = c(
                        "All" = "All",
                        setNames(as.character(1:19), paste("Chr", 1:19)),
                        "X" = "X", "Y" = "Y", "M" = "M"
                    ),
                    selected = "All",
                    width = "100%"
                )
            ),
            div(
                style = "display: flex; align-items: end; gap: 10px;",
                shiny::actionButton(
                    ns("zoom_to_chr"),
                    "🔍 Zoom to Chromosome",
                    class = "btn btn-primary",
                    style = "background: #3498db; border: none; color: white; font-weight: bold;"
                ),
                shiny::actionButton(
                    ns("reset_chr_view"),
                    "🌐 Show All",
                    class = "btn btn-secondary",
                    style = "background: #7f8c8d; border: none; color: white;"
                )
            )
        )
    )
}

#' @rdname chromosomeControlsApp
#' @export
chromosomeControlsOutput <- function(id) {
    ns <- shiny::NS(id)
    div(
        h5("Current Chromosome Selection:"),
        shiny::verbatimTextOutput(ns("debug_output"))
    )
}
