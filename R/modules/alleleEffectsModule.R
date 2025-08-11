#' Allele Effects Module
#'
#' @param id Module ID
#' @param selected_peak_reactive Reactive containing the full data row for the selected peak
#'
#' @importFrom shiny moduleServer NS renderUI uiOutput reactive reactiveVal observeEvent req renderPlot plotOutput
#' @importFrom htmltools div h5 hr tagList
#' @importFrom DT DTOutput renderDT datatable
#' @importFrom stats setNames
#'
#' @export
alleleEffectsUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::uiOutput(ns("allele_effects_section"))
}

#' @rdname alleleEffectsUI
#' @export
alleleEffectsServer <- function(id, selected_peak_reactive) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # This module now directly uses the selected peak from the scan plot

        # Reactive to get the currently selected peak data
        selected_peak_data <- shiny::reactive({
            peak <- selected_peak_reactive()
            if (is.null(peak)) {
                return(NULL)
            }
            message("alleleEffectsModule: Received selected peak:", peak$marker)
            return(peak)
        })

        # Reactive to prepare allele effects data for plotting
        allele_effects_data_for_plot <- shiny::reactive({
            peak_info <- selected_peak_data()
            if (is.null(peak_info)) {
                return(NULL)
            }

            # Use pivot_peaks helper function to reshape data for plotting
            reshaped_data <- pivot_peaks(peak_info, peak_info$marker)

            if (is.null(reshaped_data) || nrow(reshaped_data) == 0) {
                message(paste("alleleEffectsModule: No allele effects data available for marker:", peak_info$marker))
                return(NULL)
            }

            # Add trait name to the data for plot labeling
            reshaped_data$trait <- peak_info$trait
            message(paste("alleleEffectsModule: Prepared allele effects data for marker:", peak_info$marker))
            reshaped_data
        })

        # Render allele effects section using peakInfoModule
        output$allele_effects_section <- shiny::renderUI({
            htmltools::tagList(
                htmltools::hr(style = "margin: 20px 0; border-top: 2px solid #3498db;"),
                htmltools::div(
                    style = "margin-bottom: 15px;",
                    htmltools::h5("Strain Effects",
                        style = "color: #2c3e50; margin-bottom: 10px; font-weight: bold;"
                    ),
                    htmltools::p("Showing strain effects for the selected peak. Click another peak on the LOD plot to update.",
                        style = "color: #7f8c8d; margin-bottom: 10px; font-size: 12px;"
                    ),
                    htmltools::div(
                        id = ns("peak_summary_info"),
                        style = "background: #f8f9fa; padding: 10px; border-radius: 5px; margin: 10px 0; border-left: 4px solid #3498db;",
                        peakInfoUI(ns("peak_info"))
                    )
                ),
                shinycssloaders::withSpinner(
                    shiny::plotOutput(ns("allele_effects_plot_output"), height = "350px"),
                    type = 8, color = "#3498db"
                )
            )
        })

        # Mount peakInfoServer for detailed info panel
        peakInfoServer(ns("peak_info"), selected_peak_data)

        # Plotting allele effects
        output$allele_effects_plot_output <- shiny::renderPlot({
            plot_data <- allele_effects_data_for_plot()
            if (is.null(plot_data)) {
                return(NULL)
            }

            # Generate plot
            ggplot_alleles(plot_data)
        })
    })
}
