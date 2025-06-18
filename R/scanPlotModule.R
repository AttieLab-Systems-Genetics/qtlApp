#' Scan Plot Module
#'
#' @param id Module ID
#' @param trait_to_scan Reactive for the trait to scan
#' @param selected_dataset_group Reactive for the selected dataset
#' @param import_reactives Reactive containing file_directory and other import data
#' @param main_par_inputs Reactive containing main parameters (LOD_thr, selected_chr, etc.)
#' @param interaction_type_reactive Reactive for interaction type
#'
#' @importFrom shiny moduleServer NS renderUI uiOutput reactive reactiveVal observeEvent req
#' @importFrom plotly plotlyOutput renderPlotly ggplotly layout config event_data event_register
#' @importFrom shinycssloaders withSpinner
#' @importFrom htmltools div h5 tagList
#' @importFrom ggplot2 labs theme element_text
#'
#' @export
scanPlotUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::uiOutput(ns("scan_plot_ui_render"))
}

#' @rdname scanPlotUI
#' @export
scanPlotServer <- function(id, trait_to_scan, selected_dataset_group, import_reactives, main_par_inputs, interaction_type_reactive = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Create local cache for peaks data
        local_peaks_cache <- new.env(parent = emptyenv())

        plot_width_rv <- shiny::reactiveVal(1200)
        plot_height_rv <- shiny::reactiveVal(600)
        clicked_plotly_point_details_lod_scan_rv <- shiny::reactiveVal(NULL)

        # Define the %||% operator for null coalescing
        `%||%` <- function(a, b) if (!is.null(a)) a else b

        current_trait_for_scan <- shiny::reactive({
            shiny::req(trait_to_scan())
            trait_val <- trait_to_scan()
            message(paste("scanPlotModule: Received trait_to_scan value:", trait_val, "(this will be passed to trait_scan)."))
            return(trait_val)
        })

        # Source the existing scanServer function from the monolithic backup
        # and use it as the core implementation
        scan_module_outputs <- scanServer(
            id = "scan_plot_internal",
            trait_to_scan = trait_to_scan,
            selected_dataset_group = selected_dataset_group,
            import_reactives = import_reactives,
            main_par_inputs = main_par_inputs,
            interaction_type_reactive = interaction_type_reactive
        )

        # Expose the UI from the internal scan server
        output$scan_plot_ui_render <- shiny::renderUI({
            scanOutput("scan_plot_internal")
        })

        # Return the outputs from the internal scan server
        return(scan_module_outputs)
    })
}
