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

        # Render allele effects section: show plot only when a peak is selected; otherwise show nothing
        output$allele_effects_section <- shiny::renderUI({
            peak_info <- selected_peak_data()

            if (!is.null(peak_info)) {
                tagList(
                    hr(style = "margin: 20px 0; border-top: 2px solid #3498db;"),
                    shiny::div(
                        style = "display: flex; gap: 10px; align-items: flex-start;",
                        shiny::div(
                            style = "flex: 1 1 40%; min-width: 260px; background: #f8f9fa; padding: 10px; border-radius: 5px;",
                            shiny::uiOutput(ns("peak_info_display"))
                        ),
                        shiny::div(
                            style = "flex: 1 1 60%;",
                            shiny::plotOutput(ns("allele_effects_plot_output"), height = "350px") %>%
                                shinycssloaders::withSpinner(type = 8, color = "#3498db")
                        )
                    )
                )
            } else {
                return(NULL)
            }
        })

        # Info panel content next to the plot
        output$peak_info_display <- shiny::renderUI({
            peak_info <- selected_peak_data()
            if (is.null(peak_info)) {
                return(NULL)
            }

            info_elements <- list()

            # Basic info (prefer gene symbol when available)
            display_trait <- if ("gene_symbol" %in% colnames(peak_info) && !is.null(peak_info$gene_symbol) && nzchar(as.character(peak_info$gene_symbol))) {
                as.character(peak_info$gene_symbol)
            } else if ("phenotype" %in% colnames(peak_info) && !is.null(peak_info$phenotype) && nzchar(as.character(peak_info$phenotype))) {
                as.character(peak_info$phenotype)
            } else {
                as.character(peak_info$trait)
            }
            info_elements$trait <- tags$div(tags$strong("Trait:"), display_trait)
            info_elements$marker <- tags$div(tags$strong("Marker:"), peak_info$marker)
            info_elements$position <- tags$div(tags$strong("Position:"), paste0(peak_info$qtl_chr, ":", round(peak_info$qtl_pos, 2), " Mb"))
            info_elements$lod <- tags$div(tags$strong("LOD:"), round(peak_info$qtl_lod, 2))

            # Cis/Trans status
            if ("cis" %in% colnames(peak_info)) {
                cis_label <- ifelse(isTRUE(peak_info$cis), "Cis", "Trans")
                status_color <- ifelse(isTRUE(peak_info$cis), "#27ae60", "#c0392b")
                info_elements$cis <- tags$div(
                    tags$strong("Status:"),
                    tags$span(cis_label, style = paste("color: white; background-color:", status_color, "; padding: 2px 6px; border-radius: 4px; font-size: 11px;"))
                )
            }

            # Confidence interval
            if ("qtl_ci_lo" %in% colnames(peak_info) && "qtl_ci_hi" %in% colnames(peak_info)) {
                info_elements$ci <- tags$div(tags$strong("CI:"), paste0("[", round(peak_info$qtl_ci_lo, 2), " - ", round(peak_info$qtl_ci_hi, 2), "] Mb"))
            }

            # Founder allele effects A-H mapped to strain names
            allele_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
            strain_names <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
            available_alleles <- allele_cols[allele_cols %in% colnames(peak_info)]
            if (length(available_alleles) > 0) {
                allele_list <- lapply(seq_along(available_alleles), function(i) {
                    col <- available_alleles[i]
                    value <- peak_info[[col]]
                    if (!is.na(value)) paste0(strain_names[i], ": ", round(value, 3))
                })
                allele_list <- Filter(Negate(is.null), allele_list)
                if (length(allele_list) > 0) {
                    info_elements$effects_header <- tags$div(tags$strong("Founder Effects:"), style = "margin-top: 5px;")
                    info_elements$effects <- tags$div(
                        style = "font-family: monospace; font-size: 11px; margin-left: 10px;",
                        lapply(allele_list, tags$div)
                    )
                }
            }

            # Two-column arrangement if needed
            half_len <- ceiling(length(info_elements) / 2)
            col1 <- info_elements[1:half_len]
            col2 <- info_elements[(half_len + 1):length(info_elements)]
            tags$div(
                style = "display: flex; flex-wrap: wrap; gap: 8px;",
                tags$div(style = "flex: 1 1 48%; min-width: 180px;", tagList(col1)),
                tags$div(style = "flex: 1 1 48%; min-width: 180px;", tagList(col2))
            )
        })

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
