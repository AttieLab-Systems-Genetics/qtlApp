#' Allele Effects Module
#'
#' @param id Module ID
#' @param trait_for_lod_scan_reactive Reactive for the current trait being analyzed
#' @param selected_dataset_reactive Reactive for the currently selected dataset
#' @param import_reactives Reactive containing file_directory and other import data
#' @param lod_threshold_reactive Reactive for the LOD threshold
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
alleleEffectsServer <- function(id, trait_for_lod_scan_reactive, selected_dataset_reactive, import_reactives, lod_threshold_reactive) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Cache for peaks data
        peaks_cache <- new.env(parent = emptyenv())

        # Reactive value to store the selected peak marker from dropdown
        selected_peak_from_dropdown <- shiny::reactiveVal(NULL)

        # Reactive to find peak data for the selected trait
        peaks_data_for_trait <- shiny::reactive({
            trait_val <- trait_for_lod_scan_reactive()
            if (is.null(trait_val)) {
                return(NULL)
            }

            shiny::req(selected_dataset_reactive(), import_reactives())

            dataset_val <- selected_dataset_reactive()

            message(paste("alleleEffectsModule: Finding peaks for trait:", trait_val, "in dataset:", dataset_val))

            # Determine trait type for this dataset to pass to peak_finder
            # Use the base dataset name for trait type determination
            base_dataset <- selected_dataset_reactive()
            trait_type_val <- get_trait_type(import_reactives(), base_dataset)

            # Use peak_finder to get peaks for this specific trait
            tryCatch(
                {
                    peaks <- peak_finder(
                        file_dir = import_reactives()$file_directory,
                        selected_dataset = dataset_val,
                        selected_trait = trait_val,
                        trait_type = trait_type_val,
                        cache_env = peaks_cache,
                        use_cache = TRUE
                    )
                    message(paste("alleleEffectsModule: Found", if (is.null(peaks)) 0 else nrow(peaks), "peaks for trait:", trait_val))
                    peaks
                },
                error = function(e) {
                    message(paste("alleleEffectsModule: Error finding peaks for trait", trait_val, ":", e$message))
                    NULL
                }
            )
        })

        # Reactive to get ALL peaks above threshold for the selected trait
        available_peaks_for_trait <- shiny::reactive({
            peaks <- peaks_data_for_trait()
            if (is.null(peaks) || nrow(peaks) == 0) {
                return(NULL)
            }

            # Filter by LOD threshold
            lod_thr <- lod_threshold_reactive()
            shiny::req(lod_thr)

            # Check what columns are available
            message(paste("alleleEffectsModule: Peak data columns:", paste(colnames(peaks), collapse = ", ")))

            # Use qtl_lod column (from peak_finder) instead of lod
            if (!"qtl_lod" %in% colnames(peaks)) {
                message("alleleEffectsModule: No qtl_lod column found in peaks data")
                return(NULL)
            }

            # Filter by LOD threshold and sort by LOD descending
            filtered_peaks <- peaks[peaks$qtl_lod >= lod_thr, ]
            if (nrow(filtered_peaks) == 0) {
                message(paste("alleleEffectsModule: No peaks above LOD threshold", lod_thr, "for trait:", trait_for_lod_scan_reactive()))
                return(NULL)
            }

            # Sort by qtl_lod descending
            filtered_peaks <- filtered_peaks[order(-filtered_peaks$qtl_lod), ]

            message(paste("alleleEffectsModule: Found", nrow(filtered_peaks), "peaks above threshold for trait:", trait_for_lod_scan_reactive()))
            filtered_peaks
        })

        # Reactive to get the currently selected peak marker (either from dropdown or default to highest)
        selected_peak_marker <- shiny::reactive({
            available_peaks <- available_peaks_for_trait()
            if (is.null(available_peaks) || nrow(available_peaks) == 0) {
                return(NULL)
            }

            # Use dropdown selection if available, otherwise default to highest peak
            dropdown_selection <- selected_peak_from_dropdown()
            if (!is.null(dropdown_selection) && dropdown_selection %in% available_peaks$marker) {
                selected_marker <- dropdown_selection
                selected_peak_info <- available_peaks[available_peaks$marker == selected_marker, ]
                message(paste("alleleEffectsModule: Using dropdown-selected peak:", selected_marker, "with LOD:", selected_peak_info$qtl_lod[1]))
            } else {
                # Default to highest peak
                selected_marker <- available_peaks$marker[1]
                message(paste("alleleEffectsModule: Using highest peak (default):", selected_marker, "with LOD:", available_peaks$qtl_lod[1]))
            }

            selected_marker
        })

        # Reactive to prepare allele effects data
        allele_effects_data <- shiny::reactive({
            peaks <- peaks_data_for_trait()
            marker <- selected_peak_marker()

            if (is.null(peaks) || is.null(marker)) {
                return(NULL)
            }

            # Use pivot_peaks helper function to reshape data for plotting
            reshaped_data <- pivot_peaks(peaks, marker)

            if (is.null(reshaped_data) || nrow(reshaped_data) == 0) {
                message(paste("alleleEffectsModule: No allele effects data available for marker:", marker))
                return(NULL)
            }

            # Add trait name to the data for plot labeling
            reshaped_data$trait <- trait_for_lod_scan_reactive()
            message(paste("alleleEffectsModule: Prepared allele effects data for marker:", marker, "with", nrow(reshaped_data), "strain effects"))
            reshaped_data
        })

        # Observer to update selected peak when dropdown changes
        shiny::observeEvent(input$peak_selection_dropdown,
            {
                new_selection <- input$peak_selection_dropdown
                if (!is.null(new_selection)) {
                    selected_peak_from_dropdown(new_selection)
                    message(paste("alleleEffectsModule: Peak dropdown selection changed to:", new_selection))
                }
            },
            ignoreNULL = FALSE
        )

        # Observer to reset peak selection when trait changes
        shiny::observeEvent(trait_for_lod_scan_reactive(),
            {
                selected_peak_from_dropdown(NULL) # Reset selection when trait changes
                message("alleleEffectsModule: Reset peak selection due to trait change")
            },
            ignoreNULL = FALSE
        )

        # Render allele effects section conditionally
        output$allele_effects_section <- shiny::renderUI({
            # Check if we have allele effects data
            available_peaks <- available_peaks_for_trait()

            message(paste("alleleEffectsModule: Checking allele effects section. Available peaks:", !is.null(available_peaks)))

            if (!is.null(available_peaks) && nrow(available_peaks) > 0) {
                # Create choices for dropdown - just show marker names
                peak_choices <- setNames(
                    available_peaks$marker,
                    available_peaks$marker
                )

                # Set default selection to highest peak if not already set
                current_selection <- selected_peak_from_dropdown()
                if (is.null(current_selection) || !current_selection %in% available_peaks$marker) {
                    default_selection <- available_peaks$marker[1] # Highest peak
                } else {
                    default_selection <- current_selection
                }

                # Get current peak info for display
                current_peak_info <- NULL
                if (!is.null(default_selection)) {
                    current_peak_info <- available_peaks[available_peaks$marker == default_selection, ][1, ]
                }

                message(paste("alleleEffectsModule: Rendering allele effects section with", length(peak_choices), "peak choices"))

                tagList(
                    hr(style = "margin: 20px 0; border-top: 2px solid #3498db;"),
                    div(
                        style = "margin-bottom: 15px;",
                        h5("Strain Effects",
                            style = "color: #2c3e50; margin-bottom: 10px; font-weight: bold;"
                        ),
                        p(paste("Select a peak to view strain effects (", length(peak_choices), "peaks found):"),
                            style = "color: #7f8c8d; margin-bottom: 10px; font-size: 12px;"
                        ),
                        shiny::selectInput(
                            ns("peak_selection_dropdown"),
                            label = NULL,
                            choices = peak_choices,
                            selected = default_selection,
                            width = "100%"
                        ),
                        # Add peak summary info
                        if (!is.null(current_peak_info)) {
                            div(
                                id = ns("peak_summary_info"),
                                style = "background: #f8f9fa; padding: 10px; border-radius: 5px; margin: 10px 0; border-left: 4px solid #3498db;",
                                shiny::uiOutput(ns("peak_info_display"))
                            )
                        } else {
                            NULL
                        }
                    ),
                    shiny::plotOutput(ns("allele_effects_plot_output"), height = "350px") %>%
                        shinycssloaders::withSpinner(type = 8, color = "#3498db")
                )
            } else {
                message("alleleEffectsModule: No allele effects section - no peaks above threshold")
                NULL
            }
        })

        # Dynamic peak info display
        output$peak_info_display <- shiny::renderUI({
            current_peaks <- available_peaks_for_trait()
            selected_marker <- selected_peak_marker()

            if (is.null(current_peaks) || is.null(selected_marker)) {
                return(NULL)
            }

            peak_info <- current_peaks[current_peaks$marker == selected_marker, ][1, ]
            if (nrow(peak_info) == 0) {
                return(NULL)
            }

            # Build summary info
            info_elements <- list()

            # Basic info
            info_elements <- c(info_elements, list(
                tags$strong("Marker: "), peak_info$marker, tags$br(),
                tags$strong("Position: "), paste0("Chr", peak_info$qtl_chr, ":", round(peak_info$qtl_pos, 2), " Mb"), tags$br(),
                tags$strong("LOD Score: "), round(peak_info$qtl_lod, 3), tags$br()
            ))

            # Cis/Trans status
            if ("cis" %in% colnames(peak_info)) {
                cis_status <- if (is.logical(peak_info$cis)) {
                    ifelse(peak_info$cis, "Cis", "Trans")
                } else if (is.character(peak_info$cis)) {
                    ifelse(toupper(peak_info$cis) %in% c("TRUE", "1", "YES"), "Cis", "Trans")
                } else {
                    "Unknown"
                }

                cis_color <- if (cis_status == "Cis") "#27ae60" else "#e74c3c"
                info_elements <- c(info_elements, list(
                    tags$strong("Type: "),
                    tags$span(cis_status, style = paste0("color: ", cis_color, "; font-weight: bold;")),
                    tags$br()
                ))
            }

            # Confidence interval
            if ("qtl_ci_lo" %in% colnames(peak_info) && "qtl_ci_hi" %in% colnames(peak_info)) {
                if (!is.na(peak_info$qtl_ci_lo) && !is.na(peak_info$qtl_ci_hi)) {
                    info_elements <- c(info_elements, list(
                        tags$strong("95% CI: "),
                        paste0(round(peak_info$qtl_ci_lo, 2), " - ", round(peak_info$qtl_ci_hi, 2), " Mb"),
                        tags$br()
                    ))
                }
            }

            # Add founder allele effects summary
            allele_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
            available_alleles <- allele_cols[allele_cols %in% colnames(peak_info)]

            if (length(available_alleles) > 0) {
                non_na_alleles <- available_alleles[!is.na(peak_info[available_alleles])]
                if (length(non_na_alleles) > 0) {
                    info_elements <- c(info_elements, list(
                        tags$strong("Founder Effects: "),
                        paste0(length(non_na_alleles), " available (", paste(non_na_alleles, collapse = ", "), ")"),
                        tags$br()
                    ))
                }
            }

            do.call(tagList, info_elements)
        })

        # Render the allele effects plot
        output$allele_effects_plot_output <- shiny::renderPlot({
            effects_data <- allele_effects_data()

            message(paste("alleleEffectsModule: Rendering allele effects plot. Data available:", !is.null(effects_data)))

            if (is.null(effects_data)) {
                message("alleleEffectsModule: No allele effects data - showing placeholder")
                # Create a placeholder plot when no data
                ggplot2::ggplot() +
                    ggplot2::theme_void() +
                    ggplot2::labs(title = "No strain effects data available for this peak") +
                    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, color = "#7f8c8d"))
            } else {
                message(paste("alleleEffectsModule: Creating allele effects plot with", nrow(effects_data), "data points"))
                # Use the ggplot_alleles function to create the plot
                ggplot_alleles(effects_data)
            }
        })

        # Return module interface
        return(list(
            available_peaks = available_peaks_for_trait,
            selected_peak_marker = selected_peak_marker,
            allele_effects_data = allele_effects_data
        ))
    })
}
