#' Interactive Analysis Module
#'
#' @param id Module ID
#' @param selected_dataset_reactive Reactive for the currently selected dataset
#'
#' @importFrom shiny moduleServer NS selectInput renderUI req reactive reactiveVal observeEvent
#' @importFrom htmltools div p h5 hr tagList
#'
#' @export
interactiveAnalysisUI <- function(id) {
    ns <- shiny::NS(id)

    # This UI will be rendered conditionally by the server
    shiny::uiOutput(ns("interactive_analysis_section"))
}

#' @rdname interactiveAnalysisUI
#' @export
interactiveAnalysisServer <- function(id, selected_dataset_reactive) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Helper function to map HC_HF dataset names to their interactive versions
        # Only maps to datasets that actually exist in the file system
        get_interactive_dataset_name <- function(base_dataset, interaction_type) {
            if (is.null(interaction_type) || interaction_type == "none") {
                return(base_dataset)
            }

            # HC_HF Liver Genes (supports both Sex and Diet interactions)
            if (grepl("HC_HF Liver Genes", base_dataset, ignore.case = TRUE)) {
                if (interaction_type == "sex") {
                    return("HC_HF Liver Genes, interactive (Sex)")
                } else if (interaction_type == "diet") {
                    return("HC_HF Liver Genes, interactive (Diet)")
                }
            }
            # HC_HF Liver Lipids (only supports Diet interaction)
            else if (grepl("HC_HF.*Liver.*Lipid", base_dataset, ignore.case = TRUE)) {
                if (interaction_type == "diet") {
                    return("HC_HF Liver Lipids, interactive (Diet)")
                }
                # No Sex interaction available for Liver Lipids - return original
            }
            # HC_HF Clinical Traits (supports both Sex and Diet interactions)
            else if (grepl("HC_HF.*Clinical", base_dataset, ignore.case = TRUE)) {
                if (interaction_type == "sex") {
                    return("HC_HF Systemic Clinical Traits, interactive (Sex)")
                } else if (interaction_type == "diet") {
                    return("HC_HF Systemic Clinical Traits, interactive (Diet)")
                }
            }
            # HC_HF Plasma Metabolites (supports Sex interaction)
            else if (grepl("HC_HF.*Plasma.*plasma_metabolite", base_dataset, ignore.case = TRUE)) {
                if (interaction_type == "sex") {
                    return("HC_HF Plasma plasma_metabolite, interactive (Sex)")
                }
                # No Diet interaction available for Plasma Metabolites - return original
            }

            # Fallback to original dataset if no mapping found
            return(base_dataset)
        }

        # Reactive to store the selected interaction type
        interaction_type_rv <- shiny::reactive({
            interaction_value <- input$interaction_type %||% "none"
            message(paste("interactiveAnalysisModule: interaction_type_rv: Current value is:", interaction_value))
            return(interaction_value)
        })

        # Store the current interaction type to preserve across UI re-renders
        current_interaction_type_rv <- shiny::reactiveVal("none")

        # Observer to update stored interaction type when input changes
        shiny::observeEvent(input$interaction_type,
            {
                new_value <- input$interaction_type
                if (!is.null(new_value)) {
                    current_interaction_type_rv(new_value)
                    message(paste("interactiveAnalysisModule: Updated interaction type to:", new_value))
                }
            },
            ignoreNULL = TRUE
        )

        # New reactive that handles the dataset name mapping for interactive analysis
        mapped_dataset_for_interaction <- shiny::reactive({
            base_dataset <- selected_dataset_reactive()
            interaction_type <- interaction_type_rv()

            if (is.null(base_dataset)) {
                return(NULL)
            }

            # Check if this is an HC_HF dataset that supports interactive analysis
            is_hc_hf_dataset <- grepl("^HC_HF", base_dataset, ignore.case = TRUE)

            if (is_hc_hf_dataset && !is.null(interaction_type) && interaction_type != "none") {
                message(paste("interactiveAnalysisModule: HC_HF dataset detected:", base_dataset, "interaction_type is:", interaction_type))

                # Use helper function to get the appropriate dataset name
                interactive_dataset <- get_interactive_dataset_name(base_dataset, interaction_type)

                if (interactive_dataset != base_dataset) {
                    message(paste("interactiveAnalysisModule: Interactive analysis mode: Using dataset", interactive_dataset, "for interaction type:", interaction_type))
                    return(interactive_dataset)
                } else {
                    message("interactiveAnalysisModule: Interaction type is none or null, using additive dataset")
                }
            }

            return(base_dataset)
        })

        # Detect if current dataset is additive or interactive based on dataset name and interaction type
        scan_type <- shiny::reactive({
            dataset_name <- mapped_dataset_for_interaction()

            if (is.null(dataset_name) || dataset_name == "") {
                return("additive") # Default to additive
            }

            # Check if dataset name contains "interactive" or if interaction type is selected
            if (grepl("interactive", dataset_name, ignore.case = TRUE)) {
                return("interactive")
            } else {
                return("additive")
            }
        })

        # Interactive Analysis section - show for all HC_HF datasets
        output$interactive_analysis_section <- shiny::renderUI({
            dataset_group <- selected_dataset_reactive()

            # Show interactive analysis controls for all HC_HF datasets (Genes, Lipids, Clinical Traits, Metabolites)
            if (!is.null(dataset_group) && grepl("^HC_HF", dataset_group, ignore.case = TRUE)) {
                # Preserve the current selection when re-rendering
                current_selection <- current_interaction_type_rv()

                # Determine what interaction types are available for this dataset
                available_interactions <- c("None (Additive only)" = "none")

                # Check what interactions are actually available based on dataset type
                if (grepl("HC_HF Liver Genes", dataset_group, ignore.case = TRUE)) {
                    available_interactions <- c(available_interactions,
                        "Sex interaction" = "sex",
                        "Diet interaction" = "diet"
                    )
                } else if (grepl("HC_HF.*Liver.*Lipid", dataset_group, ignore.case = TRUE)) {
                    available_interactions <- c(available_interactions,
                        "Diet interaction" = "diet"
                    )
                    # No Sex interaction for Liver Lipids
                } else if (grepl("HC_HF.*Clinical", dataset_group, ignore.case = TRUE)) {
                    available_interactions <- c(available_interactions,
                        "Sex interaction" = "sex",
                        "Diet interaction" = "diet"
                    )
                } else if (grepl("HC_HF.*Plasma.*Metabol", dataset_group, ignore.case = TRUE)) {
                    available_interactions <- c(available_interactions,
                        "Sex interaction" = "sex"
                    )
                    # No Diet interaction available for Plasma Metabolites
                }

                tagList(
                    hr(style = "border-top: 2px solid #e74c3c; margin: 15px 0;"),
                    h5("🧬 Interactive Analysis", style = "color: #2c3e50; margin-bottom: 15px; font-weight: bold;"),
                    shiny::selectInput(
                        ns("interaction_type"),
                        label = "Select interaction analysis:",
                        choices = available_interactions,
                        selected = if (current_selection %in% available_interactions) current_selection else "none",
                        width = "100%"
                    ),
                    shiny::conditionalPanel(
                        condition = paste0("input['", ns("interaction_type"), "'] != 'none'"),
                        div(
                            style = "margin-top: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 5px; border-left: 4px solid #e74c3c;",
                            p("ℹ️ Interactive analysis will show stacked plots: Interactive LOD scan (top) and Difference plot (Interactive - Additive, bottom).",
                                style = "font-size: 12px; color: #6c757d; margin: 0;"
                            )
                        )
                    )
                )
            } else {
                NULL # Don't show for other datasets
            }
        })

        # Return module interface
        return(list(
            interaction_type = interaction_type_rv,
            mapped_dataset = mapped_dataset_for_interaction,
            scan_type = scan_type
        ))
    })
}
