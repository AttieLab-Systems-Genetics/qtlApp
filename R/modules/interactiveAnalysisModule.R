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

        # Use centralized mapping utility from helpers.R
        get_interactive_dataset_name <- map_interactive_dataset_name

        # Reactive to store the selected interaction type
        interaction_type_rv <- shiny::reactiveVal("none")

        # Observer to update stored interaction type when input changes
        shiny::observeEvent(input$interaction_type,
            {
                new_value <- input$interaction_type
                # Guard: only update on real user changes; ignore re-renders setting same value
                if (!is.null(new_value) && !identical(new_value, interaction_type_rv())) {
                    interaction_type_rv(new_value)
                    message(paste("interactiveAnalysisModule: Updated interaction type to:", new_value))
                }
            },
            ignoreNULL = TRUE,
            ignoreInit = TRUE # Prevent firing on initial binding
        )

        # Base mapping reactive (no debounce)
        mapped_dataset_base <- shiny::reactive({
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
                # Stabilize if base already is interactive to prevent double-switch bounce
                interactive_dataset <- if (grepl("interactive", base_dataset, ignore.case = TRUE)) base_dataset else get_interactive_dataset_name(base_dataset, interaction_type)

                if (interactive_dataset != base_dataset) {
                    message(paste("interactiveAnalysisModule: Interactive analysis mode: Using dataset", interactive_dataset, "for interaction type:", interaction_type))
                    return(interactive_dataset)
                } else {
                    message("interactiveAnalysisModule: Interaction type is none or null, using additive dataset")
                }
            }

            return(base_dataset)
        })

        # Debounced mapping for interactive types to avoid teetering
        mapped_dataset_debounced <- mapped_dataset_base %>% shiny::debounce(150)

        # Exported mapping: immediate when switching back to additive, debounced for interactive
        mapped_dataset_for_interaction <- shiny::reactive({
            if (identical(interaction_type_rv(), "none")) {
                mapped_dataset_base()
            } else {
                mapped_dataset_debounced()
            }
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

            if (is.null(dataset_group) || !grepl("^HC_HF", dataset_group, ignore.case = TRUE)) {
                return(NULL)
            }

            # Preserve current selection on re-render
            current_selection <- interaction_type_rv()

            # Determine available interactions by dataset type
            available_interactions <- c("None (Additive only)" = "none")
            if (grepl("HC_HF Liver Genes", dataset_group, ignore.case = TRUE)) {
                available_interactions <- c(available_interactions, "Sex interaction" = "sex", "Diet interaction" = "diet")
            } else if (grepl("HC_HF.*Liver.*Lipid", dataset_group, ignore.case = TRUE)) {
                available_interactions <- c(available_interactions, "Sex interaction" = "sex", "Diet interaction" = "diet", "Sex x Diet interaction" = "sex_diet")
            } else if (grepl("HC_HF.*Clinical", dataset_group, ignore.case = TRUE)) {
                available_interactions <- c(available_interactions, "Sex interaction" = "sex", "Diet interaction" = "diet", "Sex x Diet interaction" = "sex_diet")
            } else if (grepl("HC_HF.*Plasma.*plasma_metabolite|HC_HF.*Plasma.*Metabol", dataset_group, ignore.case = TRUE)) {
                available_interactions <- c(available_interactions, "Sex interaction" = "sex", "Diet interaction" = "diet", "Sex x Diet interaction" = "sex_diet")
            }

            # Compute a stable selected value: prefer current_selection if still valid
            selected_value <- if (current_selection %in% available_interactions) current_selection else {
                # If previous selection is no longer valid, prefer a stable non-none option if present
                if ("sex" %in% available_interactions) "sex" else if ("diet" %in% available_interactions) "diet" else "none"
            }

            htmltools::tagList(
                htmltools::hr(style = "border-top: 2px solid #e74c3c; margin: 15px 0;"),
                htmltools::h5("🧬 Interactive Analysis", style = "color: #2c3e50; margin-bottom: 15px; font-weight: bold;"),
                shiny::selectInput(
                    ns("interaction_type"),
                    label = "Select interaction analysis:",
                    choices = available_interactions,
                    selected = selected_value,
                    width = "100%"
                ),
                shiny::conditionalPanel(
                    condition = paste0("input['", ns("interaction_type"), "'] != 'none'"),
                    htmltools::div(
                        style = "margin-top: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 5px; border-left: 4px solid #e74c3c;",
                        htmltools::p("ℹ️ Interactive analysis will show stacked plots: Interactive LOD scan (top) and Difference plot (Interactive - Additive, bottom).",
                            style = "font-size: 12px; color: #6c757d; margin: 0;"
                        )
                    )
                )
            )
        })

        # Update choices on dataset change without forcing selection to 'none'
        shiny::observeEvent(selected_dataset_reactive(), {
            dataset_group <- selected_dataset_reactive()
            if (is.null(dataset_group) || !grepl("^HC_HF", dataset_group, ignore.case = TRUE)) return()

            available_interactions <- c("None (Additive only)" = "none")
            if (grepl("HC_HF Liver Genes", dataset_group, ignore.case = TRUE)) {
                available_interactions <- c(available_interactions, "Sex interaction" = "sex", "Diet interaction" = "diet")
            } else if (grepl("HC_HF.*Liver.*Lipid", dataset_group, ignore.case = TRUE)) {
                available_interactions <- c(available_interactions, "Sex interaction" = "sex", "Diet interaction" = "diet", "Sex x Diet interaction" = "sex_diet")
            } else if (grepl("HC_HF.*Clinical", dataset_group, ignore.case = TRUE)) {
                available_interactions <- c(available_interactions, "Sex interaction" = "sex", "Diet interaction" = "diet", "Sex x Diet interaction" = "sex_diet")
            } else if (grepl("HC_HF.*Plasma.*plasma_metabolite|HC_HF.*Plasma.*Metabol", dataset_group, ignore.case = TRUE)) {
                available_interactions <- c(available_interactions, "Sex interaction" = "sex", "Diet interaction" = "diet", "Sex x Diet interaction" = "sex_diet")
            }

            shiny::updateSelectInput(session, ns("interaction_type"), choices = available_interactions, selected = interaction_type_rv())
        }, ignoreInit = TRUE)

        # Return module interface
        return(list(
            interaction_type = interaction_type_rv,
            mapped_dataset = mapped_dataset_for_interaction,
            scan_type = scan_type
        ))
    })
}
