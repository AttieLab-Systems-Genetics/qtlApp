#' Dataset Selection Module
#'
#' @param id Module ID
#' @param import_reactives Reactive containing file_directory and other import data
#'
#' @importFrom shiny moduleServer NS selectInput updateSelectInput renderUI req observeEvent reactive
#' @importFrom htmltools div p h6 span
#' @importFrom stats setNames
#' @importFrom data.table as.data.table
#'
#' @export
datasetSelectionUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        # Top navigation bar for dataset category selection
        div(
            style = "background: linear-gradient(135deg, #2c3e50, #3498db); padding: 15px; margin-bottom: 20px; border-radius: 8px;",
            div(
                style = "display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap;",
                h3("QTL Scan Visualizer",
                    style = "color: white; margin: 0; font-weight: bold;"
                ),
                div(
                    style = "display: flex; align-items: center; gap: 15px;",
                    h5("Dataset Category:",
                        style = "color: white; margin: 0; font-weight: bold;"
                    ),
                    shiny::selectInput(ns("dataset_category_selector"),
                        NULL,
                        choices = c("Loading..." = ""),
                        width = "200px"
                    )
                )
            )
        ),

        # Dataset selection section
        h5("Dataset Selection", style = "color: #2c3e50; margin-bottom: 15px; font-weight: bold;"),

        # Conditional dataset selector - hidden for HC_HF categories, shown for others
        shiny::uiOutput(ns("dataset_selection_ui"))
    )
}

#' @rdname datasetSelectionUI
#' @export
datasetSelectionServer <- function(id, import_reactives) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Reactive for file index data table
        file_index_dt <- shiny::reactive({
            shiny::req(import_reactives()$file_directory)
            dt <- data.table::as.data.table(import_reactives()$file_directory)
            shiny::validate(
                shiny::need("dataset_category" %in% names(dt), "Error: 'dataset_category' column missing in file_index.csv."),
                shiny::need("group" %in% names(dt), "Error: 'group' column missing in file_index.csv.")
            )
            return(dt)
        })

        # Initialize dataset category choices
        shiny::observe({
            shiny::req(file_index_dt())
            categories <- unique(file_index_dt()$dataset_category)
            if (length(categories) > 0) {
                shiny::updateSelectInput(session, "dataset_category_selector",
                    choices = stats::setNames(categories, categories),
                    selected = categories[1]
                )
            } else {
                shiny::updateSelectInput(session, "dataset_category_selector",
                    choices = c("No categories found" = ""), selected = ""
                )
            }
        })

        # Update specific dataset choices when category changes
        shiny::observe({
            shiny::req(file_index_dt(), input$dataset_category_selector)
            selected_cat <- input$dataset_category_selector

            if (!is.null(selected_cat) && nzchar(selected_cat) && selected_cat != "No categories found") {
                datasets_in_category <- file_index_dt()[dataset_category == selected_cat, ]
                specific_datasets_choices <- unique(datasets_in_category$group)

                # Auto-select HC_HF datasets for all categories that have them
                hc_hf_dataset <- NULL

                if (selected_cat == "Liver Genes") {
                    # Look for HC_HF Liver Genes dataset (additive version)
                    hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Liver.*Genes", specific_datasets_choices, ignore.case = TRUE) &
                        !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
                } else if (selected_cat == "Liver Lipids") {
                    # Look for HC_HF Liver Lipids dataset (additive version)
                    hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Liver.*Lipid", specific_datasets_choices, ignore.case = TRUE) &
                        !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
                } else if (selected_cat == "Clinical Traits") {
                    # Look for HC_HF Clinical Traits dataset (additive version)
                    hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Clinical", specific_datasets_choices, ignore.case = TRUE) &
                        !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
                } else if (selected_cat == "Plasma Metabolites") {
                    # Look for HC_HF Plasma Metabolites dataset
                    hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Plasma.*Metabol", specific_datasets_choices, ignore.case = TRUE) &
                        !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
                } else if (selected_cat == "Liver Isoforms") {
                    # Look for HC_HF Liver Isoforms dataset
                    hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Liver.*Isoform", specific_datasets_choices, ignore.case = TRUE) &
                        !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
                } else if (selected_cat == "Liver Splice Junctions") {
                    # Look for HC_HF Liver Splice Junctions dataset (additive version)
                    hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Liver.*Splice.*Junc", specific_datasets_choices, ignore.case = TRUE) &
                        !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
                }

                if (!is.null(hc_hf_dataset) && length(hc_hf_dataset) > 0) {
                    message(paste("Auto-selecting HC_HF dataset for", selected_cat, ":", hc_hf_dataset[1]))
                    # Set choices to just the HC_HF dataset (hide the dropdown essentially)
                    shiny::updateSelectInput(session, "specific_dataset_selector",
                        choices = stats::setNames(hc_hf_dataset[1], hc_hf_dataset[1]),
                        selected = hc_hf_dataset[1]
                    )
                } else {
                    # Fallback to normal behavior if HC_HF dataset not found
                    message(paste("Warning: No HC_HF dataset found for", selected_cat, ", showing all options"))
                    if (length(specific_datasets_choices) > 0) {
                        shiny::updateSelectInput(session, "specific_dataset_selector",
                            choices = stats::setNames(specific_datasets_choices, specific_datasets_choices),
                            selected = specific_datasets_choices[1]
                        )
                    } else {
                        shiny::updateSelectInput(session, "specific_dataset_selector",
                            choices = c("No datasets in category" = ""), selected = ""
                        )
                    }
                }
            } else {
                shiny::updateSelectInput(session, "specific_dataset_selector",
                    choices = c("Select category first" = ""), selected = ""
                )
            }
        })

        # Conditional dataset selection UI - show info for HC_HF categories, selector for others
        output$dataset_selection_ui <- shiny::renderUI({
            selected_cat <- input$dataset_category_selector

            if (is.null(selected_cat) || !nzchar(selected_cat)) {
                return(div(
                    style = "padding: 15px; text-align: center; color: #7f8c8d; background: #f8f9fa; border-radius: 5px; border: 1px solid #e9ecef;",
                    p("Select a dataset category above", style = "margin: 0; font-style: italic;")
                ))
            }

            # Show information panel for categories that have HC_HF auto-selection
            if (selected_cat %in% c("Liver Genes", "Liver Lipids", "Clinical Traits", "Plasma Metabolites", "Liver Isoforms", "Liver Splice Junctions")) {
                # Determine the dataset name and interaction info based on category
                dataset_info <- switch(selected_cat,
                    "Liver Genes" = list(
                        name = "HC_HF Liver Genes (Additive)",
                        interaction_note = "Use interaction controls for Sex/Diet effects."
                    ),
                    "Liver Lipids" = list(
                        name = "HC_HF Liver Lipids (Additive)",
                        interaction_note = "Use interaction controls for Diet effects."
                    ),
                    "Clinical Traits" = list(
                        name = "HC_HF Clinical Traits (Additive)",
                        interaction_note = "Use interaction controls for Sex/Diet effects."
                    ),
                    "Plasma Metabolites" = list(
                        name = "HC_HF Plasma Metabolites (Additive)",
                        interaction_note = "Use interaction controls for Sex/Diet effects."
                    ),
                    "Liver Isoforms" = list(
                        name = "HC_HF Liver Isoforms",
                        interaction_note = "No interactive analysis available for this dataset."
                    ),
                    "Liver Splice Junctions" = list(
                        name = "HC_HF Liver Splice Junctions (Additive)",
                        interaction_note = "No interactive analysis available for this dataset."
                    )
                )

                return(div(
                    style = "padding: 15px; background: #e8f5e8; border-radius: 5px; border-left: 4px solid #28a745;",
                    div(
                        style = "display: flex; align-items: center; gap: 10px;",
                        span("✓", style = "color: #28a745; font-weight: bold; font-size: 16px;"),
                        div(
                            h6(dataset_info$name, style = "color: #155724; margin: 0; font-weight: bold;"),
                            p(paste("Auto-selected for streamlined analysis.", dataset_info$interaction_note),
                                style = "color: #155724; margin: 5px 0 0 0; font-size: 12px;"
                            )
                        )
                    )
                ))
            } else {
                # For other categories, show normal selector
                return(shiny::selectInput(
                    ns("specific_dataset_selector"),
                    "Select Specific Dataset:",
                    choices = c("Loading..." = ""),
                    width = "100%"
                ))
            }
        })

        # Main selected dataset group reactive
        main_selected_dataset_group <- shiny::reactive({
            selected_cat <- input$dataset_category_selector

            # Auto-select HC_HF datasets for categories that have them
            if (!is.null(selected_cat) && selected_cat %in% c("Liver Genes", "Liver Lipids", "Clinical Traits", "Plasma Metabolites", "Liver Isoforms")) {
                if (!is.null(selected_cat) && selected_cat == "Liver Splice Junctions") {
                    # handled below by general logic if not in the list; extend list to include splice junctions
                }

                shiny::req(file_index_dt())
                datasets_in_category <- file_index_dt()[dataset_category == selected_cat, ]
                specific_datasets_choices <- unique(datasets_in_category$group)

                # Find the appropriate HC_HF dataset based on category
                hc_hf_dataset <- NULL

                if (selected_cat == "Liver Genes") {
                    hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Liver.*Genes", specific_datasets_choices, ignore.case = TRUE) &
                        !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
                } else if (selected_cat == "Liver Lipids") {
                    hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Liver.*Lipid", specific_datasets_choices, ignore.case = TRUE) &
                        !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
                } else if (selected_cat == "Clinical Traits") {
                    hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Clinical", specific_datasets_choices, ignore.case = TRUE) &
                        !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
                } else if (selected_cat == "Plasma Metabolites") {
                    hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Plasma.*Metabol", specific_datasets_choices, ignore.case = TRUE) &
                        !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
                } else if (selected_cat == "Liver Isoforms") {
                    hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Liver.*Isoform", specific_datasets_choices, ignore.case = TRUE) &
                        !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
                } else if (selected_cat == "Liver Splice Junctions") {
                    hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Liver.*Splice.*Junc", specific_datasets_choices, ignore.case = TRUE) &
                        !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
                }

                if (!is.null(hc_hf_dataset) && length(hc_hf_dataset) > 0) {
                    message(paste("Auto-selected dataset for", selected_cat, "category:", hc_hf_dataset[1]))
                    return(hc_hf_dataset[1])
                } else {
                    message(paste("Warning: No HC_HF dataset found in", selected_cat, "category"))
                    return(NULL)
                }
            }

            # Normal handling for other categories
            shiny::req(input$specific_dataset_selector)
            selected_group <- input$specific_dataset_selector

            if (is.null(selected_group) || !nzchar(selected_group) ||
                selected_group %in% c("Select category first", "No datasets in category")) {
                return(NULL)
            }

            message(paste("Main selected dataset group:", selected_group))
            return(selected_group)
        })

        # Selected dataset category reactive
        selected_dataset_category_reactive <- shiny::reactive({
            shiny::req(main_selected_dataset_group(), file_index_dt())
            info <- file_index_dt()[group == main_selected_dataset_group()]
            if (nrow(info) > 0) {
                return(unique(info$dataset_category)[1])
            }
            return(NULL)
        })

        # Return module interface
        return(list(
            selected_dataset = main_selected_dataset_group,
            dataset_category = selected_dataset_category_reactive,
            file_index_dt = file_index_dt
        ))
    })
}
