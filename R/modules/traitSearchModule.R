#' Trait Search Module
#'
#' @param id Module ID
#' @param import_reactives Reactive containing file_directory and other import data
#' @param selected_dataset_reactive Reactive for the currently selected dataset
#'
#' @importFrom shiny moduleServer NS selectizeInput updateSelectizeInput observeEvent reactive reactiveVal req actionButton
#' @importFrom htmltools div p h5 hr
#'
#' @export
traitSearchUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::fluidPage(
            # fluidRow to encapsulate all content
            shiny::fluidRow(
                # Top row with back button and search input/button
                shiny::fluidRow(
                    # Main content column
                    shiny::column(
                        12,
                        shiny::fluidRow(
                            shiny::column(3, style = "padding-right: 5px;"),
                            shiny::column(
                                6,
                                # Conditional Search Input: shown only when a trait is selected
                                shiny::uiOutput(ns("trait_search_input_ui"))
                            ),
                            shiny::column(3,
                                style = "padding-left: 5px;",
                                # Search button
                                shiny::actionButton(ns("search_trait"), "Search",
                                    icon = shiny::icon("search"),
                                    class = "btn-primary", style = "width: 100%;"
                                )
                            )
                        )
                    )
                ),

                # Second row for the results table
                shiny::fluidRow(
                    shiny::column(
                        12,
                        # Table to display search results
                        DT::DTOutput(ns("trait_search_results_dt"))
                    )
                )
            )
        )
    )
}

#' @rdname traitSearchUI
#' @export
traitSearchServer <- function(id, import_reactives, selected_dataset_reactive) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Reactive value to store the trait selected for LOD scanning
        trait_for_lod_scan_rv <- shiny::reactiveVal(NULL)

        # Store a flag to prevent auto-search immediately after dataset changes
        dataset_just_changed <- shiny::reactiveVal(FALSE)

        # Reset the flag after a short delay
        observe({
            if (dataset_just_changed()) {
                invalidateLater(1000) # Wait 1 second
                dataset_just_changed(FALSE)
            }
        })

        # Observer for the trait search button
        observeEvent(input$trait_search_button, {
            # Use selected_dataset_reactive() directly (it's a string), not selected_dataset_reactive()$group
            shiny::req(selected_dataset_reactive()) # Ensure a dataset is selected

            searched_trait <- input$trait_search_input # Get selected value from selectizeInput

            if (is.null(searched_trait) || !nzchar(searched_trait)) {
                shiny::showNotification("Please select a trait/gene to search.", type = "warning", duration = 3)
                return()
            }

            # Check if the searched trait is different from the current one to avoid re-triggering for no reason
            # Or if current is NULL, then definitely update.
            if (is.null(trait_for_lod_scan_rv()) || !identical(trait_for_lod_scan_rv(), searched_trait)) {
                message(paste(
                    "traitSearchModule: Trait search triggered. Trait for LOD scan set to:", searched_trait,
                    "for dataset:", selected_dataset_reactive()
                ))
                trait_for_lod_scan_rv(searched_trait)

                # Show success notification
                shiny::showNotification(
                    paste("Searching for trait:", searched_trait),
                    type = "message",
                    duration = 2
                )
            } else {
                message(paste("traitSearchModule: Trait search for already selected trait:", searched_trait, "- no change."))
            }
        })

        # AUTO-TRIGGER SEARCH WHEN TRAIT IS SELECTED FROM DROPDOWN
        # This allows users to skip clicking the button - just select from dropdown and it searches automatically
        # But NOT immediately after dataset changes to prevent slow loading
        observeEvent(input$trait_search_input,
            {
                searched_trait <- input$trait_search_input

                # Only auto-search if:
                # 1. A valid trait is selected and dataset is available
                # 2. Dataset didn't just change (to prevent slow loading)
                # 3. This is a new selection by the user (not just preservation from dataset change)
                if (!is.null(searched_trait) && nzchar(searched_trait) &&
                    !is.null(selected_dataset_reactive()) &&
                    !dataset_just_changed()) {
                    if (is.null(trait_for_lod_scan_rv()) || !identical(trait_for_lod_scan_rv(), searched_trait)) {
                        message(paste("traitSearchModule: Auto-search triggered for trait:", searched_trait))
                        trait_for_lod_scan_rv(searched_trait)

                        # Show notification for auto-search
                        shiny::showNotification(
                            paste("Auto-searching for trait:", searched_trait),
                            type = "message",
                            duration = 2
                        )
                    }
                } else if (dataset_just_changed()) {
                    message(paste("traitSearchModule: Trait selection preserved but auto-search skipped due to recent dataset change. Use search button to trigger LOD scan."))
                }
            },
            ignoreNULL = TRUE,
            ignoreInit = TRUE
        )

        # TRAIT SEARCH DROPDOWN UPDATE LOGIC
        # Updates trait search choices when dataset changes (server-side selectize for performance)
        shiny::observeEvent(shiny::req(selected_dataset_reactive()), {
            shiny::req(import_reactives()) # Ensure import data is available
            current_ds <- selected_dataset_reactive()

            message(paste("traitSearchModule: Updating trait search choices for dataset:", current_ds))

            # SET FLAG to prevent immediate auto-search after dataset change
            dataset_just_changed(TRUE)

            # PRESERVE CURRENT SELECTION: Store the current trait selection before updating
            current_trait_selection <- input$trait_search_input
            message(paste("traitSearchModule: Current trait selection before dataset change:", current_trait_selection %||% "None"))

            # RESET TRAIT SEARCH INPUT: Prevent conflicts during update
            shiny::freezeReactiveValue(input, "trait_search_input") # Temporarily freeze reactive
            shiny::updateSelectizeInput(session, "trait_search_input",
                choices = character(0),
                selected = character(0),
                options = list(placeholder = "Loading traits...")
            )

            # GET NEW TRAIT CHOICES: Based on selected dataset
            choices <- get_trait_choices(import_reactives(), current_ds) # External helper function

            if (!is.null(choices) && length(choices) > 0) {
                message(paste("traitSearchModule: Found", length(choices), "traits for dataset:", current_ds))

                # CHECK IF PREVIOUS SELECTION EXISTS IN NEW DATASET
                trait_to_select <- NULL
                if (!is.null(current_trait_selection) && nzchar(current_trait_selection)) {
                    if (current_trait_selection %in% choices) {
                        trait_to_select <- current_trait_selection
                        message(paste("traitSearchModule: Preserving trait selection:", current_trait_selection, "- found in new dataset"))
                    } else {
                        message(paste("traitSearchModule: Previous trait selection:", current_trait_selection, "- not found in new dataset, clearing selection"))
                    }
                }

                # UPDATE TRAIT SEARCH DROPDOWN: With server-side processing for large lists
                shiny::updateSelectizeInput(session, "trait_search_input",
                    choices = choices,
                    selected = trait_to_select, # Preserve selection if trait exists in new dataset
                    options = list(
                        placeholder = "Type to search traits/genes...",
                        maxItems = 1, # Only allow single selection
                        maxOptions = 10 # Limit displayed options for performance
                    ),
                    server = TRUE # Enable server-side processing
                )

                # Notify user if trait was preserved but auto-search was disabled
                if (!is.null(trait_to_select)) {
                    shiny::showNotification(
                        paste("Trait '", trait_to_select, "' preserved. Click 'Search & Plot LOD Scan' to analyze in new dataset."),
                        type = "message",
                        duration = 4
                    )
                }
            } else {
                message(paste("traitSearchModule: No traits found for dataset:", current_ds))
                # No traits available for this dataset
                shiny::updateSelectizeInput(session, "trait_search_input",
                    choices = character(0),
                    selected = character(0),
                    options = list(placeholder = "No traits available for this dataset")
                )
            }
        })

        # When the main selected dataset changes, clear the trait_for_lod_scan_rv
        # to prevent a scan from an old selection on a new plot type.
        # BUT preserve trait selection when switching between HC_HF variants (additive <-> interactive)
        shiny::observeEvent(selected_dataset_reactive(),
            {
                current_dataset <- selected_dataset_reactive()
                message("traitSearchModule: Main dataset group changed. Clearing trait_for_lod_scan_rv.")

                # Check if this is just switching between HC_HF dataset variants (any HC_HF dataset type)
                is_hc_hf_dataset_switch <- !is.null(current_dataset) &&
                    grepl("^HC_HF", current_dataset, ignore.case = TRUE)

                if (!is_hc_hf_dataset_switch) {
                    # Only clear trait selection for real dataset changes, not interaction type changes
                    trait_for_lod_scan_rv(NULL) # Clear any active LOD scan
                    # Clear the search input selection (choices will be updated by the other observer)
                    updateSelectizeInput(session, "trait_search_input",
                        selected = character(0)
                    )
                } else {
                    message("traitSearchModule: HC_HF dataset switch detected - preserving trait selection")
                }
            },
            ignoreNULL = TRUE,
            ignoreInit = TRUE
        )

        # Return module interface
        return(list(
            trait_for_lod_scan = trait_for_lod_scan_rv
        ))
    })
}
