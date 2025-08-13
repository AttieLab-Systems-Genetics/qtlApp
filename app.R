# Load required R packages
library(shiny)
library(bslib)
library(dplyr)
library(purrr)
library(stringr)
library(qtl2)
library(ggplot2)
library(plotly)
library(shinyjs)
library(shinycssloaders)
library(data.table)
library(ggiraph)
library(writexl)
library(fontawesome)
library(DT)
library(reshape2)
library(htmltools) # Added for tags
library(stats) # Added for setNames

# Source files in specific order
source("R/helpers.R")
source("R/data_handling.R")
source("R/import_data.R")
source("R/importApp.R")
source("R/downloadApp.R")
source("R/cisTransPlotApp.R")
source("R/manhattanPlotApp.R")
source("R/ui_styles.R")
source("R/plot_enhancements.R")
source("R/plot_null.R")
source("R/peak_finder.R")
source("R/trait_scan.R")
source("R/ggplot_alleles.R")
source("R/ggplot_qtl_scan.R")
source("R/ggplotly_qtl_scan.R")
source("R/peak_info.R")
source("R/QTL_plot_visualizer.R")
source("R/fst_rows.R")
source("R/traitApp.R")
source("R/traitProcessingModule.R")
source("R/scanPlotModule.R") # Source our new module
source("R/profilePlotApp.R") # Source the new profile plot module
source("R/mainUI.R") # Source our new UI module

# Set maximum file upload size
options(shiny.maxRequestSize = 20000 * 1024^2) # 20 GB

# UI Definition
ui <- mainUI()

# Server Definition
server <- function(input, output, session) {
    ns_app_controller <- shiny::NS("app_controller")

    # Source helper for allele plots if not already available
    if (!exists("ggplot_alleles", mode = "function")) {
        source("R/ggplot_alleles.R")
    }

    # Define the %||% operator for null coalescing
    `%||%` <- function(a, b) if (!is.null(a)) a else b

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
        # HC_HF Liver Lipids (supports Diet and Sex x Diet interactions)
        else if (grepl("HC_HF.*Liver.*Lipid", base_dataset, ignore.case = TRUE)) {
            if (interaction_type == "diet") {
                return("HC_HF Liver Lipids, interactive (Diet)")
            } else if (interaction_type == "sex") {
                return("HC_HF Liver Lipids, interactive (Sex)")
            } else if (interaction_type == "sex_diet") {
                return("HC_HF Liver Lipids, interactive (Sex_Diet)")
            }
        }
        # HC_HF Clinical Traits (supports Sex, Diet, and Sex x Diet interactions)
        else if (grepl("HC_HF.*Clinical", base_dataset, ignore.case = TRUE)) {
            if (interaction_type == "sex") {
                return("HC_HF Systemic Clinical Traits, interactive (Sex)")
            } else if (interaction_type == "diet") {
                return("HC_HF Systemic Clinical Traits, interactive (Diet)")
            } else if (interaction_type == "sex_diet") {
                return("HC_HF Systemic Clinical Traits, interactive (Sex_Diet)")
            }
        } else if (grepl("HC_HF.*Plasma.*Metabol", base_dataset, ignore.case = TRUE)) {
            if (interaction_type == "sex") {
                return("HC_HF Plasma plasma_metabolite, interactive (Sex)")
            } else if (interaction_type == "diet") {
                return("HC_HF Plasma plasma_metabolite, interactive (Diet)")
            } else if (interaction_type == "sex_diet") {
                return("HC_HF Plasma plasma_metabolite, interactive (Sex_Diet)")
            }
        }

        # Fallback to original dataset if no mapping found
        return(base_dataset)
    }

    trait_cache <- new.env(parent = emptyenv())
    peaks_cache <- new.env(parent = emptyenv())

    import_reactives <- importServer("import")

    # Store the current interaction type to preserve across UI re-renders
    current_interaction_type_rv <- shiny::reactiveVal("none")

    # NEW: Reactive to determine if stacked plots should be shown
    show_stacked_plots <- shiny::reactive({
        interaction_type <- current_interaction_type_rv()
        !is.null(interaction_type) && interaction_type != "none"
    })

    # Dynamic LOD threshold slider based on scan type
    output[[ns_app_controller("lod_threshold_slider")]] <- shiny::renderUI({
        interaction_type <- sidebar_interaction_type_rv()
        # For overview plots, interaction types imply difference plots using qtlxcovar files.
        scan_info <- switch(interaction_type,
            "sex" = list(type = "Sex Difference", min = 4.1),
            "diet" = list(type = "Diet Difference", min = 4.1),
            "sex_diet" = list(type = "Sex x Diet Difference", min = 9.5),
            "none" = list(type = "Additive", min = 7.5),
            list(type = "Additive", min = 7.5) # Default
        )

        # Use the input value if it exists and is valid, otherwise default to the new min
        current_val <- input[[ns_app_controller("LOD_thr")]]
        default_val <- if (!is.null(current_val) && current_val >= scan_info$min) {
            current_val
        } else {
            scan_info$min
        }

        sliderInput(shiny::NS("app_controller", "LOD_thr"),
            label = paste0("LOD Threshold (", scan_info$type, " scan):"),
            min = scan_info$min, max = 20, value = default_val, step = 0.5,
            width = "100%"
        )
    })

    # Our own LOD threshold reactive (no longer from mainParServer)
    lod_threshold_rv <- shiny::reactive({
        current_scan_type <- scan_type()
        interaction_type <- current_interaction_type_rv()
        default_threshold <- if (current_scan_type == "interactive") {
            if (!is.null(interaction_type) && interaction_type == "sex_diet") 15.7 else 10.5
        } else {
            7.5
        }
        input[[ns_app_controller("LOD_thr")]] %||% default_threshold
    }) %>% shiny::debounce(300) # Debounce LOD threshold to prevent rapid re-firing

    # Reactive for selected chromosome (for zooming into specific chromosomes)
    selected_chromosome_rv <- shiny::reactive({
        input[[ns_app_controller("selected_chr")]] %||% "All" # Default to "All" chromosomes
    })



    file_index_dt <- shiny::reactive({
        shiny::req(import_reactives()$file_directory)
        dt <- data.table::as.data.table(import_reactives()$file_directory)
        shiny::validate(
            shiny::need("dataset_category" %in% names(dt), "Error: 'dataset_category' column missing in file_index.csv."),
            shiny::need("group" %in% names(dt), "Error: 'group' column missing in file_index.csv.")
        )
        return(dt)
    })

    shiny::observe({
        shiny::req(file_index_dt())
        categories <- unique(file_index_dt()$dataset_category)
        if (length(categories) > 0) {
            shiny::updateSelectInput(session, ns_app_controller("dataset_category_selector"),
                choices = stats::setNames(categories, categories),
                selected = categories[1]
            )
        } else {
            shiny::updateSelectInput(session, ns_app_controller("dataset_category_selector"),
                choices = c("No categories found" = ""), selected = ""
            )
        }
    })

    main_selected_dataset_group <- shiny::reactive({
        selected_cat <- input[[ns_app_controller("dataset_category_selector")]]
        shiny::req(selected_cat, file_index_dt())

        datasets_in_category <- file_index_dt()[dataset_category == selected_cat, ]
        specific_datasets_choices <- unique(datasets_in_category$group)

        # Find the appropriate HC_HF dataset (additive) or the first available dataset
        hc_hf_dataset <- NULL

        # This logic now covers all categories that have an auto-selectable HC_HF dataset
        if (selected_cat %in% c("Liver Genes", "Liver Lipids", "Clinical Traits", "Plasma Metabolites", "Liver Isoforms")) {
            pattern <- switch(selected_cat,
                "Liver Genes" = "^HC_HF.*Liver.*Genes",
                "Liver Lipids" = "^HC_HF.*Liver.*Lipid",
                "Clinical Traits" = "^HC_HF.*Clinical",
                "Plasma Metabolites" = "^HC_HF.*Plasma.*Metabol",
                "Liver Isoforms" = "^HC_HF.*Liver.*Isoform"
            )
            hc_hf_dataset <- specific_datasets_choices[grepl(pattern, specific_datasets_choices, ignore.case = TRUE) &
                !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
        }

        if (!is.null(hc_hf_dataset) && length(hc_hf_dataset) > 0) {
            message(paste("Auto-selected dataset for", selected_cat, "category:", hc_hf_dataset[1]))
            return(hc_hf_dataset[1])
        }

        # Fallback for other categories or if HC_HF is not found
        if (length(specific_datasets_choices) > 0) {
            message(paste("Defaulting to first available dataset for", selected_cat, ":", specific_datasets_choices[1]))
            return(specific_datasets_choices[1])
        }

        message(paste("Warning: No datasets found for category:", selected_cat))
        return(NULL)
    })

    # New reactive that handles the dataset name mapping for interactive analysis
    mapped_dataset_for_interaction <- shiny::reactive({
        base_dataset <- main_selected_dataset_group()
        interaction_type <- current_interaction_type_rv()

        if (is.null(base_dataset)) {
            return(NULL)
        }

        # Check if this is an HC_HF dataset that supports interactive analysis
        is_hc_hf_dataset <- grepl("^HC_HF", base_dataset, ignore.case = TRUE)

        if (is_hc_hf_dataset && !is.null(interaction_type) && interaction_type != "none") {
            message(paste("mapped_dataset_for_interaction: HC_HF dataset detected:", base_dataset, "interaction_type is:", interaction_type))

            # Use helper function to get the appropriate dataset name
            interactive_dataset <- get_interactive_dataset_name(base_dataset, interaction_type)

            if (interactive_dataset != base_dataset) {
                message(paste("Interactive analysis mode: Using dataset", interactive_dataset, "for interaction type:", interaction_type))
                return(interactive_dataset)
            } else {
                message("mapped_dataset_for_interaction: Interaction type is none or null, using additive dataset")
            }
        }

        return(base_dataset)
    })

    # Detect if current dataset is additive or interactive based on dataset name and interaction type
    scan_type <- shiny::reactive({
        dataset_name <- main_selected_dataset_group()

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

    # Reactive to find peak data for the selected trait
    peaks_data_for_trait <- shiny::reactive({
        trait_val <- trait_for_lod_scan_rv()
        if (is.null(trait_val)) {
            return(NULL)
        }

        shiny::req(mapped_dataset_for_interaction(), import_reactives())

        dataset_val <- mapped_dataset_for_interaction()

        message(paste("scanApp: Finding peaks for trait:", trait_val, "in dataset:", dataset_val))

        # Determine trait type for this dataset to pass to peak_finder
        # Use the base dataset name for trait type determination
        base_dataset <- main_selected_dataset_group()
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
                message(paste("scanApp: Found", if (is.null(peaks)) 0 else nrow(peaks), "peaks for trait:", trait_val))
                peaks
            },
            error = function(e) {
                message(paste("scanApp: Error finding peaks for trait", trait_val, ":", e$message))
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
        lod_thr <- lod_threshold_rv()
        shiny::req(lod_thr)

        # Check what columns are available
        message(paste("scanApp: Peak data columns:", paste(colnames(peaks), collapse = ", ")))

        # Use qtl_lod column (from peak_finder) instead of lod
        if (!"qtl_lod" %in% colnames(peaks)) {
            message("scanApp: No qtl_lod column found in peaks data")
            return(NULL)
        }

        # Filter by LOD threshold and sort by LOD descending
        filtered_peaks <- peaks[peaks$qtl_lod >= lod_thr, ]
        if (nrow(filtered_peaks) == 0) {
            message(paste("scanApp: No peaks above LOD threshold", lod_thr, "for trait:", trait_for_lod_scan_rv()))
            return(NULL)
        }

        # Sort by qtl_lod descending
        filtered_peaks <- filtered_peaks[order(-filtered_peaks$qtl_lod), ]

        message(paste("scanApp: Found", nrow(filtered_peaks), "peaks above threshold for trait:", trait_for_lod_scan_rv()))
        filtered_peaks
    })

    # Reactive value to store the selected peak marker from dropdown
    selected_peak_from_dropdown <- shiny::reactiveVal(NULL)

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
            message(paste("scanApp: Using dropdown-selected peak:", selected_marker, "with LOD:", selected_peak_info$qtl_lod[1]))
        } else {
            # Default to highest peak
            selected_marker <- available_peaks$marker[1]
            message(paste("scanApp: Using highest peak (default):", selected_marker, "with LOD:", available_peaks$qtl_lod[1]))
        }

        selected_marker
    })

    # Reactive to prepare allele effects data, now driven by the selected peak
    allele_effects_data <- shiny::reactive({
        # This now correctly uses the output from the scanServer module
        peak_info <- scan_module_outputs$selected_peak()
        if (is.null(peak_info)) {
            message("scanApp: allele_effects_data reactive - no peak selected")
            return(NULL)
        }

        message(paste("scanApp: *** ALLELE EFFECTS DATA REACTIVE *** Processing peak:", peak_info$marker, "LOD:", peak_info$qtl_lod))

        # Use pivot_peaks helper function to reshape data for plotting
        reshaped_data <- pivot_peaks(peak_info, peak_info$marker)

        if (is.null(reshaped_data) || nrow(reshaped_data) == 0) {
            message(paste("scanApp: No allele effects data available for marker:", peak_info$marker))
            return(NULL)
        }

        # Add trait name to the data for plot labeling
        reshaped_data$trait <- peak_info$trait
        message(paste("scanApp: Successfully prepared allele effects data for marker:", peak_info$marker, "with", nrow(reshaped_data), "data points"))
        reshaped_data
    })

    # NEW: Reactive for the first difference allele plot
    diff_allele_data_1 <- shiny::reactive({
        peak_info <- scan_module_outputs$diff_peak_1()
        req(peak_info)
        message("scanApp: Preparing data for difference allele plot 1.")

        reshaped <- pivot_peaks(peak_info, peak_info$marker)
        reshaped$trait <- peak_info$trait
        reshaped$plot_label <- peak_info$plot_label # Keep the label
        return(reshaped)
    })

    # NEW: Reactive for the second difference allele plot
    diff_allele_data_2 <- shiny::reactive({
        peak_info <- scan_module_outputs$diff_peak_2()
        req(peak_info)
        message("scanApp: Preparing data for difference allele plot 2.")

        reshaped <- pivot_peaks(peak_info, peak_info$marker)
        reshaped$trait <- peak_info$trait
        reshaped$plot_label <- peak_info$plot_label # Keep the label
        return(reshaped)
    })

    selected_dataset_category_reactive <- shiny::reactive({
        shiny::req(main_selected_dataset_group(), file_index_dt())
        info <- file_index_dt()[group == main_selected_dataset_group()]
        if (nrow(info) > 0) {
            return(unique(info$dataset_category)[1])
        }
        return(NULL)
    })





    output[[ns_app_controller("conditional_plot_ui")]] <- shiny::renderUI({
        category <- selected_dataset_category_reactive()
        shiny::req(category)

        # Ensure the module IDs are unique if using the same ns_app_controller
        if (category %in% c("Liver Lipids", "Clinical Traits", "Plasma Metabolites")) {
            tagList(
                manhattanPlotUI(ns_app_controller("manhattan_plot_module"))
            )
        } else if (category %in% c("Liver Genes", "Liver Isoforms")) {
            tagList(
                cisTransPlotInput(ns_app_controller("cistrans_plot_module")),
                cisTransPlotUI(ns_app_controller("cistrans_plot_module"))
            )
        } else {
            shiny::p(paste("No specific plot type configured for category:", category))
        }
    })

    # Create our own main_par structure (no longer using mainParServer)
    active_main_par <- shiny::reactive({
        # Create a compatible structure with the expected reactives
        list(
            selected_dataset = main_selected_dataset_group, # Our dataset selection
            LOD_thr = lod_threshold_rv, # Our LOD threshold
            selected_chr = selected_chromosome_rv, # Selected chromosome
            which_trait = shiny::reactive(trait_for_lod_scan_rv()), # Currently searched trait
            dataset_category = selected_dataset_category_reactive # Our dataset category
        )
    })

    # Instantiate plot modules and capture their outputs
    manhattan_plot_outputs <- manhattanPlotServer(ns_app_controller("manhattan_plot_module"),
        import_reactives = import_reactives,
        main_par = active_main_par,
        sidebar_interaction_type = sidebar_interaction_type_rv
    )

    cistrans_plot_outputs <- cisTransPlotServer(ns_app_controller("cistrans_plot_module"),
        import_reactives = import_reactives,
        main_par = active_main_par,
        peaks_cache = peaks_cache,
        sidebar_interaction_type = sidebar_interaction_type_rv
    )

    # Reactive value to store the trait selected from either plot for LOD scanning
    trait_for_lod_scan_rv <- shiny::reactiveVal(NULL)

    # Observe clicks from Manhattan plot
    shiny::observeEvent(manhattan_plot_outputs$clicked_phenotype_for_lod_scan(),
        {
            clicked_trait <- manhattan_plot_outputs$clicked_phenotype_for_lod_scan()
            if (!is.null(clicked_trait)) {
                trait_for_lod_scan_rv(clicked_trait)
                message(paste("scanApp: Manhattan plot click detected. Trait for LOD scan:", clicked_trait))
                # Potentially switch to LOD scan tab/view here in the future
            }
        },
        ignoreNULL = TRUE,
        ignoreInit = TRUE
    ) # ignoreNULL=FALSE if we want to clear on background click

    # Observe clicks from Cis/Trans plot
    shiny::observeEvent(cistrans_plot_outputs$clicked_phenotype_for_lod_scan(),
        {
            clicked_trait <- cistrans_plot_outputs$clicked_phenotype_for_lod_scan()
            if (!is.null(clicked_trait)) {
                trait_for_lod_scan_rv(clicked_trait)
                message(paste("scanApp: Cis/Trans plot click detected. Trait for LOD scan:", clicked_trait))
                # Potentially switch to LOD scan tab/view here in the future
            }
        },
        ignoreNULL = TRUE,
        ignoreInit = TRUE
    )

    # When the main selected dataset changes, clear the trait_for_lod_scan_rv
    # to prevent a scan from an old selection on a new plot type.
    # BUT preserve trait selection when switching between HC_HF variants (additive <-> interactive)
    shiny::observeEvent(main_selected_dataset_group(),
        {
            current_dataset <- main_selected_dataset_group()
            message("scanApp: Main dataset group changed. Clearing trait_for_lod_scan_rv.")

            # Check if this is just switching between HC_HF dataset variants (any HC_HF dataset type)
            is_hc_hf_dataset_switch <- !is.null(current_dataset) &&
                grepl("^HC_HF", current_dataset, ignore.case = TRUE)

            if (!is_hc_hf_dataset_switch) {
                # Only clear trait selection for real dataset changes, not interaction type changes
                trait_for_lod_scan_rv(NULL) # Clear any active LOD scan
                # Clear the search input selection (choices will be updated by the other observer)
                updateSelectizeInput(session, ns_app_controller("trait_search_input"),
                    selected = character(0)
                )
            } else {
                message("scanApp: HC_HF dataset switch detected - preserving trait selection")
            }
        },
        ignoreNULL = TRUE,
        ignoreInit = TRUE
    )

    # For debugging: observe the final trait selected for LOD scan
    shiny::observeEvent(trait_for_lod_scan_rv(),
        {
            message(paste("scanApp: trait_for_lod_scan_rv is now:", trait_for_lod_scan_rv()))
        },
        ignoreNULL = FALSE,
        ignoreInit = TRUE
    )

    # Clear caches when dataset changes (updated for new structure)
    shiny::observeEvent(main_selected_dataset_group(),
        {
            selected_ds_val <- main_selected_dataset_group()

            if (!is.null(selected_ds_val)) {
                message(paste("Clearing caches for dataset:", selected_ds_val))
                rm(list = ls(envir = trait_cache), envir = trait_cache)
                rm(list = ls(envir = peaks_cache), envir = peaks_cache)
            } else {
                message("No dataset selected or selected_dataset is NULL, not clearing caches.")
            }
        },
        ignoreNULL = FALSE,
        ignoreInit = TRUE
    )

    # Interactive Analysis section - show for all HC_HF datasets
    output[[ns_app_controller("interactive_analysis_section")]] <- shiny::renderUI({
        dataset_group <- main_selected_dataset_group()

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
                    "Sex interaction" = "sex",
                    "Diet interaction" = "diet",
                    "Sex x Diet interaction" = "sex_diet"
                )
            } else if (grepl("HC_HF.*Clinical", dataset_group, ignore.case = TRUE)) {
                available_interactions <- c(available_interactions,
                    "Sex interaction" = "sex",
                    "Diet interaction" = "diet",
                    "Sex x Diet interaction" = "sex_diet"
                )
            } else if (grepl("HC_HF.*Plasma.*Metabol", dataset_group, ignore.case = TRUE)) {
                available_interactions <- c(available_interactions,
                    "Sex interaction" = "sex",
                    "Diet interaction" = "diet",
                    "Sex x Diet interaction" = "sex_diet"
                )
            }

            tagList(
                hr(style = "border-top: 2px solid #e74c3c; margin: 15px 0;"),
                h5("🧬 Interactive Analysis", style = "color: #2c3e50; margin-bottom: 15px; font-weight: bold;"),
                shiny::selectInput(
                    ns_app_controller("interaction_type_selector"),
                    label = "Select interaction analysis:",
                    choices = available_interactions,
                    selected = if (current_selection %in% available_interactions) current_selection else "none",
                    width = "100%"
                ),
                shiny::conditionalPanel(
                    condition = paste0("input['", ns_app_controller("interaction_type_selector"), "'] != 'none'"),
                    div(
                        style = "margin-top: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 5px; border-left: 4px solid #e74c3c;",
                        p("ℹ️ Interactive analysis will show stacked plots: Interactive LOD scan (top) and Difference plot (Interactive - Additive, bottom).",
                            style = "font-size: 12px; color: #6c757d; margin: 0;"
                        )
                    )
                )
            )
        } else {
            NULL
        }
    })

    # Store the sidebar interaction type separately (for independent sidebar plot control)
    sidebar_interaction_type_rv <- shiny::reactiveVal("none")

    # NEW: Observer to decouple interaction type selection from downstream reactivity
    shiny::observeEvent(input[[ns_app_controller("interaction_type_selector")]],
        {
            req(input[[ns_app_controller("interaction_type_selector")]])
            current_interaction_type_rv(input[[ns_app_controller("interaction_type_selector")]])
            message(paste("Interaction type updated to:", input[[ns_app_controller("interaction_type_selector")]]))

            # Reset chromosome view to "All" when interaction type changes
            # This prevents being stuck on a chromosome view from a previous scan type
            shiny::updateSelectInput(session, ns_app_controller("selected_chr"), selected = "All")
            message("Reset chromosome view to 'All' due to interaction type change.")
        },
        ignoreNULL = TRUE,
        ignoreInit = TRUE
    )

    # Observer to update sidebar interaction type when input changes (independent of main UI)
    shiny::observeEvent(input[[ns_app_controller("sidebar_interaction_type")]],
        {
            new_value <- input[[ns_app_controller("sidebar_interaction_type")]]
            if (!is.null(new_value)) {
                sidebar_interaction_type_rv(new_value)
                message(paste("Updated sidebar interaction type to:", new_value, "(for sidebar plots only)"))
            }
        },
        ignoreNULL = TRUE
    )

    # Observer to preserve sidebar interaction type during UI updates
    shiny::observeEvent(trait_for_lod_scan_rv(),
        {
            # When a new trait is selected, update the sidebar selectInput to preserve the current selection
            current_sidebar_type <- sidebar_interaction_type_rv()
            if (!is.null(current_sidebar_type) && current_sidebar_type != "none") {
                shiny::updateSelectInput(session, ns_app_controller("sidebar_interaction_type"),
                    selected = current_sidebar_type
                )
                message(paste("Preserved sidebar interaction type:", current_sidebar_type, "during trait change"))
            }
        },
        ignoreNULL = FALSE,
        ignoreInit = TRUE
    )

    # Observer for the clear LOD scan button
    shiny::observeEvent(input[[ns_app_controller("clear_lod_scan_btn")]], {
        message("scanApp: Clear LOD scan button clicked.")
        trait_for_lod_scan_rv(NULL)
        # Also clear the search input field
        updateSelectizeInput(session, ns_app_controller("trait_search_input"),
            selected = character(0)
        )
    })

    # Call scanServer, passing the necessary reactives
    # scan_module_outputs will be a list of reactives/values returned by scanServer
    scan_module_outputs <- scanServer(
        id = ns_app_controller("scan_plot_module"),
        trait_to_scan = trait_for_lod_scan_rv, # Pass the reactive directly
        selected_dataset_group = mapped_dataset_for_interaction, # Use mapped dataset for interactive analysis
        import_reactives = import_reactives, # Pass the whole list of import reactives
        main_par_inputs = active_main_par, # Pass the combined main parameters (for LOD_thr, etc.)
        interaction_type_reactive = current_interaction_type_rv, # Pass the interaction type reactive
        overlay_diet_toggle = reactive(input[[ns_app_controller("overlay_diet")]]),
        overlay_sex_toggle = reactive(input[[ns_app_controller("overlay_sex")]]),
        overlay_sex_diet_toggle = reactive(input[[ns_app_controller("overlay_sex_diet")]])
    )

    # Removed clicked point details table. Clicking a peak now only updates allele effects.

    # UI for LOD Scan plot - refactored for clarity and dynamic content
    output[[ns_app_controller("lod_scan_plot_ui_placeholder")]] <- shiny::renderUI({
        if (!is.null(trait_for_lod_scan_rv())) {
            # Logic for interactive analysis dropdown (moved from its own renderUI)
            dataset_group <- main_selected_dataset_group()
            interaction_analysis_ui <- NULL

            if (!is.null(dataset_group) && grepl("^HC_HF", dataset_group, ignore.case = TRUE)) {
                current_selection <- current_interaction_type_rv()
                available_interactions <- c("None (Additive only)" = "none")

                if (grepl("HC_HF Liver Genes", dataset_group, ignore.case = TRUE)) {
                    available_interactions <- c(available_interactions, "Sex interaction" = "sex", "Diet interaction" = "diet")
                } else if (grepl("HC_HF.*Liver.*Lipid", dataset_group, ignore.case = TRUE)) {
                    available_interactions <- c(available_interactions, "Sex interaction" = "sex", "Diet interaction" = "diet", "Sex x Diet interaction" = "sex_diet")
                } else if (grepl("HC_HF.*Clinical", dataset_group, ignore.case = TRUE)) {
                    available_interactions <- c(available_interactions, "Sex interaction" = "sex", "Diet interaction" = "diet", "Sex x Diet interaction" = "sex_diet")
                } else if (grepl("HC_HF.*Plasma.*Metabol|HC_HF.*Plasma.*plasma_metabolite", dataset_group, ignore.case = TRUE)) {
                    available_interactions <- c(available_interactions, "Sex interaction" = "sex", "Diet interaction" = "diet", "Sex x Diet interaction" = "sex_diet")
                }

                interaction_analysis_ui <- tagList(
                    div(
                        style = "margin-bottom: 15px; background: #f8f9fa; padding: 10px 15px; border-radius: 4px; border: 1px solid #bdc3c7;",
                        div(
                            style = "display: flex; align-items: flex-end; gap: 15px; flex-wrap: wrap;",
                            div(
                                style = "flex: 1 1 180px; min-width: 180px;",
                                shiny::selectInput(
                                    ns_app_controller("interaction_type_selector"),
                                    label = "Select interaction analysis:",
                                    choices = available_interactions,
                                    selected = if (current_selection %in% available_interactions) current_selection else "none",
                                    width = "100%"
                                )
                            ),
                            div(
                                style = "flex: 1 1 120px; min-width: 120px;",
                                shiny::selectInput(
                                    ns_app_controller("selected_chr"),
                                    label = "Chromosome:",
                                    choices = c(
                                        "All" = "All",
                                        setNames(as.character(1:19), paste("Chr", 1:19)),
                                        "X" = "X", "Y" = "Y", "M" = "M"
                                    ),
                                    selected = "All",
                                    width = "100%"
                                )
                            ),
                            # Show overlay toggles only for additive scans
                            shiny::conditionalPanel(
                                condition = paste0("input['", ns_app_controller("interaction_type_selector"), "'] == 'none'"),
                                shiny::uiOutput(ns_app_controller("overlay_toggles_ui"))
                            ),
                            div(
                                style = "flex: 0 0 auto;",
                                shiny::actionButton(
                                    ns_app_controller("reset_chr_view"),
                                    "🌐 Reset Zoom",
                                    class = "btn btn-sm btn-secondary",
                                    style = "background: #7f8c8d; border: none; color: white; font-size: 11px; padding: 4px 8px;"
                                )
                            )
                        )
                    ),
                    shiny::conditionalPanel(
                        condition = paste0("input['", ns_app_controller("interaction_type_selector"), "'] != 'none'"),
                        div(
                            style = "margin-bottom: 15px; padding: 10px; background-color: #e8f4fd; border-radius: 5px; border-left: 4px solid #3498db;",
                            p("ℹ️ Interactive analysis will show stacked plots: Interactive LOD scan (top) and Difference plot (Interactive - Additive, bottom).",
                                style = "font-size: 12px; color: #2c3e50; margin: 0;"
                            )
                        )
                    )
                )
            }

            # Build the full plot UI including the interaction controls
            tagList(
                interaction_analysis_ui,
                scanOutput(ns_app_controller("scan_plot_module")),
                shiny::uiOutput(ns_app_controller("allele_effects_section"))
            )
        } else {
            div(
                style = "text-align: center; padding-top: 50px; color: #7f8c8d;",
                h5("No trait selected for LOD scan"),
                p("Select a trait from the 'Data Search' tab or click on a point in an overview plot to begin.")
            )
        }
    })

    # NEW: UI for overlay toggles
    output[[ns_app_controller("overlay_toggles_ui")]] <- shiny::renderUI({
        dataset_group <- main_selected_dataset_group()
        req(dataset_group)

        # Only show for HC_HF datasets that have interactive options
        if (!grepl("^HC_HF", dataset_group, ignore.case = TRUE)) {
            return(NULL)
        }

        toggles <- list()
        # Check for Diet interaction availability (Genes, Lipids, Clinical, Plasma Metabolites)
        if (any(grepl("Genes|Lipid|Clinical|Metabol", dataset_group, ignore.case = TRUE))) {
            toggles <- c(toggles, list(
                shiny::checkboxInput(ns_app_controller("overlay_diet"), "Overlay Diet", FALSE, width = "auto")
            ))
        }

        # Check for Sex interaction availability (Genes, Lipids, Clinical, Plasma Metabolites)
        if (any(grepl("Genes|Lipid|Clinical|Metabol", dataset_group, ignore.case = TRUE))) {
            toggles <- c(toggles, list(
                shiny::checkboxInput(ns_app_controller("overlay_sex"), "Overlay Sex", FALSE, width = "auto")
            ))
        }

        # Check for Sex x Diet interaction availability (Clinical, Lipids, Plasma Metabolites)
        if (any(grepl("Clinical|Lipid|Metabol", dataset_group, ignore.case = TRUE))) {
            toggles <- c(toggles, list(
                shiny::checkboxInput(ns_app_controller("overlay_sex_diet"), "Overlay Sex x Diet", FALSE, width = "auto")
            ))
        }

        if (length(toggles) > 0) {
            div(
                style = "display: flex; gap: 8px; align-items: center; margin-left: 15px; font-size: 12px; line-height: 1.1; white-space: nowrap; flex-wrap: nowrap;",
                toggles
            )
        } else {
            NULL
        }
    })

    # Dynamic title for the main LOD scan card
    output[[ns_app_controller("main_plot_title")]] <- shiny::renderUI({
        trait <- trait_for_lod_scan_rv()
        title_text <- if (!is.null(trait)) {
            paste("LOD Scan for Trait:", trait)
        } else {
            "LOD Scan - Detailed View"
        }
        h4(title_text, style = "font-weight: bold; margin-bottom: 0;")
    })

    # Render allele effects section conditionally
    output[[ns_app_controller("allele_effects_section")]] <- shiny::renderUI({
        additive_peak <- scan_module_outputs$selected_peak()
        diff_peak_1 <- scan_module_outputs$diff_peak_1()
        diff_peak_2 <- scan_module_outputs$diff_peak_2()

        # View 1: Additive Peak - show only allele effects plot when a peak is selected
        if (!is.null(additive_peak)) {
            message(paste("scanApp: Rendering SINGLE allele effects for peak:", additive_peak$marker))
            return(tagList(
                hr(style = "margin: 20px 0; border-top: 2px solid #3498db;"),
                shiny::plotOutput(ns_app_controller("allele_effects_plot_output"), height = "450px", width = "600px") %>%
                    shinycssloaders::withSpinner(type = 8, color = "#3498db")
            ))
        }

        # View 2: Comparative Difference Peak Details - unchanged (still plots only)
        if (!is.null(diff_peak_1) || !is.null(diff_peak_2)) {
            message("scanApp: Rendering SIDE-BY-SIDE allele effects for difference plot click.")

            ui_elements <- list()

            if (!is.null(diff_peak_1)) {
                ui_elements <- c(ui_elements, list(
                    bslib::card(
                        bslib::card_header(textOutput(ns_app_controller("diff_plot_title_1"))),
                        bslib::card_body(
                            shiny::plotOutput(ns_app_controller("diff_allele_plot_1"), height = "400px") %>%
                                shinycssloaders::withSpinner(type = 8, color = "#e74c3c")
                        )
                    )
                ))
            }

            if (!is.null(diff_peak_2)) {
                ui_elements <- c(ui_elements, list(
                    bslib::card(
                        bslib::card_header(textOutput(ns_app_controller("diff_plot_title_2"))),
                        bslib::card_body(
                            shiny::plotOutput(ns_app_controller("diff_allele_plot_2"), height = "400px") %>%
                                shinycssloaders::withSpinner(type = 8, color = "#e74c3c")
                        )
                    )
                ))
            }

            return(bslib::layout_columns(
                col_widths = c(6, 6),
                !!!ui_elements
            ))
        }

        # No peak selected -> show nothing
        return(NULL)
    })

    # Removed peak_info_display; no info table is shown, only allele effects plot on click.

    # Render the allele effects plot
    output[[ns_app_controller("allele_effects_plot_output")]] <- shiny::renderPlot({
        effects_data <- allele_effects_data()

        message(paste("scanApp: Rendering allele effects plot. Data available:", !is.null(effects_data)))

        if (is.null(effects_data)) {
            message("scanApp: No allele effects data - showing placeholder")
            # Create a placeholder plot when no data
            ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::labs(title = "No strain effects data available for this peak") +
                ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, color = "#7f8c8d"))
        } else {
            message(paste("scanApp: Creating allele effects plot with", nrow(effects_data), "data points"))
            # Use the ggplot_alleles function to create the plot
            ggplot_alleles(effects_data)
        }
    })

    # NEW: Render logic for side-by-side plots
    output[[ns_app_controller("diff_plot_title_1")]] <- renderText({
        data <- diff_allele_data_1()
        req(data)
        unique(data$plot_label)
    })

    output[[ns_app_controller("diff_allele_plot_1")]] <- shiny::renderPlot({
        data <- diff_allele_data_1()
        req(data)
        ggplot_alleles(data)
    })

    output[[ns_app_controller("diff_plot_title_2")]] <- renderText({
        data <- diff_allele_data_2()
        req(data)
        unique(data$plot_label)
    })

    output[[ns_app_controller("diff_allele_plot_2")]] <- shiny::renderPlot({
        data <- diff_allele_data_2()
        req(data)
        ggplot_alleles(data)
    })

    # --- FIXED Trait Search Logic ---
    # Store a flag to prevent auto-search immediately after dataset changes
    dataset_just_changed <- shiny::reactiveVal(FALSE)

    # Reset the flag after a short delay
    observe({
        if (dataset_just_changed()) {
            invalidateLater(1000) # Wait 1 second
            dataset_just_changed(FALSE)
        }
    })

    # --- AUTO-TRIGGER SEARCH WHEN TRAIT IS SELECTED FROM DROPDOWN ---
    # This allows users to skip clicking the button - just select from dropdown and it searches automatically
    # But NOT immediately after dataset changes to prevent slow loading
    observeEvent(input[[ns_app_controller("trait_search_input")]],
        {
            searched_trait <- input[[ns_app_controller("trait_search_input")]]

            # Only auto-search if:
            # 1. A valid trait is selected and dataset is available
            # 2. Dataset didn't just change (to prevent slow loading)
            # 3. This is a new selection by the user (not just preservation from dataset change)
            if (!is.null(searched_trait) && nzchar(searched_trait) &&
                !is.null(main_selected_dataset_group()) &&
                !dataset_just_changed()) {
                if (is.null(trait_for_lod_scan_rv()) || !identical(trait_for_lod_scan_rv(), searched_trait)) {
                    message(paste("scanApp: Auto-search triggered for trait:", searched_trait))
                    trait_for_lod_scan_rv(searched_trait)

                    # Show notification for auto-search
                    shiny::showNotification(
                        paste("Auto-searching for trait:", searched_trait),
                        type = "message",
                        duration = 2
                    )
                }
            } else if (dataset_just_changed()) {
                message(paste("scanApp: Trait selection preserved but auto-search skipped due to recent dataset change. Use search button to trigger LOD scan."))
            }
        },
        ignoreNULL = TRUE,
        ignoreInit = TRUE
    )

    # ====== TRAIT SEARCH DROPDOWN UPDATE LOGIC ======
    # Updates trait search choices when dataset changes (server-side selectize for performance)
    shiny::observeEvent(shiny::req(main_selected_dataset_group()), {
        shiny::req(import_reactives()) # Ensure import data is available
        current_ds <- main_selected_dataset_group()

        message(paste("scanApp: Updating trait search choices for dataset:", current_ds))

        # SET FLAG to prevent immediate auto-search after dataset change
        dataset_just_changed(TRUE)

        # PRESERVE CURRENT SELECTION: Store the current trait selection before updating
        current_trait_selection <- input[[ns_app_controller("trait_search_input")]]
        message(paste("scanApp: Current trait selection before dataset change:", current_trait_selection %||% "None"))

        # RESET TRAIT SEARCH INPUT: Prevent conflicts during update
        shiny::freezeReactiveValue(input, ns_app_controller("trait_search_input")) # Temporarily freeze reactive
        shiny::updateSelectizeInput(session, ns_app_controller("trait_search_input"),
            choices = character(0),
            selected = character(0),
            options = list(placeholder = "Loading traits...")
        )

        # GET NEW TRAIT CHOICES: Based on selected dataset
        choices <- get_trait_choices(import_reactives(), current_ds) # External helper function

        if (!is.null(choices) && length(choices) > 0) {
            message(paste("scanApp: Found", length(choices), "traits for dataset:", current_ds))

            # CHECK IF PREVIOUS SELECTION EXISTS IN NEW DATASET
            trait_to_select <- NULL
            if (!is.null(current_trait_selection) && nzchar(current_trait_selection)) {
                if (current_trait_selection %in% choices) {
                    trait_to_select <- current_trait_selection
                    message(paste("scanApp: Preserving trait selection:", current_trait_selection, "- found in new dataset"))
                } else {
                    message(paste("scanApp: Previous trait selection:", current_trait_selection, "- not found in new dataset, clearing selection"))
                }
            }

            # UPDATE TRAIT SEARCH DROPDOWN: With server-side processing for large lists
            shiny::updateSelectizeInput(session, ns_app_controller("trait_search_input"),
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
            message(paste("scanApp: No traits found for dataset:", current_ds))
            # No traits available for this dataset
            shiny::updateSelectizeInput(session, ns_app_controller("trait_search_input"),
                choices = character(0),
                selected = character(0),
                options = list(placeholder = "No traits available for this dataset")
            )
        }
    })

    # ====== CHROMOSOME ZOOM FUNCTIONALITY ======
    # Observer for chromosome selection dropdown
    observeEvent(input[[ns_app_controller("selected_chr")]],
        {
            selected_chr <- input[[ns_app_controller("selected_chr")]]
            if (!is.null(selected_chr) && selected_chr != "All") {
                shiny::showNotification(
                    paste("Zoomed to chromosome", selected_chr),
                    type = "message",
                    duration = 2
                )
            }
        },
        ignoreInit = TRUE
    )

    # Observer for "Show All Chromosomes" button
    observeEvent(input[[ns_app_controller("reset_chr_view")]], {
        message("scanApp: Resetting to show all chromosomes")
        shiny::updateSelectInput(session, ns_app_controller("selected_chr"),
            selected = "All"
        )
        shiny::showNotification(
            "Showing all chromosomes",
            type = "message",
            duration = 2
        )
    })

    # ====== ADDITIONAL ANALYSES PLOT OUTPUTS ======
    # Profile Plot output
    profilePlotServer(
        ns_app_controller("profile_plot_module"),
        selected_dataset_category = selected_dataset_category_reactive,
        trait_to_profile = trait_for_lod_scan_rv # Pass the selected trait
    )

    # Correlation module
    correlationServer(
        ns_app_controller("correlation_module"),
        import_reactives = import_reactives,
        main_par = active_main_par
    )

    # Peak Analysis dropdown for sidebar - separate from main UI interaction controls
    output[[ns_app_controller("peak_selection_sidebar")]] <- shiny::renderUI({
        dataset_group <- main_selected_dataset_group()
        # Only build the interaction menu; peak selection removed here
        if (!is.null(dataset_group) && grepl("^HC_HF", dataset_group, ignore.case = TRUE)) {
            available_interactions <- c("None (Additive only)" = "none")
            if (grepl("HC_HF Liver Genes", dataset_group, ignore.case = TRUE)) {
                available_interactions <- c(available_interactions,
                    "Sex interaction" = "sex",
                    "Diet interaction" = "diet"
                )
            } else if (grepl("HC_HF.*Liver.*Lipid", dataset_group, ignore.case = TRUE)) {
                available_interactions <- c(available_interactions,
                    "Sex interaction" = "sex",
                    "Diet interaction" = "diet",
                    "Sex x Diet interaction" = "sex_diet"
                )
            } else if (grepl("HC_HF.*Clinical", dataset_group, ignore.case = TRUE)) {
                available_interactions <- c(available_interactions,
                    "Sex interaction" = "sex",
                    "Diet interaction" = "diet",
                    "Sex x Diet interaction" = "sex_diet"
                )
            } else if (grepl("HC_HF.*Plasma.*Metabol", dataset_group, ignore.case = TRUE)) {
                available_interactions <- c(available_interactions,
                    "Sex interaction" = "sex",
                    "Diet interaction" = "diet",
                    "Sex x Diet interaction" = "sex_diet"
                )
            }

            tagList(
                h6("Sidebar Plot Analysis:", style = "color: #2c3e50; margin-bottom: 8px; font-weight: bold; font-size: 12px;"),
                shiny::selectInput(
                    ns_app_controller("sidebar_interaction_type"),
                    label = NULL,
                    choices = available_interactions,
                    selected = if (sidebar_interaction_type_rv() %in% available_interactions) sidebar_interaction_type_rv() else "none",
                    width = "100%"
                )
            )
        } else {
            NULL
        }
    })

    # Interactive Analysis section - show for all HC_HF datasets
}

# Launch the app
shinyApp(ui = ui, server = server)
