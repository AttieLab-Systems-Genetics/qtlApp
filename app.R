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

# Source files in specific order (updated paths)
source("R/utils/helpers.R")
source("R/data/data_handling.R")
source("R/data/import_data.R")
source("R/modules/importApp.R")
source("R/modules/downloadApp.R")
source("R/modules/cisTransPlotApp.R")
source("R/modules/manhattanPlotApp.R")
source("R/ui/ui_styles.R")
source("R/plots/plot_enhancements.R")
source("R/plots/plot_null.R")
source("R/data/peak_finder.R")
source("R/data/trait_scan.R")
source("R/plots/ggplot_alleles.R")
source("R/plots/ggplot_qtl_scan.R")
source("R/plots/ggplotly_qtl_scan.R")
source("R/plots/peak_info.R")
source("R/plots/QTL_plot_visualizer.R")
source("R/data/fst_rows.R")
source("R/modules/traitApp.R")
source("R/modules/traitProcessingModule.R")
source("R/modules/scanPlotModule.R") # Source our scan module
source("R/modules/profilePlotApp.R") # Profile plot module
source("R/modules/interactiveAnalysisModule.R") # Interactive analysis controls/module
source("R/modules/splitAlleleEffectsModule.R") # Split allele effects module
source("R/modules/overlayControlsModule.R") # Overlay controls module
source("R/modules/lodThresholdModule.R") # LOD threshold module
source("R/ui/mainUI.R") # Main UI module

# Set maximum file upload size
options(shiny.maxRequestSize = 20000 * 1024^2) # 20 GB

# UI Definition
ui <- mainUI()

# Server Definition
server <- function(input, output, session) {
    ns_app_controller <- shiny::NS("app_controller")

    # Source helper for allele plots if not already available
    if (!exists("ggplot_alleles", mode = "function")) {
        source("R/plots/ggplot_alleles.R")
    }

    # Define the %||% operator for null coalescing
    `%||%` <- function(a, b) if (!is.null(a)) a else b

    # Use interactiveAnalysisModule for dataset mapping; remove local duplicate helper

    trait_cache <- new.env(parent = emptyenv())
    peaks_cache <- new.env(parent = emptyenv())

    import_reactives <- importServer("import")

    # Placeholder; module will be instantiated after dataset reactive is defined
    interactive_analysis <- NULL
    current_interaction_type_rv <- NULL
    mapped_dataset_for_interaction <- NULL
    scan_type <- NULL
    show_stacked_plots <- NULL

    # LOD Threshold module (replaces inline slider)
    output[[ns_app_controller("lod_threshold_slider")]] <- shiny::renderUI({
        lodThresholdUI(ns_app_controller("lod_thr"))
    })
    lod_thr_module <- lodThresholdServer(
        id = ns_app_controller("lod_thr"),
        interaction_type_reactive = sidebar_interaction_type_rv
    )
    lod_threshold_rv <- lod_thr_module$lod_threshold

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

    # Initialize interactive analysis module now that dataset reactive exists (single instance)
    interactive_analysis <- interactiveAnalysisServer(
        id = ns_app_controller("interactive_analysis_module"),
        selected_dataset_reactive = main_selected_dataset_group
    )
    current_interaction_type_rv <- interactive_analysis$interaction_type
    mapped_dataset_for_interaction <- interactive_analysis$mapped_dataset
    scan_type <- interactive_analysis$scan_type
    show_stacked_plots <- shiny::reactive({
        type <- current_interaction_type_rv()
        !is.null(type) && type != "none"
    })

    # Reset chromosome view when interaction type changes to avoid stale zoom state
    shiny::observeEvent(current_interaction_type_rv(), {
        shiny::updateSelectInput(session, ns_app_controller("selected_chr"), selected = "All")
    }, ignoreInit = TRUE)

    # Ensure scanServer only runs when mapped dataset is stable/non-null
    mapped_dataset_stable <- shiny::reactive({
        ds <- mapped_dataset_for_interaction()
        shiny::req(ds)
        ds
    }) %>% shiny::debounce(100)

    # Mount overlay controls early so we can pass their reactives to scanServer
    overlay_module <- overlayControlsServer(
        id = ns_app_controller("overlay_controls"),
        selected_dataset_group_reactive = main_selected_dataset_group
    )

    # removed: duplicate mapping and scan_type logic (handled by interactiveAnalysisModule)

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


    # Store the sidebar interaction type separately (for independent sidebar plot control)
    sidebar_interaction_type_rv <- shiny::reactiveVal("none")


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

    # Call scanServer, passing the necessary reactives (after interactive_analysis is available)
    scan_module_outputs <- scanServer(
        id = ns_app_controller("scan_plot_module"),
        trait_to_scan = trait_for_lod_scan_rv,
        selected_dataset_group = mapped_dataset_stable,
        import_reactives = import_reactives,
        main_par_inputs = active_main_par,
        interaction_type_reactive = current_interaction_type_rv,
        overlay_diet_toggle = overlay_module$diet_toggle,
        overlay_sex_toggle = overlay_module$sex_toggle,
        overlay_sex_diet_toggle = overlay_module$sex_diet_toggle
    )

    # Render the LOD scan click details table
    output[[ns_app_controller("lod_scan_click_table")]] <- DT::renderDT({
        # Check if we have scan module outputs and clicked point details
        if (!is.null(scan_module_outputs) && !is.null(scan_module_outputs$clicked_point_details)) {
            clicked_details <- scan_module_outputs$clicked_point_details()

            if (!is.null(clicked_details) && nrow(clicked_details) > 0) {
                # Format the clicked point data for display
                display_data <- data.frame(
                    Chromosome = if ("chr" %in% colnames(clicked_details)) clicked_details$chr else "N/A",
                    Marker = if ("markers" %in% colnames(clicked_details)) clicked_details$markers else "N/A",
                    Position = if ("position" %in% colnames(clicked_details)) paste0(clicked_details$position, " Mb") else "N/A",
                    LOD = if ("LOD" %in% colnames(clicked_details)) clicked_details$LOD else "N/A"
                )

                message(paste("scanApp: Displaying clicked point details for marker:", display_data$Marker))

                return(DT::datatable(
                    display_data,
                    options = list(
                        dom = "t",
                        paging = FALSE,
                        searching = FALSE,
                        columnDefs = list(list(targets = "_all", className = "dt-center"))
                    ),
                    rownames = FALSE,
                    selection = "none",
                    class = "compact hover"
                ))
            }
        }

        # Default message when no point is clicked
        return(DT::datatable(
            data.frame(Info = "Click on a point in the LOD scan plot to see details"),
            options = list(dom = "t", paging = FALSE, searching = FALSE),
            rownames = FALSE,
            selection = "none"
        ))
    })

    # UI for LOD Scan plot - refactored for clarity and dynamic content
    output[[ns_app_controller("lod_scan_plot_ui_placeholder")]] <- shiny::renderUI({
        if (!is.null(trait_for_lod_scan_rv())) {
            # Build the full plot UI including the interaction controls from the module
            is_additive <- identical(interactive_analysis$interaction_type(), "none")
            tagList(
                div(
                    style = "margin-bottom: 15px; background: #f8f9fa; padding: 10px 15px; border-radius: 4px; border: 1px solid #bdc3c7;",
                    div(
                        style = "display: flex; align-items: flex-end; gap: 15px; flex-wrap: wrap;",
                        div(
                            style = "flex: 1 1 180px; min-width: 180px;",
                            interactiveAnalysisUI(ns_app_controller("interactive_analysis_module"))
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
                        if (isTRUE(is_additive)) overlayControlsUI(ns_app_controller("overlay_controls")),
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
                if (!isTRUE(is_additive)) div(
                    style = "margin-bottom: 15px; padding: 10px; background-color: #e8f4fd; border-radius: 5px; border-left: 4px solid #3498db;",
                    p("ℹ️ Interactive analysis will show stacked plots: Interactive LOD scan (top) and Difference plot (Interactive - Additive, bottom).",
                        style = "font-size: 12px; color: #2c3e50; margin: 0;"
                    )
                ),
                scanOutput(ns_app_controller("scan_plot_module")),
                div(
                    style = "margin-top: 15px;",
                    DT::DTOutput(ns_app_controller("lod_scan_click_table"))
                ),
                splitAlleleEffectsUI(ns_app_controller("split_allele_effects"))
            )
        } else {
            div(
                style = "text-align: center; padding-top: 50px; color: #7f8c8d;",
                h5("No trait selected for LOD scan"),
                p("Select a trait from the 'Data Search' tab or click on a point in an overview plot to begin.")
            )
        }
    })

    # Mount overlay controls module and wire toggle outputs to scanServer
    overlay_module <- overlayControlsServer(
        id = ns_app_controller("overlay_controls"),
        selected_dataset_group_reactive = main_selected_dataset_group
    )

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

    # Mount allele effects module server (single or side-by-side)
    splitAlleleEffectsServer(
      id = ns_app_controller("split_allele_effects"),
      selected_peak_reactive = scan_module_outputs$selected_peak,
      diff_peak_1_reactive = scan_module_outputs$diff_peak_1,
      diff_peak_2_reactive = scan_module_outputs$diff_peak_2
    )

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

    # Correlation Plot output
    output[[ns_app_controller("correlation_plot_output")]] <- plotly::renderPlotly({
        plotly::plot_ly() %>%
            plotly::add_annotations(
                text = "Correlation Analysis Coming Soon",
                x = 0.5,
                y = 0.5,
                xref = "paper",
                yref = "paper",
                showarrow = FALSE,
                font = list(size = 20, color = "#2c3e50")
            ) %>%
            plotly::layout(
                title = "Correlation Analysis",
                showlegend = FALSE,
                xaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
                yaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
            )
    })

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
