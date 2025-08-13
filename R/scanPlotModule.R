#' Scan Plot and Details Module
#'
#' This module handles the rendering of the main LOD scan plot, including
#' interactive and difference plots, overlays, and click-based peak analysis.
#'
#' @param id shiny identifier
#' @param trait_to_scan reactive containing the trait to be scanned
#' @param selected_dataset_group reactive with the name of the selected dataset group
#' @param import_reactives reactive list with imported data (file_directory, markers)
#' @param main_par_input-s reactive list with main parameters (LOD_thr, selected_chr)
#' @param interaction_type_reactive reactive specifying the interaction type ('none', 'sex', 'diet')
#' @param overlay_diet_toggle reactive boolean for diet overlay
#' @param overlay_sex_toggle reactive boolean for sex overlay
#'
#' @importFrom shiny moduleServer NS reactive reactiveVal observeEvent updateNumericInput req
#'             observe tagList h5 div uiOutput renderUI conditionalPanel checkboxInput downloadHandler
#' @importFrom plotly plotlyOutput renderPlotly ggplotly event_data layout config
#' @importFrom shinycssloaders withSpinner
#' @importFrom dplyr filter select mutate across where arrange desc inner_join everything all_of
#' @importFrom ggplot2 ggplot labs theme_void theme element_text geom_hline ggsave cairo_pdf
#' @importFrom stringr str_to_title
#' @importFrom data.table data.table as.data.table
#' @importFrom DT renderDT datatable formatStyle
#' @importFrom htmltools tags
#'
#' @return A reactiveValues object containing:
#'   - filename: A reactive for generating download filenames.
#'   - tables: A reactiveValues object with the scan table.
#'   - plots: A reactiveValues object with the ggplot object for the scan.
#'   - clicked_point_details: A reactiveVal with details of a clicked point on the plot.
#'   - selected_peak: A reactiveVal containing the full data for the selected peak.
#' @export
scanServer <- function(id, trait_to_scan, selected_dataset_group, import_reactives, main_par_inputs, interaction_type_reactive = NULL, overlay_diet_toggle = reactive({
                           FALSE
                       }), overlay_sex_toggle = reactive({
                           FALSE
                       }), overlay_sex_diet_toggle = reactive({
                           FALSE
                       })) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Create local cache for peaks data
        local_peaks_cache <- new.env(parent = emptyenv())

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
            # HC_HF Liver Lipids (supports Sex, Diet, and Sex x Diet interactions)
            else if (grepl("HC_HF.*Liver.*Lipid", base_dataset, ignore.case = TRUE)) {
                if (interaction_type == "diet") {
                    return("HC_HF Liver Lipids, interactive (Diet)")
                } else if (interaction_type == "sex") {
                    return("HC_HF Liver Lipids, interactive (Sex)")
                } else if (interaction_type == "sex_diet") {
                    return("HC_HF Liver Lipids, interactive (Sex_Diet)")
                }
                # No Sex-only interaction available for Liver Lipids - return original
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
            }
            # HC_HF Plasma Metabolites (supports Sex and Sex x Diet interactions; allow Diet if present)
            else if (grepl("HC_HF.*Plasma.*Metabol", base_dataset, ignore.case = TRUE)) {
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

        plot_width_rv <- shiny::reactiveVal(1200)
        plot_height_rv <- shiny::reactiveVal(600)
        use_alternating_colors_rv <- shiny::reactiveVal(TRUE)
        clicked_plotly_point_details_lod_scan_rv <- shiny::reactiveVal(NULL)
        selected_peak_rv <- shiny::reactiveVal(NULL) # For single additive plot
        diff_peak_1_rv <- shiny::reactiveVal(NULL) # For side-by-side allele plot 1
        diff_peak_2_rv <- shiny::reactiveVal(NULL) # For side-by-side allele plot 2
        numb_mice_rv <- shiny::reactiveVal(NULL) # To store the number of mice

        # Helper function to get paths for interactive peak files
        get_interactive_peak_filepaths <- function(dataset_group, interaction_type) {
            base_path <- "/data/dev/miniViewer_3.0/" # Assuming path from repository rules

            # Derive base name (e.g., "HC_HF Liver Genes") from full dataset name
            base_name <- gsub(",[[:space:]]*(interactive|additive).*$", "", dataset_group, ignore.case = TRUE)
            base_name <- trimws(base_name)

            dataset_component <- NULL
            if (grepl("Liver Genes", base_name, ignore.case = TRUE)) {
                dataset_component <- "liver_genes"
            } else if (grepl("Liver Lipids", base_name, ignore.case = TRUE)) {
                dataset_component <- "liver_lipids"
            } else if (grepl("Clinical Traits", base_name, ignore.case = TRUE)) {
                dataset_component <- "clinical_traits"
            } else if (grepl("Plasma Metabolites", base_name, ignore.case = TRUE)) {
                dataset_component <- "plasma_metabolites"
            } else if (grepl("Liver Isoforms", base_name, ignore.case = TRUE)) {
                dataset_component <- "liver_isoforms"
            } else {
                warning(paste("Unsupported dataset group for interactive peak file lookup:", dataset_group))
                return(NULL)
            }

            file1 <- NULL
            file2 <- NULL
            labels <- NULL

            if (interaction_type == "sex") {
                file1 <- paste0("DO1200_", dataset_component, "_qtlxsex_peaks_in_female_mice_additive.csv")
                file2 <- paste0("DO1200_", dataset_component, "_qtlxsex_peaks_in_male_mice_additive.csv")
                labels <- c("Female", "Male")
            } else if (interaction_type == "diet") {
                file1 <- paste0("DO1200_", dataset_component, "_qtlxdiet_peaks_in_HC_mice_additive.csv")
                file2 <- paste0("DO1200_", dataset_component, "_qtlxdiet_peaks_in_HF_mice_additive.csv")
                labels <- c("HC Diet", "HF Diet")
            } else if (interaction_type == "sex_diet") {
                # No dedicated peak files defined for Sex x Diet; handle in click code
                return(NULL)
            } else {
                return(NULL)
            }

            return(list(
                file1 = file.path(base_path, file1),
                file2 = file.path(base_path, file2),
                labels = labels
            ))
        }

        shiny::observeEvent(input$plot_width,
            {
                plot_width_rv(input$plot_width)
            },
            ignoreNULL = TRUE
        )
        shiny::observeEvent(input$plot_height,
            {
                plot_height_rv(input$plot_height)
            },
            ignoreNULL = TRUE
        )

        shiny::observeEvent(input$preset_1to1, {
            shiny::updateNumericInput(session, "plot_width", value = 800)
            shiny::updateNumericInput(session, "plot_height", value = 800)
        })
        shiny::observeEvent(input$preset_3to2, {
            shiny::updateNumericInput(session, "plot_width", value = 900)
            shiny::updateNumericInput(session, "plot_height", value = 600)
        })
        shiny::observeEvent(input$preset_16to9, {
            shiny::updateNumericInput(session, "plot_width", value = 1280)
            shiny::updateNumericInput(session, "plot_height", value = 720)
        })

        shiny::observeEvent(input$color_toggle, {
            if (is.logical(input$color_toggle)) {
                use_alternating_colors_rv(input$color_toggle)
            } else {
                use_alternating_colors_rv(!use_alternating_colors_rv())
            }
        })

        current_trait_for_scan <- shiny::reactive({
            shiny::req(trait_to_scan())
            trait_val <- trait_to_scan()
            message(paste("scanServer: Received trait_to_scan value:", trait_val, "(this will be passed to trait_scan)."))
            return(trait_val)
        })

        scans <- shiny::reactive({
            req(current_trait_for_scan(), selected_dataset_group())
            trait_val <- current_trait_for_scan()
            dataset_group_val <- selected_dataset_group()

            message(paste0("scanServer (within scans reactive): ABOUT TO CALL trait_scan. Trait: '", trait_val, "', Dataset Group: '", dataset_group_val, "'"))

            file_dir_val <- import_reactives()$file_directory
            req(file_dir_val)
            if (!is.data.frame(file_dir_val) || !("group" %in% names(file_dir_val))) {
                stop("scanServer: import_reactives()$file_directory is not a valid data frame or missing 'group' column before calling trait_scan.")
            }

            # DEBUG: Show what files will be used for this dataset
            debug_files <- subset(file_dir_val, group == dataset_group_val & file_type == "scans")
            message(paste0("scanServer: Files that will be used for dataset '", dataset_group_val, "':"))
            for (i in 1:min(3, nrow(debug_files))) {
                message(paste0("  - Chr ", debug_files$ID_code[i], ": ", basename(debug_files$File_path[i])))
            }

            result_list <- tryCatch(
                {
                    trait_scan(
                        file_dir = file_dir_val,
                        selected_dataset = dataset_group_val,
                        selected_trait = trait_val,
                        cache_env = NULL
                    )
                },
                error = function(e) {
                    message(paste0("scanServer (within scans reactive): ERROR DURING trait_scan CALL. Trait: '", trait_val, "', Dataset Group: '", dataset_group_val, "'. Error: ", e$message))
                    return(NULL)
                }
            )

            if (is.null(result_list) || is.null(result_list$scan_data)) {
                numb_mice_rv(NULL)
                return(data.table::data.table())
            }

            numb_mice_rv(result_list$numb_mice)

            scan_data <- result_list$scan_data
            message(paste0("scanServer (within scans reactive): RETURNED FROM trait_scan. Trait: '", trait_val, "', Dataset Group: '", dataset_group_val, "'. Result class: ", class(scan_data), ", Result nrows: ", if (!is.null(scan_data) && (is.data.frame(scan_data) || is.data.table(scan_data))) nrow(scan_data) else "N/A"))

            # Additional check: if result is NULL due to error, or if it's an empty data frame, handle appropriately.
            if (is.null(scan_data) || ((is.data.frame(scan_data) || is.data.table(scan_data)) && nrow(scan_data) == 0)) {
                message(paste0("scanServer: trait_scan returned NULL or empty for Trait: '", trait_val, "', Dataset: '", dataset_group_val, "'. Propagating as empty result."))
                # Return an empty data.table or data.frame as expected by downstream reactives to prevent crashes
                # Make sure it has the columns expected by QTL_plot_visualizer if possible, or handle this there.
                return(data.table::data.table())
            }

            scan_data
        }) %>% shiny::debounce(150) # Add debouncing to prevent rapid re-computation

        # NEW: Reactive to determine the static line threshold based on context
        static_lod_threshold_line <- reactive({
            is_lod_scan_view <- !is.null(trait_to_scan())
            interaction_type <- if (!is.null(interaction_type_reactive)) interaction_type_reactive() else "none"

            if (is_lod_scan_view) {
                # For the main LOD scan plot
                if (interaction_type == "sex") {
                    return(10.5)
                }
                if (interaction_type == "diet") {
                    return(10.5)
                }
                if (interaction_type == "sex_diet") {
                    return(15.7)
                }
                return(7.5) # Additive
            } else {
                # This part is for any other context, but this reactive is inside scanServer,
                # which is only active for LOD scans. So this is fallback.
                return(7.5)
            }
        })

        scan_table <- shiny::reactive({
            # Access the list of reactives first
            main_par_list <- main_par_inputs()
            shiny::req(
                scans(),
                current_trait_for_scan(),
                main_par_list,
                import_reactives()$markers
            )

            # Streamlined processing - reduce debug messages
            scan_data <- scans()
            trait_val <- current_trait_for_scan()
            lod_val <- static_lod_threshold_line() # Use static threshold for data processing
            markers_data <- import_reactives()$markers

            result <- QTL_plot_visualizer(scan_data, trait_val, lod_val, markers_data)
            result
        }) %>% shiny::debounce(150) # Add debouncing to prevent rapid re-computation

        scan_table_chr <- shiny::reactive({
            main_par_list <- main_par_inputs()
            shiny::req(
                scan_table(),
                main_par_list,
                main_par_list$selected_chr,
                main_par_list$selected_chr()
            )
            current_scan_table <- scan_table()
            selected_chromosome <- main_par_list$selected_chr()
            if (selected_chromosome == "All") {
                current_scan_table
            } else {
                sel_chr_num <- selected_chromosome
                if (selected_chromosome == "X") sel_chr_num <- 20
                if (selected_chromosome == "Y") sel_chr_num <- 21
                if (selected_chromosome == "M") sel_chr_num <- 22
                sel_chr_num <- as.numeric(sel_chr_num)

                dplyr::filter(current_scan_table, chr == sel_chr_num)
            }
        }) %>% shiny::debounce(200) # Debounce scan table chr filtering

        # Observer to find the highest peak and set it as the default selected peak
        # This triggers whenever a new scan is successfully completed and loaded
        shiny::observeEvent(scan_table_chr(),
            {
                # Ensure we have a trait and scan data before proceeding
                shiny::req(current_trait_for_scan(), nrow(scan_table_chr()) > 0)

                trait_val <- current_trait_for_scan()
                dataset_group_val <- selected_dataset_group()

                message(paste("scanServer: New scan data available for trait", trait_val, ". Finding highest peak to set as default."))

                # It's safer to wrap the function call in a tryCatch in case helpers aren't loaded
                trait_type_val <- tryCatch(
                    {
                        get_trait_type(import_reactives(), dataset_group_val)
                    },
                    error = function(e) {
                        warning(paste("Could not get trait type for default peak finding:", e$message))
                        "genes" # Fallback to a sensible default
                    }
                )

                # Get all peaks for the current trait
                all_peaks <- peak_finder(
                    file_dir = import_reactives()$file_directory,
                    selected_dataset = dataset_group_val,
                    selected_trait = trait_val,
                    trait_type = trait_type_val,
                    cache_env = local_peaks_cache,
                    use_cache = TRUE
                )

                if (!is.null(all_peaks) && nrow(all_peaks) > 0) {
                    # Find the highest peak and set it as the selected one
                    highest_peak <- all_peaks[which.max(all_peaks$qtl_lod), ]
                    selected_peak_rv(highest_peak)
                    message(paste("scanServer: Default peak set to highest LOD peak:", highest_peak$marker))

                    # Also set the click details to show the highest peak info by default
                    default_click_details <- data.frame(
                        markers = highest_peak$marker,
                        chr = highest_peak$qtl_chr,
                        position = round(highest_peak$qtl_pos, 3),
                        LOD = round(highest_peak$qtl_lod, 3)
                    )
                    clicked_plotly_point_details_lod_scan_rv(default_click_details)
                    message(paste("scanServer: Default click details set for highest peak:", highest_peak$marker))
                } else {
                    selected_peak_rv(NULL) # Clear selection if no peaks are found
                    clicked_plotly_point_details_lod_scan_rv(NULL) # Clear click details too
                    message("scanServer: No peaks found for this trait, clearing selected peak and click details.")
                }
            },
            ignoreNULL = TRUE,
            ignoreInit = TRUE
        )


        # Reactive to store additive plot data for difference calculations
        additive_scan_data_rv <- shiny::reactiveVal(NULL)

        # Observer to clear the additive data cache whenever the selected trait changes
        shiny::observeEvent(current_trait_for_scan(),
            {
                # This ensures we don't use stale additive data from a previous trait
                message("scanServer: New trait selected, clearing additive data cache and diff peaks.")
                additive_scan_data_rv(NULL)
                diff_peak_1_rv(NULL)
                diff_peak_2_rv(NULL)
            },
            ignoreNULL = FALSE,
            ignoreInit = TRUE
        )


        # REVISED ADDITIVE DATA HANDLING:
        # This observer handles both manual and automatic caching of additive data.
        observe({
            trait_val <- current_trait_for_scan()
            dataset_group_val <- selected_dataset_group()
            interaction_type <- if (!is.null(interaction_type_reactive)) interaction_type_reactive() else "none"

            # Scenario 1: User is viewing an additive scan. Cache it.
            if (interaction_type == "none" && !grepl("interactive", dataset_group_val, ignore.case = TRUE)) {
                current_scan <- scan_table()
                if (!is.null(current_scan) && nrow(current_scan) > 0) {
                    additive_scan_data_rv(current_scan)
                    message("scanServer: Additive data cached from direct view.")
                }

                # Scenario 2: User switches to an interactive scan. Load additive data in the background.
            } else if (interaction_type != "none" && grepl("interactive", dataset_group_val, ignore.case = TRUE)) {
                # Only run loader if cache is empty for the current trait
                if (is.null(additive_scan_data_rv())) {
                    message("scanServer: In interactive mode with empty cache, triggering background additive load.")

                    # Derive the corresponding additive dataset name
                    base_name <- gsub(",\\s*interactive\\s*\\([^)]*\\)", "", dataset_group_val)
                    additive_dataset_name <- paste0(trimws(base_name), ", additive")

                    # Load and process the additive data
                    tryCatch(
                        {
                            file_dir_val <- import_reactives()$file_directory
                            req(file_dir_val)

                            result_list <- trait_scan(
                                file_dir = file_dir_val,
                                selected_dataset = additive_dataset_name,
                                selected_trait = trait_val,
                                cache_env = NULL
                            )

                            if (!is.null(result_list) && !is.null(result_list$scan_data) && nrow(result_list$scan_data) > 0) {
                                scan_data <- result_list$scan_data
                                processed_data <- QTL_plot_visualizer(scan_data, trait_val, 7.5, import_reactives()$markers)
                                additive_scan_data_rv(processed_data)
                                message("scanServer: Additive data cache populated by background loader.")
                            } else {
                                message("scanServer: Background loader found no additive data.")
                            }
                        },
                        error = function(e) {
                            message("scanServer: Error during background additive load: ", e$message)
                            additive_scan_data_rv(NULL) # Ensure cache is cleared on error
                        }
                    )
                }
            }
        })


        # NEW: Reactive to load DIET interactive data for overlay
        overlay_diet_scan_data <- shiny::reactive({
            # Only run if the diet overlay toggle is on
            if (!overlay_diet_toggle()) {
                return(NULL)
            }

            trait_val <- current_trait_for_scan()
            dataset_group_val <- selected_dataset_group()

            # Ensure this is an additive scan context (toggles only show for additive scans)
            if (is.null(trait_val) || is.null(dataset_group_val) ||
                grepl("interactive", dataset_group_val, ignore.case = TRUE)) {
                return(NULL)
            }

            # Derive the corresponding interactive dataset name
            interactive_dataset_name <- get_interactive_dataset_name(dataset_group_val, "diet")
            if (interactive_dataset_name == dataset_group_val) {
                return(NULL)
            } # No mapping found

            message("scanServer: Loading DIET overlay data for '", trait_val, "' from dataset '", interactive_dataset_name, "'")

            tryCatch(
                {
                    file_dir_val <- import_reactives()$file_directory
                    req(file_dir_val)

                    # Load the interactive scan data
                    result_list <- trait_scan(
                        file_dir = file_dir_val,
                        selected_dataset = interactive_dataset_name,
                        selected_trait = trait_val,
                        cache_env = NULL
                    )

                    if (is.null(result_list) || is.null(result_list$scan_data) || nrow(result_list$scan_data) == 0) {
                        return(NULL)
                    }
                    scan_data <- result_list$scan_data

                    # Process it with interactive-specific LOD threshold (10.5 instead of 7.5)
                    main_par_list <- main_par_inputs()
                    req(main_par_list$LOD_thr, import_reactives()$markers)

                    # Use higher LOD threshold for interactive scans (mimicking interactive mode)
                    interactive_lod_threshold <- 10.5

                    processed_data <- QTL_plot_visualizer(scan_data, trait_val, interactive_lod_threshold, import_reactives()$markers)

                    # Apply chromosome filtering to match main plot
                    selected_chr <- main_par_list$selected_chr()
                    if (selected_chr != "All" && !is.null(processed_data) && nrow(processed_data) > 0) {
                        sel_chr_num <- selected_chr
                        if (selected_chr == "X") sel_chr_num <- 20
                        if (selected_chr == "Y") sel_chr_num <- 21
                        if (selected_chr == "M") sel_chr_num <- 22
                        sel_chr_num <- as.numeric(sel_chr_num)
                        processed_data <- dplyr::filter(processed_data, chr == sel_chr_num)
                    }

                    return(processed_data)
                },
                error = function(e) {
                    message("scanServer: Error loading DIET overlay data: ", e$message)
                    return(NULL)
                }
            )
        }) %>% shiny::debounce(250)

        # NEW: Reactive to load SEX interactive data for overlay
        overlay_sex_scan_data <- shiny::reactive({
            # Only run if the sex overlay toggle is on
            if (!overlay_sex_toggle()) {
                return(NULL)
            }

            trait_val <- current_trait_for_scan()
            dataset_group_val <- selected_dataset_group()

            # Ensure this is an additive scan context (toggles only show for additive scans)
            if (is.null(trait_val) || is.null(dataset_group_val) ||
                grepl("interactive", dataset_group_val, ignore.case = TRUE)) {
                return(NULL)
            }

            # Derive the corresponding interactive dataset name
            interactive_dataset_name <- get_interactive_dataset_name(dataset_group_val, "sex")
            if (interactive_dataset_name == dataset_group_val) {
                return(NULL)
            } # No mapping found

            message("scanServer: Loading SEX overlay data for '", trait_val, "' from dataset '", interactive_dataset_name, "'")

            tryCatch(
                {
                    file_dir_val <- import_reactives()$file_directory
                    req(file_dir_val)

                    # Load the interactive scan data
                    result_list <- trait_scan(
                        file_dir = file_dir_val,
                        selected_dataset = interactive_dataset_name,
                        selected_trait = trait_val,
                        cache_env = NULL
                    )

                    if (is.null(result_list) || is.null(result_list$scan_data) || nrow(result_list$scan_data) == 0) {
                        return(NULL)
                    }
                    scan_data <- result_list$scan_data

                    # Process it with interactive-specific LOD threshold (10.5 instead of 7.5)
                    main_par_list <- main_par_inputs()
                    req(main_par_list$LOD_thr, import_reactives()$markers)

                    # Use higher LOD threshold for interactive scans (mimicking interactive mode)
                    interactive_lod_threshold <- 10.5

                    processed_data <- QTL_plot_visualizer(scan_data, trait_val, interactive_lod_threshold, import_reactives()$markers)

                    # Apply chromosome filtering to match main plot
                    selected_chr <- main_par_list$selected_chr()
                    if (selected_chr != "All" && !is.null(processed_data) && nrow(processed_data) > 0) {
                        sel_chr_num <- selected_chr
                        if (selected_chr == "X") sel_chr_num <- 20
                        if (selected_chr == "Y") sel_chr_num <- 21
                        if (selected_chr == "M") sel_chr_num <- 22
                        sel_chr_num <- as.numeric(sel_chr_num)
                        processed_data <- dplyr::filter(processed_data, chr == sel_chr_num)
                    }

                    return(processed_data)
                },
                error = function(e) {
                    message("scanServer: Error loading SEX overlay data: ", e$message)
                    return(NULL)
                }
            )
        }) %>% shiny::debounce(250)

        # NEW: Reactive to load SEX x DIET interactive data for overlay
        overlay_sex_diet_scan_data <- shiny::reactive({
            # Only run if the Sex x Diet overlay toggle is on
            if (!overlay_sex_diet_toggle()) {
                return(NULL)
            }

            trait_val <- current_trait_for_scan()
            dataset_group_val <- selected_dataset_group()

            # Ensure this is an additive scan context (toggles only show for additive scans)
            if (is.null(trait_val) || is.null(dataset_group_val) ||
                grepl("interactive", dataset_group_val, ignore.case = TRUE)) {
                return(NULL)
            }

            # Derive the corresponding interactive dataset name
            interactive_dataset_name <- get_interactive_dataset_name(dataset_group_val, "sex_diet")
            if (interactive_dataset_name == dataset_group_val) {
                return(NULL)
            }

            message("scanServer: Loading SEXxDIET overlay data for '", trait_val, "' from dataset '", interactive_dataset_name, "'")

            tryCatch(
                {
                    file_dir_val <- import_reactives()$file_directory
                    req(file_dir_val)

                    # Load the interactive scan data
                    result_list <- trait_scan(
                        file_dir = file_dir_val,
                        selected_dataset = interactive_dataset_name,
                        selected_trait = trait_val,
                        cache_env = NULL
                    )

                    if (is.null(result_list) || is.null(result_list$scan_data) || nrow(result_list$scan_data) == 0) {
                        return(NULL)
                    }
                    scan_data <- result_list$scan_data

                    # Process it with interactive-specific LOD threshold (Sex x Diet uses 15.7)
                    main_par_list <- main_par_inputs()
                    req(main_par_list$LOD_thr, import_reactives()$markers)
                    interactive_lod_threshold <- 15.7

                    processed_data <- QTL_plot_visualizer(scan_data, trait_val, interactive_lod_threshold, import_reactives()$markers)

                    # Apply chromosome filtering to match main plot
                    selected_chr <- main_par_list$selected_chr()
                    if (selected_chr != "All" && !is.null(processed_data) && nrow(processed_data) > 0) {
                        sel_chr_num <- selected_chr
                        if (selected_chr == "X") sel_chr_num <- 20
                        if (selected_chr == "Y") sel_chr_num <- 21
                        if (selected_chr == "M") sel_chr_num <- 22
                        sel_chr_num <- as.numeric(sel_chr_num)
                        processed_data <- dplyr::filter(processed_data, chr == sel_chr_num)
                    }

                    return(processed_data)
                },
                error = function(e) {
                    message("scanServer: Error loading SEXxDIET overlay data: ", e$message)
                    return(NULL)
                }
            )
        }) %>% shiny::debounce(250)

        # Old observers for additive data are now replaced by the single observer above

        current_scan_plot_gg <- shiny::reactive({
            # Get inputs with early validation
            scan_data <- scan_table_chr()
            main_par_list <- main_par_inputs()

            # --- ROBUSTNESS FIX ---
            # Add a strong validation check to ensure all required inputs are present
            # before attempting to render the plot. This prevents the "length zero" error.
            shiny::req(
                scan_data,
                is.data.frame(scan_data),
                nrow(scan_data) > 0,
                main_par_list,
                !is.null(main_par_list$selected_chr),
                !is.null(main_par_list$selected_chr())
            )

            # Early returns for missing data (retained as a safeguard)
            if (is.null(scan_data) || nrow(scan_data) == 0 ||
                is.null(main_par_list) || is.null(main_par_list$selected_chr)) {
                return(NULL)
            }

            selected_chr <- main_par_list$selected_chr()
            static_threshold <- static_lod_threshold_line() # Get the static threshold

            # Get overlay data if toggles are active
            diet_overlay <- if (isTRUE(overlay_diet_toggle())) overlay_diet_scan_data() else NULL
            sex_overlay <- if (isTRUE(overlay_sex_toggle())) overlay_sex_scan_data() else NULL
            sex_diet_overlay <- if (isTRUE(overlay_sex_diet_toggle())) overlay_sex_diet_scan_data() else NULL


            # Streamlined plot creation
            tryCatch(
                {
                    p <- ggplot_qtl_scan(
                        scan_data,
                        -Inf, # Pass -Inf to remove the interactive threshold bar
                        selected_chr,
                        overlay_diet_data = diet_overlay,
                        overlay_sex_data = sex_overlay,
                        overlay_sex_diet_data = sex_diet_overlay
                    )
                    if (!is.null(p) && !is.null(static_threshold)) {
                        p <- p + ggplot2::geom_hline(yintercept = static_threshold, linetype = "dashed", color = "grey20")
                    }
                    p
                },
                error = function(e) {
                    message("scanServer: Error creating plot: ", e$message)
                    return(NULL)
                }
            )
        }) %>% shiny::debounce(150) # Optimized debounce timing

        # NEW: Reactive for difference plot data to be used by plot and click handler
        diff_plot_data_reactive <- shiny::reactive({
            if (is.null(interaction_type_reactive)) {
                return(NULL)
            }
            interaction_type <- interaction_type_reactive()
            if (is.null(interaction_type) || interaction_type == "none") {
                return(NULL)
            }

            interactive_data <- scan_table_chr()
            additive_data <- additive_scan_data_rv()

            if (is.null(interactive_data) || is.null(additive_data) ||
                nrow(interactive_data) == 0 || nrow(additive_data) == 0) {
                return(NULL)
            }

            main_par_list <- main_par_inputs()
            selected_chromosome <- main_par_list$selected_chr()

            additive_data_chr <- if (selected_chromosome == "All") {
                additive_data
            } else {
                sel_chr_num <- selected_chromosome
                if (selected_chromosome == "X") sel_chr_num <- 20
                if (selected_chromosome == "Y") sel_chr_num <- 21
                if (selected_chromosome == "M") sel_chr_num <- 22
                sel_chr_num <- as.numeric(sel_chr_num)
                dplyr::filter(additive_data, chr == sel_chr_num)
            }

            if (nrow(interactive_data) != nrow(additive_data_chr)) {
                aligned_data <- dplyr::inner_join(
                    interactive_data,
                    additive_data_chr,
                    by = "markers",
                    suffix = c("_int", "_add")
                )
                if (nrow(aligned_data) == 0) {
                    return(NULL)
                }

                diff_plot_data <- aligned_data %>%
                    dplyr::mutate(LOD = LOD_int - LOD_add) %>%
                    dplyr::select(markers, chr = chr_int, position = position_int, BPcum = BPcum_int, LOD)
            } else {
                diff_plot_data <- interactive_data
                diff_plot_data$LOD <- interactive_data$LOD - additive_data_chr$LOD
            }

            diff_plot_data$LOD[is.na(diff_plot_data$LOD)] <- 0
            return(diff_plot_data)
        })

        # Reactive to create difference plot (interactive - additive)
        difference_plot_gg <- shiny::reactive({
            diff_plot_data <- diff_plot_data_reactive()
            if (is.null(diff_plot_data) || nrow(diff_plot_data) == 0) {
                return(NULL)
            }

            interaction_type <- interaction_type_reactive()
            static_diff_threshold <- if (interaction_type == "sex") 4.1 else if (interaction_type == "diet") 4.1 else if (interaction_type == "sex_diet") 9.5 else NULL

            # Display label for title
            interaction_label <- if (identical(interaction_type, "sex_diet")) "Sex x Diet" else stringr::str_to_title(interaction_type)

            tryCatch(
                {
                    main_par_list <- main_par_inputs()
                    selected_chromosome <- main_par_list$selected_chr()

                    # Create plot
                    diff_plot <- ggplot_qtl_scan(diff_plot_data, -Inf, selected_chromosome)

                    if (!is.null(diff_plot)) {
                        diff_plot <- diff_plot + ggplot2::labs(title = paste("LOD Difference:", interaction_label, "- Additive"))
                        if (!is.null(static_diff_threshold)) {
                            # Add horizontal lines for positive and negative thresholds
                            diff_plot <- diff_plot +
                                ggplot2::geom_hline(yintercept = static_diff_threshold, linetype = "dashed", color = "grey20") +
                                ggplot2::geom_hline(yintercept = -static_diff_threshold, linetype = "dashed", color = "grey20")
                        }
                    }

                    return(diff_plot)
                },
                error = function(e) {
                    message("scanServer: Error creating difference plot: ", e$message)
                    return(NULL)
                }
            )
        }) %>% shiny::debounce(300)

        # Simple UI state check
        show_stacked_plots <- shiny::reactive({
            if (is.null(interaction_type_reactive)) {
                return(FALSE)
            }
            interaction_type <- interaction_type_reactive()
            !is.null(interaction_type) && interaction_type != "none"
        })

        output$scan_plot_subheader <- renderText({
            req(numb_mice_rv())
            paste("Number of mice scanned:", numb_mice_rv(), "/ 1157")
        })

        output$scan_plot_ui_render <- shiny::renderUI({
            # Minimal dependencies for faster rendering
            plot_gg <- current_scan_plot_gg()
            shiny::req(plot_gg)

            use_stacked <- show_stacked_plots()

            if (isTRUE(use_stacked)) {
                # Calculate height for each plot (split the total height)
                individual_plot_height <- plot_height_rv() / 2

                shiny::tagList(
                    # Interactive LOD plot (top)
                    shiny::div(
                        style = "margin-bottom: 10px;",
                        shiny::div("Interactive LOD Scan", style = "text-align: center; font-weight: bold; margin-bottom: 5px;"),
                        shiny::div(textOutput(ns("scan_plot_subheader")), style = "text-align: center; font-style: italic; margin-bottom: 5px;"),
                        plotly::plotlyOutput(ns("render_plotly_plot"),
                            width = paste0(plot_width_rv(), "px"),
                            height = paste0(individual_plot_height, "px")
                        ) |>
                            shinycssloaders::withSpinner(type = 8, color = "#3498db")
                    ),
                    # Difference plot (bottom)
                    shiny::div(
                        shiny::div("LOD Difference (Interactive - Additive)", style = "text-align: center; font-weight: bold; margin-bottom: 5px;"),
                        plotly::plotlyOutput(ns("render_difference_plot"),
                            width = paste0(plot_width_rv(), "px"),
                            height = paste0(individual_plot_height, "px")
                        ) |>
                            shinycssloaders::withSpinner(type = 8, color = "#e74c3c")
                    )
                )
            } else {
                # Show only the main plot (additive or regular datasets)
                shiny::tagList(
                    shiny::div("Additive LOD Scan", style = "text-align: center; font-weight: bold; margin-bottom: 5px;"),
                    shiny::div(textOutput(ns("scan_plot_subheader")), style = "text-align: center; font-style: italic; margin-bottom: 5px;"),
                    plotly::plotlyOutput(ns("render_plotly_plot"),
                        width = paste0(plot_width_rv(), "px"),
                        height = paste0(plot_height_rv(), "px")
                    ) |>
                        shinycssloaders::withSpinner(type = 8, color = "#3498db")
                )
            }
        })

        output$render_plotly_plot <- plotly::renderPlotly({
            plot_gg <- current_scan_plot_gg()
            shiny::req(plot_gg)

            # Streamlined plotly creation
            plt <- plotly::ggplotly(plot_gg,
                source = ns("qtl_scan_plotly"),
                tooltip = c("text")
            ) %>%
                plotly::layout(
                    dragmode = "zoom",
                    hovermode = "closest",
                    title = list(text = NULL),
                    xaxis = list(
                        title = plot_gg$labels$x,
                        fixedrange = FALSE
                    ),
                    yaxis = list(
                        title = plot_gg$labels$y,
                        fixedrange = TRUE
                    )
                ) %>%
                plotly::config(
                    displaylogo = FALSE,
                    modeBarButtonsToRemove = c(
                        "select2d", "lasso2d", "hoverClosestCartesian",
                        "hoverCompareCartesian", "toggleSpikelines", "sendDataToCloud"
                    )
                ) %>%
                plotly::event_register("plotly_click")

            plt
        })

        # NEW: Click handler for the difference plot
        shiny::observeEvent(plotly::event_data("plotly_click", source = ns("difference_plotly")), {
            ev_data <- plotly::event_data("plotly_click", source = ns("difference_plotly"))
            req(ev_data)

            message("scanServer: Click detected on difference plot.")

            # Get data from the difference plot
            diff_data <- diff_plot_data_reactive()
            req(diff_data, nrow(diff_data) > 0)

            # Short-circuit for Sex x Diet: no dedicated peak files defined
            interaction_type <- interaction_type_reactive()
            if (!is.null(interaction_type) && interaction_type == "sex_diet") {
                shiny::showNotification("Peak file lookup is not available for Sex x Diet interactive differences.", type = "warning", duration = 4)
                return()
            }

            # Find nearest point
            main_par_list <- main_par_inputs()
            selected_chr <- main_par_list$selected_chr()
            xvar <- if (selected_chr == "All") "BPcum" else "position"

            # Find the nearest point in the data to the click event
            distances <- sqrt((diff_data[[xvar]] - ev_data$x)^2 + (diff_data$LOD - ev_data$y)^2)
            nearest_point <- diff_data[which.min(distances), ]

            message(paste("scanServer (Diff Click): Nearest point is marker", nearest_point$markers, "on chr", nearest_point$chr, "at pos", nearest_point$position))

            # Get file paths for interactive peaks
            dataset_group <- selected_dataset_group()

            paths <- get_interactive_peak_filepaths(dataset_group, interaction_type)
            if (is.null(paths)) {
                shiny::showNotification(paste("Peak files not found for this interaction type."), type = "error")
                return()
            }

            if (!file.exists(paths$file1) || !file.exists(paths$file2)) {
                shiny::showNotification(paste("Peak files not found for this interaction type."), type = "error")
                warning(paste("Missing interactive peak files:", paths$file1, "or", paths$file2))
                return()
            }

            # Load peak data
            peak_data_1 <- data.table::fread(paths$file1)
            peak_data_2 <- data.table::fread(paths$file2)

            # Find the corresponding peak within a 4Mb window
            trait <- current_trait_for_scan()

            find_peak_in_window <- function(peak_df, target_trait, target_chr, target_pos) {
                trait_col <- if ("gene_symbol" %in% colnames(peak_df)) "gene_symbol" else "phenotype"
                peak_dt <- data.table::as.data.table(peak_df)

                # Find peaks within the window
                peak_subset <- peak_dt[get(trait_col) == target_trait &
                    qtl_chr == target_chr &
                    abs(qtl_pos - target_pos) <= 4, ]

                if (nrow(peak_subset) > 0) {
                    # Return the one with the highest LOD score in the window
                    return(as.data.frame(peak_subset[which.max(peak_subset$qtl_lod), ]))
                }
                return(NULL)
            }

            # Search for peaks in both files
            found_peak_1 <- find_peak_in_window(peak_data_1, trait, chr_XYM(nearest_point$chr), nearest_point$position)
            found_peak_2 <- find_peak_in_window(peak_data_2, trait, chr_XYM(nearest_point$chr), nearest_point$position)

            # Add labels for plotting
            if (!is.null(found_peak_1)) found_peak_1$plot_label <- paths$labels[1]
            if (!is.null(found_peak_2)) found_peak_2$plot_label <- paths$labels[2]

            # Update reactive values
            diff_peak_1_rv(found_peak_1)
            diff_peak_2_rv(found_peak_2)

            # Clear the single peak selection to avoid confusion
            selected_peak_rv(NULL)
            clicked_plotly_point_details_lod_scan_rv(NULL)
        })

        # Render function for the difference plot
        output$render_difference_plot <- plotly::renderPlotly({
            # Only render for interactive datasets
            shiny::req(show_stacked_plots())

            diff_plot_gg <- difference_plot_gg()

            if (is.null(diff_plot_gg)) {
                # Show a placeholder when no difference plot is available
                placeholder_plot <- ggplot2::ggplot() +
                    ggplot2::theme_void() +
                    ggplot2::labs(title = "No difference plot available") +
                    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, color = "#7f8c8d"))

                return(plotly::ggplotly(placeholder_plot))
            }

            # Streamlined difference plot creation
            plotly::ggplotly(diff_plot_gg,
                source = ns("difference_plotly"),
                tooltip = c("text")
            ) %>%
                plotly::layout(
                    dragmode = "zoom",
                    hovermode = "closest",
                    title = list(text = NULL),
                    xaxis = list(
                        title = diff_plot_gg$labels$x,
                        fixedrange = FALSE
                    ),
                    yaxis = list(
                        title = diff_plot_gg$labels$y,
                        fixedrange = TRUE
                    )
                ) %>%
                plotly::config(
                    displaylogo = FALSE,
                    modeBarButtonsToRemove = c(
                        "select2d", "lasso2d", "hoverClosestCartesian",
                        "hoverCompareCartesian", "toggleSpikelines", "sendDataToCloud"
                    )
                ) %>%
                plotly::event_register("plotly_click")
        })

        # Simple click handler using distance-based approach (like the old app)
        shiny::observeEvent(plotly::event_data("plotly_click", source = ns("qtl_scan_plotly")), {
            ev_data <- plotly::event_data("plotly_click", source = ns("qtl_scan_plotly"))

            if (!is.null(ev_data)) {
                current_scan_data <- scan_table_chr()

                if (is.null(current_scan_data) || nrow(current_scan_data) == 0) {
                    message("scanServer: No scan data available for click detection")
                    return()
                }

                # Get click coordinates
                x_clicked <- ev_data$x
                y_clicked <- ev_data$y

                message(paste("scanServer: Click detected at coordinates x:", x_clicked, "y:", y_clicked))

                # Determine which x variable to use based on chromosome selection
                main_par_list <- main_par_inputs()
                selected_chr <- main_par_list$selected_chr()
                xvar <- if (selected_chr == "All") "BPcum" else "position"

                # Calculate distances to find nearest point (like the old app)
                if (xvar %in% colnames(current_scan_data) && "LOD" %in% colnames(current_scan_data)) {
                    distances <- sqrt((current_scan_data[[xvar]] - x_clicked)^2 + (current_scan_data$LOD - y_clicked)^2)
                    nearest_idx <- which.min(distances)
                    nearest_point <- current_scan_data[nearest_idx, ]

                    message(paste("scanServer: *** CLICK SUCCESS *** Nearest point marker:", nearest_point$markers, "LOD:", nearest_point$LOD))

                    # Debug: Print available columns and values
                    message("scanServer: Available columns in nearest_point: ", paste(colnames(nearest_point), collapse = ", "))
                    message("scanServer: nearest_point chr value:", nearest_point$chr, "position value:", nearest_point$position)

                    # Find corresponding peak and update selected_peak_rv
                    current_trait <- current_trait_for_scan()
                    current_dataset <- selected_dataset_group()

                    if (!is.null(current_trait) && !is.null(current_dataset)) {
                        trait_type_val <- get_trait_type(import_reactives(), current_dataset)
                        all_peaks <- peak_finder(
                            file_dir = import_reactives()$file_directory,
                            selected_dataset = current_dataset,
                            selected_trait = current_trait,
                            trait_type = trait_type_val,
                            cache_env = local_peaks_cache,
                            use_cache = TRUE
                        )

                        if (!is.null(all_peaks) && nrow(all_peaks) > 0) {
                            # Find closest peak to clicked location
                            clicked_chr <- nearest_point$chr
                            clicked_pos <- nearest_point$position

                            closest_peak <- all_peaks %>%
                                dplyr::filter(qtl_chr == chr_XYM(clicked_chr)) %>%
                                dplyr::mutate(pos_diff = abs(qtl_pos - clicked_pos)) %>%
                                dplyr::filter(pos_diff == min(pos_diff))

                            if (nrow(closest_peak) > 0) {
                                # Update the selected peak reactive
                                selected_peak_rv(closest_peak[1, ])
                                message(paste(
                                    "scanServer: *** PEAK UPDATED *** Selected peak:", closest_peak[1, ]$marker,
                                    "LOD:", closest_peak[1, ]$qtl_lod
                                ))
                            }
                        }
                    }

                    # Create click details for display with raw column names (as expected by app.R)
                    click_details <- data.frame(
                        markers = if ("markers" %in% colnames(nearest_point)) nearest_point$markers else "Unknown",
                        chr = if ("chr" %in% colnames(nearest_point)) chr_XYM(nearest_point$chr) else "Unknown",
                        position = if ("position" %in% colnames(nearest_point)) round(nearest_point$position, 3) else "Unknown",
                        LOD = if ("LOD" %in% colnames(nearest_point)) round(nearest_point$LOD, 3) else "Unknown"
                    )

                    clicked_plotly_point_details_lod_scan_rv(click_details)
                } else {
                    message("scanServer: Required columns not found in scan data for click detection")
                    clicked_plotly_point_details_lod_scan_rv(data.frame(Info = "Click data not available"))
                }
            } else {
                message("scanServer: No click event data received")
            }
        })

        output$plot_click_dt <- DT::renderDT({
            details <- clicked_plotly_point_details_lod_scan_rv()
            if (is.null(details) || nrow(details) == 0) {
                return(DT::datatable(data.frame(Info = "Click on the plot to see point details."),
                    options = list(dom = "t"), rownames = FALSE
                ))
            }

            # Transpose the data for better display
            if (ncol(details) > 1) {
                # Create a two-column table with Property and Value
                transposed_details <- data.frame(
                    Property = names(details),
                    Value = as.character(unlist(details[1, ])),
                    stringsAsFactors = FALSE
                )

                # Clean up column names for better display
                transposed_details$Property <- gsub("_", " ", transposed_details$Property)
                transposed_details$Property <- gsub("Founder ", "Founder ", transposed_details$Property)
                transposed_details$Property <- ifelse(transposed_details$Property == "markers", "Marker",
                    ifelse(transposed_details$Property == "chr", "Chromosome",
                        ifelse(transposed_details$Property == "position", "Position (Mb)",
                            ifelse(transposed_details$Property == "CisOrTrans", "Cis/Trans",
                                ifelse(transposed_details$Property == "CI Range", "Confidence Interval",
                                    transposed_details$Property
                                )
                            )
                        )
                    )
                )

                return(DT::datatable(
                    transposed_details,
                    options = list(
                        dom = "t",
                        paging = FALSE,
                        searching = FALSE,
                        columnDefs = list(
                            list(targets = 0, className = "dt-left", width = "40%"),
                            list(targets = 1, className = "dt-left", width = "60%")
                        )
                    ),
                    rownames = FALSE,
                    selection = "none",
                    colnames = c("Property", "Value"),
                    class = "compact hover"
                ) %>%
                    DT::formatStyle(
                        "Property",
                        fontWeight = "bold",
                        backgroundColor = "#f8f9fa"
                    ))
            } else {
                return(DT::datatable(details, options = list(dom = "t", paging = FALSE), rownames = FALSE, selection = "none"))
            }
        })

        output$download_qtl_plot_png <- shiny::downloadHandler(
            filename = function() {
                main_par_list <- main_par_inputs()
                trait_name <- current_trait_for_scan() %||% "plot" # Use the current scanned trait
                chr_suffix <- if (main_par_list$selected_chr() != "All") paste0("_chr", main_par_list$selected_chr()) else ""
                paste0("lod_plot_", trait_name, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".png")
            },
            content = function(file) {
                shiny::req(current_scan_plot_gg())
                ggplot2::ggsave(file,
                    plot = current_scan_plot_gg(),
                    width = plot_width_rv() / 96,
                    height = plot_height_rv() / 96,
                    dpi = 300, units = "in"
                )
            }
        )

        output$download_qtl_plot_pdf <- shiny::downloadHandler(
            filename = function() {
                main_par_list <- main_par_inputs()
                trait_name <- current_trait_for_scan() %||% "plot"
                chr_suffix <- if (main_par_list$selected_chr() != "All") paste0("_chr", main_par_list$selected_chr()) else ""
                paste0("lod_plot_", trait_name, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".pdf")
            },
            content = function(file) {
                shiny::req(current_scan_plot_gg())
                ggplot2::ggsave(file,
                    plot = current_scan_plot_gg(),
                    width = plot_width_rv() / 96,
                    height = plot_height_rv() / 96,
                    device = cairo_pdf, units = "in"
                )
            }
        )

        file_name_reactive <- shiny::reactive({
            main_par_list <- main_par_inputs()
            trait_name <- current_trait_for_scan() %||% "scan"
            instanceID <- trait_name
            selected_chromosome <- main_par_list$selected_chr()
            if (!is.null(selected_chromosome) && selected_chromosome != "All") {
                instanceID <- paste0(instanceID, "_chr", selected_chromosome)
            }
            paste("scan", instanceID, sep = "_")
        })

        # Return reactive values that might be useful for other modules (e.g., download peak info)
        # This is currently not used by the main app, but good practice for a module.
        return(shiny::reactiveValues(
            filename = file_name_reactive,
            tables = shiny::reactiveValues(scan = scan_table_chr),
            plots = shiny::reactiveValues(scan = current_scan_plot_gg),
            clicked_point_details = clicked_plotly_point_details_lod_scan_rv,
            selected_peak = selected_peak_rv,
            diff_peak_1 = diff_peak_1_rv,
            diff_peak_2 = diff_peak_2_rv
        ))
    })
}

#' @rdname scanServer
#' @export
scanUI <- function(id) {
    ns <- shiny::NS(id)
    DT::DTOutput(ns("plot_click_dt"))
}

#' @rdname scanServer
#' @export
scanOutput <- function(id) {
    ns <- shiny::NS(id)
    shiny::uiOutput(ns("scan_plot_ui_render"))
}
