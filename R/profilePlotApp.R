#' Profile Plot Module
#'
#' This module displays boxplots of phenotype values, grouped by selected
#' covariates. It dynamically loads data for a specific dataset category
#' for high performance and low memory usage.
#'
#' @param id shiny identifier
#' @param selected_dataset_category A reactive string representing the currently
#'   selected dataset category (e.g., "Liver Genes").
#' @param trait_to_profile A reactive string with the name of the trait to display.
#' @importFrom plotly plotlyOutput renderPlotly ggplotly
#' @importFrom ggplot2 ggplot aes_string geom_boxplot geom_jitter labs theme_minimal theme
#' @importFrom shinycssloaders withSpinner
#' @importFrom stringr str_replace_all
#' @importFrom data.table .
#' @export
profilePlotUI <- function(id) {
    ns <- shiny::NS(id)
    # The UI is now fully rendered in the server to allow controls to be placed above the plot
    uiOutput(ns("plot_ui_wrapper"))
}

#' @rdname profilePlotUI
#' @export
profilePlotServer <- function(id, selected_dataset_category, trait_to_profile) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Reactive to store the path to the current category's data files
        current_data_paths <- reactive({
            req(selected_dataset_category())
            category <- selected_dataset_category()

            sanitized_category <- tolower(category)
            sanitized_category <- str_replace_all(sanitized_category, "[^a-z0-9]+", "_")

            list(
                fst = file.path("/data/dev/miniViewer_3.0", paste0("pheno_data_long_", sanitized_category, ".fst")),
                rds = file.path("/data/dev/miniViewer_3.0", paste0("trait_names_", sanitized_category, ".rds"))
            )
        })

        # This reactive now directly uses the trait passed from the main app
        selected_trait_data <- reactive({
            req(trait_to_profile(), nzchar(trait_to_profile()))
            paths <- current_data_paths()

            req(file.exists(paths$fst))

            # Read the FST file and filter for the specific trait
            # Note: fst::fst() doesn't preserve data.table keys, so we need to read and filter manually

            tryCatch(
                {
                    # Read the full FST file (this might be memory intensive for large files)
                    full_data <- fst::read_fst(paths$fst, as.data.table = TRUE)

                    # Filter for the specific trait
                    trait_data <- full_data[Trait_Name == trait_to_profile()]

                    req(trait_data, nrow(trait_data) > 0)

                    return(as.data.table(trait_data)) # Ensure it's a mutable data.table
                },
                error = function(e) {
                    warning("Profile Plot Data ERROR: ", e$message)
                    return(NULL)
                }
            )
        })

        # Dynamic UI for the plot area
        output$plot_ui_wrapper <- renderUI({
            # Show a message if no trait is selected yet
            if (is.null(trait_to_profile()) || !nzchar(trait_to_profile())) {
                return(div(
                    style = "display: flex; justify-content: center; align-items: center; height: 400px; color: #7f8c8d; font-size: 1.2em;",
                    "Please select a trait in the 'Trait Search & LOD Scan' tab to see its profile plot."
                ))
            }

            # Show a different message if the data file for this category doesn't exist
            if (!file.exists(current_data_paths()$fst)) {
                return(div(
                    style = "display: flex; justify-content: center; align-items: center; height: 400px; color: #7f8c8d; font-size: 1.2em;",
                    paste("No phenotype data available for the '", selected_dataset_category(), "' category.")
                ))
            }

            # If everything is okay, show the plot output with controls above it
            tagList(
                div(
                    style = "margin-bottom: 15px; max-width: 500px;",
                    shiny::selectInput(
                        ns("grouping_selector"),
                        "Group By:",
                        choices = c(
                            "Sex" = "Sex",
                            "Diet" = "Diet",
                            "Generation and Litter" = "Generation_Litter"
                        ),
                        selected = "Sex",
                        multiple = TRUE
                    )
                ),
                div(
                    style = "margin-bottom: 15px; max-width: 500px;",
                    shiny::selectInput(
                        ns("coloring_selector"),
                        "Color By:",
                        choices = c(
                            "None" = "None",
                            "Sex" = "Sex",
                            "Diet" = "Diet",
                            "Generation and Litter" = "Generation_Litter"
                        ),
                        selected = "None",
                        multiple = TRUE
                    )
                ),
                plotly::plotlyOutput(ns("profile_boxplot")) |>
                    shinycssloaders::withSpinner(type = 8, color = "#3498db")
            )
        })

        output$profile_boxplot <- plotly::renderPlotly({
            plot_data <- selected_trait_data()
            grouping_vars <- input$grouping_selector
            coloring_vars <- input$coloring_selector

            req(plot_data, nrow(plot_data) > 0, !is.null(grouping_vars), length(grouping_vars) > 0)

            # --- ROBUSTNESS FIX: Clean and validate grouping variable ---
            plot_data_clean <- data.table::copy(plot_data)

            # Derive Generation and Litter from GenLit (e.g., "G41_L2" -> Generation="G41", Litter="L2")
            if ("GenLit" %in% names(plot_data_clean)) {
                parts <- data.table::tstrsplit(as.character(plot_data_clean$GenLit), "_", fixed = TRUE)
                if (length(parts) >= 2) {
                    plot_data_clean[, Generation := parts[[1]]]
                    plot_data_clean[, Litter := parts[[2]]]
                } else {
                    # Fallback if splitting fails
                    plot_data_clean[, Generation := plot_data_clean$GenLit]
                    plot_data_clean[, Litter := NA_character_]
                }
            }

            # Map combined option to its component variables (x grouping)
            if ("Generation_Litter" %in% grouping_vars) {
                grouping_vars <- unique(c(setdiff(grouping_vars, "Generation_Litter"), "Generation", "Litter"))
            }

            # Map combined option to its component variables (color grouping)
            if (!is.null(coloring_vars) && length(coloring_vars) > 0) {
                # Remove "None" if present
                coloring_vars <- setdiff(coloring_vars, "None")
                if ("Generation_Litter" %in% coloring_vars) {
                    coloring_vars <- unique(c(setdiff(coloring_vars, "Generation_Litter"), "Generation", "Litter"))
                }
            } else {
                coloring_vars <- character(0)
            }

            # 1. Sanitize all character columns to prevent encoding issues
            char_cols <- names(which(sapply(plot_data_clean, is.character)))
            if (length(char_cols) > 0) {
                plot_data_clean[, (char_cols) := lapply(.SD, function(x) {
                    iconv(x, to = "ASCII//TRANSLIT", sub = "")
                }), .SDcols = char_cols]
            }

            # 2. Filter data based on expected groups for each selected variable
            if ("Sex" %in% grouping_vars || "Sex" %in% coloring_vars) {
                plot_data_clean <- plot_data_clean[Sex %in% c("F", "M")]
            }
            if ("Diet" %in% grouping_vars || "Diet" %in% coloring_vars) {
                plot_data_clean <- plot_data_clean[Diet %in% c("HC", "HF")]
            }
            # For Generation/Litter we keep all observed values; no explicit filtering

            # For all cases, remove rows where any of the grouping/coloring variables are NA or empty
            required_vars <- unique(c(grouping_vars, coloring_vars))
            for (var in required_vars) {
                if (var %in% names(plot_data_clean)) {
                    plot_data_clean <- plot_data_clean[!is.na(get(var)) & get(var) != ""]
                }
            }

            req(nrow(plot_data_clean) > 0)

            # Create a combined grouping variable if more than one is selected (x axis)
            if (length(grouping_vars) > 1) {
                plot_data_clean[, Combined_Group := do.call(paste, c(.SD, sep = " x ")), .SDcols = grouping_vars]
                grouping_col_name <- "Combined_Group"
                xaxis_title <- paste(grouping_vars, collapse = " x ")
            } else {
                grouping_col_name <- grouping_vars
                xaxis_title <- grouping_vars
            }

            # Ensure the grouping column is properly cleaned and formatted as a factor.
            # This robustly handles columns like 'Generation'/'Litter' by first converting to character,
            # then sanitizing with iconv, and finally converting to a factor for plotting.
            cleaned_col <- iconv(as.character(plot_data_clean[[grouping_col_name]]), to = "ASCII//TRANSLIT", sub = "")
            plot_data_clean[, (grouping_col_name) := as.factor(cleaned_col)]

            # HOTFIX: Filter out any remaining rows where the grouping column contains non-printable
            # characters, which indicates data corruption. This prevents the plot from failing.
            plot_data_clean <- plot_data_clean[!grepl("[^ -~]", get(grouping_col_name))]

            # Build color grouping column if requested
            color_col_name <- NULL
            if (length(coloring_vars) > 0) {
                # Only keep color vars that actually exist in the data
                valid_color_vars <- coloring_vars[coloring_vars %in% names(plot_data_clean)]
                if (length(valid_color_vars) > 1) {
                    plot_data_clean[, Combined_Color := do.call(paste, c(.SD, sep = " x ")), .SDcols = valid_color_vars]
                    color_col_name <- "Combined_Color"
                } else if (length(valid_color_vars) == 1) {
                    color_col_name <- valid_color_vars
                }
                if (!is.null(color_col_name)) {
                    cleaned_color <- iconv(as.character(plot_data_clean[[color_col_name]]), to = "ASCII//TRANSLIT", sub = "")
                    plot_data_clean[, (color_col_name) := as.factor(cleaned_color)]
                    # Filter out non-printable in color column as well
                    plot_data_clean <- plot_data_clean[!grepl("[^ -~]", get(color_col_name))]
                }
            }

            # Create safe trait name for title
            safe_trait_name <- iconv(trait_to_profile(), to = "ASCII//TRANSLIT", sub = "")

            # --- Direct Plotly Implementation for Boxplots ---
            if (!is.null(color_col_name) && nzchar(color_col_name)) {
                p <- plotly::plot_ly(
                    data = plot_data_clean,
                    x = ~ get(grouping_col_name),
                    y = ~Value,
                    color = ~ get(color_col_name),
                    colors = "Set2",
                    type = "box",
                    boxpoints = "all",
                    jitter = 0.4,
                    pointpos = 0,
                    marker = list(size = 5, opacity = 0.5),
                    line = list(width = 1.5)
                ) %>%
                    plotly::layout(
                        title = list(
                            text = paste("<b>Phenotype Distribution:", safe_trait_name, "</b>"),
                            x = 0.5
                        ),
                        xaxis = list(
                            title = paste("<b>", xaxis_title, "</b>"),
                            type = "category"
                        ),
                        yaxis = list(title = "<b>Phenotype Value</b>"),
                        boxmode = "group",
                        showlegend = TRUE
                    )
            } else {
                p <- plotly::plot_ly(
                    data = plot_data_clean,
                    x = ~ get(grouping_col_name),
                    y = ~Value,
                    type = "box",
                    boxpoints = "all",
                    jitter = 0.4,
                    pointpos = 0, # Position points to the left of boxes
                    marker = list(size = 5, opacity = 0.5),
                    line = list(width = 1.5)
                ) %>%
                    plotly::layout(
                        title = list(
                            text = paste("<b>Phenotype Distribution:", safe_trait_name, "</b>"),
                            x = 0.5
                        ),
                        xaxis = list(
                            title = paste("<b>", xaxis_title, "</b>"),
                            type = "category" # Force the x-axis to be treated as categorical
                        ),
                        yaxis = list(title = "<b>Phenotype Value</b>"),
                        showlegend = FALSE
                    )
            }

            return(p)
        })
    })
}
