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
                            "Genetic Litter" = "GenLit"
                        ),
                        selected = "Sex",
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

            req(plot_data, nrow(plot_data) > 0, !is.null(grouping_vars), length(grouping_vars) > 0)

            # --- ROBUSTNESS FIX: Clean and validate grouping variable ---
            plot_data_clean <- data.table::copy(plot_data)

            # 1. Sanitize all character columns to prevent encoding issues
            char_cols <- names(which(sapply(plot_data_clean, is.character)))
            if (length(char_cols) > 0) {
                plot_data_clean[, (char_cols) := lapply(.SD, function(x) {
                    iconv(x, to = "ASCII//TRANSLIT", sub = "")
                }), .SDcols = char_cols]
            }

            # 2. Filter data based on expected groups for each selected variable
            if ("Sex" %in% grouping_vars) {
                plot_data_clean <- plot_data_clean[Sex %in% c("F", "M")]
            }
            if ("Diet" %in% grouping_vars) {
                plot_data_clean <- plot_data_clean[Diet %in% c("HC", "HF")]
            }

            # For all cases, remove rows where any of the grouping variables are NA or empty
            for (var in grouping_vars) {
                if (var %in% names(plot_data_clean)) {
                    plot_data_clean <- plot_data_clean[!is.na(get(var)) & get(var) != ""]
                }
            }

            req(nrow(plot_data_clean) > 0)

            # Create a combined grouping variable if more than one is selected
            if (length(grouping_vars) > 1) {
                plot_data_clean[, Combined_Group := do.call(paste, c(.SD, sep = " x ")), .SDcols = grouping_vars]
                grouping_col_name <- "Combined_Group"
                xaxis_title <- paste(grouping_vars, collapse = " x ")
            } else {
                grouping_col_name <- grouping_vars
                xaxis_title <- grouping_vars
            }


            # Create safe trait name for title
            safe_trait_name <- iconv(trait_to_profile(), to = "ASCII//TRANSLIT", sub = "")

            # --- Direct Plotly Implementation for Boxplots ---
            p <- plotly::plot_ly(
                data = plot_data_clean,
                x = ~ get(grouping_col_name),
                y = ~Value,
                color = ~ get(grouping_col_name),
                colors = "Set2",
                type = "box",
                # Add jittered points directly, this is a robust method
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
                    xaxis = list(title = paste("<b>", xaxis_title, "</b>")),
                    yaxis = list(title = "<b>Phenotype Value</b>"),
                    showlegend = FALSE
                )

            return(p)
        })
    })
}
