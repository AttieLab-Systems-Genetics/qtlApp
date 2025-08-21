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
#' @importFrom stats t.test wilcox.test p.adjust
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
                            "Sex" = "Sex",
                            "Diet" = "Diet",
                            "Generation and Litter" = "Generation_Litter"
                        ),
                        selected = "Sex",
                        multiple = FALSE
                    )
                ),
                div(
                    style = "display: flex; gap: 12px; flex-wrap: wrap; margin-bottom: 10px; max-width: 700px;",
                    shiny::selectInput(
                        ns("test_method"),
                        "Group-vs-Rest Test:",
                        choices = c(
                            "Wilcoxon (rank-sum)" = "wilcox",
                            "t-test (Welch)" = "ttest"
                        ),
                        selected = "wilcox"
                    ),
                    shiny::selectInput(
                        ns("p_adjust_method"),
                        "P-adjust:",
                        choices = c(
                            "BH (FDR)" = "BH",
                            "Bonferroni" = "bonferroni",
                            "None" = "none"
                        ),
                        selected = "BH"
                    ),
                    shiny::checkboxInput(ns("show_p_annotations"), "Show p on plot", value = TRUE)
                ),
                plotly::plotlyOutput(ns("profile_boxplot")) |>
                    shinycssloaders::withSpinner(type = 8, color = "#3498db")
            )
        })

        output$profile_boxplot <- plotly::renderPlotly({
            plot_data <- selected_trait_data()
            grouping_vars <- input$grouping_selector
            coloring_var <- input$coloring_selector

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

            # Determine single color variable
            color_col_name <- NULL
            if (!is.null(coloring_var) && nzchar(coloring_var)) {
                if (identical(coloring_var, "Generation_Litter")) {
                    # For color-by Generation and Litter, use original GenLit combined code
                    color_col_name <- if ("GenLit" %in% names(plot_data_clean)) "GenLit" else "Generation"
                } else {
                    color_col_name <- coloring_var
                }
            }

            # 1. Sanitize all character columns to prevent encoding issues
            char_cols <- names(which(sapply(plot_data_clean, is.character)))
            if (length(char_cols) > 0) {
                plot_data_clean[, (char_cols) := lapply(.SD, function(x) {
                    iconv(x, to = "ASCII//TRANSLIT", sub = "")
                }), .SDcols = char_cols]
            }

            # 2. Filter data based on expected groups for each selected variable
            if ("Sex" %in% grouping_vars || "Sex" %in% coloring_var) {
                plot_data_clean <- plot_data_clean[Sex %in% c("F", "M")]
            }
            if ("Diet" %in% grouping_vars || "Diet" %in% coloring_var) {
                plot_data_clean <- plot_data_clean[Diet %in% c("HC", "HF")]
            }
            # For Generation/Litter we keep all observed values; no explicit filtering

            # For all cases, remove rows where any of the grouping/coloring variables are NA or empty
            required_vars <- unique(c(grouping_vars, if (!is.null(color_col_name)) color_col_name else character(0)))
            for (var in required_vars) {
                if (var %in% names(plot_data_clean)) {
                    plot_data_clean <- plot_data_clean[!is.na(get(var)) & get(var) != ""]
                }
            }

            req(nrow(plot_data_clean) > 0)

            # Create a combined grouping variable if more than one is selected (x axis)
            if (length(grouping_vars) > 1) {
                # Special handling whenever both Generation and Litter are present (alone or with others)
                if (all(c("Generation", "Litter") %in% grouping_vars)) {
                    plot_data_clean[, GL_combined := paste0(Generation, "_", Litter)]
                    remaining <- setdiff(grouping_vars, c("Generation", "Litter"))
                    cols <- c("GL_combined", remaining)
                    if (length(cols) > 1) {
                        plot_data_clean[, Combined_Group := do.call(paste, c(.SD, sep = " x ")), .SDcols = cols]
                    } else {
                        plot_data_clean[, Combined_Group := GL_combined]
                    }
                    title_parts <- c("Generation and Litter", remaining)
                    xaxis_title <- paste(title_parts, collapse = " x ")
                } else {
                    plot_data_clean[, Combined_Group := do.call(paste, c(.SD, sep = " x ")), .SDcols = grouping_vars]
                    xaxis_title <- paste(grouping_vars, collapse = " x ")
                }
                grouping_col_name <- "Combined_Group"
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

            # Clean and validate color column if present
            if (!is.null(color_col_name) && (color_col_name %in% names(plot_data_clean))) {
                cleaned_color <- iconv(as.character(plot_data_clean[[color_col_name]]), to = "ASCII//TRANSLIT", sub = "")
                plot_data_clean[, (color_col_name) := as.factor(cleaned_color)]
                # Ensure deterministic level order for consistent color mapping
                if (identical(color_col_name, "Sex")) {
                    if (all(c("F", "M") %in% levels(plot_data_clean$Sex))) {
                        plot_data_clean[, Sex := factor(Sex, levels = c("F", "M"))]
                    }
                } else if (identical(color_col_name, "Diet")) {
                    if (all(c("HC", "HF") %in% levels(plot_data_clean$Diet))) {
                        plot_data_clean[, Diet := factor(Diet, levels = c("HC", "HF"))]
                    }
                }
                plot_data_clean <- plot_data_clean[!grepl("[^ -~]", get(color_col_name))]
            } else {
                color_col_name <- NULL
            }

            # Create safe trait name for title
            safe_trait_name <- iconv(trait_to_profile(), to = "ASCII//TRANSLIT", sub = "")

            # Prepare numeric positions for categorical x with jittered markers
            plot_data_clean[, x_index := as.numeric(get(grouping_col_name))]
            plot_data_clean[, x_jitter := jitter(x_index, amount = 0.2)]
            tick_vals <- sort(unique(plot_data_clean$x_index))
            tick_text <- levels(plot_data_clean[[grouping_col_name]])
            # Dynamic x-axis label formatting to prevent clutter
            has_gen_lit <- all(c("Generation", "Litter") %in% grouping_vars)
            has_sex <- "Sex" %in% grouping_vars
            has_diet <- "Diet" %in% grouping_vars
            is_heavy_combo <- has_gen_lit && has_sex && has_diet
            tick_count <- length(tick_text)
            tick_angle <- if (is_heavy_combo || tick_count > 10) -45 else 0
            tick_font_size <- if (is_heavy_combo || tick_count > 14) 8 else if (tick_count > 8) 10 else 12

            # Compute per-group vs all-others p-values (on-the-fly)
            test_method <- if (!is.null(input$test_method) && nzchar(input$test_method)) input$test_method else "wilcox"
            padj_method <- if (!is.null(input$p_adjust_method) && nzchar(input$p_adjust_method)) input$p_adjust_method else "BH"
            group_levels <- levels(plot_data_clean[[grouping_col_name]])
            pvals <- rep(NA_real_, length(group_levels))
            for (i in seq_along(group_levels)) {
                lev <- group_levels[i]
                in_group <- plot_data_clean[[grouping_col_name]] == lev
                n_in <- sum(in_group, na.rm = TRUE)
                n_out <- sum(!in_group, na.rm = TRUE)
                if (is.finite(n_in) && is.finite(n_out) && n_in >= 3 && n_out >= 3) {
                    pvals[i] <- tryCatch(
                        {
                            if (identical(test_method, "wilcox")) {
                                stats::wilcox.test(plot_data_clean$Value[in_group], plot_data_clean$Value[!in_group], exact = FALSE)$p.value
                            } else {
                                stats::t.test(plot_data_clean$Value[in_group], plot_data_clean$Value[!in_group], var.equal = FALSE)$p.value
                            }
                        },
                        error = function(e) NA_real_
                    )
                } else {
                    pvals[i] <- NA_real_
                }
            }
            pvals_adj <- if (!identical(padj_method, "none")) stats::p.adjust(pvals, method = padj_method) else pvals

            # Compute y positions for annotations above each group's box
            y_by_group <- plot_data_clean[, .(y_max = max(Value, na.rm = TRUE)), by = c(grouping_col_name)]
            # Reorder to match level order
            y_by_group <- y_by_group[match(group_levels, y_by_group[[grouping_col_name]]), ]
            y_range <- max(plot_data_clean$Value, na.rm = TRUE) - min(plot_data_clean$Value, na.rm = TRUE)
            y_margin <- if (is.finite(y_range) && y_range > 0) 0.05 * y_range else 0.1
            ann_y <- y_by_group$y_max + y_margin
            ann_text <- vapply(seq_along(group_levels), function(i) {
                pv <- pvals_adj[i]
                if (is.na(pv)) "p: n/a" else paste0("p=", formatC(pv, format = "g", digits = 3))
            }, character(1))
            annotations <- lapply(seq_along(group_levels), function(i) {
                list(
                    x = tick_vals[i],
                    y = ann_y[i],
                    text = ann_text[i],
                    xref = "x",
                    yref = "y",
                    showarrow = FALSE,
                    align = "center",
                    font = list(size = 10, color = "#2c3e50")
                )
            })

            # Build uncolored boxplots (one per x group)
            p <- plotly::plot_ly()
            p <- plotly::add_trace(
                p,
                data = plot_data_clean,
                x = ~x_index,
                y = ~Value,
                type = "box",
                boxpoints = FALSE,
                line = list(color = "#000000", width = 1.5),
                fillcolor = "rgba(0,0,0,0)",
                marker = list(size = 5),
                showlegend = FALSE
            )

            # Choose palette for markers when coloring is enabled
            colors_map <- NULL
            if (!is.null(color_col_name)) {
                if (identical(color_col_name, "Sex")) {
                    # Force two explicit traces to guarantee colors
                    if ("F" %in% plot_data_clean$Sex) {
                        p <- plotly::add_markers(
                            p,
                            data = plot_data_clean[Sex == "F"],
                            x = ~x_jitter,
                            y = ~Value,
                            marker = list(size = 6, opacity = 0.8, color = "#e41a1c"),
                            name = "Female",
                            legendgroup = "Sex",
                            showlegend = TRUE
                        )
                    }
                    if ("M" %in% plot_data_clean$Sex) {
                        p <- plotly::add_markers(
                            p,
                            data = plot_data_clean[Sex == "M"],
                            x = ~x_jitter,
                            y = ~Value,
                            marker = list(size = 6, opacity = 0.8, color = "#377eb8"),
                            name = "Male",
                            legendgroup = "Sex",
                            showlegend = TRUE
                        )
                    }
                } else if (identical(color_col_name, "Diet")) {
                    # Force two explicit traces to guarantee colors
                    if ("HC" %in% plot_data_clean$Diet) {
                        p <- plotly::add_markers(
                            p,
                            data = plot_data_clean[Diet == "HC"],
                            x = ~x_jitter,
                            y = ~Value,
                            marker = list(size = 6, opacity = 0.8, color = "#ff7f00"),
                            name = "HC",
                            legendgroup = "Diet",
                            showlegend = TRUE
                        )
                    }
                    if ("HF" %in% plot_data_clean$Diet) {
                        p <- plotly::add_markers(
                            p,
                            data = plot_data_clean[Diet == "HF"],
                            x = ~x_jitter,
                            y = ~Value,
                            marker = list(size = 6, opacity = 0.8, color = "#ffd92f"),
                            name = "HF",
                            legendgroup = "Diet",
                            showlegend = TRUE
                        )
                    }
                } else {
                    # Fallback to a vibrant palette distinct from Sex/Diet
                    colors_map <- "Dark2"
                    p <- plotly::add_markers(
                        p,
                        data = plot_data_clean,
                        x = ~x_jitter,
                        y = ~Value,
                        color = ~ get(color_col_name),
                        colors = colors_map,
                        marker = list(size = 6, opacity = 0.8),
                        showlegend = TRUE
                    )
                }
            } else {
                p <- plotly::add_markers(
                    p,
                    data = plot_data_clean,
                    x = ~x_jitter,
                    y = ~Value,
                    marker = list(size = 5, opacity = 0.6, color = "rgba(127, 140, 141, 0.8)"),
                    showlegend = FALSE
                )
            }

            p <- plotly::layout(
                p,
                title = list(
                    text = paste("<b>Phenotype Distribution:", safe_trait_name, "</b>"),
                    x = 0.5
                ),
                xaxis = list(
                    title = paste("<b>", xaxis_title, "</b>"),
                    tickmode = "array",
                    tickvals = tick_vals,
                    ticktext = tick_text,
                    tickangle = tick_angle,
                    tickfont = list(size = tick_font_size),
                    automargin = TRUE
                ),
                yaxis = list(title = "<b>Phenotype Value</b>")
            )

            if (isTRUE(input$show_p_annotations)) {
                p <- plotly::layout(p, annotations = annotations)
            }

            return(p)
        })
    })
}
