# R/manhattanPlotApp.R

#' Helper function to map dataset and interaction type to qtlxcovar file path for Manhattan plots
get_qtlxcovar_file_path_manhattan <- function(base_dataset, interaction_type) {
  qtlxcovar_dir <- "/data/dev/miniViewer_3.0"

  # Map dataset categories to file prefixes
  if (grepl("Clinical", base_dataset, ignore.case = TRUE)) {
    file_prefix <- "DO1200_clinical_traits_all_mice"
  } else if (grepl("Liver.*Lipid", base_dataset, ignore.case = TRUE)) {
    file_prefix <- "DO1200_liver_lipids_all_mice"
  } else if (grepl("Plasma.*Metabol", base_dataset, ignore.case = TRUE)) {
    file_prefix <- "DO1200_plasma_metabolites_all_mice"
  } else {
    message("get_qtlxcovar_file_path_manhattan: Unknown dataset category for:", base_dataset)
    return(NULL)
  }

  # Map interaction type to file suffix
  if (interaction_type == "sex") {
    file_suffix <- "qtlxsex_peaks.csv"
  } else if (interaction_type == "diet") {
    file_suffix <- "qtlxdiet_peaks.csv"
  } else {
    message("get_qtlxcovar_file_path_manhattan: Unknown interaction type:", interaction_type)
    return(NULL)
  }

  file_path <- file.path(qtlxcovar_dir, paste0(file_prefix, "_", file_suffix))

  if (file.exists(file_path)) {
    message(paste("get_qtlxcovar_file_path_manhattan: Found qtlxcovar file:", file_path))
    return(file_path)
  } else {
    message(paste("get_qtlxcovar_file_path_manhattan: File not found:", file_path))
    return(NULL)
  }
}

#' Manhattan Plot Module UI
#'
#' @param id Module ID.
#' @export
manhattanPlotInput <- function(id) {
  # This input is no longer needed as the plot will be determined by the main app's selection.
  # ns <- shiny::NS(id)
  # shiny::tagList(
  #   if (exists("create_select_input", mode = "function")) {
  #     create_select_input(ns("phenotype_class_selector"),
  #                         label = "Filter by Phenotype Class:",
  #                         choices = c("Clinical Traits" = "clinical_trait",
  #                                       "Liver Lipids" = "Lipids"),
  #                         selected = "clinical_trait")
  #   } else {
  #     shiny::selectInput(ns("phenotype_class_selector"),
  #                        label = "Filter by Phenotype Class:",
  #                        choices = c("Clinical Traits" = "clinical_trait",
  #                                      "Liver Lipids" = "Lipids"),
  #                        selected = "clinical_trait")
  #   }
  # )
  return(NULL) # Return NULL as there's no UI content for this input anymore
}

#' Manhattan Plot Module UI Output
#'
#' @param id Module ID.
#' @export
manhattanPlotUI <- function(id) {
  ns <- shiny::NS(id)
  shinycssloaders::withSpinner(plotly::plotlyOutput(ns("manhattan_plot_output"), height = "calc(65vh - 40px)"))
}

#' Manhattan Plot Module Server
#'
#' @param id Module ID.
#' @param import_reactives Reactive that returns a list containing file_directory.
#' @param main_par Reactive that re turns a list containing selected_dataset (which is a reactive group name).
#'
#' @importFrom dplyr %>% filter select mutate arrange distinct group_by summarise
#' @importFrom data.table fread as.data.table setkey
#' @import ggplot2
#' @import plotly
#' @export
manhattanPlotServer <- function(id, import_reactives, main_par, sidebar_interaction_type = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    selected_peaks_file <- shiny::reactive({
      shiny::req(
        import_reactives(),
        import_reactives()$file_directory,
        main_par(), # main_par is a reactive returning a list
        main_par()$selected_dataset, # This element is the reactive for the group name
        main_par()$selected_dataset() # This is the actual selected group name value
      )

      file_dir_dt <- data.table::as.data.table(import_reactives()$file_directory)
      # selected_dataset_group_name is the 'group' from file_index.csv, passed by scanApp
      selected_dataset_group_name <- main_par()$selected_dataset()

      # Check if we should use qtlxcovar data
      interaction_type <- if (!is.null(sidebar_interaction_type)) sidebar_interaction_type() else "none"
      base_dataset <- selected_dataset_group_name

      message(paste("ManhattanPlot: Received selected_dataset_group_name:", selected_dataset_group_name, "interaction_type:", interaction_type))

      if (is.null(selected_dataset_group_name) || !nzchar(selected_dataset_group_name)) {
        shiny::showNotification("ManhattanPlot: No dataset group selected from main app.", type = "warning")
        return(NULL)
      }

      # Check if we should use qtlxcovar data instead of regular peaks
      if (!is.null(interaction_type) && interaction_type != "none") {
        qtlxcovar_file <- get_qtlxcovar_file_path_manhattan(base_dataset, interaction_type)
        if (!is.null(qtlxcovar_file)) {
          message(paste("ManhattanPlot: Using qtlxcovar file for", interaction_type, "interaction:", qtlxcovar_file))
          return(qtlxcovar_file)
        } else {
          message(paste("ManhattanPlot: No qtlxcovar file found for", interaction_type, "interaction, falling back to regular peaks"))
        }
      }

      # Find the row for this specific group that is a 'peaks' file.
      target_peaks_row <- file_dir_dt[group == selected_dataset_group_name & file_type == "peaks"]

      if (nrow(target_peaks_row) == 0) {
        shiny::showNotification(paste("ManhattanPlot: No 'peaks' file found in file_index.csv for group:", selected_dataset_group_name), type = "warning", duration = NULL)
        return(NULL)
      }

      if (nrow(target_peaks_row) > 1) {
        shiny::showNotification(paste("ManhattanPlot: Multiple 'peaks' files found for group:", selected_dataset_group_name, ". Selecting the first one."), type = "info")
        target_peaks_row <- target_peaks_row[1, ]
      }

      peaks_file_path <- target_peaks_row$File_path[[1]]

      if (!is.na(peaks_file_path) && nzchar(peaks_file_path)) {
        if (file.exists(peaks_file_path)) {
          message(paste("ManhattanPlot: Loading peaks file for group ", selected_dataset_group_name, ":", peaks_file_path))
          return(peaks_file_path)
        } else {
          shiny::showNotification(paste("ManhattanPlot: Peaks file specified for group does not exist:", peaks_file_path), type = "error", duration = NULL)
          return(NULL)
        }
      } else {
        shiny::showNotification(paste("ManhattanPlot: 'File_path' for peaks is empty or NA for group:", selected_dataset_group_name), type = "warning", duration = NULL)
        return(NULL)
      }
    })

    # Add caching for raw peaks data
    raw_peaks_data <- shiny::reactive({
      shiny::req(selected_peaks_file())
      file_path <- selected_peaks_file()

      # Create cache key from file path and modification time
      cache_key <- paste0(file_path, "_", file.info(file_path)$mtime)

      # Check if data is already cached
      if (exists("manhattan_cache", envir = .GlobalEnv)) {
        if (cache_key %in% names(.GlobalEnv$manhattan_cache)) {
          message("ManhattanPlot: Using cached peaks data")
          return(.GlobalEnv$manhattan_cache[[cache_key]])
        }
      } else {
        .GlobalEnv$manhattan_cache <- list()
      }

      tryCatch(
        {
          message(paste("ManhattanPlot: Loading and caching peaks data from:", basename(file_path)))
          dt <- data.table::fread(file_path, showProgress = FALSE)
          req_cols <- c("phenotype_class", "phenotype", "qtl_lod", "qtl_chr", "qtl_pos")
          missing_cols <- req_cols[!req_cols %in% colnames(dt)]
          if (length(missing_cols) > 0) {
            shiny::showNotification(paste("Missing required columns in peaks file (", basename(file_path), "):", paste(missing_cols, collapse = ", ")),
              type = "error", duration = NULL
            )
            return(NULL)
          }

          # Convert to data.table and optimize column types upfront
          dt <- data.table::as.data.table(dt)
          if ("qtl_lod" %in% colnames(dt)) dt[, qtl_lod := as.numeric(qtl_lod)]
          if ("qtl_pos" %in% colnames(dt)) dt[, qtl_pos := as.numeric(qtl_pos)]
          if ("lod_diff" %in% colnames(dt)) dt[, lod_diff := as.numeric(lod_diff)]

          # Cache the processed data
          .GlobalEnv$manhattan_cache[[cache_key]] <- dt

          # Keep cache size manageable (max 5 datasets)
          if (length(.GlobalEnv$manhattan_cache) > 5) {
            .GlobalEnv$manhattan_cache <- .GlobalEnv$manhattan_cache[-(1:2)]
          }

          return(dt)
        },
        error = function(e) {
          shiny::showNotification(paste("Error loading peaks file (", basename(file_path), "):", e$message), type = "error", duration = NULL)
          return(NULL)
        }
      )
    }) %>% shiny::debounce(100)

    plot_data_prep <- shiny::reactive({
      shiny::req(raw_peaks_data(), main_par(), main_par()$selected_dataset, main_par()$LOD_thr, import_reactives(), import_reactives()$file_directory)

      # Get reactive values once to avoid multiple evaluations
      peaks_dt <- raw_peaks_data() # Already a data.table, no need to copy unless modifying
      selected_group_name <- main_par()$selected_dataset()
      lod_threshold <- main_par()$LOD_thr()
      file_index <- data.table::as.data.table(import_reactives()$file_directory)
      interaction_type <- if (!is.null(sidebar_interaction_type)) sidebar_interaction_type() else "none"

      # Early exit checks
      if (is.null(peaks_dt) || nrow(peaks_dt) == 0) {
        return(data.table::data.table())
      }

      # Check if this is qtlxcovar data (already determined during loading)
      is_qtlxcovar_data <- "lod_diff" %in% colnames(peaks_dt)

      # Fast dataset category lookup (memoize this lookup)
      current_dataset_info <- file_index[group == selected_group_name]
      if (nrow(current_dataset_info) == 0) {
        shiny::showNotification(paste("ManhattanPlot: No data found for group:", selected_group_name), type = "warning")
        return(data.table::data.table())
      }
      current_dataset_category <- current_dataset_info$dataset_category[1]

      # Streamlined phenotype class mapping
      target_phenotype_class_values <- switch(current_dataset_category,
        "Clinical Traits" = "clinical_trait",
        "Liver Lipids" = "liver_lipid",
        "Plasma Metabolites" = c("plasma_13C_metabolite", "plasma_2H_metabolite"),
        NULL
      )

      if (is.null(target_phenotype_class_values)) {
        shiny::showNotification(paste("ManhattanPlot: Unsupported category:", current_dataset_category), type = "warning")
        return(data.table::data.table())
      }

      # Efficient filtering chain
      peaks_dt <- peaks_dt[phenotype_class %in% target_phenotype_class_values]

      if (nrow(peaks_dt) == 0) {
        return(data.table::data.table())
      }

      # Streamlined LOD column setup (columns already numeric from caching)
      if (is_qtlxcovar_data && interaction_type != "none") {
        peaks_dt <- peaks_dt[!is.na(lod_diff)]
        peaks_dt[, `:=`(plot_lod = abs(lod_diff), plot_lod_signed = lod_diff)]
        lod_col_name <- "plot_lod"
      } else {
        peaks_dt <- peaks_dt[!is.na(qtl_lod)]
        peaks_dt[, `:=`(plot_lod = qtl_lod, plot_lod_signed = qtl_lod)]
        lod_col_name <- "plot_lod"
      }

      # Apply LOD threshold filter
      peaks_dt <- peaks_dt[get(lod_col_name) >= lod_threshold]

      if (nrow(peaks_dt) == 0) {
        return(data.table::data.table())
      }

      # Efficient chromosome factor creation
      peaks_dt[, chr_factor := factor(qtl_chr, levels = c(as.character(1:19), "X", "Y", "M"), ordered = TRUE)]
      peaks_dt <- peaks_dt[!is.na(chr_factor) & !is.na(qtl_pos)]

      if (nrow(peaks_dt) == 0) {
        return(data.table::data.table())
      }

      # Set key for faster operations
      data.table::setkey(peaks_dt, chr_factor, qtl_pos)

      # Add metadata for plot creation
      attr(peaks_dt, "is_qtlxcovar_data") <- is_qtlxcovar_data
      attr(peaks_dt, "interaction_type") <- interaction_type

      return(peaks_dt)
    }) %>% shiny::debounce(150)

    output$manhattan_plot_output <- plotly::renderPlotly({
      df_to_plot <- plot_data_prep()
      shiny::req(df_to_plot, nrow(df_to_plot) > 0)

      # Get metadata from plot data
      is_qtlxcovar_data <- attr(df_to_plot, "is_qtlxcovar_data") %||% FALSE
      interaction_type <- attr(df_to_plot, "interaction_type") %||% "none"
      lod_display_col <- "plot_lod_signed"

      # Optimize for large datasets - sample points if too many
      if (nrow(df_to_plot) > 5000) {
        # Keep top peaks and random sample of others
        top_peaks <- df_to_plot[order(-get(lod_display_col))][1:min(2000, nrow(df_to_plot))]
        if (nrow(df_to_plot) > 2000) {
          # Use anti-join to get remaining rows - ensure marker column exists
          if ("marker" %in% colnames(df_to_plot)) {
            remaining <- df_to_plot[!top_peaks, on = "marker"]
          } else {
            # Fallback: use row indexing if marker column doesn't exist
            top_indices <- df_to_plot[order(-get(lod_display_col))][1:min(2000, nrow(df_to_plot)), which = TRUE]
            remaining <- df_to_plot[-top_indices]
          }
          sampled <- remaining[sample(.N, min(3000, .N))]
          df_to_plot <- rbind(top_peaks, sampled)
        } else {
          df_to_plot <- top_peaks
        }
        data.table::setkey(df_to_plot, chr_factor, qtl_pos)
      }

      # Streamlined hover text creation
      df_to_plot[, hover_text := paste0(
        "Phenotype: ", phenotype,
        "<br>Marker: ", marker,
        "<br>", if (is_qtlxcovar_data && interaction_type != "none") "LOD Difference" else "LOD Score", ": ", round(get(lod_display_col), 2),
        "<br>Chr: ", qtl_chr,
        "<br>Pos (Mbp): ", round(qtl_pos, 2)
      )]

      # Pre-calculate y-axis limits
      y_axis_limits <- if (is_qtlxcovar_data && interaction_type != "none") {
        c(4, max(abs(df_to_plot[[lod_display_col]]), na.rm = TRUE) * 1.1)
      } else {
        NULL
      }

      # Streamlined plot creation
      p <- ggplot2::ggplot(df_to_plot, ggplot2::aes(x = qtl_pos, y = get(lod_display_col), text = hover_text, key = marker, customdata = phenotype)) +
        ggplot2::geom_point(alpha = 0.6, size = 1.2, color = "#2c3e50") +
        {
          if (is.null(y_axis_limits)) {
            ggplot2::scale_y_continuous(
              name = if (is_qtlxcovar_data && interaction_type != "none") "LOD Difference" else "LOD Score",
              expand = ggplot2::expansion(mult = c(0.05, 0.05))
            )
          } else {
            ggplot2::scale_y_continuous(
              name = "LOD Difference", limits = y_axis_limits,
              expand = ggplot2::expansion(mult = c(0.02, 0.05))
            )
          }
        } +
        ggplot2::labs(
          title = paste0(
            if (is_qtlxcovar_data && interaction_type != "none") {
              paste("Manhattan Plot -", stringr::str_to_title(interaction_type), "Interaction Difference")
            } else {
              paste("Manhattan Plot for", main_par()$selected_dataset())
            },
            " (|LOD| ≥ ", main_par()$LOD_thr(), ")"
          ),
          x = "Position (Mbp)"
        ) +
        ggplot2::facet_grid(. ~ chr_factor, scales = "free_x", space = "free_x", switch = "x") +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
          panel.spacing.x = ggplot2::unit(0.1, "lines"),
          strip.background = ggplot2::element_blank(),
          strip.placement = "outside",
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11)
        )

      # Add horizontal line for difference plots
      if (is_qtlxcovar_data && interaction_type != "none") {
        p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7)
      }

      # Optimized plotly conversion
      plotly_p <- plotly::ggplotly(p, tooltip = "text", source = ns("manhattan_plotly")) %>%
        plotly::layout(
          dragmode = "pan",
          margin = list(l = 60, r = 40, b = 40, t = 60, pad = 4), # Add plot margins
          xaxis = list(
            title = "Chromosome",
            showgrid = FALSE,
            tickangle = -45,
            tickfont = list(size = 10),
            tickmode = "array",
            tickvals = 1:20,
            ticktext = c(1:19, "X", "Y", "M"),
            range = c(0.5, 20.5)
          ),
          yaxis = list(
            title = if (is_qtlxcovar_data && interaction_type != "none") "LOD Difference" else "LOD Score",
            showgrid = TRUE,
            gridcolor = "lightgray",
            gridwidth = 0.5,
            zeroline = TRUE,
            zerolinecolor = "black",
            zerolinewidth = 1,
            tickfont = list(size = 10),
            range = y_axis_limits
          ),
          hovermode = "closest",
          hoverlabel = list(
            bgcolor = "white",
            bordercolor = "gray",
            font = list(size = 12)
          ),
          legend = list(
            orientation = "h",
            x = 0.5,
            y = -0.1,
            font = list(size = 10)
          ),
          title = list(
            text = paste0(
              if (is_qtlxcovar_data && interaction_type != "none") {
                paste("Manhattan Plot -", stringr::str_to_title(interaction_type), "Interaction Difference")
              } else {
                paste("Manhattan Plot for", main_par()$selected_dataset())
              },
              " (|LOD| ≥ ", main_par()$LOD_thr(), ")"
            ),
            x = 0.5,
            y = 0.95,
            font = list(size = 12, face = "bold")
          ),
          showlegend = FALSE,
          autosize = TRUE
        ) %>%
        plotly::config(
          displaylogo = FALSE,
          modeBarButtonsToRemove = c("select2d", "lasso2d", "hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines"),
          responsive = TRUE
        ) %>%
        plotly::event_register("plotly_click")

      return(plotly_p)
    }) %>% shiny::debounce(200)

    clicked_phenotype_for_lod_scan_rv <- shiny::reactiveVal(NULL)

    shiny::observeEvent(plotly::event_data("plotly_click", source = ns("manhattan_plotly")), {
      click_data <- plotly::event_data("plotly_click", source = ns("manhattan_plotly"))
      if (!is.null(click_data) && !is.null(click_data$customdata) && length(click_data$customdata) > 0) {
        # customdata directly gives the value from the aesthetic mapping
        clicked_phenotype <- click_data$customdata[[1]]
        clicked_phenotype_for_lod_scan_rv(clicked_phenotype)
        message(paste("ManhattanPlot Clicked! Phenotype for LOD scan:", clicked_phenotype))
      } else {
        # Clear if click is not on a point with customdata (e.g., background)
        clicked_phenotype_for_lod_scan_rv(NULL)
      }
    })

    return(list(
      clicked_phenotype_for_lod_scan = clicked_phenotype_for_lod_scan_rv
    ))
  })
}
