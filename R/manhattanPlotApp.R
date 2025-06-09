# R/manhattanPlotApp.R

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
manhattanPlotServer <- function(id, import_reactives, main_par) {
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

      message(paste("ManhattanPlot: Received selected_dataset_group_name:", selected_dataset_group_name))

      if (is.null(selected_dataset_group_name) || !nzchar(selected_dataset_group_name)) {
        shiny::showNotification("ManhattanPlot: No dataset group selected from main app.", type = "warning")
        return(NULL)
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

    raw_peaks_data <- shiny::reactive({
      shiny::req(selected_peaks_file())
      file_path <- selected_peaks_file()
      tryCatch(
        {
          dt <- data.table::fread(file_path)
          req_cols <- c("phenotype_class", "phenotype", "qtl_lod", "qtl_chr", "qtl_pos")
          missing_cols <- req_cols[!req_cols %in% colnames(dt)]
          if (length(missing_cols) > 0) {
            shiny::showNotification(paste("Missing required columns in peaks file (", basename(file_path), "):", paste(missing_cols, collapse = ", ")),
              type = "error", duration = NULL
            )
            return(NULL)
          }
          return(data.table::as.data.table(dt))
        },
        error = function(e) {
          shiny::showNotification(paste("Error loading peaks file (", basename(file_path), "):", e$message), type = "error", duration = NULL)
          return(NULL)
        }
      )
    })

    plot_data_prep <- shiny::reactive({
      shiny::req(raw_peaks_data(), main_par(), main_par()$selected_dataset, main_par()$LOD_thr, import_reactives(), import_reactives()$file_directory)

      peaks_dt <- data.table::copy(raw_peaks_data())
      selected_group_name <- main_par()$selected_dataset() # This is the 'group' name from file_index
      lod_threshold <- main_par()$LOD_thr() # Get the LOD threshold from the slider
      file_index <- data.table::as.data.table(import_reactives()$file_directory)

      # Determine the dataset_category for the selected_group_name
      current_dataset_info <- file_index[group == selected_group_name]
      current_dataset_category <- NULL
      if (nrow(current_dataset_info) > 0) {
        current_dataset_category <- current_dataset_info$dataset_category[1]
      } else {
        shiny::showNotification(paste("ManhattanPlot: Could not find 'group'", selected_group_name, "in file_index.csv to determine category."), type = "warning")
        return(NULL) # Or an empty data.table
      }

      shiny::validate(
        shiny::need(current_dataset_category, "ManhattanPlot: Dataset category could not be determined.")
      )

      # Determine the phenotype_class string to filter by, based on the dataset category
      target_phenotype_class_value <- NULL
      if (current_dataset_category == "Clinical Traits") {
        target_phenotype_class_value <- "clinical_trait"
      } else if (current_dataset_category == "Liver Lipids") {
        target_phenotype_class_value <- "liver_lipid"
      } else if (current_dataset_category == "Plasma 2H Metabolites") {
        target_phenotype_class_value <- "plasma_2H_metabolite"
      } else {
        # Fallback or error if the category is unexpected for Manhattan plots
        shiny::showNotification(paste("ManhattanPlot: Unexpected dataset category '", current_dataset_category, "' for Manhattan plot."), type = "warning")
        return(NULL) # Or an empty data.table
      }

      message(paste0(
        "ManhattanPlot DEBUG: For group '", selected_group_name, "' (category '", current_dataset_category, "'), ",
        "filtering by phenotype_class = '", target_phenotype_class_value, "'."
      ))

      if (nrow(peaks_dt) > 0 && "phenotype_class" %in% colnames(peaks_dt)) {
        message(paste0("ManhattanPlot DEBUG: Unique 'phenotype_class' values BEFORE filtering: ", paste(unique(peaks_dt$phenotype_class), collapse = ", ")))
      } else if (nrow(peaks_dt) == 0) {
        message("ManhattanPlot DEBUG: raw_peaks_data for plot_data_prep is empty BEFORE filtering.")
      } else {
        message("ManhattanPlot DEBUG: 'phenotype_class' column not found in raw_peaks_data for plot_data_prep.")
        return(NULL) # Critical column missing
      }

      # Actual filtering based on the determined target_phenotype_class_value
      peaks_dt <- peaks_dt[phenotype_class == target_phenotype_class_value]

      if (nrow(peaks_dt) == 0) {
        shiny::showNotification(
          paste0(
            "No data after filtering for phenotype class: ", target_phenotype_class_value,
            " (for dataset '", selected_group_name, "'). Check if the peaks file contains this phenotype class and the correct string."
          ),
          type = "warning", duration = 10
        )

        return(data.table::data.table()) # Return empty table to avoid plot errors
      } else {
        message(paste0("ManhattanPlot DEBUG: Filtering for '", target_phenotype_class_value, "' resulted in ", nrow(peaks_dt), " rows."))
      }

      # Ensure all subsequent code uses peaks_dt_filtered
      if (!is.numeric(peaks_dt$qtl_lod)) {
        peaks_dt[, qtl_lod := as.numeric(as.character(qtl_lod))]
      }
      peaks_dt <- peaks_dt[!is.na(qtl_lod)]

      # FILTER BY LOD THRESHOLD - Only show points at or above the threshold
      message(paste0("ManhattanPlot DEBUG: Applying LOD threshold filter. Before: ", nrow(peaks_dt), " rows, threshold: ", lod_threshold))
      peaks_dt <- peaks_dt[qtl_lod >= lod_threshold]
      message(paste0("ManhattanPlot DEBUG: After LOD threshold filter: ", nrow(peaks_dt), " rows"))

      peaks_dt[, chr_factor := factor(qtl_chr,
        levels = c(as.character(1:19), "X", "Y", "M"),
        ordered = TRUE
      )]
      peaks_dt <- peaks_dt[!is.na(chr_factor)]

      if (!is.numeric(peaks_dt$qtl_pos)) {
        peaks_dt[, qtl_pos := as.numeric(as.character(qtl_pos))]
      }
      peaks_dt <- peaks_dt[!is.na(qtl_pos)]

      if (nrow(peaks_dt) == 0) {
        message("ManhattanPlot DEBUG: No data remaining after all filters")
        return(data.table::data.table())
      }

      data.table::setkey(peaks_dt, chr_factor, qtl_pos)

      return(peaks_dt)
    })

    output$manhattan_plot_output <- plotly::renderPlotly({
      df_to_plot <- plot_data_prep()
      shiny::req(df_to_plot, nrow(df_to_plot) > 0)

      df_to_plot[, hover_text := paste(
        "Phenotype:", phenotype, # Original phenotype for display
        "<br>Marker:", marker,
        "<br>LOD:", round(qtl_lod, 2),
        "<br>Chr:", qtl_chr,
        "<br>Pos (Mbp):", round(qtl_pos, 2)
      )]

      p <- ggplot(df_to_plot, aes(x = qtl_pos, y = qtl_lod, text = hover_text, key = marker, customdata = phenotype)) +
        geom_point(alpha = 0.7, size = 1.5, color = "#2c3e50") +
        scale_y_continuous(name = "LOD Score", expand = expansion(mult = c(0, 0.05))) +
        labs(title = paste("Manhattan Plot for", main_par()$selected_dataset(), "(LOD ≥", main_par()$LOD_thr(), ")"), x = "Position (Mbp)") +
        facet_grid(. ~ chr_factor, scales = "free_x", space = "free_x", switch = "x") + # Facet by chromosome
        theme_minimal(base_size = 11) +
        theme(
          panel.spacing.x = unit(0.1, "lines"), # Reduced space between facets
          strip.background = element_blank(), # Remove facet label background
          strip.placement = "outside", # Place facet labels outside plot
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8), # Angled X-axis text for readability
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold")
        )

      plotly_p <- plotly::ggplotly(p, tooltip = "text", source = ns("manhattan_plotly")) %>%
        plotly::layout(dragmode = "pan") %>%
        plotly::config(displaylogo = FALSE)

      # Register the plotly_click event
      plotly_p <- plotly::event_register(plotly_p, "plotly_click")

      return(plotly_p)
    })

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
