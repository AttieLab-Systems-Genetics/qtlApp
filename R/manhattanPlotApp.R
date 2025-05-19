# R/manhattanPlotApp.R

#' Manhattan Plot Module UI
#'
#' @param id Module ID.
#' @export
manhattanPlotInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    if (exists("create_select_input", mode = "function")) {
      create_select_input(ns("phenotype_class_selector"),
                          label = "Select Phenotype Class:",
                          choices = c("Clinical Traits" = "clinical_trait", 
                                        "Liver Lipids" = "Lipids"), # Changed to "Lipids" (capital L)
                          selected = "clinical_trait")
    } else {
      shiny::selectInput(ns("phenotype_class_selector"),
                         label = "Select Phenotype Class:",
                         choices = c("Clinical Traits" = "clinical_trait", 
                                       "Liver Lipids" = "Lipids"), # Changed to "Lipids" (capital L)
                         selected = "clinical_trait")
    }
  )
}

#' Manhattan Plot Module UI Output
#'
#' @param id Module ID.
#' @export
manhattanPlotUI <- function(id) {
  ns <- shiny::NS(id)
  shinycssloaders::withSpinner(plotly::plotlyOutput(ns("manhattan_plot_output"), height = "600px"))
}

#' Manhattan Plot Module Server
#'
#' @param id Module ID.
#' @param import_reactives Reactive list containing file_directory.
#' @param main_par Reactive list containing selected_dataset.
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
        input$phenotype_class_selector, # Now reacts to this
        import_reactives(),
        import_reactives()$file_directory
      )
      file_dir_dt <- data.table::as.data.table(import_reactives()$file_directory)
      selected_phen_class_value <- input$phenotype_class_selector # "clinical_trait" or "Lipids"

      # Determine target trait_type in file_index.csv based on selector
      # The values of input$phenotype_class_selector are "clinical_trait" and "Lipids"
      # In file_index.csv, the trait_type column has "Clinical Traits" and "Lipids"
      # So, we need a mapping if they are not identical.
      # For "Lipids", it's a direct match.
      # For "clinical_trait", the file_index.csv uses "Clinical Traits" (with a space and capital T)
      
      target_trait_type <- ""
      if (selected_phen_class_value == "clinical_trait") {
        target_trait_type <- "Clinical Traits"
      } else if (selected_phen_class_value == "Lipids") {
        target_trait_type <- "Lipids"
      } else {
        shiny::showNotification(paste("ManhattanPlot: Unknown phenotype_class_selector value:", selected_phen_class_value), type = "error")
        return(NULL)
      }

      # Filter for rows that are 'peaks' files and match the target_trait_type
      # Also, let's add a preference for scan_type == "additive" if multiple exist
      
      if (!("trait_type" %in% names(file_dir_dt))) {
          shiny::showNotification("ManhattanPlot: Column 'trait_type' not found in file_index.csv. Cannot select peaks file based on phenotype class.", type="error", duration=NULL)
          return(NULL)
      }
      if (!("file_type" %in% names(file_dir_dt))) {
          shiny::showNotification("ManhattanPlot: Column 'file_type' not found in file_index.csv.", type="error", duration=NULL)
          return(NULL)
      }
      if (!("File_path" %in% names(file_dir_dt))) {
          shiny::showNotification("ManhattanPlot: Column 'File_path' not found in file_index.csv.", type="error", duration=NULL)
          return(NULL)
      }
       if (!("scan_type" %in% names(file_dir_dt))) { # Ensure scan_type column exists for preference
          shiny::showNotification("ManhattanPlot: Column 'scan_type' not found in file_index.csv. Cannot prefer additive scan.", type="warning", duration=NULL)
          # Proceed without scan_type preference if column is missing
          peaks_rows <- file_dir_dt[file_type == "peaks" & trait_type == target_trait_type]
      } else {
          peaks_rows <- file_dir_dt[file_type == "peaks" & trait_type == target_trait_type]
      }


      if (nrow(peaks_rows) == 0) {
        shiny::showNotification(paste("ManhattanPlot: No 'peaks' file found in file_index.csv for trait_type:", target_trait_type), type = "warning", duration = NULL)
        return(NULL)
      }
      
      selected_row <- NULL
      if (nrow(peaks_rows) > 1 && "scan_type" %in% names(peaks_rows)) {
        # Try to find an "additive" scan_type
        additive_peaks_rows <- peaks_rows[scan_type == "additive"]
        if (nrow(additive_peaks_rows) > 0) {
          selected_row <- additive_peaks_rows[1, ] # Take the first additive one
          if (nrow(additive_peaks_rows) > 1) {
            shiny::showNotification(paste("ManhattanPlot: Multiple 'additive' peaks files found for trait_type:", target_trait_type, ". Selecting the first one."), type = "info")
          }
        } else {
          # No additive found, take the first one from the original multi-row list
          selected_row <- peaks_rows[1, ]
          shiny::showNotification(paste("ManhattanPlot: Multiple peaks files found for trait_type:", target_trait_type, "and none are 'additive'. Selecting the first one found."), type = "info")
        }
      } else {
        # Only one row found, or scan_type column was missing
        selected_row <- peaks_rows[1, ]
      }

      if (!is.null(selected_row)) {
        peaks_file_path <- selected_row$File_path[[1]]
        if (!is.na(peaks_file_path) && nzchar(peaks_file_path)) {
          if (file.exists(peaks_file_path)) {
            message(paste("ManhattanPlot: Loading peaks file for", target_trait_type, ":", peaks_file_path)) # Added message
            return(peaks_file_path)
          } else {
            shiny::showNotification(paste("ManhattanPlot: Peaks file specified in file_index.csv does not exist:", peaks_file_path), type = "error", duration = NULL)
            return(NULL)
          }
        } else {
          shiny::showNotification(paste("ManhattanPlot: 'File_path' for selected peaks file is empty or NA. Trait_type:", target_trait_type), type = "warning", duration = NULL)
          return(NULL)
        }
      } else {
        # This case should ideally not be reached if nrow(peaks_rows) > 0
        shiny::showNotification(paste("ManhattanPlot: Could not select a specific peaks file for trait_type:", target_trait_type), type = "error")
        return(NULL)
      }
    })

    raw_peaks_data <- shiny::reactive({
      shiny::req(selected_peaks_file())
      file_path <- selected_peaks_file()
      tryCatch({
        dt <- data.table::fread(file_path)
        # Ensure required columns are present
        req_cols <- c("phenotype_class", "phenotype", "qtl_lod", "qtl_chr", "qtl_pos")
        missing_cols <- req_cols[!req_cols %in% colnames(dt)]
        if (length(missing_cols) > 0) {
          shiny::showNotification(paste("Missing required columns in peaks file:", paste(missing_cols, collapse = ", ")),
                                type = "error", duration = NULL)
          return(NULL)
        }
        return(data.table::as.data.table(dt))
      }, error = function(e) {
        shiny::showNotification(paste("Error loading peaks file:", file_path, e$message), type = "error", duration = NULL)
        return(NULL)
      })
    })

    plot_data_prep <- shiny::reactive({
      shiny::req(raw_peaks_data(), input$phenotype_class_selector)
      peaks_dt <- data.table::copy(raw_peaks_data())
      
      peaks_dt <- peaks_dt[phenotype_class == input$phenotype_class_selector]
      
      if (nrow(peaks_dt) == 0) {
        shiny::showNotification(paste("No data for phenotype class:", input$phenotype_class_selector), type = "warning")
        return(data.table::data.table()) 
      }

      if (!is.numeric(peaks_dt$qtl_lod)) {
          peaks_dt[, qtl_lod := as.numeric(as.character(qtl_lod))]
      }
      peaks_dt <- peaks_dt[!is.na(qtl_lod)]

      # Convert qtl_chr to a factor for ordered faceting
      # Keep original qtl_chr for display, make a factor for ordering
      peaks_dt[, chr_factor := factor(qtl_chr, 
                                      levels = c(as.character(1:19), "X", "Y", "M"), 
                                      ordered = TRUE)]
      peaks_dt <- peaks_dt[!is.na(chr_factor)] # Remove rows where chr is not in 1-19,X,Y,M
      
      if (!is.numeric(peaks_dt$qtl_pos)) {
          peaks_dt[, qtl_pos := as.numeric(as.character(qtl_pos))]
      }
      peaks_dt <- peaks_dt[!is.na(qtl_pos)]
            
      if (nrow(peaks_dt) == 0) {
        shiny::showNotification(paste("No valid data after processing for phenotype class:", input$phenotype_class_selector), type = "warning")
        return(data.table::data.table())
      }
      
      # Sort by chromosome factor and position for consistent plotting (though facet handles separation)
      data.table::setkey(peaks_dt, chr_factor, qtl_pos)
      
      return(peaks_dt)
    })

    output$manhattan_plot_output <- plotly::renderPlotly({
      df_to_plot <- plot_data_prep()
      shiny::req(df_to_plot, nrow(df_to_plot) > 0)
      
      # Prepare hover text (marker column is assumed to exist from your CSV header)
      df_to_plot[, hover_text := paste(
        "Phenotype:", phenotype,
        "<br>Marker:", marker, # Using 'marker' column for tooltip
        "<br>LOD:", round(qtl_lod, 2),
        "<br>Chr:", qtl_chr,
        "<br>Pos (Mbp):", round(qtl_pos, 2)
      )]
      
      # Chromosome colors for alternating pattern in facets (optional, but common)
      # For now, all points are one color, faceting provides separation.

      p <- ggplot(df_to_plot, aes(x = qtl_pos, y = qtl_lod, text = hover_text, key = marker)) + # Added key = marker for plotly
        geom_point(alpha = 0.7, size = 1.5, color = "#2c3e50") + 
        scale_y_continuous(name = "LOD Score", expand = expansion(mult = c(0, 0.05))) +
        labs(title = paste("Manhattan Plot for", input$phenotype_class_selector), x = "Position (Mbp)") +
        facet_grid(. ~ chr_factor, scales = "free_x", space = "free_x", switch = "x") + # Facet by chromosome
        theme_minimal(base_size = 11) +
        theme(
          panel.spacing.x = unit(0.1, "lines"), # Reduced space between facets
          strip.background = element_blank(), # Remove facet label background
          strip.placement = "outside", # Place facet labels outside plot
          axis.text.x = element_text(angle = 45, hjust = 1, size=8), # Angled X-axis text for readability
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
      
      # Convert to plotly
      # Setting a source for event_data if needed later for interactivity, e.g. clicking a point
      plotly_p <- plotly::ggplotly(p, tooltip = "text", source = ns("manhattan_plotly")) %>%
        plotly::layout(dragmode = "pan") %>%
        plotly::config(displaylogo = FALSE)
      
      return(plotly_p)
    })
  })
} 