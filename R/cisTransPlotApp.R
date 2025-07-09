# R/cisTransPlotApp.R

#' Cis-Trans Plot Module UI components
#'
#' @param id Module ID
#' @export
library(ggplot2)
library(dplyr)

#' Helper function to map dataset and interaction type to qtlxcovar file path
get_qtlxcovar_file_path <- function(base_dataset, interaction_type) {
  qtlxcovar_dir <- "/data/dev/miniViewer_3.0"

  # Map dataset categories to file prefixes
  if (grepl("Liver.*Genes", base_dataset, ignore.case = TRUE)) {
    file_prefix <- "DO1200_liver_genes_all_mice"
  } else if (grepl("Clinical", base_dataset, ignore.case = TRUE)) {
    file_prefix <- "DO1200_clinical_traits_all_mice"
  } else {
    message("get_qtlxcovar_file_path: Unknown dataset category for:", base_dataset)
    return(NULL)
  }

  # Map interaction type to file suffix
  if (interaction_type == "sex") {
    file_suffix <- "qtlxsex_peaks.csv"
  } else if (interaction_type == "diet") {
    file_suffix <- "qtlxdiet_peaks.csv"
  } else {
    message("get_qtlxcovar_file_path: Unknown interaction type:", interaction_type)
    return(NULL)
  }

  file_path <- file.path(qtlxcovar_dir, paste0(file_prefix, "_", file_suffix))
  return(file_path)
}

#' @param id Module ID
#' @export
cisTransPlotInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    # The UI for selecting dataset within cisTransPlotApp is removed,
    # as this will be driven by the main app's selection.
    # uiOutput(ns("dataset_selector_ui"))
  )
}

#' Cis-Trans Plot Module Output UI
#'
#' @param id Module ID
#' @export
cisTransPlotUI <- function(id) {
  ns <- shiny::NS(id)
  shinycssloaders::withSpinner(
    plotly::plotlyOutput(
      ns("cis_trans_plot_output"),
      height = "100%",
      width = "100%"
    )
  )
}

#' Cis-Trans Plot Module Server
#'
#' @param id Module ID
#' @param import_reactives Reactive list containing file_directory and annotation_list
#' @param main_par A REACTIVE that returns a list, expected to contain selected_dataset (which is a reactive group name)
#' @param peaks_cache Environment for caching peak finder results.
#' @export
cisTransPlotServer <- function(id, import_reactives, main_par, peaks_cache, sidebar_interaction_type = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Define the %||% operator for null coalescing
    `%||%` <- function(a, b) if (!is.null(a)) a else b

    # Removed internal dataset selection UI and logic
    # output$dataset_selector_ui <- renderUI({...})
    # cistrans_dataset_choices <- shiny::reactive({...})

    selected_dataset <- shiny::reactive({
      shiny::req(main_par(), main_par()$selected_dataset, main_par()$selected_dataset())
      dataset_group_name <- main_par()$selected_dataset()
      message("cisTransPlotServer: Received selected_dataset from main_par: ", dataset_group_name)
      return(dataset_group_name)
    })

    trait_type <- shiny::reactive({
      shiny::req(import_reactives(), selected_dataset())
      # Ensure import_reactives() is called to get the list
      current_trait_type <- get_trait_type(import_reactives(), selected_dataset())
      message("cisTransPlotServer: Determined Trait Type: ", current_trait_type)
      current_trait_type
    })


    # Create a debounced reactive for interaction type to prevent excessive reloads
    interaction_type_debounced <- shiny::reactive({
      if (!is.null(sidebar_interaction_type)) sidebar_interaction_type() else "none"
    }) %>% shiny::debounce(500) # Debounce to prevent rapid firing

    peaks_data <- shiny::reactive({
      shiny::req(import_reactives(), selected_dataset(), trait_type())

      # Check if we should use qtlxcovar data
      interaction_type <- interaction_type_debounced()
      base_dataset <- selected_dataset()

      message("cisTransPlotServer: peaks_data reactive. Dataset:", base_dataset, "Interaction:", interaction_type)

      # Use qtlxcovar data for HC_HF datasets with interaction types
      use_qtlxcovar <- !is.null(interaction_type) && interaction_type != "none" &&
        !is.null(base_dataset) && grepl("^HC_HF", base_dataset, ignore.case = TRUE)

      if (use_qtlxcovar) {
        message("cisTransPlotServer: Using qtlxcovar peaks for interaction:", interaction_type)

        # Map dataset to qtlxcovar file
        qtlxcovar_file <- get_qtlxcovar_file_path(base_dataset, interaction_type)

        if (!is.null(qtlxcovar_file) && file.exists(qtlxcovar_file)) {
          message("cisTransPlotServer: Loading qtlxcovar file:", qtlxcovar_file)
          found_peaks <- tryCatch(
            {
              data.table::fread(qtlxcovar_file)
            },
            error = function(e) {
              message("cisTransPlotServer: Error reading qtlxcovar file:", e$message)
              return(NULL)
            }
          )

          if (!is.null(found_peaks) && nrow(found_peaks) > 0) {
            message("cisTransPlotServer: Loaded", nrow(found_peaks), "qtlxcovar peaks")
          } else {
            message("cisTransPlotServer: No qtlxcovar peaks loaded")
          }
        } else {
          message("cisTransPlotServer: QTLxcovar file not found or invalid path:", qtlxcovar_file)
          found_peaks <- NULL
        }
      } else {
        message("cisTransPlotServer: Using regular peaks via peak_finder")

        current_peaks_cache <- if (is.environment(peaks_cache)) peaks_cache else NULL

        # Ensure import_reactives() is called to get the list
        found_peaks <- peak_finder(import_reactives()$file_directory,
          selected_dataset(),
          trait_type = trait_type(),
          cache_env = current_peaks_cache,
          use_cache = TRUE
        )
      }
      message("cisTransPlotServer: Data returned by peak_finder:")
      if (is.data.frame(found_peaks)) {
        message("peak_finder returned a data.frame with ", nrow(found_peaks), " rows and columns: ", paste(colnames(found_peaks), collapse = ", "))
        if (nrow(found_peaks) > 0) {
          message("Head of peak_finder output:")
          print(utils::head(found_peaks))
        } else {
          message("peak_finder returned an empty data.frame.")
        }
      } else {
        message("peak_finder did not return a data.frame. Object type: ", class(found_peaks))
        print(found_peaks)
      }
      found_peaks
    })


    plot_data <- shiny::reactive({
      df <- peaks_data()
      shiny::req(main_par(), main_par()$LOD_thr) # Ensure LOD threshold is available
      lod_threshold <- main_par()$LOD_thr() # Get the LOD threshold from the slider
      interaction_type <- interaction_type_debounced()

      message("cisTransPlotServer: Preparing plot_data. Initial df from peaks_data() has ", if (is.null(df)) "NULL" else nrow(df), " rows.") # DEBUG
      if (is.null(df) || nrow(df) == 0) {
        message("cisTransPlotServer: plot_data - peaks_data() is NULL or empty, returning NULL.") # DEBUG
        return(NULL)
      }

      # Check if this is qtlxcovar data (interactive analysis)
      is_qtlxcovar_data <- "lod_diff" %in% colnames(df)
      use_lod_diff <- is_qtlxcovar_data && !is.null(interaction_type) && interaction_type != "none"

      # Determine which LOD column to use
      lod_column <- if (use_lod_diff) "lod_diff" else "qtl_lod"

      # Required columns depend on whether we're using qtlxcovar data
      req_cols <- c("qtl_chr", "qtl_pos", "cis", "gene_chr", "gene_start", lod_column)
      message("cisTransPlotServer: Required columns for plot: ", paste(req_cols, collapse = ", ")) # DEBUG
      message("cisTransPlotServer: Columns in df from peaks_data(): ", paste(colnames(df), collapse = ", ")) # DEBUG
      message("cisTransPlotServer: Using LOD column: ", lod_column, " (qtlxcovar: ", is_qtlxcovar_data, ")") # DEBUG

      # Check for gene_symbol column specifically
      if ("gene_symbol" %in% colnames(df)) {
        message("cisTransPlotServer: gene_symbol column found with ", sum(!is.na(df$gene_symbol)), " non-NA values out of ", nrow(df), " rows")
        # Show sample of gene symbols
        sample_symbols <- head(unique(df$gene_symbol[!is.na(df$gene_symbol)]), 5)
        message("cisTransPlotServer: Sample gene symbols: ", paste(sample_symbols, collapse = ", "))
      } else {
        message("cisTransPlotServer: WARNING - gene_symbol column NOT found in peaks data!")
        message("cisTransPlotServer: Available columns: ", paste(colnames(df), collapse = ", "))
      }

      if (!all(req_cols %in% colnames(df))) {
        missing_cols <- req_cols[!req_cols %in% colnames(df)]
        message("cisTransPlotServer: plot_data - Missing required columns: ", paste(missing_cols, collapse = ", "), ". Returning NULL.") # DEBUG
        return(NULL)
      }

      # Filter out NA values - use the appropriate LOD column
      if (use_lod_diff) {
        df_filtered <- dplyr::filter(df, !is.na(qtl_chr) & !is.na(qtl_pos) & !is.na(cis) & !is.na(gene_chr) & !is.na(gene_start) & !is.na(lod_diff))
      } else {
        df_filtered <- dplyr::filter(df, !is.na(qtl_chr) & !is.na(qtl_pos) & !is.na(cis) & !is.na(gene_chr) & !is.na(gene_start) & !is.na(qtl_lod))
      }
      message("cisTransPlotServer: plot_data - After filtering NAs, df_filtered has ", nrow(df_filtered), " rows.") # DEBUG

      if (nrow(df_filtered) == 0) {
        message("cisTransPlotServer: plot_data - After NA filtering, data frame is empty. Returning NULL.") # DEBUG
        return(NULL)
      }

      # FILTER BY LOD THRESHOLD - Use appropriate column and handle absolute values for difference data
      if (use_lod_diff) {
        # For difference data, use absolute value for threshold filtering
        message(paste0("cisTransPlotServer: Applying LOD difference threshold filter (abs value). Before: ", nrow(df_filtered), " rows, threshold: ", lod_threshold))
        df_filtered <- dplyr::filter(df_filtered, abs(lod_diff) >= lod_threshold)

        # Add a standardized qtl_lod column for plotting (use absolute value of lod_diff)
        df_filtered$qtl_lod <- abs(df_filtered$lod_diff)

        # Store original lod_diff for hover text
        df_filtered$lod_diff_original <- df_filtered$lod_diff
      } else {
        message(paste0("cisTransPlotServer: Applying LOD threshold filter. Before: ", nrow(df_filtered), " rows, threshold: ", lod_threshold))
        df_filtered <- dplyr::filter(df_filtered, qtl_lod >= lod_threshold)
      }
      message(paste0("cisTransPlotServer: After LOD threshold filter: ", nrow(df_filtered), " rows"))

      if (nrow(df_filtered) == 0) {
        message("cisTransPlotServer: plot_data - After LOD threshold filtering, data frame is empty. Returning NULL.") # DEBUG
        return(NULL)
      }

      chrom_levels <- c(as.character(1:19), "X")
      df_filtered$qtl_chr <- factor(chr_XYM(df_filtered$qtl_chr), levels = chrom_levels, ordered = TRUE)
      df_filtered$gene_chr <- factor(chr_XYM(df_filtered$gene_chr), levels = chrom_levels, ordered = TRUE)
      # Ensure cis is logical or character
      if (is.logical(df_filtered$cis)) {
        df_filtered$cis <- as.character(df_filtered$cis)
      } else if (!is.character(df_filtered$cis)) {
        df_filtered$cis <- as.character(as.logical(df_filtered$cis))
      }

      # Add metadata for plot rendering
      attr(df_filtered, "is_qtlxcovar_data") <- is_qtlxcovar_data
      attr(df_filtered, "interaction_type") <- interaction_type

      df_filtered
    })

    output$cis_trans_plot_output <- plotly::renderPlotly({
      plot_data_filtered <- plot_data()

      if (is.null(plot_data_filtered) || nrow(plot_data_filtered) == 0) {
        message("cisTransPlotServer: No data to plot, rendering null plot.")
        return(plotly::ggplotly(plot_null("No cis/trans data available for this dataset.")))
      }

      current_selected_dataset_name <- selected_dataset()

      plot_data_filtered_dt <- data.table::as.data.table(plot_data_filtered)

      # Filter out rows with NA gene_symbol
      plot_data_filtered_dt <- plot_data_filtered_dt[!is.na(gene_symbol)]

      # Ensure 'cis' column is character "TRUE" or "FALSE" for scale_color_manual
      if (is.logical(plot_data_filtered_dt$cis)) {
        plot_data_filtered_dt[, cis_char := ifelse(cis, "TRUE", "FALSE")]
      } else {
        plot_data_filtered_dt[, cis_char := as.character(cis)] # Assume it's already proper if not logical
      }

      # Get metadata from plot data
      is_qtlxcovar_data <- attr(plot_data_filtered, "is_qtlxcovar_data") %||% FALSE
      interaction_type <- attr(plot_data_filtered, "interaction_type") %||% "none"

      # Create appropriate hover text based on data type
      if (is_qtlxcovar_data && interaction_type != "none") {
        # For difference data, show both original difference and absolute value
        plot_data_filtered_dt[, hover_text := paste0(
          gene_symbol, "<br>",
          "LOD Difference: ", round(lod_diff_original, 2), "<br>",
          "|LOD Difference|: ", round(qtl_lod, 2)
        )]
      } else {
        # For regular data, show standard LOD
        plot_data_filtered_dt[, hover_text := paste0(
          gene_symbol, "<br>",
          "LOD: ", round(qtl_lod, 2)
        )]
      }

      cis.colors <- c("FALSE" = "#E41A1C", "TRUE" = "blue")
      names(cis.colors) <- c("FALSE", "TRUE")

      # Reverse the order of gene_chr factor levels
      plot_data_filtered_dt[, gene_chr := factor(gene_chr, levels = rev(levels(gene_chr)))]

      # Try native plotly approach for better hover performance with dense data
      # Create separate traces for cis and trans points
      cis_points <- plot_data_filtered_dt[cis_char == "TRUE"]
      trans_points <- plot_data_filtered_dt[cis_char == "FALSE"]

      fig <- plotly::plot_ly(source = ns("cistrans_plotly")) %>%
        # Add cis points (blue)
        plotly::add_markers(
          data = cis_points,
          x = ~qtl_pos, y = ~gene_start,
          color = I("#0066CC"),
          size = I(8),
          alpha = 0.7,
          text = ~hover_text,
          customdata = ~gene_symbol,
          hovertemplate = "%{text}<extra></extra>",
          name = "Cis",
          showlegend = TRUE
        ) %>%
        # Add trans points (red)
        plotly::add_markers(
          data = trans_points,
          x = ~qtl_pos, y = ~gene_start,
          color = I("#E41A1C"),
          size = I(8),
          alpha = 0.7,
          text = ~hover_text,
          customdata = ~gene_symbol,
          hovertemplate = "%{text}<extra></extra>",
          name = "Trans",
          showlegend = TRUE
        ) %>%
        plotly::layout(
          title = if (is_qtlxcovar_data && interaction_type != "none") {
            paste("Cis/Trans QTL Plot -", stringr::str_to_title(interaction_type), "Interaction Differences")
          } else {
            "Cis/Trans QTL Plot"
          },
          xaxis = list(title = "QTL Position"),
          yaxis = list(title = "Gene Position"),
          hovermode = "closest",
          dragmode = "pan",
          showlegend = TRUE,
          margin = list(l = 50, r = 50, t = 50, b = 50)
        ) %>%
        plotly::config(
          displaylogo = FALSE,
          modeBarButtonsToRemove = c("select2d", "lasso2d", "toggleSpikelines", "sendDataToCloud"),
          scrollZoom = TRUE,
          doubleClick = "reset",
          responsive = TRUE
        )

      fig <- plotly::event_register(fig, "plotly_click")

      return(fig)
    })

    clicked_phenotype_for_lod_scan_rv <- shiny::reactiveVal(NULL)

    shiny::observeEvent(plotly::event_data("plotly_click", source = ns("cistrans_plotly")), {
      click_data <- plotly::event_data("plotly_click", source = ns("cistrans_plotly"))
      if (!is.null(click_data) && !is.null(click_data$customdata) && length(click_data$customdata) > 0) {
        clicked_identifier <- click_data$customdata[[1]]
        clicked_phenotype_for_lod_scan_rv(clicked_identifier)
        message(paste("CisTransPlot Clicked! Identifier for LOD scan (should be symbol):", clicked_identifier))
      } else {
        clicked_phenotype_for_lod_scan_rv(NULL)
      }
    })

    return(list(
      clicked_phenotype_for_lod_scan = clicked_phenotype_for_lod_scan_rv
    ))
  })
}
