# R/cisTransPlotApp.R

#' Cis-Trans Plot Module UI components
#'
#' @param id Module ID
#' @export
library(ggplot2)
library(dplyr)

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
  shinycssloaders::withSpinner(plotly::plotlyOutput(ns("cis_trans_plot_output"), height = "calc(65vh - 40px)"))
}

#' Cis-Trans Plot Module Server
#'
#' @param id Module ID
#' @param import_reactives Reactive list containing file_directory and annotation_list
#' @param main_par A REACTIVE that returns a list, expected to contain selected_dataset (which is a reactive group name)
#' @param peaks_cache Environment for caching peak finder results.
#' @export
cisTransPlotServer <- function(id, import_reactives, main_par, peaks_cache) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

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


    peaks_data <- shiny::reactive({
      shiny::req(import_reactives(), selected_dataset(), trait_type())
      message("cisTransPlotServer: Calling peak_finder with dataset: ", selected_dataset(), ", trait_type: ", trait_type())

      current_peaks_cache <- if (is.environment(peaks_cache)) peaks_cache else NULL

      # Ensure import_reactives() is called to get the list
      found_peaks <- peak_finder(import_reactives()$file_directory,
        selected_dataset(),
        trait_type = trait_type(),
        cache_env = current_peaks_cache,
        use_cache = TRUE
      )
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

      message("cisTransPlotServer: Preparing plot_data. Initial df from peaks_data() has ", if (is.null(df)) "NULL" else nrow(df), " rows.") # DEBUG
      if (is.null(df) || nrow(df) == 0) {
        message("cisTransPlotServer: plot_data - peaks_data() is NULL or empty, returning NULL.") # DEBUG
        return(NULL)
      }

      req_cols <- c("qtl_chr", "qtl_pos", "cis", "gene_chr", "gene_start", "qtl_lod")
      message("cisTransPlotServer: Required columns for plot: ", paste(req_cols, collapse = ", ")) # DEBUG
      message("cisTransPlotServer: Columns in df from peaks_data(): ", paste(colnames(df), collapse = ", ")) # DEBUG

      if (!all(req_cols %in% colnames(df))) {
        missing_cols <- req_cols[!req_cols %in% colnames(df)]
        message("cisTransPlotServer: plot_data - Missing required columns: ", paste(missing_cols, collapse = ", "), ". Returning NULL.") # DEBUG
        return(NULL)
      }

      df_filtered <- dplyr::filter(df, !is.na(qtl_chr) & !is.na(qtl_pos) & !is.na(cis) & !is.na(gene_chr) & !is.na(gene_start) & !is.na(qtl_lod))
      message("cisTransPlotServer: plot_data - After filtering NAs, df_filtered has ", nrow(df_filtered), " rows.") # DEBUG

      if (nrow(df_filtered) == 0) {
        message("cisTransPlotServer: plot_data - After NA filtering, data frame is empty. Returning NULL.") # DEBUG
        return(NULL)
      }

      # FILTER BY LOD THRESHOLD - Only show points at or above the threshold
      message(paste0("cisTransPlotServer: Applying LOD threshold filter. Before: ", nrow(df_filtered), " rows, threshold: ", lod_threshold))
      df_filtered <- dplyr::filter(df_filtered, qtl_lod >= lod_threshold)
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

      # Ensure 'cis' column is character "TRUE" or "FALSE" for scale_color_manual
      if (is.logical(plot_data_filtered_dt$cis)) {
        plot_data_filtered_dt[, cis_char := ifelse(cis, "TRUE", "FALSE")]
      } else {
        plot_data_filtered_dt[, cis_char := as.character(cis)] # Assume it's already proper if not logical
      }

      plot_data_filtered_dt[, hover_text := paste(
        "Symbol:", gene_symbol,
        "<br>ID:", gene_id,
        "<br>Gene Chr:", gene_chr,
        "<br>Gene Pos (Mbp):", round(gene_start / 1e6, 2),
        "<br>QTL Chr:", qtl_chr,
        "<br>QTL Pos (Mbp):", round(qtl_pos / 1e6, 2),
        "<br>LOD:", round(qtl_lod, 2),
        "<br>QTL Type:", ifelse(cis_char == "TRUE", "Cis", "Trans"),
        "<br>Dataset:", current_selected_dataset_name
      )]

      cis.colors <- c("FALSE" = "#E41A1C", "TRUE" = "blue") # From your old code

      p <- ggplot(
        plot_data_filtered_dt,
        aes(
          x = qtl_pos, # Using raw qtl_pos as in old example
          y = gene_start, # Using raw gene_start as in old example
          color = cis_char, # Use the character version of 'cis'
          key = gene_id,
          text = hover_text,
          customdata = gene_symbol
        )
      ) +
        geom_point(alpha = 0.5, size = 1.2) + # Old point aesthetics
        scale_color_manual(
          values = cis.colors,
          labels = c("FALSE" = "Trans", "TRUE" = "Cis"),
          name = "QTL Type", # Legend title
          drop = FALSE
        ) +
        labs(
          title = paste("Cis/Trans Plot for", current_selected_dataset_name, "(LOD ≥", main_par()$LOD_thr(), ")"),
          x = NULL, # Old axis labs
          y = NULL # Old axis labs
        ) +
        theme_minimal(base_size = 11) +
        theme( # Theme elements from your old code
          panel.grid = element_blank(),
          panel.border = element_rect(color = "grey70", fill = NA),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.spacing = unit(0.1, "lines"),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "top" # Or "right" or remove for default
        )

      fig <- plotly::ggplotly(p, tooltip = "text", source = ns("cistrans_plotly")) %>%
        plotly::layout(dragmode = "pan", hovermode = "closest") %>%
        plotly::config(displaylogo = FALSE, modeBarButtonsToRemove = c("select2d", "lasso2d", "hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines", "sendDataToCloud"))

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
