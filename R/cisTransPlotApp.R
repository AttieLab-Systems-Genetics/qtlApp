# R/cisTransPlotApp.R

#' Cis-Trans Plot Module UI components
#'
#' @param id Module ID
#' @export
cisTransPlotInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::selectInput(ns("dataset_select"), "Select Dataset:", choices = NULL), # Choices populated by server
    shiny::numericInput(ns("min_lod_filter"), "Minimum LOD to Display:", value = 7.5, min = 0, step = 0.5)
    # Add other controls as needed (e.g., highlight specific genes)
  )
}

#' @rdname cisTransPlotInput
#' @export
cisTransPlotUI <- function(id) {
  ns <- shiny::NS(id)
  shinycssloaders::withSpinner(plotly::plotlyOutput(ns("cis_trans_plot_output"), height = "800px"))
}

#' Cis-Trans Plot Module Server
#'
#' @param id Module ID
#' @param import_reactives Reactive list containing file_directory and annotation_list
#' @param peaks_cache Environment for caching peak finder results.
#'
#' @importFrom shiny moduleServer NS reactive req selectInput numericInput observe updateSelectInput
#' @importFrom dplyr filter select left_join mutate case_when rename distinct
#' @importFrom ggplot2 ggplot aes geom_point facet_grid labs theme_minimal geom_abline scale_color_manual element_text theme
#' @importFrom plotly plotlyOutput renderPlotly ggplotly event_data layout config
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @export
cisTransPlotServer <- function(id, import_reactives, peaks_cache) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Populate dataset choices
    shiny::observe({
      imp_data <- import_reactives()
      shiny::req(imp_data, imp_data$file_directory)
      dataset_choices <- unique(imp_data$file_directory$group)
      shiny::updateSelectInput(session, "dataset_select", choices = dataset_choices, selected = if(!is.null(dataset_choices) && length(dataset_choices)>0) dataset_choices[1] else NULL)
    })

    # Reactive: Get all peaks for the selected dataset
    all_peaks_for_dataset <- shiny::reactive({
      shiny::req(input$dataset_select, import_reactives()$file_directory)
      message(paste("cisTransPlot: Fetching all peaks for dataset:", input$dataset_select))
      peaks <- peak_finder(import_reactives()$file_directory, input$dataset_select, 
                           selected_trait = NULL, cache_env = peaks_cache)
      message(paste("cisTransPlot: peak_finder returned", nrow(peaks), "peaks for dataset."))
      peaks
    })

    # Reactive: Prepare data for plotting (join peaks with gene annotations)
    plot_data <- shiny::reactive({
      shiny::req(all_peaks_for_dataset(), import_reactives()$annotation_list)
      
      peaks_orig <- all_peaks_for_dataset()
      if (nrow(peaks_orig) == 0) return(data.frame())
      
      if (!("lod" %in% colnames(peaks_orig))){
          message("cisTransPlot: 'lod' column not found in peaks data.")
          return(data.frame())
      }
      peaks_lod_filtered <- dplyr::filter(peaks_orig, .data$lod >= input$min_lod_filter)
      if (nrow(peaks_lod_filtered) == 0) return(data.frame())
      message(paste("cisTransPlot: After LOD filter, have", nrow(peaks_lod_filtered), "peaks."))

      annotations <- import_reactives()$annotation_list
      gene_annots <- annotations$genes
      if (is.null(gene_annots) || !all(c("symbol", "chr", "start") %in% colnames(gene_annots))) {
          message("cisTransPlot: Gene annotations or required columns (symbol, chr, start) missing.")
          return(data.frame())
      }
      
      if (!("trait" %in% colnames(peaks_lod_filtered))){
          message("cisTransPlot: 'trait' column (expected gene/isoform ID) not found in peaks data.")
          return(data.frame())
      }

      # Ensure A-H columns are present for tooltips, rename other necessary columns
      peaks_renamed <- dplyr::rename(peaks_lod_filtered, 
                                   qtl_chr_raw = .data$chr, 
                                   qtl_pos = .data$pos, 
                                   gene_id_from_peak = .data$trait,
                                   qtl_marker = .data$marker) 
      # Keep A-H if they exist
      allele_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
      existing_allele_cols <- intersect(allele_cols, colnames(peaks_renamed))
      
      gene_annots_simplified <- dplyr::select(gene_annots, 
                                            gene_id_for_join = .data$symbol, 
                                            gene_chr_raw = .data$chr, 
                                            gene_pos_y = .data$start)

      plot_df <- dplyr::left_join(peaks_renamed, gene_annots_simplified, 
                                  by = c("gene_id_from_peak" = "gene_id_for_join"),
                                  relationship = "many-to-many")
      
      plot_df <- dplyr::filter(plot_df, !is.na(.data$gene_chr_raw) & !is.na(.data$gene_pos_y) & !is.na(.data$qtl_chr_raw))
      message(paste("cisTransPlot: After joining, have", nrow(plot_df), "peaks with gene positions."))
      if (nrow(plot_df) == 0) return(data.frame())
      
      # Add a column for cis/trans status (using original chr values)
      plot_df <- dplyr::mutate(plot_df, type = dplyr::case_when(
        .data$qtl_chr_raw == .data$gene_chr_raw ~ "cis",
        TRUE ~ "trans"
      ))
      
      message(paste("cisTransPlot: Final plot_df for SINGLE PANEL plot has", nrow(plot_df), "rows."))
      if (nrow(plot_df) == 0) return(data.frame())
      
      # Select necessary columns plus A-H
      # Use original chr names (_raw) for tooltips if needed
      cols_to_keep <- c("qtl_marker", "qtl_chr_raw", "qtl_pos", "gene_id_from_peak", 
                        "gene_chr_raw", "gene_pos_y", "lod", "type", existing_allele_cols)
      plot_df <- dplyr::select(plot_df, dplyr::all_of(intersect(cols_to_keep, colnames(plot_df))))
      
      # Construct tooltip text dynamically
      tooltip_text_parts <- list(
        paste0("QTL Marker: ", plot_df$qtl_marker),
        paste0("QTL: Chr", plot_df$qtl_chr_raw, "@", round(plot_df$qtl_pos/1e6,1), " Mb"),
        paste0("Gene: ", plot_df$gene_id_from_peak, " (Chr", plot_df$gene_chr_raw, "@", round(plot_df$gene_pos_y/1e6,1), " Mb)"),
        paste0("LOD: ", plot_df$lod)
      )
      allele_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
      for(col in allele_cols){
          if(col %in% colnames(plot_df)){
              tooltip_text_parts[[length(tooltip_text_parts) + 1]] <- paste0(col, ": ", round(plot_df[[col]], 3))
          }
      }
      plot_df$tooltip <- sapply(1:nrow(plot_df), function(i) paste0(lapply(tooltip_text_parts, `[[`, i), collapse = "<br>"))
      
      plot_df
    })

    # Render the plot
    output$cis_trans_plot_output <- plotly::renderPlotly({
      p_data <- plot_data()
      shiny::req(p_data, nrow(p_data) > 0)

      # Create the ggplot object - NO FACETING
      p <- ggplot2::ggplot(p_data, ggplot2::aes(x = .data$qtl_pos / 1e6, 
                                                y = .data$gene_pos_y / 1e6, 
                                                color = .data$type,
                                                text = .data$tooltip)) + 
        ggplot2::geom_point(alpha = 0.6, size = 2) + 
        ggplot2::scale_color_manual(values = c("cis" = "red", "trans" = "blue")) +
        ggplot2::scale_x_continuous(name = "QTL Position (Mb)", expand = ggplot2::expansion(mult = 0.05)) +
        ggplot2::scale_y_continuous(name = "Gene Start Position (Mb)", expand = ggplot2::expansion(mult = 0.05)) +
        ggplot2::labs(
          title = paste("Cis/Trans Plot for Dataset:", input$dataset_select),
          subtitle = paste("LOD >=", input$min_lod_filter),
          color = "QTL Type"
        ) +
        ggplot2::theme_minimal(base_size = 11) + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), 
                       legend.position = "top")
      
      # Convert to plotly and configure
      plotly::ggplotly(p, tooltip = "text", source = ns("cis_trans_plotly")) %>%
        plotly::layout(
          dragmode = "zoom",
          hovermode = "closest",
          xaxis = list(title = "QTL Position (Mb)"), # Ensure axis titles are set in layout too
          yaxis = list(title = "Gene Start Position (Mb)")
        ) %>%
        plotly::config(
          displaylogo = FALSE,
          modeBarButtonsToRemove = c("select2d", "lasso2d", "hoverClosestCartesian", 
                                     "hoverCompareCartesian", "toggleSpikelines", "sendDataToCloud")
        )
    })

  })
} 