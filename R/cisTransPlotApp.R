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

    shiny::observe({
      imp_data <- import_reactives()
      shiny::req(imp_data, imp_data$file_directory)
      dataset_choices <- unique(imp_data$file_directory$group)
      shiny::updateSelectInput(session, "dataset_select", choices = dataset_choices, selected = if(!is.null(dataset_choices) && length(dataset_choices)>0) dataset_choices[1] else NULL)
    })

    current_trait_type <- shiny::reactive({
        shiny::req(input$dataset_select, import_reactives())
        get_trait_type(import_reactives(), input$dataset_select)
    })

    all_peaks_for_dataset <- shiny::reactive({
      shiny::req(input$dataset_select, import_reactives()$file_directory, import_reactives()$annotation_list, current_trait_type())
      peaks <- peak_finder(import_reactives()$file_directory, 
                           input$dataset_select, 
                           selected_trait = NULL, 
                           trait_type = current_trait_type(),
                           cache_env = peaks_cache, 
                           use_cache = TRUE)
      peaks
    })

    plot_data <- shiny::reactive({
      shiny::req(all_peaks_for_dataset(), import_reactives()$annotation_list, current_trait_type())
      
      trait_type_val <- current_trait_type()
      if (!(trait_type_val %in% c("genes", "isoforms"))) {
        return(data.frame())
      }

      peaks_orig <- all_peaks_for_dataset()
      if (is.null(peaks_orig) || nrow(peaks_orig) == 0) return(data.frame())
      
      if (!("lod" %in% colnames(peaks_orig))){
          warning("cisTransPlot: 'lod' column not found in peaks data.")
          return(data.frame())
      }
      peaks_lod_filtered <- dplyr::filter(peaks_orig, .data$lod >= input$min_lod_filter)
      if (nrow(peaks_lod_filtered) == 0) return(data.frame())

      annotations <- import_reactives()$annotation_list
      gene_annots <- annotations$genes
      if (is.null(gene_annots)) {
          warning("cisTransPlot: Gene annotations (annotations$genes) are NULL.")
          return(data.frame())
      }
      if (!all(c("gene.id", "symbol", "chr", "start") %in% colnames(gene_annots))) {
          warning("cisTransPlot: Gene annotations missing one of required columns: gene.id, symbol, chr, start. Actual colnames: ", paste(colnames(gene_annots), collapse=", "))
          return(data.frame())
      }
      
      required_peak_cols <- c("chr", "pos", "lod", "marker", "gene_id", "gene_symbol")
      if (!all(required_peak_cols %in% colnames(peaks_lod_filtered))){
          missing_p_cols <- setdiff(required_peak_cols, colnames(peaks_lod_filtered))
          warning(paste("cisTransPlot: Peaks data missing required columns for cis/trans plot:", paste(missing_p_cols, collapse=", ")))
          return(data.frame())
      }

      gene_annots_simplified <- dplyr::select(gene_annots, 
                                            gene_id_for_join = .data$`gene.id`,
                                            gene_symbol_annot = .data$symbol,
                                            gene_chr_raw = .data$chr, 
                                            gene_pos_y = .data$start)

      plot_df <- dplyr::left_join(peaks_lod_filtered, gene_annots_simplified, 
                                  by = c("gene_id" = "gene_id_for_join"))
      
      plot_df <- dplyr::filter(plot_df, !is.na(.data$gene_chr_raw) & !is.na(.data$gene_pos_y) & !is.na(.data$chr))
      if (nrow(plot_df) == 0) return(data.frame())
      
      plot_df <- dplyr::mutate(plot_df, 
                               qtl_chr_raw_char = as.character(.data$chr),
                               gene_chr_raw_char = as.character(.data$gene_chr_raw))
      
      plot_df <- dplyr::mutate(plot_df, type = dplyr::case_when(
        .data$qtl_chr_raw_char == .data$gene_chr_raw_char ~ "cis",
        TRUE ~ "trans"
      ))
      
      allele_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
      existing_allele_cols <- intersect(allele_cols, colnames(plot_df))
      
      plot_df$tooltip <- sapply(1:nrow(plot_df), function(i) {
        parts <- list(
            paste0("QTL Marker: ", plot_df$marker[i]),
            paste0("QTL: Chr", plot_df$qtl_chr_raw_char[i], "@", round(plot_df$pos[i]/1e6,1), " Mb"),
            paste0("Gene: ", plot_df$gene_symbol_annot[i], " (Chr", plot_df$gene_chr_raw_char[i], "@", round(plot_df$gene_pos_y[i]/1e6,1), " Mb)"),
            paste0("LOD: ", round(plot_df$lod[i], 2)),
            paste0("Type: ", plot_df$type[i])
        )
        for(al_col in existing_allele_cols){
            parts[[length(parts) + 1]] <- paste0(al_col, ": ", round(plot_df[[al_col]][i], 3))
        }
        paste0(parts, collapse = "<br>")
      })
      
      # Ensure factor levels for chromosomes for consistent faceting
      all_chrs_qtl <- unique(plot_df$qtl_chr_raw_char)
      all_chrs_gene <- unique(plot_df$gene_chr_raw_char)
      all_chrs_combined <- unique(c(all_chrs_qtl, all_chrs_gene)) # Get unique chr values present in data
      
      # Define the complete, desired order for all possible chromosomes
      standard_chr_order <- c(as.character(1:19), "X", "Y", "M")
      
      # Convert numeric chromosome names in our data to X, Y, M for matching with standard order
      data_chr_labels <- chr_XYM(all_chrs_combined)
      
      # Determine the factor levels to use: only those present in the data, but in the standard order
      actual_levels_in_standard_order <- intersect(standard_chr_order, data_chr_labels)
      # If some chr_labels from data aren't in standard_chr_order (e.g. non-canonical), append them at the end
      # This ensures they are still part of the levels and don't get turned into NA by factor()
      custom_chrs_present <- setdiff(data_chr_labels, standard_chr_order)
      final_factor_levels <- c(actual_levels_in_standard_order, custom_chrs_present)
      
      plot_df$qtl_chr_fac <- factor(chr_XYM(plot_df$qtl_chr_raw_char), levels = final_factor_levels)
      plot_df$gene_chr_fac <- factor(chr_XYM(plot_df$gene_chr_raw_char), levels = final_factor_levels)

      # Filter out rows where faceting variables became NA (if a chr in data wasn't in final_factor_levels - should not happen with above logic)
      plot_df <- dplyr::filter(plot_df, !is.na(qtl_chr_fac) & !is.na(gene_chr_fac))
      
      cols_for_plot <- c("pos", "gene_pos_y", "type", "tooltip", "qtl_chr_fac", "gene_chr_fac", "lod")
      # Use intersect to select only existing columns from cols_for_plot
      actual_cols_for_plot <- intersect(cols_for_plot, colnames(plot_df))
      plot_df_final <- plot_df[, actual_cols_for_plot, drop = FALSE]
      
      return(plot_df_final)
    })

    output$cis_trans_plot_output <- plotly::renderPlotly({
      p_data <- plot_data()
      if (is.null(p_data) || nrow(p_data) == 0) {
        return(plotly::plotly_empty(type = "scatter", mode = "text") %>%
                 plotly::layout(title = "Cis/trans plot is only available for gene/isoform datasets with annotations, or no data meets LOD filter.")) # Updated message
      }

      cis_trans_colors <- c("cis" = "red", "trans" = "blue")
      
      p <- ggplot2::ggplot(p_data, ggplot2::aes(x = .data$pos / 1e6,
                                                y = .data$gene_pos_y / 1e6,
                                                color = .data$type,
                                                text = .data$tooltip)) + 
        ggplot2::geom_point(alpha = 0.6, size = 1.2) + # Slightly smaller points for dense plots
        ggplot2::scale_color_manual(values = cis_trans_colors) +
        # Remove space = "free" to get uniformly sized panels
        ggplot2::facet_grid(gene_chr_fac ~ qtl_chr_fac, scales = "free", drop = TRUE) + 
        ggplot2::labs(
          x = "QTL Position (Mb)", 
          y = "Gene Start Position (Mb)",
          title = paste("Cis/Trans Plot for Dataset:", input$dataset_select),
          subtitle = paste("LOD >=", input$min_lod_filter),
          color = "Type"
        ) +
        # Suggest max 5 breaks for y-axis in each panel, adjust as needed
        ggplot2::scale_y_continuous(n.breaks = 5) +
        ggplot2::scale_x_continuous(n.breaks = 3) + # Suggest fewer breaks for x-axis too
        ggplot2::theme_bw(base_size = 9) + # Slightly smaller base size for text
        ggplot2::theme( 
          panel.spacing = ggplot2::unit(0.1, "lines"),
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = rel(0.8)), # Relative size
          axis.text.y = ggplot2::element_text(size = rel(0.8)), # Relative size
          # Adjust strip text size for facet labels
          strip.text = ggplot2::element_text(size = rel(0.85), margin = ggplot2::margin(t=1, b=1, unit="pt")),
          strip.background = ggplot2::element_rect(fill="grey92", color = "grey70", linewidth = 0.5),
          legend.position = "top",
          plot.title = ggplot2::element_text(size = rel(1.2)),
          plot.subtitle = ggplot2::element_text(size = rel(1.0))
        )
      
      # Convert to plotly
      plotly::ggplotly(p, tooltip = "text", source = ns("cis_trans_plotly")) %>% 
        plotly::layout(
          dragmode = "pan",
          hovermode = "closest",
          margin = list(l = 70, r = 20, b = 70, t = 90) 
        ) %>% 
        plotly::config(
          displaylogo = FALSE,
          modeBarButtonsToRemove = c("select2d", "lasso2d", "hoverClosestCartesian", 
                                     "hoverCompareCartesian", "toggleSpikelines", "sendDataToCloud")
        )
    })

  })
} 