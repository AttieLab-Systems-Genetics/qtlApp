# R/cisTransPlotApp.R

#' Cis-Trans Plot Module UI components
#'
#' @param id Module ID
#' @export
cisTransPlotInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::selectInput(ns("dataset_select"), "Select Dataset:", choices = NULL), # Choices populated by server
    shiny::numericInput(ns("min_lod_filter"), "Minimum LOD to Display:", value = 0, min = 0, step = 0.5)
    # Add other controls as needed (e.g., highlight specific genes)
  )
}

#' @rdname cisTransPlotInput
#' @export
cisTransPlotUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::plotOutput(ns("cis_trans_plot_output"), height = "800px") # Plot area
}

#' Cis-Trans Plot Module Server
#'
#' @param id Module ID
#' @param import_reactives Reactive list containing file_directory and annotation_list
#'
#' @importFrom shiny moduleServer NS reactive req selectInput numericInput observe eventReactive
#' @importFrom dplyr filter select left_join mutate case_when rename
#' @importFrom ggplot2 ggplot aes geom_point facet_grid labs theme_minimal geom_abline scale_color_manual
#' @importFrom rlang .data
#' @export
cisTransPlotServer <- function(id, import_reactives) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Populate dataset choices
    shiny::observe({
      imp_data <- import_reactives()
      shiny::req(imp_data, imp_data$file_directory)
      dataset_choices <- unique(imp_data$file_directory$group)
      shiny::updateSelectInput(session, "dataset_select", choices = dataset_choices, selected = dataset_choices[1])
    })

    # Reactive: Get all peaks for the selected dataset
    all_peaks_for_dataset <- shiny::reactive({
      shiny::req(input$dataset_select, import_reactives()$file_directory)
      message(paste("cisTransPlot: Fetching all peaks for dataset:", input$dataset_select))
      # Call peak_finder with selected_trait = NULL to get all peaks
      peaks <- peak_finder(import_reactives()$file_directory, input$dataset_select, selected_trait = NULL)
      message(paste("cisTransPlot: peak_finder returned", nrow(peaks), "peaks for dataset."))
      peaks
    })

    # Reactive: Prepare data for plotting (join peaks with gene annotations)
    plot_data <- shiny::reactive({
      shiny::req(all_peaks_for_dataset(), import_reactives()$annotation_list)
      
      peaks <- all_peaks_for_dataset()
      if (nrow(peaks) == 0) return(NULL)
      
      # Filter by LOD threshold
      peaks <- dplyr::filter(peaks, .data$lod >= input$min_lod_filter)
      if (nrow(peaks) == 0) return(NULL)
      message(paste("cisTransPlot: After LOD filter, have", nrow(peaks), "peaks."))

      annotations <- import_reactives()$annotation_list
      # Assuming 'lodcolumn' in peaks contains the gene/isoform ID
      # and annotation_list has $genes and $isoforms data frames
      # This part needs to be robust to match IDs from peaks$lodcolumn to annotation_list
      
      # For simplicity, let's assume lodcolumn contains gene symbols for now
      # and gene_annotations has 'symbol', 'chr', 'start' (gene_pos_y)
      # We'll need to standardize this based on your actual annotation_list structure
      
      # Attempt to join with genes first, then isoforms if genes are not primary
      # This needs careful handling based on how 'lodcolumn' is structured (e.g. Ensembl IDs vs Symbols)
      # Let's assume 'lodcolumn' can be directly matched to 'symbol' in annotation_list$genes or $isoforms
      # And that your annotation_list$genes and $isoforms have 'symbol', 'chr', 'start' columns
      
      gene_annots <- annotations$genes # Or $isoforms, or a combined version
      if (is.null(gene_annots) || !("symbol" %in% colnames(gene_annots))) {
          message("cisTransPlot: Gene annotations or 'symbol' column missing.")
          return(NULL)
      }
      
      # Standardize column name in peaks to join (e.g. if it's 'lodcolumn' or 'trait')
      # For now, assume peak_finder already created a 'trait' column which holds the gene/isoform identifier
      if (!("trait" %in% colnames(peaks))){
          message("cisTransPlot: 'trait' column (expected to hold gene/isoform ID) not found in peaks data.")
          return(NULL)
      }

      # Rename for clarity before join
      peaks <- dplyr::rename(peaks, qtl_chr = .data$chr, qtl_pos = .data$pos, gene_id_from_peak = .data$trait)

      # Prepare gene annotations (select relevant columns, rename for clarity)
      # This needs to be adapted to your actual annotation_list structure for genes/isoforms
      # Assuming 'symbol' is the join key from peaks$gene_id_from_peak
      # Assuming 'chr' and 'start' exist in gene_annots for gene positions
      if(!("chr" %in% colnames(gene_annots)) || !("start" %in% colnames(gene_annots))){
          message("cisTransPlot: Annotation list missing 'chr' or 'start' columns for gene positions.")
          return(NULL)
      }
      
      gene_annots_simplified <- dplyr::select(gene_annots, 
                                            gene_id_for_join = .data$symbol, # Or gene.id / transcript.id
                                            gene_chr = .data$chr, 
                                            gene_pos_y = .data$start) # Using start position for y-axis

      # Join peaks with gene annotations
      # The join key will be tricky: peaks$gene_id_from_peak needs to match gene_annots_simplified$gene_id_for_join
      # This will depend on whether peaks$trait column contains symbols or Ensembl IDs, etc.
      # For this example, assuming it's symbols.
      
      plot_df <- dplyr::left_join(peaks, gene_annots_simplified, by = c("gene_id_from_peak" = "gene_id_for_join"))
      
      # Remove rows where gene annotation was not found
      plot_df <- dplyr::filter(plot_df, !is.na(.data$gene_chr) & !is.na(.data$gene_pos_y))
      message(paste("cisTransPlot: After joining with annotations, have", nrow(plot_df), "peaks with gene positions."))
      if (nrow(plot_df) == 0) return(NULL)

      # Convert chromosome columns to factor for plotting order (ensure X, Y, M are last)
      chr_levels <- c(as.character(1:19), "X", "Y", "M")
      # Apply chr_XYM element-wise using sapply
      plot_df$qtl_chr_f <- factor(sapply(plot_df$qtl_chr, chr_XYM), levels = chr_levels)
      plot_df$gene_chr_f <- factor(sapply(plot_df$gene_chr, chr_XYM), levels = chr_levels)
      
      # Add a column for cis/trans status
      plot_df <- dplyr::mutate(plot_df, type = dplyr::case_when(
        .data$qtl_chr == .data$gene_chr ~ "cis",
        TRUE ~ "trans"
      ))
      
      plot_df
    })

    # Render the plot
    output$cis_trans_plot_output <- shiny::renderPlot({
      p_data <- plot_data()
      shiny::req(p_data)
      if (nrow(p_data) == 0) {
        # Ensure plot_null is available (it should be if sourced by app.R)
        if (exists("plot_null", mode = "function")) {
            return(plot_null(msg = "No peaks to plot with selected filters."))
        } else {
            return(ggplot2::ggplot() + ggplot2::labs(title="No peaks to plot with selected filters.") + ggplot2::theme_void())
        }
      }

      # Create the plot - NO FACETING
      p <- ggplot2::ggplot(p_data, ggplot2::aes(x = .data$qtl_pos / 1e6, 
                                                y = .data$gene_pos_y / 1e6, 
                                                color = .data$type,
                                                text = paste("QTL Marker:", .data$marker,
                                                             "<br>QTL Chr:", .data$qtl_chr_f, "@", round(.data$qtl_pos/1e6,1), "Mb",
                                                             "<br>Gene:", .data$gene_id_from_peak, 
                                                             "<br>Gene Chr:", .data$gene_chr_f, "@", round(.data$gene_pos_y/1e6,1), "Mb",
                                                             "<br>LOD:", .data$lod))) + 
        ggplot2::geom_point(alpha = 0.6, size = 2) + 
        ggplot2::scale_color_manual(values = c("cis" = "red", "trans" = "blue")) +
        ggplot2::labs(
          title = paste("Cis/Trans eQTL Plot for Dataset:", input$dataset_select),
          x = "QTL Position (Mb on its chromosome)",
          y = "Gene Start Position (Mb on its chromosome)",
          color = "QTL Type"
        ) +
        ggplot2::theme_minimal(base_size = 11) + # Set a base font size
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), 
                       legend.position = "top")
      
      # For interactivity later, you would convert to plotly here:
      # return(plotly::ggplotly(p, tooltip="text"))
      p
    }, res = 96)

  })
} 