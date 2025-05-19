#' Plot the QTL if doing a new scan
#' 
#' @param scan_table dataframe with data to plot
#' @param LOD_thr threshold for LOD score
#' @param selected_chr chromosome to plot (default is "All")
#' 
#' @importFrom dplyr group_by summarise
#' @importFrom ggplot2 aes aes_string element_blank element_line element_text
#'             expansion geom_hline geom_line ggplot margin scale_color_manual
#'             scale_x_continuous theme theme_bw xlab ylab ylim
#' @importFrom rlang .data
#' @export
ggplot_qtl_scan <- function(scan_table, LOD_thr = NULL, selected_chr = "All") {
  if (!exists("create_modern_theme", mode = "function")) {
    source("R/plot_enhancements.R")
  }
  
  if(is.null(scan_table) || nrow(scan_table) == 0){
    warning("ggplot_qtl_scan: scan_table is NULL or empty. Returning empty plot.")
    # Return a blank ggplot object or use plot_null()
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title="No scan data to display"))
  }
  
  # Ensure scan_table has 'markers', 'order' (numeric chr), 'LOD', 'position', 'BPcum'
  # QTL_plot_visualizer is responsible for providing these with correct types.
  required_cols <- c("markers", "order", "LOD", "position", "BPcum")
  if(!all(required_cols %in% colnames(scan_table))){
      missing_cols <- setdiff(required_cols, colnames(scan_table))
      stop(paste("ggplot_qtl_scan: scan_table is missing required columns:", paste(missing_cols, collapse=", ")))
  }
  
  # Convert numeric 'order' column to character chromosome labels (1-19, X, Y, M)
  # This new column will be used for factor levels and color mapping.
  scan_table$chr_labels <- chr_XYM(scan_table$order)
  
  # Define the desired order of chromosome labels for factors and scales
  # This ensures consistent ordering in legends and on axes.
  all_possible_chr_labels <- c(as.character(1:19), "X", "Y", "M")
  # Filter to include only those present in the current scan_table
  present_chr_labels <- intersect(all_possible_chr_labels, unique(scan_table$chr_labels))
  # Make sure it's a factor with these levels
  scan_table$chr_factor <- factor(scan_table$chr_labels, levels = present_chr_labels)
  
  # Set x variable based on chromosome selection
  xvar <- if (selected_chr == "All") "BPcum" else "position"
  
  # Create axis labels data frame (using numeric 'order' for calculation, character 'chr_labels' for display)
  axisdf <- scan_table %>%
    dplyr::group_by(chr_factor, order) %>% # Group by both factor and numeric order
    dplyr::summarise(center = (max(.data[[xvar]], na.rm = TRUE) + min(.data[[xvar]], na.rm = TRUE))/2, .groups = "drop") %>%
    dplyr::arrange(order) # Arrange by numeric order to ensure correct sequence
  
  # Create color palette named by the factor levels
  num_colors_needed <- length(present_chr_labels)
  plot_colors <- if (exists("create_modern_palette", mode = "function")) {
    rep(create_modern_palette(2), length.out = num_colors_needed)
  } else {
    rep(c("#3498db", "#2c3e50"), length.out = num_colors_needed)
  }
  names(plot_colors) <- present_chr_labels
  
  # Create plot object
  p <- ggplot2::ggplot(scan_table, ggplot2::aes_string(x = xvar, y = "LOD", customdata = "markers")) +
    ggplot2::geom_line(ggplot2::aes(color = chr_factor, group = order), alpha = 0.8, linewidth = 1) +
    ggplot2::scale_color_manual(values = plot_colors, name = "Chromosome", breaks = present_chr_labels) +
    ggplot2::scale_x_continuous(
      label = axisdf$chr_factor, # Use the factor for labels
      breaks = axisdf$center,
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::ylim(0, max(scan_table$LOD, na.rm = TRUE) * 1.25)
  
  # Apply theme
  if (exists("create_modern_theme", mode = "function")) {
    p <- p + create_modern_theme() +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          angle = 0, hjust = 0.5, vjust = 0.5,
          margin = ggplot2::margin(t = 10)
        ),
        plot.margin = ggplot2::margin(b = 40, l = 20, r = 20, t = 20, unit = "pt")
      )
  } else { # Fallback theme
    p <- p + ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "#2c3e50", size = 0.5),
      axis.text = ggplot2::element_text(size = 12, color = "#2c3e50"),
      axis.title = ggplot2::element_text(size = 14, face = "bold", color = "#2c3e50"),
      axis.text.x = ggplot2::element_text(
        angle = 0, hjust = 0.5, vjust = 0.5,
        margin = ggplot2::margin(t = 10)
      ),
      plot.margin = ggplot2::margin(b = 40, l = 20, r = 20, t = 20, unit = "pt")
    )
  }
  
  # Add axis labels
  p <- p + ggplot2::xlab(if(selected_chr == "All") "Chromosome" else paste("Position on Chr", selected_chr, "(Mb)")) +
       ggplot2::ylab("LOD Score")
  
  # Add LOD threshold line
  if(!is.null(LOD_thr) && is.numeric(LOD_thr) && LOD_thr > 0) { # Ensure LOD_thr is valid
    p <- p + ggplot2::geom_hline(yintercept = LOD_thr, color = "#e74c3c", linetype = "dashed", linewidth = 0.8)
  }
  
  return(p)
}
