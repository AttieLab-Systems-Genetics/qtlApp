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
  # Source plot enhancement functions if not already loaded
  if (!exists("create_modern_theme", mode = "function")) {
    source("R/plot_enhancements.R")
  }
  
  # If no data, return a null plot
  if(is.null(scan_table) || !nrow(scan_table)) {
    # Use the null plot function if available
    if (exists("create_null_plot", mode = "function")) {
      source("R/plot_null.R")
      return(create_null_scan_plot())
    } else {
      return(ggplot2::ggplot())
    }
  }
  
  # Set x variable based on chromosome selection
  if (selected_chr == "All") {
    xvar = "BPcum"
  } else {
    xvar = "position"
    # `scan_table` is probably already filtered by chromosome, but make sure.
    scan_table <- dplyr::filter(scan_table, .data$chr == selected_chr)
  }
  
  # Create axis labels
  axisdf = dplyr::group_by(scan_table, .data$order) |>
    dplyr::summarize(center = (max(.data[[xvar]]) + min(.data[[xvar]]))/2)
    
  # Convert chromosome labels
  axisdf$order[axisdf$order == 20] <- "X"
  axisdf$order[axisdf$order == 21] <- "Y"
  axisdf$order[axisdf$order == 22] <- "M"
  axisdf$order <- factor(axisdf$order, levels=c(as.character(1:19), "X", "Y", "M"))
  
  # Create modern color palette - use alternating colors from modern palette
  if (exists("create_modern_palette", mode = "function")) {
    colors <- rep(create_modern_palette(2), 22)
  } else {
    colors <- rep(c("#3498db", "#2c3e50"), 22)  # Default modern colors
  }
  
  # Create plot object with modern styling
  p <- ggplot2::ggplot(scan_table, ggplot2::aes_string(x=xvar, y="LOD")) +
    ggplot2::geom_line(ggplot2::aes(color=as.factor(chr)), alpha=0.8, linewidth=1) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_x_continuous(
      label = axisdf$order, 
      breaks = axisdf$center,
      expand = ggplot2::expansion(mult = c(0.01, 0.01))  # Reduced padding to 1% on each side
    ) +
    ggplot2::ylim(0, max(scan_table$LOD, na.rm=TRUE) * 1.25)
  
  # Apply modern theme if available, otherwise use enhanced classic theme
  if (exists("create_modern_theme", mode = "function")) {
    p <- p + create_modern_theme() +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          angle = 0,           # Keep labels straight for better readability
          hjust = 0.5,         # Center labels
          vjust = 0.5,         # Center vertically
          margin = ggplot2::margin(t = 10)  # Add margin at top of labels
        ),
        # Add more bottom margin for labels
        plot.margin = ggplot2::margin(b = 40, l = 20, r = 20, t = 20, unit = "pt")
      )
  } else {
    # Fallback to enhanced classic theme
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
          angle = 0,           # Keep labels straight for better readability
          hjust = 0.5,         # Center labels
          vjust = 0.5,         # Center vertically
        margin = ggplot2::margin(t = 10)  # Add margin at top of labels
      ),
        # Add more bottom margin to accommodate labels
      plot.margin = ggplot2::margin(b = 40, l = 20, r = 20, t = 20, unit = "pt")
      )
  }
  
  # Add axis labels
  p <- p + ggplot2::xlab(if(selected_chr == "All") "Chromosome" else paste("Position on Chr", selected_chr, "(Mb)")) +
       ggplot2::ylab("LOD Score")
  
  # Add LOD threshold line with modern styling
  if(!is.null(LOD_thr)) {
    p <- p + ggplot2::geom_hline(ggplot2::aes(yintercept = LOD_thr),
      color = "#e74c3c",  # Modern red color
      linetype = "dashed",
      size = 0.8)  # Slightly thicker line
  }
  
  return(p)
}
