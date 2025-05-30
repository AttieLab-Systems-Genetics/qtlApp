#' Plot the Alleles
#'
#' @param peak object with data to plot
#' @param colors character vector of colors for the plot
#'
#' @importFrom ggplot2 aes element_blank element_line element_text
#'             geom_hline geom_point ggplot labs scale_color_manual
#'             theme theme_bw theme_void
#' @importFrom rlang .data
#' @export
ggplot_alleles <- function(peak, colors = NULL) {
  # Source plot enhancement functions if not already loaded
  if (!exists("create_modern_theme", mode = "function")) {
    source("R/plot_enhancements.R")
  }
  
  # If no data or input is NULL, return a null plot with message
  # Explicitly check for NULL first
  if(is.null(peak) || !is.data.frame(peak) || !nrow(peak)) {
    # Use the null plot function if available
    if (exists("create_null_plot", mode = "function")) {
      source("R/plot_null.R")
      return(create_null_allele_plot())
    } else {
      return(ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "No allele effects data available for this peak.") + # More informative message
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  }
  
  # Define modern strain colors if not provided
  if (is.null(colors)) {
    colors <- c(
    "AJ" = "#000000",
    "B6" = "#96989A",
    "129" = "#E69F00",
    "NOD" = "#0072B2",
    "NZO" = "#619BFF",
    "CAST" = "#009E73",
    "PWK" = "#D55E00",
    "WSB" = "#CC79A7")
  }
  
  # Create plot with modern styling
  p <- ggplot2::ggplot(data = peak) +
    ggplot2::aes(x = .data$marker, y = .data$value, color = .data$variable) +
    ggplot2::geom_point(size = 5, alpha = 0.8) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::geom_hline(yintercept = 0, color = "#7f8c8d", linetype = "dashed", size = 0.5)
  
  # Apply modern theme if available, otherwise use enhanced classic theme
  if (exists("create_modern_theme", mode = "function")) {
    p <- p + create_modern_theme() +
      ggplot2::theme(
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 12, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)
      )
  } else {
    # Fallback to enhanced classic theme
    p <- p + ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "#2c3e50", size = 0.5),
        axis.text = ggplot2::element_text(size = 12, color = "#2c3e50"),
        axis.title = ggplot2::element_text(size = 14, face = "bold", color = "#2c3e50"),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 12, face = "bold", color = "#2c3e50"),
        legend.text = ggplot2::element_text(size = 11, color = "#2c3e50"),
        plot.title = ggplot2::element_text(size = 16, face = "bold", color = "#2c3e50"),
        plot.subtitle = ggplot2::element_text(size = 12, color = "#7f8c8d")
      )
  }
  
  # Add labels
  p <- p + ggplot2::labs(
    x = "Marker ID", 
    y = "Founder Allele Effect", 
    color = "Strain",
    title = paste0("Strain Effects at ", peak$marker[1]),
    subtitle = if("trait" %in% colnames(peak)) paste0("Trait: ", peak$trait[1]) else NULL
  )
  
  return(p)
}
