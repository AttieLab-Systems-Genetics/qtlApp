#' Plot the Alleles
#'
#' @param peak object with data to plot (long format from pivot_peaks)
#' @param colors character vector of colors for the plot
#' @param trait_name character string for the plot subtitle
#' @param peak_marker character string for the plot title
#'
#' @importFrom ggplot2 aes element_blank element_line element_text geom_hline 
#'             geom_point ggplot labs scale_color_manual theme theme_minimal ggtitle
#' @importFrom rlang .data
#' @export
ggplot_alleles <- function(peak_long, colors = NULL, trait_name = NULL, peak_marker = NULL) {
  if(is.null(peak_long) || !nrow(peak_long)) {
    return(ggplot2::ggplot() + 
             ggplot2::theme_void() + 
             ggplot2::labs(title = "No allele effects data available for this peak.") + 
             ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  
  # Define default colors if not provided
  if (is.null(colors)) {
      colors <- c(
          "AJ" = "#F0E442", "B6" = "#000000", "129" = "#E69F00",
          "NOD" = "#56B4E9", "NZO" = "#0072B2", "CAST" = "#009E73",
          "PWK" = "#D55E00", "WSB" = "#CC79A7"
      )
  }
  
  # Ensure the variable column is a factor with levels in the desired order
  ordered_strains <- intersect(names(colors), unique(peak_long$variable))
  peak_long$variable <- factor(peak_long$variable, levels = ordered_strains)
  
  plot_title <- if(!is.null(peak_marker)) paste0("Strain Effects at ", peak_marker) else "Strain Effects"
  plot_subtitle <- if(!is.null(trait_name)) paste0("Trait: ", trait_name) else NULL

  p <- ggplot2::ggplot(data = peak_long, ggplot2::aes(x = variable, y = value, color = variable)) +
    ggplot2::geom_point(size = 4, alpha = 0.8) +
    ggplot2::scale_color_manual(values = colors, breaks = ordered_strains, name = "Strain") + # Use breaks to ensure order
    ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 12, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size=10), # Rotate labels
        axis.title.x = ggplot2::element_blank(), # Remove x-axis title
        axis.title.y = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 11)
    ) +
    ggplot2::labs(y = "Founder Allele Effect", title = plot_title, subtitle = plot_subtitle)
    
    return(p)
} 