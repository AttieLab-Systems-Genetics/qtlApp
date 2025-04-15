#' Plot the Alleles
#'
#' @param peak object with data to plot
#' @param colors character vector of colors for the plot
#'
#' @importFrom ggplot2 aes element_blank element_line element_text
#'             geom_hline geom_point ggplot labs scale_color_manual theme theme_bw
#' @importFrom rlang .data
#' @export
ggplot_alleles <- function(peak, colors = newClrs) {
  # Define colors
  newClrs <- c(
    "AJ" = "#000000",
    "B6" = "#96989A",
    "129" = "#E69F00",
    "NOD" = "#0072B2",
    "NZO" = "#619BFF",
    "CAST" = "#009E73",
    "PWK" = "#D55E00",
    "WSB" = "#CC79A7")
  ggplot2::ggplot(data = peak) +
    ggplot2::aes(x = .data$marker, y = .data$value, color = .data$variable)) +
    ggplot2::geom_point(size = 10) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 18),
      panel.border = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      axis.text = ggplot2::element_text(size = 18),
      axis.title = ggplot2::element_text(size = 20)) +
    ggplot2::labs(x = "Marker ID", y = "Founder allele effect", color = "Strain") +
    ggplot2::geom_hline(yintercept = 0, color = "black")
}
