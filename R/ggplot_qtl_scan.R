#' Plot the QTL if doing a new scan
#' 
#' @param qtl_plot_obj object with data to plot
#' @param LOD_thr threshold for LOD score
#' @param xvar variable for x-axis (default: "BPcum")
#' 
#' @importFrom dplyr group_by summarise
#' @importFrom ggplot2 aes aes_string element_blank element_line element_text
#'             expansion geom_hline geom_line ggplot margin scale_color_manual
#'             scale_x_continuous theme theme_bw xlab ylab ylim
#' @importFrom rlang .data
#' @export
ggplot_qtl_scan <- function(qtl_plot_obj, LOD_thr = NULL, xvar = "BPcum") {
  # Create axis labels
  axisdf = dplyr::group_by(qtl_plot_obj, .data$order) |>
    dplyr::summarize(center = (max(.data[[xvar]]) + min(.data[[xvar]]))/2)
  # Convert chromosome labels
  axisdf$order[axisdf$order == 20] <- "X"
  axisdf$order[axisdf$order == 21] <- "Y"
  axisdf$order[axisdf$order == 22] <- "M"
  axisdf$order <- factor(axisdf$order, levels=c(as.character(1:19), "X", "Y", "M"))
  # Create plot object.
  p <- ggplot2::ggplot(qtl_plot_obj, ggplot2::aes_string(x=xvar, y="LOD")) +
    ggplot2::geom_line(ggplot2::aes(color=as.factor(chr)), alpha=0.8, linewidth=.5) +
    ggplot2::scale_color_manual(values = rep(c("black", "darkgrey"), 22)) +
    ggplot2::scale_x_continuous(
      label = axisdf$order, 
      breaks = axisdf$center,
      expand = ggplot2::expansion(mult = 0.15)  # Increased padding
    ) +
    ggplot2::ylim(0, max(qtl_plot_obj$LOD, na.rm=TRUE) * 1.25) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      axis.text = ggplot2::element_text(size = 18),
      axis.title = ggplot2::element_text(size = 20),
      axis.text.x = ggplot2::element_text(
        angle = 45,           # Rotate labels 45 degrees
        hjust = 1,           # Adjust horizontal position
        vjust = 1,           # Adjust vertical position
        margin = ggplot2::margin(t = 10)  # Add margin at top of labels
      ),
      # Add more bottom margin to accommodate rotated labels
      plot.margin = ggplot2::margin(b = 40, l = 20, r = 20, t = 20, unit = "pt")
    ) +
    ggplot2::xlab("Chromosome") +
    ggplot2::ylab("LOD")
  if(!is.null(LOD_thr)) {
    p <- p + ggplot2::geom_hline(ggplot2::aes(yintercept = LOD_thr),
      color = "black", linetype = "dashed")
  }
  return(p)
}
