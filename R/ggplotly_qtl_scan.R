#' Plot the QTL if doing a new scan
#' 
#' @param peak_table data frame with peak information
#' @param scan_table data frame with scan data
#' @param selected_chr selected chromosome
#' @param plot_width width of the plot
#' @param plot_height height of the plot
#'
#' @importFrom plotly config ggplotly layout
#' @export
ggplotly_qtl_scan <- function(scan_table, peak_table,
                              selected_chr = "2", plot_width = 900, plot_height = 600) {
  # Create formatted trait text for plot title using the official gene symbol
  trait_text <- paste0("<b style='font-size: 24px;'>", official_gene_symbol(), "</b>")
  # Show the highest LOD peak
  if (nrow(peak_table) > 0) {
    peak_table <- peak_table |>
      dplyr::arrange(desc(lod)) |>
      dplyr::slice(1)
    peak_point <- scan_table |>
      dplyr::filter(markers == peak_table$marker)
    if (nrow(peak_point) > 0) {
      # Add red diamond at the peak
      if (input$selected_chr == "All") {
        xvar = "BPcum"
      } else {
        xvar = "position"
      }
      scan_plot <- ggplot_qtl_scan(scan_table) + 
        ggplot2::geom_point(data = peak_point,
          ggplot2::aes(x = .data[[xvar]], y = .data$LOD),
          color = "red",
          size = 3,
          shape = 20)
    }
  }
  # Create subtitle with peak information
  subtitle <-
    if (nrow(peak_table) > 0) {
      chr_label <- 
        if (peak_table$chr %in% c(20,21,22)) {
          c("X","Y","M")[peak_table$chr-19]
        } else {
          peak_table$chr
        }
      paste0(
        "<span style='font-size: 16px;'>",
        "<b>Peak Marker:</b> ", peak_table$marker,
        " (Chr", chr_label, ":", round(peak_table$pos, 2), " Mb) | ",
        "<b>LOD:</b> ", round(peak_table$lod, 2),
        "</span>"
      )
    } else {
      "<span style='font-size: 16px; color: #7f8c8d;'>No significant peaks</span>"
    }
  # Convert ggplot to plotly with custom dimensions and removed features
  plt <- plotly::ggplotly(scan_plot, source = "scan_plot",
    width = plot_width, height = plot_height, tooltip = c("x", "y", "chr")) |>
    plotly::layout(
      title = list(
        text = paste0(trait_text, '<br>', subtitle),
        font = list(family = "Arial"),
        x = 0,
        xanchor = "left",
        y = 0.95,  
        yanchor = "top",
        pad = list(b = 20)
      ),
      margin = list(t = 80),
      hoverlabel = list(
        bgcolor = "white",
        font = list(family = "Arial", size = 12, color = "#2c3e50"),
        bordercolor = "#95a5a6"
      ),
      hovermode = "closest",
      # Add double click event to reset to full view
      doubleclick = if (selected_chr != "All") TRUE else FALSE
    ) |>
    # Remove unwanted modebar buttons
    plotly::config(
      displaylogo = FALSE,
      modeBarButtonsToRemove = c(
        "select2d", "lasso2d", "autoScale2d", "hoverClosestCartesian",
        "hoverCompareCartesian", "toggleSpikelines"
      )
    )
  return(plt)
}