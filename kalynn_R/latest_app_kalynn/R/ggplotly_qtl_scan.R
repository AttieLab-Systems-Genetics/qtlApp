#' Convert ggplot QTL Scan to Plotly
#' 
#' Takes a ggplot object representing a QTL scan and converts it to an interactive Plotly object,
#' adding peak information and configuring layout.
#' 
#' @param scan_plot ggplot object created by ggplot_qtl_scan or similar.
#' @param plot_data The data frame used to generate the ggplot object.
#' @param peak_table Data frame with peak information for the selected trait.
#' @param official_trait_symbol The official trait symbol for the title.
#' @param selected_chr The selected chromosome view ("All" or specific chr).
#' @param plot_width Width of the plot in pixels.
#' @param plot_height Height of the plot in pixels.
#' @param source Source ID for plotly events.
#' 
#' @importFrom plotly ggplotly layout config event_register
#' @importFrom ggplot2 geom_point aes
#' @importFrom dplyr filter arrange slice
#' @importFrom rlang .data
#' @export
ggplotly_qtl_scan <- function(scan_plot, plot_data, peak_table,
                              official_trait_symbol,
                              selected_chr = "All",
                              plot_width = 1000, plot_height = 600,
                              source = "scanly_plot") {
  
  if (is.null(scan_plot) || is.null(plot_data)) return(NULL)
  
  # Determine x-axis variable based on view
  xvar <- if (selected_chr == "All") "BPcum" else "position"

  # Find the highest peak in the provided peak_table
  highest_peak_info <- NULL
  if (!is.null(peak_table) && nrow(peak_table) > 0 && "lod" %in% colnames(peak_table)) {
    highest_peak_info <- peak_table %>%
      dplyr::arrange(dplyr::desc(.data$lod)) %>%
      dplyr::slice(1)
  }

  # Add marker for the highest peak to the ggplot object
  if (!is.null(highest_peak_info) && nrow(highest_peak_info) > 0) {
    peak_point <- plot_data %>%
      dplyr::filter(.data$markers == highest_peak_info$marker)
    
    if (nrow(peak_point) > 0) {
      scan_plot <- scan_plot +
        ggplot2::geom_point(data = peak_point,
                            ggplot2::aes(x = .data[[xvar]], y = .data$LOD),
                            color = "red",
                            size = 3,
                            shape = 20) # Use filled circle shape
    }
  }
  
  # Create title and subtitle
  trait_text <- paste0("<b style='font-size: 24px;'>", 
                       if(!is.null(official_trait_symbol)) official_trait_symbol else "Trait", 
                       "</b>")
  
  subtitle <- if (!is.null(highest_peak_info) && nrow(highest_peak_info) > 0) {
    chr_label <- chr_XYM(highest_peak_info$chr)
    paste0(
      "<span style='font-size: 16px;'>",
      "<b>Peak Marker:</b> ", highest_peak_info$marker,
      " (Chr", chr_label, ":", round(highest_peak_info$pos, 2), " Mb) | ",
      "<b>LOD:</b> ", round(highest_peak_info$lod, 2),
      "</span>"
    )
  } else {
    "<span style='font-size: 16px; color: #7f8c8d;'>No significant peaks found</span>"
  }

  # Convert ggplot to plotly
  plt <- plotly::ggplotly(scan_plot, source = source, 
                          width = plot_width, height = plot_height, 
                          tooltip = c("x", "y")) %>% # Keep tooltip simple initially
    # Add layout configuration
    plotly::layout(
      title = list(
        text = paste0(trait_text, '<br>', subtitle),
        font = list(family = "Arial"),
        x = 0,
        xanchor = "left",
        y = 0.98, # Adjust y position slightly
        yanchor = "top",
        pad = list(t = 20, b = 20) # Adjust padding
      ),
      margin = list(t = 100), # Increase top margin for title
      hoverlabel = list(
        bgcolor = "white",
        font = list(family = "Arial", size = 12, color = "#2c3e50"),
        bordercolor = "#95a5a6"
      ),
      hovermode = "closest",
      xaxis = list(title = scan_plot$labels$x), # Get x label from ggplot
      yaxis = list(title = scan_plot$labels$y) # Get y label from ggplot
    ) %>% 
    # Register plotly events
    plotly::event_register('plotly_click') %>%
    plotly::event_register('plotly_doubleclick') %>% 
    # Configure mode bar
    plotly::config(
      displaylogo = FALSE,
      modeBarButtonsToRemove = c(
        "select2d", "lasso2d", "autoScale2d", "hoverClosestCartesian",
        "hoverCompareCartesian", "toggleSpikelines"
      )
    )
    
   
   plt <- plt %>% 
    plotly::style(text = ~paste("Marker: ", markers, 
                                "<br>Chr: ", chr_XYM(chr),
                                "<br>Pos: ", round(position, 2), " Mb",
                                "<br>LOD: ", round(LOD, 2)), 
                    traces = 1) # Assuming the main line is the first trace

  return(plt)
} 