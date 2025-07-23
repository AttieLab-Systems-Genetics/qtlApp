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
ggplot_qtl_scan <- function(scan_table, LOD_thr = NULL, selected_chr = "All",
                            overlay_diet_data = NULL, overlay_sex_data = NULL) {
  if (!exists("create_modern_theme", mode = "function")) {
    source("R/plot_enhancements.R")
  }

  if (is.null(scan_table) || nrow(scan_table) == 0) {
    # Return a blank ggplot object if no base data
    return(ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "No scan data to display"))
  }

  xvar <- if (selected_chr == "All") "BPcum" else "position"

  # --- Data Preparation ---
  # Label the base data
  scan_table$type <- "Additive"

  # Combine all data into one data frame for robust plotting
  plot_data <- scan_table

  if (!is.null(overlay_diet_data) && nrow(overlay_diet_data) > 0) {
    overlay_diet_data$type <- "Diet Interactive"
    plot_data <- rbind(plot_data, overlay_diet_data)
    message("ggplot_qtl_scan: Merged DIET overlay data.")
  }

  if (!is.null(overlay_sex_data) && nrow(overlay_sex_data) > 0) {
    overlay_sex_data$type <- "Sex Interactive"
    plot_data <- rbind(plot_data, overlay_sex_data)
    message("ggplot_qtl_scan: Merged SEX overlay data.")
  }

  plot_data$type <- factor(plot_data$type, levels = c("Additive", "Diet Interactive", "Sex Interactive"))

  # --- Color and Theme Setup ---
  # Define colors for each type of scan
  color_map <- c(
    "Additive" = "#3498db", # Blue
    "Diet Interactive" = "#2c3e50", # Dark Blue
    "Sex Interactive" = "#e74c3c" # Light Red
  )

  axisdf <- scan_table %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(center = (max(.data[[xvar]], na.rm = TRUE) + min(.data[[xvar]], na.rm = TRUE)) / 2, .groups = "drop") %>%
    dplyr::arrange(chr)

  axisdf$chr <- chr_XYM(axisdf$chr)

  # --- Plot Construction ---
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[xvar]], y = .data$LOD, color = .data$type, group = interaction(.data$type, .data$chr))) +
    ggplot2::geom_line(aes(linewidth = .data$type, alpha = .data$type)) +

    # Manual scales for aesthetics
    ggplot2::scale_color_manual(values = color_map, name = "Scan Type") +
    ggplot2::scale_linewidth_manual(values = c("Additive" = 0.7, "Diet Interactive" = 1.2, "Sex Interactive" = 1.2), guide = "none") +
    ggplot2::scale_alpha_manual(values = c("Additive" = 0.8, "Diet Interactive" = 0.85, "Sex Interactive" = 0.85), guide = "none") +

    # Axis and theme setup
    ggplot2::scale_x_continuous(label = axisdf$chr, breaks = axisdf$center, expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::ylim(0, max(plot_data$LOD, na.rm = TRUE) * 1.25) +
    create_modern_theme() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = if (selected_chr == "All") "Chromosome" else paste("Position on Chr", selected_chr, "(Mb)"),
      y = "LOD Score"
    )

  # Add LOD threshold line only for the additive scan
  if (!is.null(LOD_thr) && is.numeric(LOD_thr) && LOD_thr > 0) {
    p <- p + ggplot2::geom_hline(yintercept = LOD_thr, color = "#e74c3c", linetype = "dashed", linewidth = 0.8)
  }

  # Conditionally hide legend if only one scan type is present (no overlays)
  if (length(unique(plot_data$type)) <= 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  return(p)
}
