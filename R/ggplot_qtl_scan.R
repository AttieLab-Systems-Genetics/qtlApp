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
                            overlay_diet_data = NULL, overlay_sex_data = NULL,
                            overlay_sex_diet_data = NULL,
                            show_thresholds = TRUE,
                            thresholds_by_type = c(
                              "Additive" = 7.5,
                              "Diet Interactive" = 10.5,
                              "Sex Interactive" = 10.5,
                              "Sex x Diet Interactive" = 15.7
                            )) {
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
  # Label the base data (respect pre-labeled type if present)
  if (!"type" %in% colnames(scan_table)) {
    scan_table$type <- "Additive"
  }

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

  if (!is.null(overlay_sex_diet_data) && nrow(overlay_sex_diet_data) > 0) {
    overlay_sex_diet_data$type <- "Sex x Diet Interactive"
    plot_data <- rbind(plot_data, overlay_sex_diet_data)
    message("ggplot_qtl_scan: Merged SEXxDIET overlay data.")
  }

  plot_data$type <- factor(plot_data$type, levels = c("Additive", "Diet Interactive", "Sex Interactive", "Sex x Diet Interactive"))

  # Build formatted hover text: LOD and Position as Chr#:##.##Mb
  plot_data$chr_char <- chr_XYM(plot_data$chr)
  plot_data$hover_text <- paste0(
    "LOD: ", round(plot_data$LOD, 2),
    "<br>Chr", plot_data$chr_char, ":", round(plot_data$position, 1), "Mb"
  )

  # --- Color and Theme Setup ---
  # Define colors for each type of scan
  color_map <- c(
    "Additive" = "#3498db", # Blue
    "Diet Interactive" = "#2c3e50", # Dark Blue
    "Sex Interactive" = "#e74c3c", # Light Red
    "Sex x Diet Interactive" = "#f6ae2d" # Dark burnt yellow
  )

  axisdf <- scan_table %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(center = (max(.data[[xvar]], na.rm = TRUE) + min(.data[[xvar]], na.rm = TRUE)) / 2, .groups = "drop") %>%
    dplyr::arrange(chr)

  axisdf$chr <- chr_XYM(axisdf$chr)

  # --- Plot Construction ---
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[xvar]], y = .data$LOD, color = .data$type, group = interaction(.data$type, .data$chr), text = .data$hover_text)) +
    ggplot2::geom_line(aes(linewidth = .data$type, alpha = .data$type)) +

    # Manual scales for aesthetics
    ggplot2::scale_color_manual(values = color_map, name = "Scan Type") +
    ggplot2::scale_linewidth_manual(values = c("Additive" = 0.7, "Diet Interactive" = 0.7, "Sex Interactive" = 0.7, "Sex x Diet Interactive" = 0.7), guide = "none") +
    ggplot2::scale_alpha_manual(values = c("Additive" = 0.8, "Diet Interactive" = 0.85, "Sex Interactive" = 0.85, "Sex x Diet Interactive" = 0.85), guide = "none") +

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

  # Add series-specific threshold lines if requested
  if (isTRUE(show_thresholds) && length(thresholds_by_type) > 0) {
    present_types <- levels(plot_data$type)[levels(plot_data$type) %in% unique(as.character(plot_data$type))]
    xmin <- suppressWarnings(min(plot_data[[xvar]], na.rm = TRUE))
    xmax <- suppressWarnings(max(plot_data[[xvar]], na.rm = TRUE))

    # Always draw additive if present
    if ("Additive" %in% present_types && !is.na(thresholds_by_type["Additive"])) {
      p <- p + ggplot2::geom_hline(yintercept = thresholds_by_type["Additive"], color = color_map["Additive"], linetype = "dashed", linewidth = 0.6)
    }

    # Handle Diet and Sex interactive thresholds
    diet_present <- "Diet Interactive" %in% present_types && !is.na(thresholds_by_type["Diet Interactive"]) && !is.null(color_map["Diet Interactive"]) && !is.infinite(xmin) && !is.infinite(xmax)
    sex_present <- "Sex Interactive" %in% present_types && !is.na(thresholds_by_type["Sex Interactive"]) && !is.null(color_map["Sex Interactive"]) && !is.infinite(xmin) && !is.infinite(xmax)

    same_threshold <- diet_present && sex_present && identical(as.numeric(thresholds_by_type["Diet Interactive"]), as.numeric(thresholds_by_type["Sex Interactive"]))

    if (same_threshold) {
      # Build alternating color segments along the full x-range
      n_segments <- 120
      xs <- seq(xmin, xmax, length.out = n_segments + 1)
      seg_df <- data.frame(
        x0 = xs[-length(xs)],
        x1 = xs[-1],
        idx = seq_len(n_segments)
      )
      # Diet draws odd segments, Sex draws even (or vice versa)
      seg_diet <- seg_df[seg_df$idx %% 2 == 1, ]
      seg_sex <- seg_df[seg_df$idx %% 2 == 0, ]
      ythr <- as.numeric(thresholds_by_type["Diet Interactive"]) # same as Sex Interactive

      if (nrow(seg_diet) > 0) {
        p <- p + ggplot2::geom_segment(data = seg_diet, ggplot2::aes(x = x0, xend = x1, y = ythr, yend = ythr), inherit.aes = FALSE, color = color_map["Diet Interactive"], linewidth = 0.6)
      }
      if (nrow(seg_sex) > 0) {
        p <- p + ggplot2::geom_segment(data = seg_sex, ggplot2::aes(x = x0, xend = x1, y = ythr, yend = ythr), inherit.aes = FALSE, color = color_map["Sex Interactive"], linewidth = 0.6)
      }
    } else {
      # Draw them separately if only one present or thresholds differ
      if (diet_present) {
        p <- p + ggplot2::geom_hline(yintercept = thresholds_by_type["Diet Interactive"], color = color_map["Diet Interactive"], linetype = "dashed", linewidth = 0.6)
      }
      if (sex_present) {
        p <- p + ggplot2::geom_hline(yintercept = thresholds_by_type["Sex Interactive"], color = color_map["Sex Interactive"], linetype = "dashed", linewidth = 0.6)
      }
    }

    # Sex x Diet interactive
    if ("Sex x Diet Interactive" %in% present_types && !is.na(thresholds_by_type["Sex x Diet Interactive"])) {
      p <- p + ggplot2::geom_hline(yintercept = thresholds_by_type["Sex x Diet Interactive"], color = color_map["Sex x Diet Interactive"], linetype = "dashed", linewidth = 0.6)
    }
  }

  # Add LOD threshold line only for the additive scan (kept for backward compatibility when show_thresholds is FALSE)
  if (isFALSE(show_thresholds) && !is.null(LOD_thr) && is.numeric(LOD_thr) && LOD_thr > 0) {
    p <- p + ggplot2::geom_hline(yintercept = LOD_thr, color = "#e74c3c", linetype = "dashed", linewidth = 0.8)
  }

  # Conditionally hide legend if only one scan type is present (no overlays)
  if (length(unique(plot_data$type)) <= 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  return(p)
}
