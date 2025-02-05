# plot the QTL if doing a new scan
QTL_plot_visualizer <- function(qtl.temp, phenotype, LOD_thr, mrkrs) {
  # The LOD values are already in the "LOD" column
  qtl.temp <- qtl.temp[, c("marker", "LOD")]
  colnames(qtl.temp) <- c("markers", "LOD")
  
  mrkrs2 <- mrkrs[c("marker", "chr", "bp_grcm39")]
  colnames(mrkrs2) <- c("markers", "chr", "position")
  mrkrs2$position <- as.numeric(mrkrs2$position)/(10^6)
  
  # Join the data
  qtl.temp <- merge(qtl.temp, mrkrs2, by.x="markers", by.y="markers", all.x=FALSE)
  qtl.temp <- na.omit(qtl.temp)
  
  # Convert chr to numeric
  qtl.temp$chr <- as.character(qtl.temp$chr)
  qtl.temp$order <- sapply(qtl.temp$chr, function(x) {
    if (x == "X") return(20)
    else if (x == "Y") return(21)
    else if (x == "M") return(22)
    else return(as.numeric(x))
  })
  
  # Remove any rows where chromosome conversion failed
  qtl.temp <- qtl.temp[!is.na(qtl.temp$order), ]
  qtl.temp$chr <- qtl.temp$order
  
  # Create plot object
  qtl_plot_obj <- qtl.temp |>
    dplyr::group_by(chr) |>
    dplyr::summarise(chr_len = max(position)) |>
    dplyr::mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) |>
    dplyr::select(-chr_len) |>
    dplyr::left_join(qtl.temp, ., by = c("chr" = "chr")) |>
    dplyr::arrange(order, position) |>
    dplyr::mutate(BPcum = position + tot)
  
  # Create axis labels
  axisdf = qtl_plot_obj |>
    dplyr::group_by(order) |>
    dplyr::summarize(center = (max(BPcum) + min(BPcum))/2)
  
  # Convert chromosome labels
  axisdf$order[axisdf$order == 20] <- "X"
  axisdf$order[axisdf$order == 21] <- "Y"
  axisdf$order[axisdf$order == 22] <- "M"
  axisdf$order <- factor(axisdf$order, levels=c(as.character(1:19),"X","Y","M"))
  
  # Create plot with x-axis
  plot_QTL <- ggplot2::ggplot(qtl_plot_obj, ggplot2::aes(x=BPcum, y=LOD)) +
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
    ggplot2::ylab("LOD") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=LOD_thr), color="black", linetype="dashed")
  
  plots <- list(plot_QTL, qtl_plot_obj)
  return(plots)
}
