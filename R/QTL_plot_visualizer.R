#' Plot the QTL if doing a new scan
#' 
#' @param qtl.temp data frame of QTL scans
#' @param phenotype character string
#' @param LOD_thr value of LOD threshold
#' @param mrkrs list of markers
#' 
#' @importFrom dplyr arrange as_tibble group_by left_join mutate select summarise
#' @importFrom data.table as.data.table setorder
#' @export
QTL_plot_visualizer <- function(qtl.temp, phenotype, LOD_thr, mrkrs) {
  # Create a data.table for faster processing
  qtl.temp <- as.data.table(qtl.temp)
  # Select only the columns we need
  if ("marker" %in% colnames(qtl.temp) && "LOD" %in% colnames(qtl.temp)) {
    # If data already has the right columns, just select them
    qtl.temp <- qtl.temp[, .(markers = marker, LOD)]
  } else {
    # Otherwise, try to find the right columns
    qtl.temp <- qtl.temp[, .(markers = marker, LOD)]
  }
  # Convert markers to data.table for faster joins
  mrkrs2 <- data.table::as.data.table(mrkrs)[, .(markers = marker, chr, position = bp_grcm39/1e6)]
  # Use data.table join instead of merge (much faster)
  qtl.temp <- mrkrs2[qtl.temp, on = "markers", nomatch = NULL]
  # Convert chr to numeric efficiently
  qtl.temp[, chr := as.character(chr)]
  qtl.temp[chr == "X", chr := "20"]
  qtl.temp[chr == "Y", chr := "21"]
  qtl.temp[chr == "M", chr := "22"]
  qtl.temp[, chr := as.numeric(chr)]
  # Remove rows with NA in chr
  qtl.temp <- qtl.temp[!is.na(chr)]
  # Add order column (same as chr for numeric processing)
  qtl.temp[, order := chr]
  # Calculate chromosome lengths and cumulative positions
  chr_info <- qtl.temp[, .(chr_len = max(position)), by = chr]
  chr_info[, tot := cumsum(chr_len) - chr_len]
  # Join back to main data
  qtl.temp <- chr_info[qtl.temp, on = "chr", nomatch = NULL]
  # Calculate cumulative position
  qtl.temp[, BPcum := position + tot]
  # Sort by chromosome and position
  data.table::setorder(qtl.temp, chr, position)
  # Convert back to tibble for ggplot
  qtl_plot_obj <- dplyr::as_tibble(qtl.temp)
  return(qtl_plot_obj)
}
