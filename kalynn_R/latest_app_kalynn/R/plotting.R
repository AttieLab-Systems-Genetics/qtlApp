# Plotting helper functions

# Prepares data for the main LOD scan plot
# Converts chromosome names, calculates cumulative positions
# Returns a data.table suitable for ggplot
QTL_plot_visualizer <- function(qtl_scan_data, mrkrs) {
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("Package 'data.table' needed for QTL_plot_visualizer. Please install it.", call. = FALSE)
  }
  
  # Ensure input is a data.table
  if (!is.data.table(qtl_scan_data)) {
     qtl_scan_data <- data.table::as.data.table(qtl_scan_data)
  }
  
  # Standardize input column names (make copies first)
  qtl_plot_dt <- data.table::copy(qtl_scan_data)
  required_cols <- c("marker", "LOD")
  if (!all(required_cols %in% colnames(qtl_plot_dt))) {
      # Attempt renaming based on common alternatives
      if (!"marker" %in% colnames(qtl_plot_dt) && "markers" %in% colnames(qtl_plot_dt)) {
          data.table::setnames(qtl_plot_dt, "markers", "marker")
      } else if (!"marker" %in% colnames(qtl_plot_dt)){
          stop("Input data must have a 'marker' (or 'markers') column.")
      }
      # Note: LOD standardization happens within trait_scan
      if (!"LOD" %in% colnames(qtl_plot_dt)) {
          stop("Input data must have a 'LOD' column.")
      }
  }
  
  # Select only necessary columns for plotting
  qtl_plot_dt <- qtl_plot_dt[, .(marker, LOD)]
  
  # Prepare marker data (convert to data.table, select columns, calculate Mb)
  if (!is.data.table(mrkrs)) mrkrs <- data.table::as.data.table(mrkrs)
  # Ensure bp_grcm39 column exists before attempting division
  if (!"bp_grcm39" %in% colnames(mrkrs)) stop("Markers data must contain 'bp_grcm39' column.")
  mrkrs_dt <- mrkrs[, .(marker = marker, chr, position = bp_grcm39 / 1e6)]
  
  # Perform efficient join using data.table
  qtl_plot_dt <- mrkrs_dt[qtl_plot_dt, on = "marker", nomatch = NULL] # Use nomatch=NULL to drop markers not found
  
  # Handle chromosome conversion (X, Y, M to numeric) using data.table syntax
  qtl_plot_dt[, chr_original := chr] # Keep original chr name
  qtl_plot_dt[, chr_numeric := data.table::fcase(
    chr == "X", 20L, # Use integer literal L
    chr == "Y", 21L,
    chr == "M", 22L,
    # Try converting others to integer, default to NA if fails
    default = suppressWarnings(as.integer(as.character(chr)))
  )]
  
  # Remove rows with NA chromosome numbers after conversion
  qtl_plot_dt <- qtl_plot_dt[!is.na(chr_numeric)]
  if (nrow(qtl_plot_dt) == 0) {
      warning("No valid chromosome data remaining after merging with markers and converting chromosome names.")
      return(qtl_plot_dt) # Return empty data.table
  }
  
  
  data.table::setorder(qtl_plot_dt, chr_numeric, position)
  # Calculate max position per chromosome efficiently
  chr_info <- qtl_plot_dt[, .(chr_len = max(position)), keyby = .(chr_numeric)] # Use keyby for sorted result
  # Calculate cumulative start position for each chromosome
  chr_info[, tot := cumsum(shift(chr_len, fill = 0))] # More direct way using shift
  
  # Join cumulative start positions back to main data using efficient join
  qtl_plot_dt <- chr_info[, .(chr_numeric, tot)][qtl_plot_dt, on = "chr_numeric"] 
  
  # Calculate cumulative bp position (BPcum)
  qtl_plot_dt[, BPcum := position + tot]
  
  # Ensure final data is sorted correctly (already sorted by keyby and join order)
  # data.table::setorder(qtl_plot_dt, chr_numeric, position)
  
  # Return the processed data.table (ggplot2 generally handles data.table well)
  return(qtl_plot_dt)
} 