#' Plot the QTL if doing a new scan
#' 
#' @param scan_data data frame of QTL scans
#' @param phenotype_name character string
#' @param lod_threshold value of LOD threshold
#' @param markers_info list of markers
#' 
#' @importFrom dplyr arrange as_tibble group_by left_join mutate select summarise
#' @importFrom data.table as.data.table setorder
#' @export
QTL_plot_visualizer <- function(scan_data, phenotype_name, lod_threshold, markers_info) {
  # Ensure scan_data is a data.table
  if (!is.data.table(scan_data)) {
    scan_data <- data.table::as.data.table(scan_data)
  }

  # Check for essential columns from trait_scan output
  # trait_scan should return 'marker' (original marker ID) and 'LOD'
  if (!all(c("marker", "LOD") %in% colnames(scan_data))) {
    stop("QTL_plot_visualizer: Input scan_data must contain 'marker' and 'LOD' columns.")
  }
  
  # Select and rename marker column for clarity if needed, keep LOD
  # Use a copy to avoid modifying the original reactive if scan_data is passed directly
  plot_dt <- data.table::copy(scan_data)
  if ("marker" %in% colnames(plot_dt) && !("markers" %in% colnames(plot_dt))) {
      data.table::setnames(plot_dt, "marker", "markers")
  } else if (!("markers" %in% colnames(plot_dt))){
      stop("QTL_plot_visualizer: No 'marker' or 'markers' column in scan_data.")
  }
  
  # Keep only necessary columns: markers and LOD from scan_data
  plot_dt <- plot_dt[, .SD, .SDcols = c("markers", "LOD")]

  # Prepare marker information
  if (!is.data.table(markers_info)) {
    markers_dt <- data.table::as.data.table(markers_info)
  } else {
    markers_dt <- data.table::copy(markers_info)
  }
  
  # Ensure markers_dt has the correct columns ('marker', 'chr', 'bp_grcm39')
  # Renaming 'marker' to 'markers' in markers_dt for the join
  if ("marker" %in% colnames(markers_dt) && !("markers" %in% colnames(markers_dt))){
      data.table::setnames(markers_dt, "marker", "markers")
  }
  if (!all(c("markers", "chr", "bp_grcm39") %in% colnames(markers_dt))){
      stop("QTL_plot_visualizer: markers_info must contain 'marker', 'chr', and 'bp_grcm39' columns.")
  }
  markers_dt <- markers_dt[, .(markers, chr_orig = chr, position = bp_grcm39 / 1e6)] # Keep original chr for X,Y,M

  # Join scan data with marker information
  # This keeps LOD from scan_data (plot_dt)
  plot_dt <- markers_dt[plot_dt, on = "markers", nomatch = 0L] # Use nomatch=0L to drop non-matching rows
  
  if(nrow(plot_dt) == 0){
      warning("QTL_plot_visualizer: No matching markers found after join. Check marker consistency.")
      return(dplyr::as_tibble(data.frame())) # Return empty tibble
  }

  # Convert chr to numeric factor for ordering (X=20, Y=21, M=22)
  plot_dt[, chr_num := as.character(chr_orig)]
  plot_dt[chr_orig == "X", chr_num := "20"]
  plot_dt[chr_orig == "Y", chr_num := "21"]
  plot_dt[chr_orig == "M", chr_num := "22"]
  plot_dt[, chr_num := as.numeric(chr_num)]
  plot_dt <- plot_dt[!is.na(chr_num)]
  data.table::setnames(plot_dt, "chr_num", "chr") # Use 'chr' for the numeric version internally
  plot_dt[, order := chr] # 'order' will be numeric for sorting and x-axis breaks

  # Calculate cumulative positions
  data.table::setkey(plot_dt, chr, position)
  chr_info <- plot_dt[, .(chr_len = max(position, na.rm = TRUE)), by = chr]
  chr_info[, tot := cumsum(as.numeric(chr_len)) - as.numeric(chr_len)] # Ensure numeric for cumsum
  
  plot_dt <- chr_info[plot_dt, on = "chr"]
  plot_dt[, BPcum := position + tot]
  
  # Keep original chr_orig for axis labels if needed, or use the numeric chr for calculations
  # For ggplot_qtl_scan, it expects 'order' (numeric) and uses it to make chr labels
  # Final required columns for ggplot_qtl_scan: BPcum, LOD, order (numeric chr), position
  
  # Convert back to tibble for ggplot compatibility
  return(dplyr::as_tibble(plot_dt))
}
