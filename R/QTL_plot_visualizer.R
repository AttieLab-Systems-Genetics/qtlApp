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
  message("--- QTL_plot_visualizer: Start ---")

  # --- CRITICAL INPUT VALIDATION ---
  if (is.null(scan_data)) {
    message("QTL_plot_visualizer ERROR: scan_data input is NULL. Returning empty plot.")
    return(dplyr::as_tibble(data.frame()))
  }
  if (!is.data.frame(scan_data) && !is.data.table(scan_data)) {
    message("QTL_plot_visualizer ERROR: scan_data input is not a data.frame or data.table. Returning empty plot.")
    str(scan_data) # Print structure of whatever it is
    return(dplyr::as_tibble(data.frame()))
  }
  if (nrow(scan_data) == 0) {
    message("QTL_plot_visualizer WARNING: scan_data input has 0 rows. Returning empty plot.")
    return(dplyr::as_tibble(data.frame()))
  }
  if (is.null(markers_info)) {
    message("QTL_plot_visualizer ERROR: markers_info input is NULL. Returning empty plot.")
    return(dplyr::as_tibble(data.frame()))
  }
  if (!is.data.frame(markers_info) && !is.data.table(markers_info)) {
    message("QTL_plot_visualizer ERROR: markers_info input is not a data.frame or data.table. Returning empty plot.")
    str(markers_info)
    return(dplyr::as_tibble(data.frame()))
  }
  if (nrow(markers_info) == 0) {
    message("QTL_plot_visualizer WARNING: markers_info input has 0 rows. Returning empty plot.")
    return(dplyr::as_tibble(data.frame()))
  }
  # --- END CRITICAL INPUT VALIDATION ---

  message("QTL_plot_visualizer DEBUG: Initial scan_data (head & str) after passing critical validation:")
  print(head(scan_data))
  str(scan_data)
  message("QTL_plot_visualizer DEBUG: Initial markers_info (head & str):")
  print(head(markers_info))
  str(markers_info)

  if (!is.data.table(scan_data)) {
    scan_data <- data.table::as.data.table(scan_data)
  }
  if (!all(c("marker", "LOD") %in% colnames(scan_data))) {
    stop("QTL_plot_visualizer: Input scan_data must contain 'marker' and 'LOD' columns.")
  }
  
  plot_dt <- data.table::copy(scan_data)
  if ("marker" %in% colnames(plot_dt) && !("markers" %in% colnames(plot_dt))) {
      data.table::setnames(plot_dt, "marker", "markers")
  } else if (!("markers" %in% colnames(plot_dt))){
      stop("QTL_plot_visualizer: No 'marker' or 'markers' column in scan_data.")
  }
  
  plot_dt <- plot_dt[, .SD, .SDcols = c("markers", "LOD")]
  message("QTL_plot_visualizer DEBUG: plot_dt after selecting columns (head & str):")
  print(head(plot_dt))
  str(plot_dt)

  if (is.data.table(plot_dt) && nrow(plot_dt) > 0 && "markers" %in% colnames(plot_dt)) {
    plot_dt[, markers := as.character(markers)]
  } else {
    message("QTL_plot_visualizer WARNING: plot_dt is not a valid data.table with a markers column before character conversion, or is empty.")
  }

  if (!is.data.table(markers_info)) {
    markers_dt <- data.table::as.data.table(markers_info)
  } else {
    markers_dt <- data.table::copy(markers_info)
  }
  
  if ("marker" %in% colnames(markers_dt) && !("markers" %in% colnames(markers_dt))){
      data.table::setnames(markers_dt, "marker", "markers")
  }
  if (is.data.table(markers_dt) && nrow(markers_dt) > 0 && "markers" %in% colnames(markers_dt)) {
    markers_dt[, markers := as.character(markers)]
  } else {
    message("QTL_plot_visualizer WARNING: markers_dt is not a valid data.table with a markers column before character conversion, or is empty.")
  }

  if (!all(c("markers", "chr", "bp_grcm39") %in% colnames(markers_dt))){
      stop("QTL_plot_visualizer: markers_info must contain 'markers', 'chr', and 'bp_grcm39' columns.")
  }
  markers_dt <- markers_dt[, .(markers, chr_orig = chr, position = bp_grcm39 / 1e6)]
  message("QTL_plot_visualizer DEBUG: markers_dt after selecting columns (head & str):")
  print(head(markers_dt))
  str(markers_dt)

  plot_dt <- data.table::setDF(plot_dt)
  markers_dt <- data.table::setDF(markers_dt)

  # --- DEBUGGING MARKER JOIN --- 
  message("QTL_plot_visualizer DEBUG: Structure of plot_dt$markers before join:")
  if (nrow(plot_dt) > 0) str(plot_dt$markers) else message("plot_dt is empty or plot_dt$markers does not exist")
  message("QTL_plot_visualizer DEBUG: Head of plot_dt$markers (dput):")
  if (nrow(plot_dt) > 0) dput(head(plot_dt$markers)) else message("plot_dt is empty or plot_dt$markers does not exist")
  
  message("QTL_plot_visualizer DEBUG: Structure of markers_dt$marker before join:")
  # Note: after renaming, it's markers_dt$markers, not markers_dt$marker for the column name
  if (nrow(markers_dt) > 0 && "markers" %in% colnames(markers_dt)) str(markers_dt$markers) else message("markers_dt is empty or markers_dt$markers does not exist")
  message("QTL_plot_visualizer DEBUG: Head of markers_dt$markers (dput):")
  if (nrow(markers_dt) > 0 && "markers" %in% colnames(markers_dt)) dput(head(markers_dt$markers)) else message("markers_dt is empty or markers_dt$markers does not exist")
  
  # Ensure marker columns are character for robust joining
  # This conversion is now done earlier and conditionally
  
  message("QTL_plot_visualizer DEBUG: Structure of plot_dt$markers AFTER initial as.character():")
  if (nrow(plot_dt) > 0) str(plot_dt$markers) else message("plot_dt is empty or plot_dt$markers does not exist")
  message("QTL_plot_visualizer DEBUG: Structure of markers_dt$markers AFTER initial as.character():")
  if (nrow(markers_dt) > 0 && "markers" %in% colnames(markers_dt)) str(markers_dt$markers) else message("markers_dt is empty or markers_dt$markers does not exist")
  # --- END DEBUGGING MARKER JOIN ---

  # Perform the join only if both tables have rows and the join column
  if (nrow(plot_dt) > 0 && nrow(markers_dt) > 0 && "markers" %in% colnames(plot_dt) && "markers" %in% colnames(markers_dt)) {
      plot_dt <- dplyr::inner_join(plot_dt, markers_dt, by = "markers") # Simplified by = "markers" if both are named "markers"
  } else {
      message("QTL_plot_visualizer: Skipping join because one or both tables are empty or missing 'markers' column.")
      plot_dt <- data.frame() # Ensure plot_dt is an empty data.frame to prevent downstream errors
  }

  if (nrow(plot_dt) == 0) {
      warning("QTL_plot_visualizer: No matching markers found after join or input tables were empty. Check marker consistency.")
      return(dplyr::as_tibble(data.frame())) # Return empty tibble
  }

  # Convert chr to numeric factor for ordering (X=20, Y=21, M=22)
  # Ensure plot_dt is a data.table for := operations
  if (!is.data.table(plot_dt)) plot_dt <- data.table::as.data.table(plot_dt)
  
  # Check if chr_orig exists before trying to use it (it should, if join was successful)
  if ("chr_orig" %in% colnames(plot_dt)) {
    plot_dt[, chr_num := as.character(chr_orig)]
    plot_dt[chr_orig == "X", chr_num := "20"]
    plot_dt[chr_orig == "Y", chr_num := "21"]
    plot_dt[chr_orig == "M", chr_num := "22"]
    plot_dt[, chr_num := as.numeric(chr_num)]
    plot_dt <- plot_dt[!is.na(chr_num)]
    data.table::setnames(plot_dt, "chr_num", "chr") 
    plot_dt[, order := chr]
  } else {
    warning("QTL_plot_visualizer: chr_orig column not found after join. Cannot calculate cumulative positions.")
    return(dplyr::as_tibble(data.frame()))
  }

  # Calculate cumulative positions
  data.table::setkey(plot_dt, chr, position)
  # Ensure 'position' column exists and is numeric before this step
  if (!("position" %in% colnames(plot_dt)) || !is.numeric(plot_dt$position)){
      warning("QTL_plot_visualizer: 'position' column missing or not numeric in plot_dt before cumulative calculation.")
      return(dplyr::as_tibble(data.frame()))
  }
  chr_info <- plot_dt[, .(chr_len = max(position, na.rm = TRUE)), by = chr]
  chr_info[, tot := cumsum(as.numeric(chr_len)) - as.numeric(chr_len)]
  
  plot_dt <- chr_info[plot_dt, on = "chr"]
  # Ensure 'position' and 'tot' exist before this addition
  if (!all(c("position", "tot") %in% colnames(plot_dt))){
      warning("QTL_plot_visualizer: 'position' or 'tot' column missing before BPcum calculation.")
      return(dplyr::as_tibble(data.frame()))
  }
  plot_dt[, BPcum := position + tot]
  
  message("--- QTL_plot_visualizer: End ---")
  return(dplyr::as_tibble(plot_dt))
}
