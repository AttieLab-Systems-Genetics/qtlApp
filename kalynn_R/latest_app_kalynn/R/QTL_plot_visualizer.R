#' Plot the QTL if doing a new scan
#' 
#' @param qtl.temp data frame of QTL scans
#' @param phenotype character string
#' @param LOD_thr value of LOD threshold
#' @param mrkrs list of markers
#' 
#' @importFrom dplyr arrange as_tibble group_by left_join mutate select summarise
#' @importFrom data.table as.data.table setorder `:=`
#' @export
QTL_plot_visualizer <- function(qtl.temp, phenotype, LOD_thr, mrkrs) {
  # Create a data.table for faster processing
  qtl.temp <- data.table::as.data.table(qtl.temp)
 
  # Select only the columns we need
  if ("marker" %in% colnames(qtl.temp) && "LOD" %in% colnames(qtl.temp)) {
      # If data already has the right columns, just select them
      qtl.temp <- qtl.temp[, .(markers = marker, LOD)]
  } else {
      # Otherwise, try to find the right columns
      # Assuming 'marker' and 'LOD' exist with potentially different casing
      marker_col <- grep("^marker$", colnames(qtl.temp), ignore.case = TRUE, value = TRUE)[1]
      lod_col <- grep("^LOD$", colnames(qtl.temp), ignore.case = TRUE, value = TRUE)[1]
      if(is.na(marker_col) || is.na(lod_col)) {
          stop("Could not find required 'marker' and 'LOD' columns in qtl.temp")
      }
      qtl.temp <- qtl.temp[, .(markers = get(marker_col), LOD = get(lod_col))]
  }
 
  # Convert markers to data.table for faster joins
  # Ensure mrkrs is a data.table or data.frame
  if(!is.data.table(mrkrs) && !is.data.frame(mrkrs)) stop("'mrkrs' must be a data.table or data.frame")
  mrkrs2 <- data.table::as.data.table(mrkrs)[, .(markers = marker, chr, position = bp_grcm39/1e6)]
 
  # Use data.table join instead of merge (much faster)
  # Make sure marker columns are compatible type
  qtl.temp[, markers := as.character(markers)]
  mrkrs2[, markers := as.character(markers)]
  qtl.temp <- mrkrs2[qtl.temp, on = "markers", nomatch = NULL]
 
  # Convert chr to numeric 
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
  # Ensure 'position' column exists and is numeric
  if(!("position" %in% colnames(qtl.temp)) || !is.numeric(qtl.temp$position)) {
      stop("'position' column is missing or not numeric in qtl.temp after join")
  }
  chr_info <- qtl.temp[, .(chr_len = max(position, na.rm = TRUE)), by = chr] # Add na.rm=TRUE
  # Check for issues with max if all positions are NA for a chr
  if(anyNA(chr_info$chr_len)) warning("NA values encountered when calculating max chromosome length.")
  chr_info[, tot := cumsum(chr_len) - chr_len]
 
  # Join back to main data
  qtl.temp <- chr_info[qtl.temp, on = "chr", nomatch = NULL]
 
  # Calculate cumulative position
  qtl.temp[, BPcum := position + tot]
 
  # Sort by chromosome and position
  data.table::setorder(qtl.temp, chr, position)
 
  qtl_plot_obj <- dplyr::as_tibble(qtl.temp)
 
  return(qtl_plot_obj)
} 