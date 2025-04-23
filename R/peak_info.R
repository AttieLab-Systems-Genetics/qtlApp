peak_info <- function(peak_table, scan_table, which_peak = NULL, click_table = NULL) {
  if (is.null(click_table)) {
    # If no click data, try to get information from the selected peak instead
    if(is.null(which_peak)) return(NULL)
      return(peak_info_table(peak_table, which_peak))
  }
  # If click data is available, process it with scan_table.
  if (is.null(scan_table)) return(NULL)
  # Find nearest point if any.
  distances <- sqrt(
    (scan_table$BPcum - click_table$x)^2 + (scan_table$LOD - click_table$y)^2)
  nearest_point <- scan_table[which.min(distances), ]
  if (nrow(nearest_point) == 0) return(NULL)
  # Return the peak closest to the clicked point.
  return(peak_info_table(peak_table, which_peak, nearest_point))  
}
peak_info_table <- function(peak_table,
                            which_peak = nearest_point$markers,
                            nearest_point = NULL) {
  # Find peak by marker, considering possible column name differences
  peak_match <- NULL
  if ("marker" %in% colnames(peak_table) & !is.null(which_peak)) {
    peak_match <- dplyr::filter(peak_table, .data$marker == which_peak)
  }
  if(is.null(peak_match) || !nrow(peak_match)) {
    peak_match <- data.frame(marker = NA, chr = NA, pos = NA, lod = NA,
      trait = NA, cis = NA, ci_lo = NA, ci_hi = NA,
      A = NA, B = NA, C = NA, D = NA, E = NA, F = NA, G = NA, H = NA)
  }
  # Create peak info table with safer column access
  out <- data.frame(
    Marker = peak_match$marker[1],
    Chromosome = {
      chr_val <- peak_match$chr[1]
      if (!is.na(chr_val) && is.numeric(chr_val) && chr_val > 19) {
        chr_val <- c("X", "Y", "M")[chr_val - 19]
      }
      chr_val
    },
    Position = round(peak_match$pos[1], 3),
    LOD = round(peak_match$lod[1], 3),
    Trait = peak_match$trait[1],
    Cis = peak_match$cis[1],
    CI_Low = round(peak_match$ci_lo[1], 3),
    CI_High = round(peak_match$ci_hi[1], 3),
    # Add strain effects if available
    AJ = round(peak_match$A[1], 3),
    B6 = round(peak_match$B[1], 3),
    `129` = round(peak_match$C[1], 3),
    NOD = round(peak_match$D[1], 3),
    NZO = round(peak_match$E[1], 3),
    CAST = round(peak_match$F[1], 3),
    PWK = round(peak_match$G[1], 3),
    WSB = round(peak_match$H[1], 3))
  if(!is.null(nearest_point)) {
    if(is.na(out$Marker)) out$Marker <- nearest_point$markers
    if(is.na(out$Chromosome)) {
      out$Chromosome <- 
        if(nearest_point$chr %in% c(20,21,22)) c("X","Y","M")[nearest_point$chr-19]
        else nearest_point$chr
    }
    if(is.na(out$Position)) out$Position <- round(nearest_point$position, 3)
    if(is.na(out$LOD)) out$LOD <- round(nearest_point$LOD, 3)
  }
  out
}