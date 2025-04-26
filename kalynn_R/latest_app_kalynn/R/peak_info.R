#' Get Peak/Point Information
#'
#' Determines the information to display in the info table based on either 
#' a Plotly click event or the currently selected peak.
#'
#' @param peak_table Data frame of peaks for the current trait.
#' @param plot_data Data frame used for the QTL plot.
#' @param selected_peak_marker Marker ID of the peak selected in the dropdown (from peak module).
#' @param click_event Plotly click event data (from `event_data("plotly_click", ...)`).
#' @param selected_chr Currently selected chromosome view ("All" or specific chr).
#'
#' @return A data frame with one row containing information about the selected/
#'         clicked point, or NULL if no relevant point is found.
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @export
peak_info <- function(peak_table, plot_data, selected_peak_marker = NULL, 
                      click_event = NULL, selected_chr = "All") {
  
  # --- Case 1: Handle Plotly Click Event --- 
  if (!is.null(click_event)) {
    message("peak_info: Handling plotly click event.")
    if (is.null(plot_data)) return(NULL)
    
    # Find nearest point in plot_data based on click
    x_clicked <- click_event$x
    y_clicked <- click_event$y
    xvar <- if(selected_chr == "All") "BPcum" else "position"
    
    if(!(xvar %in% colnames(plot_data))){
         warning("Column ", xvar, " not found in plot data for click detection.")
         return(NULL)
    }
        
    distances <- sqrt((plot_data[[xvar]] - x_clicked)^2 + (plot_data$LOD - y_clicked)^2)
    nearest_point <- plot_data[which.min(distances), ]
       
    if (nrow(nearest_point) == 0) return(NULL)
    
    message("peak_info: Nearest point to click: ", nearest_point$markers)
    
    # Check if this nearest point is in the peak_table
    peak_match <- NULL
    if (!is.null(peak_table) && "marker" %in% colnames(peak_table) && nrow(peak_table) > 0) {
        peak_match <- peak_table %>% dplyr::filter(.data$marker == nearest_point$markers)
    } 
        
    if (!is.null(peak_match) && nrow(peak_match) > 0) {
        message("peak_info: Clicked point matches a peak in peak_table.")
        # Use the full info from the matched peak
        return(create_peak_info_df(peak_match[1,])) # Use helper function
    } else {
        message("peak_info: Clicked point not in peak_table. Displaying basic info.")
        # If no match in peak table, display basic info from the nearest plot point
        point_info <- data.frame(
            Marker = nearest_point$markers,
            Chromosome = chr_XYM(nearest_point$chr),
            Position = round(nearest_point$position, 3),
            LOD = round(nearest_point$LOD, 3)
        )
        return(point_info)
    }
  }
  
  # --- Case 2: Handle Selected Peak (No Click Event) --- 
  if (!is.null(selected_peak_marker) && selected_peak_marker != "") {
      message("peak_info: No click event, using selected peak: ", selected_peak_marker)
      if (is.null(peak_table) || !("marker" %in% colnames(peak_table))) { 
           warning("peak_info: peak_table is NULL or missing 'marker' column.")
           return(NULL)
      }
      
      peak_match <- peak_table %>% dplyr::filter(.data$marker == selected_peak_marker)
      
      if(nrow(peak_match) > 0) {
          message("peak_info: Found selected peak in peak_table.")
          return(create_peak_info_df(peak_match[1,])) # Use helper function
      } else {
          warning("peak_info: Selected peak ", selected_peak_marker, " not found in peak_table.")
          return(NULL)
      }
  }
  
  # --- Case 3: No Click, No Selection --- 
  message("peak_info: No click event or selected peak.")
  return(NULL)
}

#' Create Peak Info Data Frame Helper
#' 
#' Internal helper function to construct the standard info data frame from a peak row.
#' 
#' @param peak_row A single row data frame or list corresponding to a peak.
#' @return A single-row data frame with formatted peak information.
create_peak_info_df <- function(peak_row) {
    if(is.null(peak_row) || nrow(peak_row) == 0) return(NULL)
    
    # Helper to safely get and round values
    get_val <- function(col_name, digits = 3) {
        val <- peak_row[[col_name]]
        if(is.null(val) || length(val) == 0 || is.na(val[1])) return(NA)
        if(is.numeric(val) && !is.na(digits)) round(val[1], digits) else val[1]
    }

    # Build the data frame column by column
    info_df <- data.frame(row.names = 1) # Ensure it's a data frame
    info_df$Marker <- get_val("marker", NA)
    info_df$Chromosome <- chr_XYM(get_val("chr", NA))
    info_df$Position <- get_val("pos")
    info_df$LOD <- get_val("lod")
    if ("trait" %in% names(peak_row)) info_df$Trait <- get_val("trait", NA)
    if ("cis" %in% names(peak_row)) info_df$Cis <- get_val("cis", NA)
    if ("ci_lo" %in% names(peak_row)) info_df$CI_Low <- get_val("ci_lo")
    if ("ci_hi" %in% names(peak_row)) info_df$CI_High <- get_val("ci_hi")
    if ("A" %in% names(peak_row)) info_df$AJ <- get_val("A")
    if ("B" %in% names(peak_row)) info_df$B6 <- get_val("B")
    if ("C" %in% names(peak_row)) info_df$`129` <- get_val("C")
    if ("D" %in% names(peak_row)) info_df$NOD <- get_val("D")
    if ("E" %in% names(peak_row)) info_df$NZO <- get_val("E")
    if ("F" %in% names(peak_row)) info_df$CAST <- get_val("F")
    if ("G" %in% names(peak_row)) info_df$PWK <- get_val("G")
    if ("H" %in% names(peak_row)) info_df$WSB <- get_val("H")
    
    return(info_df)
} 