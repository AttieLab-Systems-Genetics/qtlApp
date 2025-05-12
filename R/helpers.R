#' @export
create_cache <- function(caches = c("peaks", "trait")) {
  
}
#' @export
get_trait_choices <- function(import, selected_dataset = NULL) {
  trait_type <- get_trait_type(import, selected_dataset)
  if(is.null(trait_type)) {
    return(NULL)
  }
  trait_id <- get_trait_id(trait_type)
  trait_list = get_trait_list(import, trait_type)
  
  annotation_list <- import$annotation_list
  if(!(trait_type %in% c("genes", "isoforms"))) {
    choices <- annotation_list[[trait_type]][[trait_id]]
  } else { # "genes", "isoforms"
    choices <- join_symbol_id(annotation_list, trait_type, trait_id)
  }
  choices
}
#' @importFrom stringr str_remove
#' @export
join_symbol_id <- function(annotation_list,
                           trait_type = c("genes","isoforms"),
                           trait_id = get_trait_id(trait_type)) {
  
  trait_type <- match.arg(trait_type)
  
  id <- stringr::str_remove(annotation_list[[trait_type]][[trait_id]],
                            pattern = "ENSMUS[GT]0+")
  paste(annotation_list[[trait_type]]$symbol, id, sep = "_")
}
#' @importFrom dplyr arrange desc filter
#' @importFrom rlang .data
#' @export
highest_peaks <- function(peak_table, LOD_thr) {
  filtered_peaks <- dplyr::filter(peak_table, .data$lod >= LOD_thr) |>
    dplyr::arrange(dplyr::desc(.data$lod))
  if (nrow(filtered_peaks) == 0) return(NULL)
  filtered_peaks
}
#' @importFrom stringr str_split
#' @export
get_selected_trait <- function(import, which_trait, selected_dataset) {
 
  trait_type <- get_trait_type(import, selected_dataset)
  if(trait_type %in% c("genes", "isoforms")) {
    which_trait <- stringr::str_split(which_trait, pattern="_")[[1]][1]
  }
  which_trait
}
#' @importFrom reshape2 melt
#' @export
pivot_peaks <- function(peaks, which_peak) {
  # Initialize peak as NULL to handle cases where no valid peak data is found
  peak_reshaped <- NULL
  
  # Ensure 'peaks' is a data.frame and has rows
  if (is.null(peaks) || !is.data.frame(peaks) || nrow(peaks) == 0) {
    warning("pivot_peaks: Input 'peaks' data is NULL, not a data.frame, or has 0 rows.")
    return(peak_reshaped)
  }
  
  # Ensure 'which_peak' is provided
  if (is.null(which_peak) || which_peak == "" || is.na(which_peak)){
    warning("pivot_peaks: 'which_peak' (marker for filtering) is NULL, empty or NA.")
    return(peak_reshaped)
  }

  # Filter for the specific peak based on the 'marker' column
  # Ensure 'marker' column exists
  if (!("marker" %in% colnames(peaks))) {
    warning("pivot_peaks: 'marker' column not found in peaks data.")
    return(peak_reshaped)
  }
  peak_subset <- subset(peaks, marker == which_peak)
  
  # Check if we have data after filtering for the specific peak
  if (nrow(peak_subset) > 0) {
    # Define the allele columns we expect
    allele_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
    strain_names <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
    
    # Check if all required allele columns exist in the peak_subset
    if (all(allele_cols %in% colnames(peak_subset))) {
      # Select the marker and allele columns
      peak_for_melt <- peak_subset[, c("marker", allele_cols), drop = FALSE]
      
      # Rename allele columns to strain names for melting
      # Ensure the number of strain_names matches allele_cols
      if (length(allele_cols) == length(strain_names)) {
          colnames(peak_for_melt)[(which(colnames(peak_for_melt) %in% allele_cols))] <- strain_names
      } else {
          warning("pivot_peaks: Mismatch between allele_cols and strain_names length. Using original allele column names.")
          # If mismatch, melt with original A-H names (less ideal for plot labels)
          strain_names <- allele_cols # Fallback
      }
      
      # Reshape data using the (potentially renamed) strain_names as measure.vars
      # Ensure marker is the id.var
      peak_reshaped <- reshape2::melt(peak_for_melt, 
                                      id.vars = "marker", 
                                      measure.vars = strain_names, # Use the actual column names present after renaming
                                      variable.name = "variable", 
                                      value.name = "value")
      message(paste("pivot_peaks: Successfully reshaped data for marker:", which_peak))
    } else {
      warning(paste("pivot_peaks: Not all allele columns (A-H) found for marker:", which_peak, ". Cannot reshape."))
    }
  } else {
    warning(paste("pivot_peaks: No peak data found for marker:", which_peak, "after subsetting."))
  }
  
  # The decision based on 'intcovar' might still be needed if the *source* of allele effects 
  # (A-H columns) fundamentally differs for additive vs. interactive, 
  # or if additional interaction-specific columns need to be pivoted.
  # For now, this uses A-H for both, assuming they are appropriate.
  # If interactive scans have different columns for effects (e.g. dietA_eff, dietB_eff), 
  # then the column selection and melting logic would need to be conditional on 'intcovar'.

  return(peak_reshaped)
}
# internal helper functions
get_trait_type <- function(import, selected_dataset = NULL) {
  if(is.null(selected_dataset)) {
    return(NULL)
  }
  file_directory <- import$file_directory
  # Filter for the selected dataset
  file_directory_subset <- subset(file_directory, group == selected_dataset)
  
  if(nrow(file_directory_subset) == 0) {
    warning(paste("No entry found for dataset:", selected_dataset, "in file directory. Cannot determine trait_type."))
    return(NULL)
  }
  
  # Get trait_type from the first entry of the subset and convert to lowercase
  trait_type_from_file <- tolower(file_directory_subset$trait_type[1])
  
  # Standardize "clinical traits" to "clinical"
  if (trait_type_from_file == "clinical traits") {
    return("clinical")
  }
  
  return(trait_type_from_file)
}
get_trait_list <- function(import, trait_type) {
  annotation_list <- import$annotation_list
  annotation_list[[trait_type]]
}
get_trait_id <- function(trait_type) {
  switch(trait_type,
    genes    = "gene.id",
    isoforms = "transcript.id",
    "data_name")
}
chr_XYM <- function(chr) {
  if(is.null(chr)) return(NA)
  
  if (chr %in% c(20, 21, 22)) {
    c("X", "Y", "M")[chr - 19]
  } else {
    chr
  }
}
