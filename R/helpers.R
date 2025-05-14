#' Creates cache objects for storing peaks and traits
#' @param caches A character vector specifying which caches to create (default: c("peaks", "trait")).
#' @export
create_cache <- function(caches = c("peaks", "trait")) {
  
}
#' @export
get_trait_choices <- function(import_data, selected_dataset = NULL) {
  shiny::req(import_data, selected_dataset)
  trait_type <- get_trait_type(import_data, selected_dataset)
  if(is.null(trait_type)) {
    return(character(0))
  }

  trait_list_df <- get_trait_list(import_data, trait_type)
  if(is.null(trait_list_df)){
      return(character(0))
  }

  choices <- character(0)
  if(trait_type %in% c("genes", "isoforms")) {
    if("symbol" %in% colnames(trait_list_df)) {
      choices <- unique(trait_list_df$symbol)
    } else {
      return(character(0))
    }
  } else {
    trait_id <- get_trait_id(trait_type)
    if(is.null(trait_list_df) || !is.data.frame(trait_list_df)){
        warning("Annotation data (trait_list_df) for type '", trait_type, "' is NULL or not a data frame.")
        return(character(0))
    }
    if(!(trait_id %in% colnames(trait_list_df))){
        warning("Trait ID column '", trait_id, "' not found in annotation data frame for type: ", trait_type, ". Available columns: ", paste(colnames(trait_list_df), collapse=", "))
        return(character(0))
    }
    choices <- unique(trait_list_df[[trait_id]])
    valid_choices_count <- sum(!is.na(choices) & choices != "")
    if(valid_choices_count == 0 && length(choices) > 0){
        warning("Column '", trait_id, "' exists but contains only NAs or empty strings.")
    }
  }
  valid_choices <- choices[!is.na(choices) & choices != ""]
  sort(unique(valid_choices))
}
#' @importFrom dplyr arrange desc filter
#' @importFrom rlang .data
#' @export
highest_peaks <- function(peak_table, LOD_thr) {
  # Filter peaks that meet or exceed the LOD threshold
  filtered_peaks <- dplyr::filter(peak_table, .data$lod >= LOD_thr) |>
    dplyr::arrange(dplyr::desc(.data$lod))
  if (nrow(filtered_peaks) == 0) return(NULL)
  filtered_peaks
}
#' @importFrom stringr str_split
#' @export
get_selected_trait <- function(import_data, which_trait, selected_dataset) {
  shiny::req(which_trait)
  trimmed_trait <- trimws(which_trait)
  return(trimmed_trait)
}
#' @importFrom reshape2 melt
#' @export
pivot_peaks <- function(peaks, which_peak) {
  peak_reshaped <- NULL
  if (is.null(peaks) || !is.data.frame(peaks) || nrow(peaks) == 0) {
    warning("pivot_peaks: Input 'peaks' data is NULL, not a data.frame, or has 0 rows.")
    return(peak_reshaped)
  }
  if (is.null(which_peak) || which_peak == "" || is.na(which_peak)){
    warning("pivot_peaks: 'which_peak' (marker for filtering) is NULL, empty or NA.")
    return(peak_reshaped)
  }
  if (!("marker" %in% colnames(peaks))) {
    warning("pivot_peaks: 'marker' column not found in peaks data.")
    return(peak_reshaped)
  }
  peak_subset <- subset(peaks, marker == which_peak)
  if (nrow(peak_subset) > 0) {
    allele_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
    strain_names <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
    if (all(allele_cols %in% colnames(peak_subset))) {
      peak_for_melt <- peak_subset[, c("marker", allele_cols), drop = FALSE]
      if (length(allele_cols) == length(strain_names)) {
          colnames(peak_for_melt)[(which(colnames(peak_for_melt) %in% allele_cols))] <- strain_names
      } else {
          warning("pivot_peaks: Mismatch between allele_cols and strain_names length. Using original allele column names.")
          strain_names <- allele_cols
      }
      peak_reshaped <- reshape2::melt(peak_for_melt, 
                                      id.vars = "marker", 
                                      measure.vars = strain_names, 
                                      variable.name = "variable", 
                                      value.name = "value")
    } else {
      warning(paste("pivot_peaks: Not all allele columns (A-H) found for marker:", which_peak, ". Cannot reshape."))
    }
  } else {
    warning(paste("pivot_peaks: No peak data found for marker:", which_peak, "after subsetting."))
  }
  return(peak_reshaped)
}
# internal helper functions
get_trait_type <- function(import_data, selected_dataset = NULL) {
  if(is.null(import_data) || is.null(import_data$file_directory)){
      warning("get_trait_type: import_data or import_data$file_directory is NULL.")
      return(NULL)
  }
  file_directory <- import_data$file_directory
  
  if(is.null(selected_dataset) || selected_dataset == "") {
    return(NULL)
  }
  
  # Filter directory for the selected dataset
  file_dir_subset <- subset(file_directory, group == selected_dataset)
  if(nrow(file_dir_subset) == 0) {
    warning("No entry found for dataset: ", selected_dataset, " in file directory.")
    return(NULL)
  }
  
  # Assuming trait_type column exists
  if(!("trait_type" %in% colnames(file_dir_subset))) {
     warning("'trait_type' column missing from file directory.")
     return(NULL)
  }
  trait_type_raw <- tolower(file_dir_subset$trait_type[1])
  
  # Standardize trait types
  if (grepl("clinical", trait_type_raw, ignore.case = TRUE)) { # Using grepl for consistency
      return("clinical")
  }
  if (grepl("gene", trait_type_raw, ignore.case = TRUE)) { # Check if "gene" is present anywhere
      return("genes")
  }
  if (grepl("isoform", trait_type_raw, ignore.case = TRUE)) { # Check if "isoform" is present anywhere
      return("isoforms")
  }
  
  return(trait_type_raw) # Return the original (lowercased) if no specific keyword matched
}
get_trait_list <- function(import_data, trait_type) {
  if(is.null(import_data) || is.null(import_data$annotation_list)){
      warning("get_trait_list: import_data or import_data$annotation_list is NULL.")
      return(NULL)
  }
  annotation_list <- import_data$annotation_list
  
  if(is.null(trait_type) || !(trait_type %in% names(annotation_list))) {
     warning("Trait type ", trait_type, " not found in annotation list.")
     return(NULL)
  }
  annotation_list[[trait_type]]
}
get_trait_id <- function(trait_type) {
  switch(trait_type,
    genes    = "gene.id",
    isoforms = "transcript.id",
    "data_name")
}
chr_XYM <- function(chr_vector) {
  if(is.null(chr_vector)) return(NA_character_) # Return NA of appropriate type
  
  # Ensure it's character for processing
  chr_vector_char <- as.character(chr_vector)
  
  # Replace values
  # Use a more robust way to handle vectors for replacement
  new_chr_vector <- chr_vector_char # Initialize with original character values
  new_chr_vector[chr_vector_char == "20"] <- "X"
  new_chr_vector[chr_vector_char == "21"] <- "Y"
  new_chr_vector[chr_vector_char == "22"] <- "M"
  
  # Handle cases where input might not be 1-22 (e.g. already X, Y, M, or other strings)
  # The above replacements are safe. If it was already "X", it remains "X".
  
  return(new_chr_vector)
}
