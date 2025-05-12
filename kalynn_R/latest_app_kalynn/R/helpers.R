#' Create Cache Environments
#' 
#' Creates global environments for caching peaks and trait data if they don't exist.
#' 
#' @param caches A character vector specifying which caches to create (default: c("peaks", "trait")).
#' @export
create_cache <- function(caches = c("peaks", "trait")) {
  message("--- Running create_cache() ---")
  for(cache_name in caches) {
    cache_obj_name <- paste0(cache_name, "_cache")
    message("  Assigning ", cache_obj_name, " to globalenv()")
    assign(cache_obj_name, new.env(parent = emptyenv()),
      pos = globalenv())
    # Verify assignment
    if(exists(cache_obj_name, envir = globalenv())){
        message("  Verified: ", cache_obj_name, " exists in globalenv()")
    } else {
        warning("  Verification FAILED: ", cache_obj_name, " NOT found in globalenv() after assign()")
    }
  }
  message("--- Finished create_cache() ---")
}

#' Get Trait Type from Selected Dataset
#' 
#' Extracts the trait type (e.g., "genes", "clinical") from the dataset group name.
#' 
#' @param file_directory The file directory data frame.
#' @param selected_dataset The currently selected dataset string.
#' @return The trait type string in lowercase, or NULL if not found.
get_trait_type <- function(file_directory, selected_dataset = NULL) {
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
  trait_type <- tolower(file_dir_subset$trait_type[1])
  return(trait_type)
}

#' Get Trait ID Column Name
#' 
#' Returns the standard ID column name based on the trait type.
#' 
#' @param trait_type The trait type string (e.g., "genes", "isoforms", "clinical").
#' @return The corresponding ID column name string.
get_trait_id <- function(trait_type) {
  switch(trait_type,
    genes    = "gene.id",
    isoforms = "transcript.id",
    # Add other types as needed
    "data_name") # Default for clinical, etc.
}

#' Get Trait List from Annotations
#' 
#' Retrieves the specific annotation data frame based on the trait type.
#' 
#' @param annotation_list The list containing all annotation data frames.
#' @param trait_type The trait type string.
#' @return The corresponding annotation data frame, or NULL if not found.
get_trait_list <- function(annotation_list, trait_type) {
  if(is.null(trait_type) || !(trait_type %in% names(annotation_list))) {
     warning("Trait type ", trait_type, " not found in annotation list.")
     return(NULL)
  }
  annotation_list[[trait_type]]
}

#' Join Symbol and ID for Gene/Isoform Choices
#' 
#' Creates formatted choices (e.g., "Symbol_ID") for gene/isoform selection.
#' 
#' @param annotation_list The list containing all annotation data frames.
#' @param trait_type The trait type ("genes" or "isoforms").
#' @param trait_id The corresponding ID column name.
#' @return A character vector of formatted choices.
#' @importFrom stringr str_remove
join_symbol_id <- function(annotation_list,
                           trait_type = c("genes","isoforms"),
                           trait_id = get_trait_id(trait_type)) {
  trait_type <- match.arg(trait_type)
  trait_data <- get_trait_list(annotation_list, trait_type)
  if(is.null(trait_data) || !(trait_id %in% colnames(trait_data)) || !("symbol" %in% colnames(trait_data))){
      warning("Required columns (", trait_id, ", symbol) not found for trait type: ", trait_type)
      return(character(0))
  }
  
  # Remove "ENSMUSG" (genes) or "ENSMUST" (transcripts) and leading zeros
  id <- stringr::str_remove(trait_data[[trait_id]], pattern = "ENSMUS[GT]0+")
  paste(trait_data$symbol, id, sep = "_")
}

#' Get Trait Choices for Selectize Input
#' 
#' Generates the choices list for the trait selection input based on dataset.
#' 
#' @param file_directory The file directory data frame.
#' @param annotation_list The list containing all annotation data frames.
#' @param selected_dataset The currently selected dataset string.
#' @return A character vector of choices suitable for selectizeInput.
get_trait_choices <- function(file_directory, annotation_list, selected_dataset = NULL) {
  trait_type <- get_trait_type(file_directory, selected_dataset)
  if(is.null(trait_type)) {
    return(NULL)
  }
  
  trait_list_df <- get_trait_list(annotation_list, trait_type)
  if(is.null(trait_list_df)){
      return(NULL)
  }
  
  trait_id <- get_trait_id(trait_type)
  
  if(!(trait_type %in% c("genes", "isoforms"))) {
      # For non-gene/isoform types, use the ID directly
      if(!(trait_id %in% colnames(trait_list_df))){
          warning("Trait ID column ", trait_id, " not found in annotation for type: ", trait_type)
          return(NULL)
      }
    choices <- unique(trait_list_df[[trait_id]])
    # Removed explicit list structure
  } else { # "genes", "isoforms"
    # For genes/isoforms, create Symbol_ID for value using helper
    choices <- join_symbol_id(annotation_list, trait_type, trait_id)
    # Removed explicit list structure
  }
  
  sort(unique(choices)) # Return sorted unique vector
}

#' Get Selected Trait Symbol/Name
#' 
#' Extracts the actual trait name (symbol or full name) from the selectize input value,
#' handling the special "Symbol_ID" format for genes/isoforms.
#' 
#' @param file_directory The file directory data frame.
#' @param which_trait The value from the trait selectizeInput.
#' @param selected_dataset The currently selected dataset string.
#' @return The extracted trait name string.
#' @importFrom stringr str_split str_remove
get_selected_trait <- function(file_directory, which_trait, selected_dataset) {
  if(is.null(which_trait) || which_trait == "") return(NULL)
  trait_type <- get_trait_type(file_directory, selected_dataset)
  if(is.null(trait_type)) return(which_trait) # Return original if type unknown
  
  # If gene or isoform, split Symbol_ID to get Symbol
  if(trait_type %in% c("genes", "isoforms")) {
      # Use str_remove to handle cases where symbol itself might contain underscore
      selected_trait_name <- stringr::str_remove(which_trait, "_[^_]+$") 
  } else {
      selected_trait_name <- which_trait
  }
  selected_trait_name
}

#' Get Highest Peaks Above Threshold
#' 
#' Filters a peak table for peaks above a LOD threshold and orders them by LOD score.
#' 
#' @param peak_table A data frame containing peak information with 'lod' and 'marker' columns.
#' @param LOD_thr The LOD threshold value.
#' @return A data frame of filtered and sorted peaks, or NULL if no peaks meet the criteria.
#' @importFrom dplyr arrange desc filter
#' @importFrom rlang .data
highest_peaks <- function(peak_table, LOD_thr) {
  if(is.null(peak_table) || nrow(peak_table) == 0 || !("lod" %in% colnames(peak_table))){
      return(NULL)
  }
  # Ensure LOD_thr is numeric
  if(!is.numeric(LOD_thr)) {
      warning("LOD_thr must be numeric.")
      return(peak_table) # Return unfiltered if threshold invalid
  }
  filtered_peaks <- dplyr::filter(peak_table, .data$lod >= LOD_thr) |>
    dplyr::arrange(dplyr::desc(.data$lod))
    
  if (nrow(filtered_peaks) == 0) return(NULL)
  
  filtered_peaks
}

#' Pivot Peak Data for Allele Effects Plot
#' 
#' Selects the chosen peak and reshapes the founder allele effect columns (A-H)
#' into long format suitable for ggplot.
#' 
#' @param peaks A data frame of peaks, potentially including A-H allele effect columns.
#' @param which_peak The marker ID string of the selected peak.
#' @return A data frame in long format with 'marker', 'variable' (Strain), 
#'         and 'value' (Effect) columns, or NULL if data is missing or not additive.
#' @importFrom reshape2 melt
pivot_peaks <- function(peaks, which_peak) {
  if(is.null(peaks) || nrow(peaks) == 0 || is.null(which_peak) || which_peak == ""){
      return(NULL)
  }
  # Check if the scan type allows for allele effects (assuming additive has A-H columns)
  required_cols <- c("marker", "A","B","C","D","E","F","G","H")
  if (!all(required_cols %in% colnames(peaks))) {
      message("Allele effect columns (A-H) not found. Cannot create strain effects plot.")
      return(NULL)
  }
  
  # Filter for the selected peak
  peak <- subset(peaks, marker == which_peak)
  
  # Check if we have data after filtering
  if (nrow(peak) > 0) {
    # Select and rename columns
    peak_effects <- peak[1, required_cols] # Take only the first row if multiple matches
    colnames(peak_effects)[2:9] <- c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")
    # Reshape data
    peak_long <- reshape2::melt(peak_effects, id.vars = "marker", variable.name = "variable", value.name = "value")
    return(peak_long)
  } else {
    warning("Selected peak marker (", which_peak, ") not found in peaks table.")
    return(NULL)
  }
}

#' Convert Chromosome Number to Character (X, Y, M)
#' 
#' Converts numeric chromosome identifiers 20, 21, 22 to "X", "Y", "M".
#' 
#' @param chr A numeric or character vector of chromosome identifiers.
#' @return A character vector with 20, 21, 22 replaced by "X", "Y", "M".
chr_XYM <- function(chr) {
  if(is.null(chr)) return(NA)
  chr_char <- as.character(chr)
  chr_char[chr_char == "20"] <- "X"
  chr_char[chr_char == "21"] <- "Y"
  chr_char[chr_char == "22"] <- "M"
  return(chr_char)
} 