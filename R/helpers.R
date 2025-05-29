# =============================================================================
# HELPERS.R - Core Utility Functions
# =============================================================================
# Contains general-purpose helper functions used throughout the QTL mapping app.
# Functions are organized into categories: caching, trait selection, peak analysis,
# data reshaping, and chromosome handling.

# =============================================================================
# CACHING UTILITIES
# =============================================================================

#' Creates cache objects for storing peaks and traits
#' 
#' Sets up caching infrastructure to improve app performance by storing
#' frequently accessed data like peak results and trait information.
#' 
#' @param caches A character vector specifying which caches to create 
#'   (default: c("peaks", "trait")). Options include:
#'   - "peaks": Cache for QTL peak analysis results
#'   - "trait": Cache for trait-related data and metadata
#' @return Cache objects (implementation depends on caching strategy)
#' @export
#' @examples
#' # Create default caches for peaks and traits
#' create_cache()
#' 
#' # Create only peaks cache
#' create_cache("peaks")
create_cache <- function(caches = c("peaks", "trait")) {
  # Implementation needed - likely creates environment or list-based cache
}

# =============================================================================
# TRAIT SELECTION AND VALIDATION
# =============================================================================

#' Get available trait choices for a selected dataset
#' 
#' Extracts and validates trait options based on the dataset type.
#' Handles different trait types (genes, isoforms, clinical, lipids) and
#' ensures only valid, non-empty choices are returned.
#' 
#' @param import_data List containing imported data structure with:
#'   - file_directory: DataFrame with dataset metadata
#'   - annotation_list: Named list of annotation data by trait type
#' @param selected_dataset String identifying which dataset to get traits for
#' @return Character vector of available trait choices, sorted alphabetically.
#'   Returns empty vector if no valid choices found.
#' @export
#' @examples
#' # Get trait choices for a clinical dataset
#' traits <- get_trait_choices(import_data, "clinical_study_1")
get_trait_choices <- function(import_data, selected_dataset = NULL) {
  # Validate required inputs using Shiny's req() for reactive contexts
  shiny::req(import_data, selected_dataset)
  
  # Determine what type of traits this dataset contains
  trait_type <- get_trait_type(import_data, selected_dataset)
  if(is.null(trait_type)) {
    return(character(0))
  }

  # Get the annotation data for this trait type
  trait_list_df <- get_trait_list(import_data, trait_type)
  if(is.null(trait_list_df)){
      return(character(0))
  }

  choices <- character(0)
  
  # Handle gene/isoform data - use 'symbol' column for user-friendly names
  if(trait_type %in% c("genes", "isoforms")) {
    if("symbol" %in% colnames(trait_list_df)) {
      choices <- unique(trait_list_df$symbol)
    } else {
      return(character(0))
    }
  } else {
    # Handle other trait types (clinical, lipids, etc.)
    trait_id <- get_trait_id(trait_type)
    
    # Validate annotation data structure
    if(is.null(trait_list_df) || !is.data.frame(trait_list_df)){
        warning("Annotation data (trait_list_df) for type '", trait_type, "' is NULL or not a data frame.")
        return(character(0))
    }
    
    # Check if required ID column exists
    if(!(trait_id %in% colnames(trait_list_df))){
        warning("Trait ID column '", trait_id, "' not found in annotation data frame for type: ", trait_type, ". Available columns: ", paste(colnames(trait_list_df), collapse=", "))
        return(character(0))
    }
    
    # Extract unique trait identifiers
    choices <- unique(trait_list_df[[trait_id]])
    
    # Warn if column exists but contains no valid data
    valid_choices_count <- sum(!is.na(choices) & choices != "")
    if(valid_choices_count == 0 && length(choices) > 0){
        warning("Column '", trait_id, "' exists but contains only NAs or empty strings.")
    }
  }
  
  # Filter out invalid choices and return sorted results
  valid_choices <- choices[!is.na(choices) & choices != ""]
  sort(unique(valid_choices))
}

#' Extract and clean selected trait name
#' 
#' Simple utility to clean trait names by removing leading/trailing whitespace.
#' Used when processing user-selected traits from UI components.
#' 
#' @param import_data Imported data structure (required by shiny::req but not used)
#' @param which_trait String containing the selected trait name
#' @param selected_dataset String identifying the dataset (required by shiny::req but not used)
#' @return Trimmed trait name string
#' @export
#' @examples
#' clean_trait <- get_selected_trait(import_data, "  gene_name  ", "dataset1")
#' # Returns: "gene_name"
get_selected_trait <- function(import_data, which_trait, selected_dataset) {
  shiny::req(which_trait)
  trimmed_trait <- trimws(which_trait)
  return(trimmed_trait)
}

# =============================================================================
# PEAK ANALYSIS UTILITIES  
# =============================================================================

#' Filter and sort peaks by LOD threshold
#' 
#' Identifies significant QTL peaks that meet or exceed a specified LOD score
#' threshold. Results are sorted by LOD score in descending order.
#' 
#' @param peak_table DataFrame containing QTL peak results with required columns:
#'   - lod: LOD score for each peak
#'   - (other columns preserved in output)
#' @param LOD_thr Numeric threshold for minimum LOD score
#' @return DataFrame of filtered peaks sorted by descending LOD score,
#'   or NULL if no peaks meet threshold
#' @importFrom dplyr arrange desc filter
#' @importFrom rlang .data
#' @export
#' @examples
#' # Get peaks with LOD >= 3.0
#' significant_peaks <- highest_peaks(peak_data, LOD_thr = 3.0)
highest_peaks <- function(peak_table, LOD_thr) {
  # Filter peaks that meet or exceed the LOD threshold
  filtered_peaks <- dplyr::filter(peak_table, .data$lod >= LOD_thr) |>
    dplyr::arrange(dplyr::desc(.data$lod))
  
  # Return NULL if no significant peaks found
  if (nrow(filtered_peaks) == 0) return(NULL)
  filtered_peaks
}

# =============================================================================
# DATA RESHAPING UTILITIES
# =============================================================================

#' Reshape peak data for visualization
#' 
#' Converts wide-format allele effect data to long format suitable for plotting.
#' Transforms founder strain allele columns (A-H) into strain-specific rows
#' with corresponding effect values.
#' 
#' @param peaks DataFrame containing peak data with:
#'   - marker: Marker identifier column
#'   - A, B, C, D, E, F, G, H: Allele effect columns for 8 founder strains
#' @param which_peak String specifying which marker to reshape
#' @return Long-format DataFrame with columns:
#'   - marker: Original marker identifier  
#'   - variable: Strain name (AJ, B6, 129, NOD, NZO, CAST, PWK, WSB)
#'   - value: Allele effect value
#'   Returns NULL if input validation fails.
#' @importFrom reshape2 melt
#' @export
#' @examples
#' # Reshape peak data for marker "rs123456"
#' long_data <- pivot_peaks(peak_results, "rs123456")
pivot_peaks <- function(peaks, which_peak) {
  peak_reshaped <- NULL
  
  # Validate input data structure
  if (is.null(peaks) || !is.data.frame(peaks) || nrow(peaks) == 0) {
    warning("pivot_peaks: Input 'peaks' data is NULL, not a data.frame, or has 0 rows.")
    return(peak_reshaped)
  }
  
  # Validate peak identifier
  if (is.null(which_peak) || which_peak == "" || is.na(which_peak)){
    warning("pivot_peaks: 'which_peak' (marker for filtering) is NULL, empty or NA.")
    return(peak_reshaped)
  }
  
  # Check for required marker column
  if (!("marker" %in% colnames(peaks))) {
    warning("pivot_peaks: 'marker' column not found in peaks data.")
    return(peak_reshaped)
  }
  
  # Filter to specific peak/marker
  peak_subset <- subset(peaks, marker == which_peak)
  
  if (nrow(peak_subset) > 0) {
    # Define mapping between allele columns and strain names
    allele_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
    strain_names <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
    
    # Verify all required allele columns are present
    if (all(allele_cols %in% colnames(peak_subset))) {
      # Select relevant columns for reshaping
      peak_for_melt <- peak_subset[, c("marker", allele_cols), drop = FALSE]
      
      # Rename allele columns to strain names for clearer output
      if (length(allele_cols) == length(strain_names)) {
          colnames(peak_for_melt)[(which(colnames(peak_for_melt) %in% allele_cols))] <- strain_names
      } else {
          warning("pivot_peaks: Mismatch between allele_cols and strain_names length. Using original allele column names.")
          strain_names <- allele_cols
      }
      
      # Reshape from wide to long format
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

# =============================================================================
# CHROMOSOME HANDLING UTILITIES
# =============================================================================
# Functions for converting between chromosome labels and numeric values
# Handles standard autosomes (1-19) plus sex chromosomes (X,Y) and mitochondrial (M)

#' Convert chromosome label to numeric value
#' 
#' Standardizes chromosome identifiers by converting character labels
#' to numeric values. Handles sex chromosomes and mitochondrial DNA.
#' 
#' @param chr Character or numeric chromosome label(s). Can be:
#'   - Numeric: 1-19 (autosomes)  
#'   - Character: "X", "Y", "M" (sex chromosomes, mitochondrial)
#' @return Numeric vector where:
#'   - 1-19: Autosomal chromosomes
#'   - 20: X chromosome
#'   - 21: Y chromosome  
#'   - 22: Mitochondrial chromosome
#'   - NA: Invalid input
#' @export
#' @examples
#' chr_to_numeric(c("1", "X", "Y", "M"))
#' # Returns: c(1, 20, 21, 22)
chr_to_numeric <- function(chr) {
  if (is.null(chr)) return(NA_real_)
  chr_char <- as.character(chr)
  
  # Handle vectors properly - initialize result vector
  result <- rep(NA_real_, length(chr_char))
  
  # Map special chromosomes to numeric codes
  result[chr_char == "X"] <- 20
  result[chr_char == "Y"] <- 21
  result[chr_char == "M"] <- 22
  
  # For standard autosomes, convert directly to numeric
  numeric_mask <- !chr_char %in% c("X", "Y", "M")
  result[numeric_mask] <- as.numeric(chr_char[numeric_mask])
  
  return(result)
}

#' Convert numeric chromosome value to label
#' 
#' Reverse of chr_to_numeric - converts numeric chromosome codes back
#' to standard chromosome labels for display purposes.
#' 
#' @param chr_num Numeric chromosome value(s):
#'   - 1-19: Autosomal chromosomes
#'   - 20: X chromosome  
#'   - 21: Y chromosome
#'   - 22: Mitochondrial chromosome
#' @return Character vector of chromosome labels
#' @export
#' @examples
#' numeric_to_chr(c(1, 20, 21, 22))
#' # Returns: c("1", "X", "Y", "M")
numeric_to_chr <- function(chr_num) {
  if (is.null(chr_num)) return(NA_character_)
  
  # Handle vectors properly - start with character conversion
  result <- as.character(chr_num)
  
  # Map special numeric codes back to chromosome labels
  result[chr_num == 20] <- "X"
  result[chr_num == 21] <- "Y"
  result[chr_num == 22] <- "M"
  
  return(result)
}

#' Legacy function for chromosome label conversion
#' 
#' Maintained for backward compatibility. Use numeric_to_chr() instead.
#' 
#' @param chr_vector Numeric chromosome values to convert
#' @return Character vector of chromosome labels
#' @export
chr_XYM <- function(chr_vector) {
  numeric_to_chr(chr_vector)
}

# =============================================================================
# INTERNAL HELPER FUNCTIONS
# =============================================================================
# These functions support the main exported functions above and handle
# data structure navigation and validation.

#' Determine trait type for a selected dataset
#' 
#' Internal function that examines the file directory to determine what
#' type of traits are available in the selected dataset.
#' 
#' @param import_data List containing file_directory DataFrame
#' @param selected_dataset String identifying the dataset
#' @return String indicating trait type: "clinical", "lipids", "genes", 
#'   "isoforms", or original value if no pattern matches
get_trait_type <- function(import_data, selected_dataset = NULL) {
  # Validate input structure
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
  
  # Check for trait_type column
  if(!("trait_type" %in% colnames(file_dir_subset))) {
     warning("'trait_type' column missing from file directory.")
     return(NULL)
  }
  trait_type_raw <- tolower(file_dir_subset$trait_type[1])
  
  # Standardize trait type names using pattern matching
  if (grepl("clinical", trait_type_raw, ignore.case = TRUE)) { 
      return("clinical")
  }
  if (grepl("lipid", trait_type_raw, ignore.case = TRUE)) { 
      return("lipids")
  }
  if (grepl("gene", trait_type_raw, ignore.case = TRUE)) { 
      return("genes")
  }
  if (grepl("isoform", trait_type_raw, ignore.case = TRUE)) { 
      return("isoforms")
  }
  
  # Return original (lowercased) if no specific pattern matched
  return(trait_type_raw)
}

#' Retrieve annotation data for a specific trait type
#' 
#' Internal function that extracts the appropriate annotation DataFrame
#' from the imported data structure based on trait type.
#' 
#' @param import_data List containing annotation_list
#' @param trait_type String specifying which annotation set to retrieve
#' @return DataFrame containing annotation data, or NULL if not found
get_trait_list <- function(import_data, trait_type) {
  # Validate input structure
  if(is.null(import_data) || is.null(import_data$annotation_list)){
      warning("get_trait_list: import_data or import_data$annotation_list is NULL.")
      return(NULL)
  }
  annotation_list <- import_data$annotation_list
  
  # Check if requested trait type exists in annotations
  if(is.null(trait_type) || !(trait_type %in% names(annotation_list))) {
     warning("Trait type ", trait_type, " not found in annotation list.")
     return(NULL)
  }
  
  annotation_list[[trait_type]]
}

#' Get the appropriate ID column name for a trait type
#' 
#' Internal function that maps trait types to their corresponding
#' identifier column names in annotation data.
#' 
#' @param trait_type String specifying the trait type
#' @return String containing the appropriate column name for trait IDs
get_trait_id <- function(trait_type) {
  switch(trait_type,
    genes    = "gene.id",      # Gene-based studies use gene.id
    isoforms = "transcript_id", # Isoform studies use transcript_id  
    "data_name")               # Default for other types (clinical, lipids, etc.)
}

# =============================================================================
# END OF HELPERS.R
# =============================================================================