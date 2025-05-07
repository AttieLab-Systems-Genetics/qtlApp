#' Find the peaks
#'
#' @param file_dir data frame with file directory information
#' @param selected_dataset character string
#' @param selected_trait character string
#' @param cache_env environment to cache results (optional)
#'
#' @importFrom dplyr across mutate relocate where
#' @export
peak_finder <- function(file_dir, selected_dataset, selected_trait = NULL, cache_env = NULL) {
  cache_key <- if (is.null(selected_trait)) {
    selected_dataset  # Key for all peaks in dataset
  } else {
    paste(selected_dataset, tolower(selected_trait), sep = "_") # Key for specific trait within dataset
  }

  # Check cache only if cache_env is provided
  if (!is.null(cache_env) && !is.null(cache_env[[cache_key]])) {
     message("Using cached peaks data for ", 
            if (is.null(selected_trait)) selected_dataset else paste(selected_trait, "in", selected_dataset))
     return(cache_env[[cache_key]])
  }
  
  message("Loading peaks data for dataset: ", selected_dataset, 
          if(!is.null(selected_trait)) paste("(for trait:", selected_trait, ")") else "(all traits)")
  empty_peaks <- create_empty_peaks()
  
  # Filter file_dir for the peaks file for the selected dataset
  # IMPORTANT: peak_finder should read the WHOLE dataset's peaks file, 
  #            filtering by selected_trait happens LATER if needed (e.g., in peakServer)
  #            OR if selected_trait is passed explicitly to peak_finder for that purpose.
  #            The current logic reads the whole file if selected_trait=NULL (correct for cisTransPlot) 
  #            or if caching is off/misses.
  #            If selected_trait is provided AND cache is off/misses, it reads the whole file then filters.
  
  file_dir_filtered <- subset(file_dir, group == selected_dataset & file_type == "peaks")
  
  if (nrow(file_dir_filtered) == 0) {
    warning("No peaks files found for dataset: ", selected_dataset)
    # Cache the empty result if cache_env is provided
    if (!is.null(cache_env)) cache_env[[cache_key]] <- empty_peaks 
    return(empty_peaks)
  }
  
  csv_path <- file_dir_filtered$File_path[1]
  if (!file.exists(csv_path)) {
    warning("Peaks file does not exist: ", csv_path)
    if (!is.null(cache_env)) cache_env[[cache_key]] <- empty_peaks
    return(empty_peaks)
  }
  
  peaks_data_to_return <- empty_peaks
  tryCatch({
      message("Reading peaks file: ", basename(csv_path), " for dataset: ", selected_dataset)
      peaks_data <- data.table::fread(csv_path)
      message("Original peaks file columns: ", paste(colnames(peaks_data), collapse=", "))
      
      # Check for required columns and standardize names
      if ("lodcolumn" %in% colnames(peaks_data)) {
        peaks_data$trait <- peaks_data$lodcolumn
        message("Added trait column based on lodcolumn")
        message(paste("peak_finder: First few unique values in original 'lodcolumn':", paste(head(unique(peaks_data$lodcolumn), 10), collapse=", ")))
      }
      if (any(grepl("phenotype", tolower(colnames(peaks_data))))) {
        phenotype_col <- grep("phenotype", tolower(colnames(peaks_data)), value = TRUE)[1]
        peaks_data$trait <- peaks_data[[phenotype_col]]
        message("Added trait column based on phenotype column")
        message(paste("peak_finder: First few unique values in 'trait' after phenotype check:", paste(head(unique(peaks_data$trait), 10), collapse=", ")))
      }
      
      # Add a message here to show the state of the 'trait' column before it might be further filtered by selected_trait
      if ("trait" %in% colnames(peaks_data)) {
        message(paste("peak_finder: Final unique values in 'trait' column before trait-specific filtering (if any):", paste(head(unique(peaks_data$trait), 20), collapse=", ")))
      } else {
        message("peak_finder: 'trait' column was NOT created from lodcolumn or phenotype.")
      }

      # Filter for the specific trait ONLY if selected_trait was provided to this function call
      if (!is.null(selected_trait)) {
          message("Filtering for trait: ", selected_trait, " in dataset: ", selected_dataset)
          rows_to_keep <- rep(FALSE, nrow(peaks_data))
          if ("trait" %in% colnames(peaks_data)) { rows_to_keep <- rows_to_keep | (tolower(peaks_data$trait) == tolower(selected_trait)) }
          if ("lodcolumn" %in% colnames(peaks_data)) { rows_to_keep <- rows_to_keep | (tolower(peaks_data$lodcolumn) == tolower(selected_trait)) }
          phenotype_cols <- grep("phenotype", tolower(colnames(peaks_data)), value = TRUE)
          for (col in phenotype_cols) { rows_to_keep <- rows_to_keep | (tolower(peaks_data[[col]]) == tolower(selected_trait)) }
          peaks_data <- peaks_data[rows_to_keep]
          message("Found ", nrow(peaks_data), " peaks for specific trait: ", selected_trait)
      } 
      # Ensure we have all required columns
      # Check and rename columns if needed
      col_mapping <- list(
        marker = c("marker", "markers", "id", "snp", "SNP"),
        chr = c("chr", "chrom", "chromosome"),
        pos = c("pos", "position", "bp", "location"),
        lod = c("lod", "LOD", "score", "pvalue")
      )
      required_cols <- names(col_mapping)
      # Try to map columns
      for (req_col in required_cols) {
        if (!(req_col %in% colnames(peaks_data))) {
          # Try to find a matching column
          for (alt_name in col_mapping[[req_col]]) {
            if (alt_name %in% colnames(peaks_data)) {
              # Rename the column
              message("Renaming column '", alt_name, "' to '", req_col, "'")
              data.table::setnames(peaks_data, alt_name, req_col)
              break
            }
          }
        }
      }
      # Check if we have all required columns
      if (!all(required_cols %in% colnames(peaks_data))) {
        missing_cols <- required_cols[!(required_cols %in% colnames(peaks_data))]
        warning("Missing required columns in peaks file: ", paste(missing_cols, collapse=", "))
        warning("Missing required columns... returning empty peaks.")
        if (!is.null(cache_env)) cache_env[[cache_key]] <- empty_peaks
        return(empty_peaks)
      }
      # Print final column names for debug
      message("Final peaks data columns: ", paste(colnames(peaks_data), collapse=", "))
      # Convert to data.frame for compatibility
      peaks_data <- as.data.frame(peaks_data)
      # Remove unnecessary column `Which_mice` if present
      peaks_data[["Which_mice"]] <- NULL
      # Round off numeric columns to 5 significant digits
      peaks_data <- dplyr::mutate(peaks_data,
        dplyr::across(dplyr::where(is.numeric), function(x) signif (x, 5)))
      # Assign to the variable to be returned
      peaks_data_to_return <- peaks_data
          
    }, error = function(e) {
      warning("Error reading/processing peaks file ", csv_path, ": ", e$message)
      # Return empty frame on error (already initialized)
    })

  # Cache the result (either processed data or empty frame) if cache_env is provided
  if (!is.null(cache_env)) {
    cache_env[[cache_key]] <- peaks_data_to_return
  }
  
  return(peaks_data_to_return)
}
create_empty_peaks <- function() {
  data.frame(
    marker = character(0),
    trait = character(0),
    chr = character(0),
    pos = numeric(0),
    lod = numeric(0),
    A = numeric(0),
    B = numeric(0),
    C = numeric(0),
    D = numeric(0),
    E = numeric(0),
    F = numeric(0),
    G = numeric(0),
    H = numeric(0)
  )
}
