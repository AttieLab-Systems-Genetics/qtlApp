#' Find the peaks
#'
#' @param file_dir data frame with file directory information
#' @param selected_dataset character string
#' @param selected_trait character string
#'
#' @importFrom dplyr across mutate relocate where
#' @export
peak_finder <- function(file_dir, selected_dataset, selected_trait = NULL) {
  # Create a unique cache key that includes the dataset
  cache_key <- if(is.null(selected_trait)) {
    selected_dataset  # Include full dataset name
  } else {
    paste(selected_dataset, tolower(selected_trait), sep = "_")  # Include full dataset name
  }
  # Check if we already have this data in cache
  if(is.null(peaks_cache[[cache_key]])) {
    message("Loading peaks data for dataset: ", selected_dataset)
    empty_peaks <- create_empty_peaks()
    # Filter the file directory for the selected dataset and peaks files
    file_dir <- subset(file_dir, group == selected_dataset & file_type == "peaks")
    # Check if we found any matching files
    if(nrow(file_dir) == 0) {
      warning("No peaks files found for dataset: ", selected_dataset)
      # Return empty data frame with correct structure
      peaks_cache[[cache_key]] <- empty_peaks
      return(empty_peaks)
    }
    # We now have a single consolidated peaks file
    message("Reading consolidated peaks file for dataset: ", selected_dataset)
    csv_path <- file_dir$File_path[1]
    if(!file.exists(csv_path)) {
      warning("Peaks file does not exist: ", csv_path)
      peaks_cache[[cache_key]] <- empty_peaks
      return(empty_peaks)
    }
    tryCatch({
      # Read the consolidated CSV file - this might be large, so use data.table for efficiency
      message("Reading peaks file: ", basename(csv_path), " for dataset: ", selected_dataset)
      peaks_data <- data.table::fread(csv_path)
      # Print out column names for debug
      message("Original peaks file columns: ", paste(colnames(peaks_data), collapse=", "))
      # Check for required columns and standardize names
      if("lodcolumn" %in% colnames(peaks_data)) {
        # We'll keep lodcolumn as is but also add a trait column for compatibility
        peaks_data$trait <- peaks_data$lodcolumn
        message("Added trait column based on lodcolumn")
      }
      if(any(grepl("phenotype", tolower(colnames(peaks_data))))) {
        phenotype_col <- grep("phenotype", tolower(colnames(peaks_data)), value = TRUE)[1]
        # We'll keep original column and add a trait column
        peaks_data$trait <- peaks_data[[phenotype_col]]
        message("Added trait column based on phenotype column")
      }
      # Filter for the specific trait if provided (case-insensitive)
      if(!is.null(selected_trait)) {
        message("Filtering for trait: ", selected_trait, " in dataset: ", selected_dataset)
        # Try various columns that might contain the trait
        rows_to_keep <- rep(FALSE, nrow(peaks_data))
        # Check trait column
        if("trait" %in% colnames(peaks_data)) {
          message("Checking trait column")
          rows_to_keep <- rows_to_keep | (tolower(peaks_data$trait) == tolower(selected_trait))
        }
        # Check lodcolumn
        if("lodcolumn" %in% colnames(peaks_data)) {
          message("Checking lodcolumn column")
          rows_to_keep <- rows_to_keep | (tolower(peaks_data$lodcolumn) == tolower(selected_trait))
        }
        # Check phenotype columns
        phenotype_cols <- grep("phenotype", tolower(colnames(peaks_data)), value = TRUE)
        for(col in phenotype_cols) {
          message("Checking phenotype column: ", col)
          rows_to_keep <- rows_to_keep | (tolower(peaks_data[[col]]) == tolower(selected_trait))
        }
        # Apply the filter
        peaks_data <- peaks_data[rows_to_keep]
        message("Found ", nrow(peaks_data), " peaks for trait: ", selected_trait,
          " in dataset: ", selected_dataset)
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
      for(req_col in required_cols) {
        if(!(req_col %in% colnames(peaks_data))) {
          # Try to find a matching column
          for(alt_name in col_mapping[[req_col]]) {
            if(alt_name %in% colnames(peaks_data)) {
              # Rename the column
              message("Renaming column '", alt_name, "' to '", req_col, "'")
              data.table::setnames(peaks_data, alt_name, req_col)
              break
            }
          }
        }
      }
      # Check if we have all required columns
      if(!all(required_cols %in% colnames(peaks_data))) {
        missing_cols <- required_cols[!(required_cols %in% colnames(peaks_data))]
        warning("Missing required columns in peaks file: ", paste(missing_cols, collapse=", "))
        # Return empty data frame
        peaks_cache[[cache_key]] <- empty_peaks
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
        dplyr::across(dplyr::where(is.numeric), function(x) signif(x, 5)))
      # Cache the result
      peaks_cache[[cache_key]] <- peaks_data
          
    }, error = function(e) {
      warning("Error reading peaks file ", csv_path, ": ", e$message)
      # Return empty data frame
      peaks_cache[[cache_key]] <- empty_peaks
    })
  } else {
    message("Using cached peaks data for ",
      if(is.null(selected_trait)) selected_dataset else paste(selected_trait, "in", selected_dataset))
  }
  return(peaks_cache[[cache_key]])
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