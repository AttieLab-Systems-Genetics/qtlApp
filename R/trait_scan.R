#' Find trait scan
#' 
#' @param file_dir data frame with file directory information
#' @param selected_dataset character string
#' @param selected_trait character string
#' @param cache_env environment to store cached results
#' 
#' @importFrom fst read_fst
#' @importFrom stringr str_detect str_replace
#' @importFrom data.table rbindlist setnames
#' @export
trait_scan <- function(file_dir, selected_dataset, selected_trait, cache_env = NULL) {
  # Trim selected_trait at the very beginning of the function
  selected_trait <- trimws(selected_trait)
  message("trait_scan: Processing trait '", selected_trait, "' for dataset '", selected_dataset, "'")
  
  # Check cache first
  cache_key <- paste(selected_dataset, tolower(selected_trait), sep = "_")
  if (!is.null(cache_env) && !is.null(cache_env[[cache_key]])) {
    message("Using cached data for trait: ", selected_trait, " in dataset: ", selected_dataset)
    return(cache_env[[cache_key]])
  }
  
  # Filter for scan files in the selected dataset
  file_dir <- subset(file_dir, group == selected_dataset & file_type == "scans")
  if (nrow(file_dir) == 0) {
    stop("No matching files found for the selected dataset: ", selected_dataset)
  }
  
  message("Processing ", nrow(file_dir), " scan files for dataset: ", selected_dataset)
  all_data <- list()
  
  for (i in 1:nrow(file_dir)) {
    chr_num <- file_dir$ID_code[i]
    original_fst_path <- file_dir$File_path[i]
    trait_type <- tolower(file_dir$trait_type[i])
    
    # Correct file path based on trait type
    corrected_fst_path <- correct_file_path(original_fst_path, trait_type)
    
    if (!file.exists(corrected_fst_path)) {
      warning("File not found, skipping: ", corrected_fst_path)
      next
    }
    
    # Ensure file is FST format
    fst_path <- ensure_fst_format(corrected_fst_path)
    if (is.null(fst_path)) {
      warning("Could not process file format: ", corrected_fst_path)
      next
    }
    
    # Find or create row index file
    row_index_path <- get_or_create_row_index(fst_path)
    if (is.null(row_index_path)) {
      message("No valid index file found for: ", basename(fst_path), ". Skipping.")
      next
    }
    
    # Process the trait data
    trait_data <- process_trait_from_file(fst_path, row_index_path, selected_trait, chr_num)
    if (!is.null(trait_data) && nrow(trait_data) > 0) {
      all_data[[length(all_data) + 1]] <- trait_data
    }
  }
  
  if (length(all_data) == 0) {
    stop("Trait '", selected_trait, "' not found in any chromosome for dataset: ", selected_dataset)
  }
  
  combined_data <- data.table::rbindlist(all_data, fill = TRUE)
  message("Combined data: ", nrow(combined_data), " rows for trait: ", selected_trait)
  
  # Cache the result
  if (!is.null(cache_env)) {
    cache_env[[cache_key]] <- combined_data
  }
  
  return(combined_data)
}

# Helper function to correct file paths based on trait type
correct_file_path <- function(original_path, trait_type) {
  if (is.na(trait_type) || !nzchar(trait_type)) {
    return(original_path)
  }
  
  processed_trait_type <- trait_type
  if (trait_type == "clinical traits") {
    processed_trait_type <- "clinical"
  }
  
  corrected_path <- original_path
  
  if (processed_trait_type == "clinical" || processed_trait_type == "liver_lipids") {
    if (grepl("_with_symbols\\.fst$", original_path)) {
      corrected_path <- sub("_with_symbols\\.fst$", "_processed.fst", original_path)
    } else if (!grepl("_processed\\.fst$", original_path)) {
      corrected_path <- paste0(tools::file_path_sans_ext(original_path), "_processed.fst")
    }
  } else if (processed_trait_type %in% c("genes", "isoforms")) {
    if (grepl("_processed\\.fst$", original_path)) {
      corrected_path <- sub("_processed\\.fst$", "_with_symbols.fst", original_path)
    } else if (!grepl("_with_symbols\\.fst$", original_path)) {
      corrected_path <- paste0(tools::file_path_sans_ext(original_path), "_with_symbols.fst")
    }
  }
  
  return(corrected_path)
}

# Helper function to ensure FST format
ensure_fst_format <- function(file_path) {
  if (stringr::str_detect(file_path, "fst$")) {
    return(file_path)
  }
  
  if (stringr::str_detect(file_path, "csv$")) {
    fst_path <- stringr::str_replace(file_path, "csv$", "fst")
    if (file.exists(fst_path)) {
      message("Switched from CSV to FST: ", basename(fst_path))
      return(fst_path)
    }
  }
  
  return(NULL)
}

# Helper function to get or create row index
get_or_create_row_index <- function(fst_path) {
  index_path_new <- sub("\\.fst$", "_rows.fst", fst_path)
  index_path_legacy <- sub("\\.fst$", "_row.fst", fst_path)
  
  if (file.exists(index_path_new)) {
    return(index_path_new)
  } else if (file.exists(index_path_legacy)) {
    return(index_path_legacy)
  } else {
    # Try to generate index on-the-fly
    tryCatch({
      row_index_path <- fst_rows(fst_path)
      if (file.exists(row_index_path)) {
        message("Generated index file: ", basename(row_index_path))
        return(row_index_path)
      }
    }, error = function(e) {
      warning("Error creating index file for ", basename(fst_path), ": ", e$message)
    })
  }
  
  return(NULL)
}

# Helper function to process trait data from a file
process_trait_from_file <- function(fst_path, row_index_path, selected_trait, chr_num) {
  tryCatch({
    # Read the row index to find the trait
    trait_index <- fst::read_fst(row_index_path, as.data.table = TRUE)
    trait_index[, Phenotype := tolower(Phenotype)]
    trait_rows <- trait_index[Phenotype == tolower(selected_trait), ]
    
    if (nrow(trait_rows) == 0) {
      return(NULL)
    }
    
    # Handle both old (from/to) and new (.row_min/.row_max) column naming
    if ("from" %in% colnames(trait_rows) && "to" %in% colnames(trait_rows)) {
      from_row <- trait_rows$from
      to_row <- trait_rows$to
    } else if (".row_min" %in% colnames(trait_rows) && ".row_max" %in% colnames(trait_rows)) {
      from_row <- trait_rows$.row_min
      to_row <- trait_rows$.row_max
    } else {
      warning("Row index file has unexpected column names for chromosome ", chr_num)
      return(NULL)
    }
    
    message("Found trait in chromosome ", chr_num, " at rows ", from_row, "-", to_row)
    
    # Read only the rows for this trait
    data <- fst::read_fst(fst_path,
      from = from_row,
      to = to_row,
      as.data.table = TRUE)
    
    # Ensure required columns are present
    data <- ensure_required_columns(data, fst_path)
    if (is.null(data)) return(NULL)
    
    # Filter by phenotype if column exists
    if ("Phenotype" %in% colnames(data)) {
      data <- data[tolower(Phenotype) == tolower(selected_trait)]
    }
    
    if (nrow(data) > 0) {
      message("Adding ", nrow(data), " rows from chromosome ", chr_num)
      return(data)
    }
    
  }, error = function(e) {
    warning("Error processing chromosome ", chr_num, ": ", e$message)
  })
  
  return(NULL)
}

# Helper function to ensure required columns exist
ensure_required_columns <- function(data, file_path) {
  # Check for LOD column
  if (!"LOD" %in% colnames(data)) {
    possible_lod_cols <- grep("lod|LOD|score", colnames(data), ignore.case = TRUE, value = TRUE)
    if (length(possible_lod_cols) > 0) {
      data.table::setnames(data, possible_lod_cols[1], "LOD")
    } else {
      warning("LOD column not found in file: ", file_path)
      return(NULL)
    }
  }
  
  # Check for marker column
  if (!"marker" %in% colnames(data)) {
    possible_marker_cols <- grep("marker|id|snp", colnames(data), ignore.case = TRUE, value = TRUE)
    if (length(possible_marker_cols) > 0) {
      data.table::setnames(data, possible_marker_cols[1], "marker")
    } else {
      warning("marker column not found in file: ", file_path)
      return(NULL)
    }
  }
  
  return(data)
}
