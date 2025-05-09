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
  
  cache_key <- paste(selected_dataset, tolower(selected_trait), sep = "_")
  # Check cache only if cache_env is provided
  if (!is.null(cache_env) && !is.null(cache_env[[cache_key]])) {
    message("Using cached data for trait: ", selected_trait, " in dataset: ", selected_dataset)
    return(cache_env[[cache_key]])
  }
  message("Searching for trait: ", selected_trait, " in dataset: ", selected_dataset)
  
  file_dir <- subset(file_dir, group == selected_dataset & file_type == "scans")
  # Check if we found any matching files
  if (nrow(file_dir) == 0) {
    stop("No matching files found for the selected dataset: ", selected_dataset)
  }
  
  
  all_data <- list()
  # Process each FST file (one per chromosome)
  for (i in 1:nrow(file_dir)) {
   
    chr_num <- file_dir$ID_code[i]
    fst_path <- file_dir$File_path[i]

    # No longer strictly necessary, but keep as a safeguard
    if (!file.exists(fst_path)) {
      warning("File check consistency error, skipping: ", fst_path)
      next
    }
    
    if (!stringr::str_detect(fst_path, "fst$")) {
      fst_path <- stringr::str_replace(fst_path, "csv$", "fst")
      if (!file.exists(fst_path)) {
        warning("FST file not found: ", fst_path)
        next
      }
    }
    message("Checking chromosome ", chr_num, " for trait: ", selected_trait, " in dataset: ", selected_dataset)
    # Create row index if it doesn't exist
    row_index_path <- fst_rows(fst_path)
    tryCatch({
      # Read the row index to find the trait
      trait_index <- fst::read_fst(row_index_path, as.data.table = TRUE)
      # Convert Phenotype column to lowercase for case-insensitive matching
      trait_index[, Phenotype := tolower(Phenotype)]
      # Check if the trait is present in this chromosome (case-insensitive)
      # ** this will have problems with duplicate genes or isoforms **
      trait_rows <- trait_index[Phenotype == tolower(selected_trait), ]
      if (nrow(trait_rows) > 0) {
        message("Found trait in chromosome ", chr_num, " at rows ", trait_rows$from, "-", trait_rows$to)
        # Read only the rows for this trait
        data <- fst::read_fst(fst_path,
          from = trait_rows$from,
          to = trait_rows$to,
          as.data.table = TRUE)
        # Ensure required columns are present
        if (!"LOD" %in% colnames(data)) {
          possible_lod_cols <- grep("lod|LOD|score", colnames(data), ignore.case = TRUE, value = TRUE)
          if (length(possible_lod_cols) > 0) {
            data.table::setnames(data, possible_lod_cols[1], "LOD")
          } else {
            warning("LOD column not found in file: ", fst_path)
            next
          }
        }
        if (!"marker" %in% colnames(data)) {
          possible_marker_cols <- grep("marker|id|snp", colnames(data), ignore.case = TRUE, value = TRUE)
          if (length(possible_marker_cols) > 0) {
            data.table::setnames(data, possible_marker_cols[1], "marker")
          } else {
            warning("marker column not found in file: ", fst_path)
            next
          }
        }
        # Verify that we have the correct trait data (case-insensitive)
        if ("Phenotype" %in% colnames(data)) {
            # Double-check that all rows are for the requested trait
            data <- data[tolower(Phenotype) == tolower(selected_trait)]
            message("Verified ", nrow(data), " rows for trait: ", selected_trait)
        }
        if (nrow(data) > 0) {
          message("Adding ", nrow(data), " rows from chromosome ", chr_num)
          all_data[[length(all_data) + 1]] <- data
        }
      } else {
        message("Trait not found in chromosome ", chr_num)
      }
    }, error = function(e) {
      warning("Error processing chromosome ", chr_num, ": ", e$message)
    })
  }
  
  if (length(all_data) == 0) {
    stop("Trait '", selected_trait, "' not found in any chromosome for dataset: ", selected_dataset)
  }
  
  combined_data <- data.table::rbindlist(all_data, fill = TRUE)
  message("Total rows in combined data: ", nrow(combined_data))

  
  if (!is.null(cache_env)) {
    cache_env[[cache_key]] <- combined_data
  }
  return(combined_data)
}
