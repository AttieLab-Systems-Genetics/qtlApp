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
    original_fst_path <- file_dir$File_path[i] # Path from file_index.csv
    
    # Get trait_type for this specific file row and convert to lowercase
    file_specific_trait_type_from_file <- tolower(file_dir$trait_type[i]) 

    # Standardize "clinical traits" to "clinical" for consistent logic
    processed_trait_type <- file_specific_trait_type_from_file
    if (file_specific_trait_type_from_file == "clinical traits") {
      processed_trait_type <- "clinical"
    }

    corrected_fst_path <- original_fst_path

    if (!is.na(processed_trait_type) && processed_trait_type != "") {
        if (processed_trait_type == "clinical") {
            if (grepl("_with_symbols\\.fst$", original_fst_path)) {
                corrected_fst_path <- sub("_with_symbols\\.fst$", "_processed.fst", original_fst_path)
                message(paste("Clinical type for file '", basename(original_fst_path), "', adjusted path to: '", basename(corrected_fst_path), "'"))
            } else if (!grepl("_processed\\.fst$", original_fst_path)) {
                # If it's clinical and not _processed and not _with_symbols, assume it's a base name and add _processed.fst
                corrected_fst_path <- paste0(tools::file_path_sans_ext(original_fst_path), "_processed.fst")
                message(paste("Clinical type for file '", basename(original_fst_path), "', assumed base, adjusted path to: '", basename(corrected_fst_path), "'"))
            }
        } else if (processed_trait_type %in% c("genes", "isoforms")) {
            if (grepl("_processed\\.fst$", original_fst_path)) {
                corrected_fst_path <- sub("_processed\\.fst$", "_with_symbols.fst", original_fst_path)
                message(paste("Gene/Isoform type for file '", basename(original_fst_path), "', adjusted path to: '", basename(corrected_fst_path), "'"))
            } else if (!grepl("_with_symbols\\.fst$", original_fst_path)) {
                # If it's gene/isoform and not _with_symbols and not _processed, assume it's a base name and add _with_symbols.fst
                corrected_fst_path <- paste0(tools::file_path_sans_ext(original_fst_path), "_with_symbols.fst")
                message(paste("Gene/Isoform type for file '", basename(original_fst_path), "', assumed base, adjusted path to: '", basename(corrected_fst_path), "'"))
            }
        } else if (processed_trait_type == "liver_lipids") {
            if (grepl("_with_symbols\\.fst$", original_fst_path)) {
                corrected_fst_path <- sub("_with_symbols\\.fst$", "_processed.fst", original_fst_path)
                message(paste("Liver Lipids type for file '", basename(original_fst_path), "', adjusted path from _with_symbols to: '", basename(corrected_fst_path), "'"))
            } else if (!grepl("_processed\\.fst$", original_fst_path)) {
                # If it's liver_lipids and not _processed and not _with_symbols, assume it's a base name and add _processed.fst
                corrected_fst_path <- paste0(tools::file_path_sans_ext(original_fst_path), "_processed.fst")
                message(paste("Liver Lipids type for file '", basename(original_fst_path), "', assumed base, adjusted path to: '", basename(corrected_fst_path), "'"))
            }
        } else {
            message(paste("Unknown or unspecified trait_type ('", processed_trait_type, "') for file '", basename(original_fst_path), "'. Using original path."))
        }
    } else {
        warning(paste("Missing or empty trait_type for file '", basename(original_fst_path), "'. Using original path. Please check file_index.csv."))
    }

    fst_path <- corrected_fst_path # Use the corrected (or original if no correction applied) path

    # --- ADD DEBUGGING --- 
    message(paste("DEBUG trait_scan: About to check file.exists() for fst_path:", fst_path))
    message(paste("DEBUG trait_scan: dput(fst_path) is:", paste(capture.output(dput(fst_path)), collapse=""))) # Shows exact R string
    parent_dir <- dirname(fst_path)
    message(paste("DEBUG trait_scan: Checking parent directory:", parent_dir))
    if(dir.exists(parent_dir)){
        message(paste("DEBUG trait_scan: Parent directory exists. Listing some files in parent_dir:"))
        tryCatch({
            print(utils::head(list.files(parent_dir)))
        }, error = function(e) { message(paste("DEBUG trait_scan: Error listing files:", e$message))})
    } else {
        message(paste("DEBUG trait_scan: PARENT DIRECTORY DOES NOT EXIST:", parent_dir))
    }
    # --- END DEBUGGING ---

    # Safeguard: Check if the (potentially corrected) file exists
    if (!file.exists(fst_path)) {
      warning("File check consistency error or corrected path invalid, skipping: ", fst_path, " (Original was: ", original_fst_path, ")")
      next
    }
    
    # This check for .fst extension might be redundant if corrections ensure .fst, but kept for safety
    if (!stringr::str_detect(fst_path, "fst$")) {
      # Attempt to replace .csv with .fst only if it's a .csv, otherwise log a warning
      if (stringr::str_detect(fst_path, "csv$")){
        fst_path_csv_replaced <- stringr::str_replace(fst_path, "csv$", "fst")
        if (!file.exists(fst_path_csv_replaced)) {
            warning("Original path was not FST, and replacing .csv with .fst also not found: ", fst_path_csv_replaced, " (Original non-FST path was: ", fst_path, ")")
            next
        } else {
            message("Path was CSV, switched to FST: ", fst_path_csv_replaced)
            fst_path <- fst_path_csv_replaced
        }
      } else {
        warning("File path does not end with .fst and is not .csv, skipping: ", fst_path)
        next
      }
    }
    message("Checking chromosome ", chr_num, " for trait: ", selected_trait, " in dataset: ", selected_dataset)
    
    # Determine potential index file paths
    # New convention: _rows.fst
    index_path_new <- sub("\\\\.fst$", "_rows.fst", fst_path)
    if (index_path_new == fst_path) { # Safety for paths not ending in .fst or if sub() fails
        index_path_new <- paste0(fst_path, "_rows.fst")
    }

    # Legacy convention: _row.fst (this is also what fst_rows() creates)
    index_path_legacy <- sub("\\\\.fst$", "_row.fst", fst_path)
    if (index_path_legacy == fst_path) { # Safety for paths not ending in .fst or if sub() fails
        index_path_legacy <- paste0(fst_path, "_row.fst")
    }

    row_index_path <- NULL # Initialize

    if (file.exists(index_path_new)) {
        message(paste("INFO: Using pre-existing NEWER convention index file:", basename(index_path_new), "(for FST:", basename(fst_path), ")"))
        row_index_path <- index_path_new
    } else if (file.exists(index_path_legacy)) {
        message(paste("INFO: Using pre-existing LEGACY convention index file:", basename(index_path_legacy), "(for FST:", basename(fst_path), ")"))
        row_index_path <- index_path_legacy
    } else {
        message(paste("INFO: Neither '_rows.fst' nor '_row.fst' index found for", basename(fst_path), ". Attempting to generate '_row.fst' on-the-fly using fst_rows()."))
        # Attempt to create the legacy _row.fst index file on the fly
        # fst_rows() will create it if it doesn't exist and return its path.
        tryCatch({
            row_index_path <- fst_rows(fst_path) # This will generate `_row.fst`
            if (file.exists(row_index_path)) {
                 message(paste("INFO: Successfully generated and using index file:", basename(row_index_path)))
            } else {
                 warning(paste("CRITICAL: fst_rows() was called but failed to create or return a valid path for index file:", basename(row_index_path), "(for FST:", basename(fst_path), "). Skipping this file."))
                 row_index_path <- NULL 
            }
        }, error = function(e_create) {
            warning(paste("CRITICAL: Error during on-the-fly creation of index file for", basename(fst_path), "using fst_rows():", e_create$message, ". Skipping this file."))
            row_index_path <- NULL 
        })
    }

    # If after all checks and potential creation, row_index_path is still NULL or doesn't exist, skip.
    if (is.null(row_index_path) || !file.exists(row_index_path)) {
        if(!is.null(row_index_path)) { 
             warning(paste("CRITICAL: Index file path was determined as '", basename(row_index_path), "' but the file does not exist. Skipping FST file:", basename(fst_path)))
        } else {
             message(paste("INFO: No valid index file found or generated for FST file:", basename(fst_path), ". Skipping processing for this file."))
        }
        next # Skip to the next file in the loop
    }

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

  message(paste("trait_scan DEBUG: Final combined_data for trait:", selected_trait))
  message("trait_scan DEBUG: str(combined_data):")
  str(combined_data)
  message("trait_scan DEBUG: print(head(combined_data)):")
  print(head(combined_data))
  message("trait_scan DEBUG: colnames(combined_data):")
  print(colnames(combined_data))
  if (nrow(combined_data) > 0 && "marker" %in% colnames(combined_data)) {
    message("trait_scan DEBUG: dput(head(combined_data$marker)):")
    dput(head(combined_data$marker))
    message("trait_scan DEBUG: Any NAs in combined_data$marker? ", any(is.na(combined_data$marker)))
    message("trait_scan DEBUG: Number of unique markers in combined_data: ", length(unique(combined_data$marker)))
  } else {
    message("trait_scan DEBUG: combined_data has 0 rows or no 'marker' column.")
  }

  # Save RDS for specific traits/datasets for external debugging
  # Ensure selected_dataset is available in this scope or passed if necessary
  # Assuming selected_dataset variable is accessible here from the function arguments or a higher scope
  if (exists("selected_dataset", inherits = FALSE)) { # Check if selected_dataset is in the local function scope
    if (tolower(selected_trait) == "bmp_18_2_22_6" || grepl("lipid", tolower(selected_dataset), ignore.case = TRUE)) {
        save_path <- paste0("/tmp/debug_trait_scan_output_", gsub("[^a-zA-Z0-9_.-]", "_", selected_dataset), "_", gsub("[^a-zA-Z0-9_.-]", "_", selected_trait), ".rds")
        message(paste("trait_scan DEBUG: SAVING combined_data for trait:", selected_trait, "from dataset:", selected_dataset, "TO:", save_path))
        tryCatch({
            saveRDS(combined_data, file = save_path)
            message("trait_scan DEBUG: Successfully saved to ", save_path)
        }, error = function(e) {
            message(paste("trait_scan DEBUG: FAILED to save RDS:", e$message))
        })
    }
  } else {
    message("trait_scan DEBUG: selected_dataset variable not found in trait_scan scope for RDS saving condition.")
  }

  return(combined_data)
}
