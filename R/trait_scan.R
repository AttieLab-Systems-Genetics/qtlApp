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
  selected_trait_trimmed <- trimws(selected_trait)
  message(paste0("--- trait_scan ENTERED --- Original Trait: '", selected_trait, "', Trimmed Trait: '", selected_trait_trimmed, "', Dataset: '", selected_dataset, "', Cache Used: ", !is.null(cache_env)))
  # Use selected_trait_trimmed throughout the function
  selected_trait <- selected_trait_trimmed 

  message("Initial file_dir nrow: ", if(!is.null(file_dir) && is.data.frame(file_dir)) nrow(file_dir) else "NULL or not a dataframe")
  if(!is.null(file_dir) && is.data.frame(file_dir) && nrow(file_dir) > 0) {
    message("Initial file_dir head:")
    print(utils::head(file_dir))
  }
  
  cache_key <- paste(selected_dataset, tolower(selected_trait), sep = "_")
  if (!is.null(cache_env) && !is.null(cache_env[[cache_key]])) {
    message("Using cached data for trait: ", selected_trait, " in dataset: ", selected_dataset)
    return(cache_env[[cache_key]])
  }
  message("Searching for trait: ", selected_trait, " in dataset: ", selected_dataset)
  
  file_dir <- subset(file_dir, group == selected_dataset & file_type == "scans")
  if (nrow(file_dir) == 0) {
    stop("No matching files found for the selected dataset: ", selected_dataset)
  }
  
  message(paste("DEBUG trait_scan: About to loop through", nrow(file_dir), "scan file entries for dataset:", selected_dataset))
  all_data <- list()
  for (i in 1:nrow(file_dir)) {
    # REMOVED: message(paste0("DEBUG trait_scan: Processing loop iteration i = ", i, " for dataset: ", selected_dataset))
    
    # REMOVED: chr_num <- NA
    # REMOVED: original_fst_path <- NA
    
    # REMOVED: tryCatch block start
    chr_num <- file_dir$ID_code[i]
    original_fst_path <- file_dir$File_path[i] # Path from file_index.csv
    # REMOVED: message(paste0("  Successfully assigned chr_num: ", chr_num, ", original_fst_path: ", original_fst_path))
    # REMOVED: tryCatch block end and error function
    # REMOVED: if (is.na(chr_num) || is.na(original_fst_path)) check

    file_specific_trait_type_from_file <- tolower(file_dir$trait_type[i]) 

    processed_trait_type <- file_specific_trait_type_from_file
    if (file_specific_trait_type_from_file == "clinical traits") {
      processed_trait_type <- "clinical"
    }

    corrected_fst_path <- original_fst_path

    # --- START PATH CORRECTION LOGIC (with initial message for original_fst_path) ---
    message(paste0("DEBUG trait_scan (Loop Iteration for Chr: ", chr_num, ", Dataset: ", selected_dataset, "):"))
    message(paste0("  Original FST path from file_index: ", original_fst_path))
    message(paste0("  Processed trait_type for path correction: ", processed_trait_type))

    if (!is.na(processed_trait_type) && nzchar(processed_trait_type)) {
        if (processed_trait_type == "clinical") {
            if (grepl("_with_symbols\\.fst$", original_fst_path)) {
                corrected_fst_path <- sub("_with_symbols\\.fst$", "_processed.fst", original_fst_path)
            } else if (!grepl("_processed\\.fst$", original_fst_path)) {
                corrected_fst_path <- paste0(tools::file_path_sans_ext(original_fst_path), "_processed.fst")
            }
        } else if (processed_trait_type %in% c("genes", "isoforms")) {
            if (grepl("_processed\\.fst$", original_fst_path)) {
                corrected_fst_path <- sub("_processed\\.fst$", "_with_symbols.fst", original_fst_path)
            } else if (!grepl("_with_symbols\\.fst$", original_fst_path)) {
                corrected_fst_path <- paste0(tools::file_path_sans_ext(original_fst_path), "_with_symbols.fst")
            }
        } else if (processed_trait_type == "liver_lipids") {
            if (grepl("_with_symbols\\.fst$", original_fst_path)) {
                corrected_fst_path <- sub("_with_symbols\\.fst$", "_processed.fst", original_fst_path)
            } else if (!grepl("_processed\\.fst$", original_fst_path)) {
                corrected_fst_path <- paste0(tools::file_path_sans_ext(original_fst_path), "_processed.fst")
            }
        } else { # ADDED ELSE to match working version's structure
            message(paste("  Unknown or unspecified trait_type ('", processed_trait_type, "') for file '", basename(original_fst_path), "'. Using original path."))
        }

        if (corrected_fst_path != original_fst_path) {
            message(paste0("  Path corrected based on trait_type. Original: ", basename(original_fst_path), ", Corrected: ", basename(corrected_fst_path)))
        } else {
            # message(paste0("  Path not changed by trait_type logic. Using: ", basename(corrected_fst_path))) # Message is slightly different now due to else above
        }
    } else {
        warning(paste("Missing or empty trait_type for file '", basename(original_fst_path), "'. Using original path. Please check file_index.csv."))
        # message(paste0("  Path not changed (missing trait_type). Using: ", basename(corrected_fst_path)))
    }
    # --- END PATH CORRECTION LOGIC ---

    fst_path <- corrected_fst_path 

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

    if (!file.exists(fst_path)) {
      warning("File check consistency error or corrected path invalid, skipping: ", fst_path, " (Original was: ", original_fst_path, ")")
      next
    }
    
    if (!stringr::str_detect(fst_path, "fst$")) {
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
    
    index_path_new <- sub("\\\\.fst$", "_rows.fst", fst_path)
    if (index_path_new == fst_path) { # Safety for paths not ending in .fst or if sub() fails
        index_path_new <- paste0(fst_path, "_rows.fst")
    }

    index_path_legacy <- sub("\\\\.fst$", "_row.fst", fst_path)
    if (index_path_legacy == fst_path) {
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
        tryCatch({
            row_index_path <- fst_rows(fst_path)
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
      trait_index[, Phenotype := tolower(Phenotype)]
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
        if ("Phenotype" %in% colnames(data)) {
            data <- data[tolower(Phenotype) == tolower(selected_trait)]
            message("Verified ", nrow(data), " rows for trait: ", selected_trait)
        }
        if (nrow(data) > 0) {
          message("Adding ", nrow(data), " rows from chromosome ", chr_num)
          all_data[[length(all_data) + 1]] <- data
        }
      } else {
        # More detailed logging when trait is not found in the current chromosome's index
        searched_trait_lower <- tolower(selected_trait)
        message(paste("DIAGNOSTIC: Trait '", selected_trait, "' (searched as '", searched_trait_lower, "') not found in chromosome ", chr_num, ". Index file: ", basename(row_index_path)))
        if (nrow(trait_index) > 0 && "Phenotype" %in% colnames(trait_index)) {
            # The Phenotype column in trait_index is already lowercased at this point by the line above.
            message(paste("DIAGNOSTIC: First 5 unique (lowercased) Phenotypes in this index are: ", paste(utils::head(unique(trait_index$Phenotype), 5), collapse=", ")))
        } else if (nrow(trait_index) == 0) {
            message(paste("DIAGNOSTIC: Index file", basename(row_index_path), "is empty."))
        } else { 
             message(paste("DIAGNOSTIC: Index file", basename(row_index_path), "is problematic (e.g. missing 'Phenotype' column). Columns: ", paste(colnames(trait_index), collapse=", ")))
        }
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
  if (exists("selected_dataset", inherits = FALSE)) {
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
