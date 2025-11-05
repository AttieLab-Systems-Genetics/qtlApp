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
  numb_mice <- NA # Initialize to NA

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
      # Capture Numb_mice from the first valid trait_data
      if (is.na(numb_mice) && "Numb_mice" %in% colnames(trait_data)) {
        numb_mice <- trait_data$Numb_mice[1]
      }
    }
  }

  if (length(all_data) == 0) {
    stop("Trait '", selected_trait, "' not found in any chromosome for dataset: ", selected_dataset)
  }

  combined_data <- data.table::rbindlist(all_data, fill = TRUE)
  # Deduplicate potential overlaps from multiple slices by marker/chr/position
  dedup_keys <- intersect(c("marker", "chr", "position"), names(combined_data))
  if (length(dedup_keys) >= 1) {
    data.table::setkeyv(combined_data, dedup_keys)
    combined_data <- unique(combined_data)
  }
  message("Combined data: ", nrow(combined_data), " rows for trait: ", selected_trait)

  # Prepare the result list
  result <- list(
    scan_data = combined_data,
    numb_mice = numb_mice
  )

  # Cache the result
  if (!is.null(cache_env)) {
    cache_env[[cache_key]] <- result
  }

  return(result)
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

  # If the corrected path does not exist, add robust fallbacks for splice junctions
  if (!file.exists(corrected_path)) {
    tt <- tolower(processed_trait_type)
    if (grepl("splice|junction", tt)) {
      # Try toggling between 'splice_juncs' and 'splice_junctions'
      alt1 <- sub("splice_juncs", "splice_junctions", corrected_path, ignore.case = TRUE)
      if (!identical(alt1, corrected_path) && file.exists(alt1)) {
        return(alt1)
      }
      alt2 <- sub("splice_junctions", "splice_juncs", corrected_path, ignore.case = TRUE)
      if (!identical(alt2, corrected_path) && file.exists(alt2)) {
        return(alt2)
      }

      # Directory scan fallback: search for a file matching chromosome token and expected keywords
      dirpath <- dirname(corrected_path)
      # Extract chromosome token from original path (e.g., 'chromosome1', 'chromosomeX')
      chr_token_match <- regexpr("chromosome[0-9XYxy]+", basename(original_path), perl = TRUE)
      chr_token <- if (chr_token_match[1] > 0) substr(basename(original_path), chr_token_match[1], chr_token_match[1] + attr(chr_token_match, "match.length") - 1) else NULL
      if (!is.null(dirpath) && nzchar(dirpath) && !is.null(chr_token) && nzchar(chr_token) && dir.exists(dirpath)) {
        # Build a permissive pattern: chromosomeN .* splice .* junc .* diet .* interactive .* (fst|csv)
        # Prefer fst
        patt_fst <- paste0("^", chr_token, ".*splice.*junc.*diet.*interactive.*\\.fst$")
        candidates <- tryCatch(list.files(dirpath, pattern = patt_fst, ignore.case = TRUE, full.names = TRUE), error = function(e) character(0))
        # Exclude index files ("_rows.fst" or "_row.fst"); we want the data file
        if (length(candidates) > 0) {
          candidates <- candidates[!grepl("(_rows|_row)\\.fst$", basename(candidates), ignore.case = TRUE)]
        }
        if (length(candidates) == 0) {
          patt_csv <- paste0("^", chr_token, ".*splice.*junc.*diet.*interactive.*\\.csv$")
          candidates <- tryCatch(list.files(dirpath, pattern = patt_csv, ignore.case = TRUE, full.names = TRUE), error = function(e) character(0))
        }
        if (length(candidates) > 0) {
          # Prefer files that contain 'data_with_symbols' or 'data_processed'
          pref <- candidates[grepl("data_with_symbols\\.fst$", candidates, ignore.case = TRUE)]
          if (length(pref) == 0) pref <- candidates[grepl("data_processed\\.fst$", candidates, ignore.case = TRUE)]
          chosen <- if (length(pref) > 0) pref[[1]] else candidates[[1]]
          return(chosen)
        }
      }
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
    tryCatch(
      {
        row_index_path <- fst_rows(fst_path)
        if (file.exists(row_index_path)) {
          message("Generated index file: ", basename(row_index_path))
          return(row_index_path)
        }
      },
      error = function(e) {
        warning("Error creating index file for ", basename(fst_path), ": ", e$message)
      }
    )
  }

  return(NULL)
}

# Helper function to process trait data from a file
process_trait_from_file <- function(fst_path, row_index_path, selected_trait, chr_num) {
  tryCatch(
    {
      # Read the row index to find the trait
      trait_index <- fst::read_fst(row_index_path, as.data.table = TRUE)
      trait_index[, Phenotype := tolower(trimws(as.character(Phenotype)))]
      sel_trait <- tolower(trimws(as.character(selected_trait)))

      # First: exact lower-case match
      trait_rows <- trait_index[Phenotype == sel_trait, ]

      # Fallback: normalized match removing non-alphanumeric characters
      if (nrow(trait_rows) == 0) {
        sel_norm <- gsub("[^a-z0-9]+", "", sel_trait)
        trait_index[, phen_norm := gsub("[^a-z0-9]+", "", Phenotype)]
        trait_rows <- trait_index[phen_norm == sel_norm, ]
        if (nrow(trait_rows) == 0) {
          # Last resort: substring search on normalized keys
          trait_rows <- trait_index[grepl(sel_norm, phen_norm, fixed = TRUE), ]
          if (nrow(trait_rows) == 0) {
            # Debug: show a few available keys to help diagnose mismatches
            sample_keys <- paste(utils::head(unique(trait_index$Phenotype), 5), collapse = "; ")
            message(sprintf(
              "process_trait_from_file: No Phenotype match for '%s' (norm='%s') in %s chr %s. Sample keys: %s",
              sel_trait, sel_norm, basename(fst_path), as.character(chr_num), sample_keys
            ))
            return(NULL)
          } else {
            message(sprintf(
              "process_trait_from_file: Using normalized substring match for '%s' (norm='%s') in %s chr %s",
              sel_trait, sel_norm, basename(fst_path), as.character(chr_num)
            ))
          }
        } else {
          message(sprintf(
            "process_trait_from_file: Using normalized exact match for '%s' (norm='%s') in %s chr %s",
            sel_trait, sel_norm, basename(fst_path), as.character(chr_num)
          ))
        }
      }

      # Handle both old (from/to) and new (.row_min/.row_max) column naming
      if ("from" %in% colnames(trait_rows) && "to" %in% colnames(trait_rows)) {
        from_row <- as.integer(trait_rows$from)
        to_row <- as.integer(trait_rows$to)
      } else if (".row_min" %in% colnames(trait_rows) && ".row_max" %in% colnames(trait_rows)) {
        from_row <- as.integer(trait_rows$.row_min)
        to_row <- as.integer(trait_rows$.row_max)
      } else {
        warning("Row index file has unexpected column names for chromosome ", chr_num)
        return(NULL)
      }

      # Defensive bounds: ensure vectors are same length and valid scalars per slice
      n_slices <- min(length(from_row), length(to_row))
      if (n_slices <= 0) {
        return(NULL)
      }
      from_row <- from_row[seq_len(n_slices)]
      to_row <- to_row[seq_len(n_slices)]

      # Read one or more ranges; rbind if multiple slices matched
      message("Found trait in chromosome ", chr_num, " at rows count=", n_slices)
      slice_list <- vector("list", n_slices)
      for (k in seq_len(n_slices)) {
        fr <- from_row[k]
        tr <- to_row[k]
        if (!is.finite(fr) || !is.finite(tr) || tr < fr) next
        slice_list[[k]] <- tryCatch(
          fst::read_fst(
            fst_path,
            from = fr,
            to = tr,
            as.data.table = TRUE
          ),
          error = function(e) {
            warning("Failed reading slice ", k, " for chr ", chr_num, ": ", e$message)
            NULL
          }
        )
      }
      slice_list <- Filter(Negate(is.null), slice_list)
      if (length(slice_list) == 0) {
        return(NULL)
      }
      data <- data.table::rbindlist(slice_list, fill = TRUE)

      # Ensure required columns are present
      data <- ensure_required_columns(data, fst_path)
      if (is.null(data)) {
        return(NULL)
      }

      # Filter by phenotype if column exists, but do not drop the slice if no match
      if ("Phenotype" %in% colnames(data)) {
        data[, Phenotype := tolower(trimws(as.character(Phenotype)))]
        filtered <- data[Phenotype == sel_trait]
        if (nrow(filtered) == 0) {
          # Try normalized equality inside slice
          data[, phen_norm := gsub("[^a-z0-9]+", "", Phenotype)]
          sel_norm <- gsub("[^a-z0-9]+", "", sel_trait)
          filtered <- data[phen_norm == sel_norm]
        }
        # If still zero after attempts, keep original 'data' (slice corresponds to target trait)
        if (nrow(filtered) > 0) {
          data <- filtered
        }
      }

      if (nrow(data) > 0) {
        message("Adding ", nrow(data), " rows from chromosome ", chr_num)
        return(data)
      }
    },
    error = function(e) {
      warning("Error processing chromosome ", chr_num, ": ", e$message)
    }
  )

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
