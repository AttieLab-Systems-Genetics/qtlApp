#' Find the peaks
#'
#' @param file_dir Data frame with file directory information
#' @param selected_dataset Character string identifying the dataset group
#' @param selected_trait Optional character string for a specific trait (gene symbol) to subset. If NULL, return all peaks for the dataset.
#' @param trait_type Character string indicating the type of traits (e.g., "genes", "clinical")
#' @param cache_env Environment for caching results
#' @param use_cache Logical, whether to use caching (default TRUE)
#'
#' @importFrom data.table fread rbindlist setnames
#' @importFrom dplyr select rename filter distinct any_of all_of
#' @importFrom tools file_path_sans_ext
#' @export
peak_finder <- function(file_dir, selected_dataset, selected_trait = NULL, trait_type = NULL, cache_env = NULL, use_cache = TRUE) {
  
  # Input validation
  if (is.null(file_dir) || !is.data.frame(file_dir) || nrow(file_dir) == 0) {
    warning("peak_finder: file_dir is missing or empty.")
    return(data.frame())
  }
  if (is.null(selected_dataset) || selected_dataset == "") {
    warning("peak_finder: selected_dataset is missing.")
    return(data.frame())
  }
  
  file_dir_subset_for_path <- dplyr::filter(file_dir, 
                                     .data$group == selected_dataset & 
                                     .data$file_type == "peaks")
  
  peaks_file_path_for_warning <- "[Unknown CSV - file_dir_subset empty or path not found]"
  if (nrow(file_dir_subset_for_path) > 0) {
      temp_path <- file_dir_subset_for_path$File_path[1]
      if (file.exists(temp_path)) {
          peaks_file_path_for_warning <- temp_path
      } else {
          potential_paths_for_warning <- c(
              temp_path,
              paste0(tools::file_path_sans_ext(temp_path), ".csv"),
              paste0(tools::file_path_sans_ext(temp_path), ".fst")
          )
          found_for_warning <- FALSE
          for(p_warn in potential_paths_for_warning){
              if(file.exists(p_warn)){
                  peaks_file_path_for_warning <- p_warn
                  found_for_warning <- TRUE
                  break
              }
          }
          if(!found_for_warning){
              peaks_file_path_for_warning <- paste("[File not found for warning: ", temp_path, "]")
          }
      }
  } else {
      peaks_file_path_for_warning <- "[No 'peaks' entry in file_index for this dataset]"
  }

  cache_key <- selected_dataset
  
  if (use_cache && !is.null(cache_env) && !is.null(cache_env[[cache_key]])) {
    all_peaks_for_dataset <- cache_env[[cache_key]]
  } else {
    
    file_dir_subset <- file_dir_subset_for_path 
                                     
    if (nrow(file_dir_subset) == 0) {
      warning("peak_finder: No 'peaks' file found for dataset: ", selected_dataset)
      return(data.frame())
    }
    
    peaks_file_path_to_load <- file_dir_subset$File_path[1]
    actual_file_to_load <- NULL
    if (file.exists(peaks_file_path_to_load)) {
        actual_file_to_load <- peaks_file_path_to_load
    } else {
        potential_paths_to_load <- c(
            peaks_file_path_to_load,
            paste0(tools::file_path_sans_ext(peaks_file_path_to_load), ".csv"),
            paste0(tools::file_path_sans_ext(peaks_file_path_to_load), ".fst")
        )
        for(p_load in potential_paths_to_load){
            if(file.exists(p_load)){
                actual_file_to_load <- p_load
                break
            }
        }
    }
    if(is.null(actual_file_to_load)){
         warning("peak_finder: Could not find peaks file to load after checking common extensions for: ", peaks_file_path_to_load)
         return(data.frame())
    }
    
    all_peaks_for_dataset <- data.frame()

    tryCatch({
      all_peaks_for_dataset_raw <- data.table::fread(actual_file_to_load, stringsAsFactors = FALSE)
      
      current_colnames <- colnames(all_peaks_for_dataset_raw)
      
      column_map <- c(
        qtl_chr = "qtl_chr", qtl_pos = "qtl_pos", qtl_lod = "qtl_lod", marker = "marker", cis = "cis",
        qtl_ci_lo = "qtl_ci_lo", qtl_ci_hi = "qtl_ci_hi", A = "A", B = "B", C = "C", D = "D", 
        E = "E", F = "F", G = "G", H = "H"
      )
      essential_old_names <- c("qtl_chr", "qtl_pos", "qtl_lod", "marker")
      if (!is.null(trait_type) && trait_type %in% c("genes", "isoforms")) {
        column_map["gene_symbol"] <- "gene_symbol"
        column_map["gene_id"]     <- "gene_id"
        column_map["trait"]       <- "phenotype"
        column_map["gene_chr"]    <- "gene_chr"
        column_map["gene_start"]  <- "gene_start"
        essential_old_names <- c(essential_old_names, "phenotype", "gene_symbol", "gene_id", "gene_chr", "gene_start", "cis")
      } else { 
        column_map["trait"] <- "phenotype"
        essential_old_names <- c(essential_old_names, "phenotype", "cis")
      }
      rename_vec <- c()
      for (new_name_app in names(column_map)) {
        old_name_csv <- column_map[new_name_app]
        if (old_name_csv %in% current_colnames) {
          rename_vec[new_name_app] <- old_name_csv
        } else {
          if (old_name_csv %in% essential_old_names) { 
            warning("peak_finder: Essential original column '", old_name_csv, "' for mapping to '", new_name_app, "' not found in CSV.") 
          }
        }
      }
      successfully_mapped_essential_old_names <- intersect(essential_old_names, unname(rename_vec))
      missing_essential_mappings <- setdiff(essential_old_names, successfully_mapped_essential_old_names)
      if(length(missing_essential_mappings) > 0){
          warning("peak_finder: Missing or unmappable ESSENTIAL original columns in CSV for type '", trait_type, "': ", paste(missing_essential_mappings, collapse=", "), ". Cannot proceed.")
          stop("Essential columns missing, stopping processing in tryCatch.")
      }
      all_peaks_for_dataset <- dplyr::select(as.data.frame(all_peaks_for_dataset_raw), dplyr::all_of(unname(rename_vec))) %>%
        dplyr::rename(!!!rename_vec)

      if("qtl_chr" %in% colnames(all_peaks_for_dataset)){ all_peaks_for_dataset$qtl_chr <- as.character(all_peaks_for_dataset$qtl_chr) }
      if("qtl_pos" %in% colnames(all_peaks_for_dataset) && !is.numeric(all_peaks_for_dataset$qtl_pos)){ all_peaks_for_dataset$qtl_pos <- suppressWarnings(as.numeric(all_peaks_for_dataset$qtl_pos)) }
      if("qtl_lod" %in% colnames(all_peaks_for_dataset) && !is.numeric(all_peaks_for_dataset$qtl_lod)){ all_peaks_for_dataset$qtl_lod <- suppressWarnings(as.numeric(all_peaks_for_dataset$qtl_lod)) }
      allele_cols_to_check <- c("A", "B", "C", "D", "E", "F", "G", "H")
      for(ac in allele_cols_to_check){
          if(ac %in% colnames(all_peaks_for_dataset) && !is.numeric(all_peaks_for_dataset[[ac]])){
              all_peaks_for_dataset[[ac]] <- suppressWarnings(as.numeric(all_peaks_for_dataset[[ac]]))
          }
      }
      if (use_cache && !is.null(cache_env)) {
        cache_env[[cache_key]] <- all_peaks_for_dataset
      }
    }, error = function(e) {
      warning("peak_finder: Error reading or processing peaks file '", actual_file_to_load, "': ", e$message)
      all_peaks_for_dataset <<- data.frame()
    })
    if(!exists("all_peaks_for_dataset") || is.null(all_peaks_for_dataset) || nrow(all_peaks_for_dataset) == 0) { 
        return(data.frame())
    }
  } 

  if (is.null(all_peaks_for_dataset) || nrow(all_peaks_for_dataset) == 0) {
    return(data.frame())
  }

  if (!is.null(selected_trait) && selected_trait != "") {
    filter_col <- NULL
    if (!is.null(trait_type) && trait_type %in% c("genes", "isoforms")) {
        filter_col <- "gene_symbol"
    } else { 
        filter_col <- "trait"
    }
    if (!(filter_col %in% colnames(all_peaks_for_dataset))) {
      warning("peak_finder: Column '", filter_col, "' needed for filtering type '", trait_type, "' not found after renaming. Returning all peaks.")
      return(all_peaks_for_dataset)
    }
    filtered_peaks <- dplyr::filter(all_peaks_for_dataset, .data[[filter_col]] == selected_trait)
    
    if(nrow(filtered_peaks) == 0 && trait_type %in% c("genes", "isoforms")){
        warning("peak_finder: No peaks found for gene symbol '", selected_trait, "'. Check if this symbol exists in the 'gene_symbol' column of the CSV: ", basename(peaks_file_path_for_warning))
    }
    if(nrow(filtered_peaks) > 0 && "marker" %in% colnames(filtered_peaks) && "qtl_lod" %in% colnames(filtered_peaks)){
        filtered_peaks <- filtered_peaks %>%
           dplyr::group_by(marker) %>%
           dplyr::filter(qtl_lod == max(qtl_lod, na.rm = TRUE)) %>%
           dplyr::slice(1) %>% 
           dplyr::ungroup()
    }
    return(filtered_peaks)
  } else {
    return(all_peaks_for_dataset)
  }
}
