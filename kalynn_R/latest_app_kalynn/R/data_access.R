# Functions for accessing and processing scan and peak data

# Function to retrieve scan data for a specific trait from FST files
trait_scan <- function(file_dir, selected_dataset, selected_trait, trait_cache) {
  # Ensure required packages are loaded (can be moved to global if preferred)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' required.")
  if (!requireNamespace("fst", quietly = TRUE)) stop("Package 'fst' required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Package 'stringr' required.")
  
  # Create a more unique cache key that includes the full dataset name
  cache_key <- paste(selected_dataset, tolower(selected_trait), sep = "_")
  
  # Check cache first
  if (!is.null(trait_cache[[cache_key]])) {
      message("Using cached data for trait: ", selected_trait, " in dataset: ", selected_dataset)
      return(trait_cache[[cache_key]])
  }
  
  message("Searching for trait: ", selected_trait, " in dataset: ", selected_dataset)
  
  # Filter the file directory for the selected dataset and scan type
  # Assuming file_dir is a data frame/table passed from global.R
  dataset_files <- base::subset(file_dir, group == selected_dataset & file_type == "scans")
  
  # Check if we found any matching files
  if (nrow(dataset_files) == 0) {
      stop("No matching scan files found for the selected dataset: ", selected_dataset)
  }
  
  # Initialize list to store data from each relevant chromosome file
  all_data <- list()
  
  # Process each FST file (one per chromosome for the selected dataset)
  for (i in 1:nrow(dataset_files)) {
      chr_num <- dataset_files$ID_code[i]
      fst_path <- dataset_files$File_path[i]
      
      if (!file.exists(fst_path)) {
          warning("Scan file does not exist: ", fst_path)
          next
      }
      
      # Ensure it's an FST file (or try to find it if only CSV path given)
      if (!stringr::str_detect(fst_path, "\\.fst$")) {
          potential_fst_path <- stringr::str_replace(fst_path, "\\.csv$", ".fst")
          if (file.exists(potential_fst_path)) {
              fst_path <- potential_fst_path
          } else {
              warning("Neither CSV nor FST file found for: ", fst_path)
              next
          }
      }
      
      message("Checking chromosome ", chr_num, " (file: ", basename(fst_path), ") for trait: ", selected_trait)
      
      # Create/get row index path using the helper function
      tryCatch({
          row_index_path <- fst_rows(fst_path) # Assumes fst_rows is available (sourced in global.R)
          trait_index <- fst::read_fst(row_index_path, as.data.table = TRUE)
          
          # Convert Phenotype column to lowercase for case-insensitive matching
          trait_index[, Phenotype := tolower(Phenotype)]
          
          # Find the rows for the selected trait (case-insensitive)
          trait_rows <- trait_index[Phenotype == tolower(selected_trait), ]
          
          if (nrow(trait_rows) > 0) {
              message("Found trait in chromosome ", chr_num, " at rows ", trait_rows$from, "-", trait_rows$to)
              
              # Read only the rows for this trait
              data <- fst::read_fst(fst_path,
                                   from = trait_rows$from,
                                   to = trait_rows$to,
                                   as.data.table = TRUE)
              
              # --- Standardize required columns --- 
              # LOD Column (using data.table::setnames for efficiency)
              if (!"LOD" %in% colnames(data)) {
                  possible_lod_cols <- grep("^(lod|score)$", colnames(data), ignore.case = TRUE, value = TRUE)
                  if (length(possible_lod_cols) > 0) {
                      data.table::setnames(data, possible_lod_cols[1], "LOD")
                  } else {
                      warning("LOD/score column not found in file: ", fst_path)
                      next # Skip this file/chromosome if LOD is missing
                  }
              }
              
              # Marker Column
              if (!"marker" %in% colnames(data)) {
                  possible_marker_cols <- grep("^(marker|id|snp)$", colnames(data), ignore.case = TRUE, value = TRUE)
                  if (length(possible_marker_cols) > 0) {
                      data.table::setnames(data, possible_marker_cols[1], "marker")
                  } else {
                      warning("marker/id/snp column not found in file: ", fst_path)
                      next # Skip this file/chromosome if marker is missing
                  }
              }
              
              # Phenotype column (crucial for verification)
              pheno_col_name <- grep("^Phenotype$", colnames(data), ignore.case = TRUE, value = TRUE)
              if (length(pheno_col_name) == 0) {
                   warning("Phenotype column not found in file: ", fst_path, ". Cannot verify trait data.")
                   # Proceed with caution - maybe add a check later?
              } else {
                  # Ensure it's exactly "Phenotype" if found with different case
                  if(pheno_col_name != "Phenotype") data.table::setnames(data, pheno_col_name, "Phenotype")
                  
                  # Verify that we have the correct trait data (case-insensitive)
                  original_rows <- nrow(data)
                  # Use := for efficient subsetting in place (though creates a copy here)
                  data <- data[tolower(Phenotype) == tolower(selected_trait)]
                  message("Verified ", nrow(data), " out of ", original_rows, " rows for trait: ", selected_trait)
              }
              # ------
              
              if (nrow(data) > 0) {
                  message("Adding ", nrow(data), " rows from chromosome ", chr_num)
                  all_data[[length(all_data) + 1]] <- data
              }
          } else {
              message("Trait not found in chromosome ", chr_num)
          }
      }, error = function(e) {
          warning("Error processing scan data for chromosome ", chr_num, " (", basename(fst_path), "): ", e$message)
      })
  }
  
  # Check if any data was collected
  if (length(all_data) == 0) {
      stop("Trait ", sQuote(selected_trait), " not found in any accessible chromosome scan file for dataset: ", selected_dataset)
  }
  
  # Combine all collected data using data.table::rbindlist
  combined_data <- data.table::rbindlist(all_data, fill = TRUE)
  message("Total rows in combined scan data for trait ", sQuote(selected_trait), ": ", nrow(combined_data))
  rm(all_data) # Clean up list
  gc()         # Suggest garbage collection
  
  # Cache the result
  trait_cache[[cache_key]] <- combined_data
  return(combined_data)
}

# Function to retrieve peak data for a dataset (and optionally a specific trait)
peak_finder <- function(file_dir, selected_dataset, selected_trait = NULL, peaks_cache) {
  # Ensure required packages are loaded
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' required.")
  
  # Create a unique cache key that includes the dataset
  cache_key <- if(is.null(selected_trait)) {
      selected_dataset  # Key for all peaks in the dataset
  } else {
      paste(selected_dataset, tolower(selected_trait), sep = "_") # Key for specific trait peaks
  }
  
  # Check cache first
  if (!is.null(peaks_cache[[cache_key]])) {
      message("Using cached peaks data for ", if(is.null(selected_trait)) selected_dataset else paste(sQuote(selected_trait), "in", selected_dataset))
      return(peaks_cache[[cache_key]])
  }
   
  message("Loading peaks data for dataset: ", selected_dataset)
   
  # Filter the file directory for the selected dataset and peaks files
  peaks_files <- base::subset(file_dir, group == selected_dataset & file_type == "peaks")
  
  # Define the expected structure for an empty peaks data.table
  empty_peaks_dt <- data.table::data.table(
       marker = character(0),
       trait = character(0),
       chr = character(0),
       pos = numeric(0),
       lod = numeric(0),
       A = numeric(0), B = numeric(0), C = numeric(0), D = numeric(0),
       E = numeric(0), F = numeric(0), G = numeric(0), H = numeric(0)
       # Add other potential columns like ci_lo, ci_hi, cis if expected
   )
   
  # Check if exactly one peaks file was found
  if (nrow(peaks_files) != 1) {
      warning("Expected exactly one peaks file for dataset: ", selected_dataset, ", found ", nrow(peaks_files), ". Please check file_index.csv.")
      peaks_cache[[cache_key]] <- empty_peaks_dt
      return(empty_peaks_dt)
  }
  
  csv_path <- peaks_files$File_path[1]
  
  if (!file.exists(csv_path)) {
      warning("Peaks file does not exist: ", csv_path)
      peaks_cache[[cache_key]] <- empty_peaks_dt
      return(empty_peaks_dt)
  }
  
  # Read the consolidated CSV file using data.table::fread
  tryCatch({
      message("Reading consolidated peaks file: ", basename(csv_path), " for dataset: ", selected_dataset)
      peaks_data <- data.table::fread(csv_path)
       
      message("Original peaks file columns: ", paste(colnames(peaks_data), collapse=", "))
      
      # --- Column Standardization using data.table::setnames --- 
      # Define mappings from possible names to standard names
      col_mapping <- list(
          trait = c("trait", "lodcolumn", "phenotype"), # Order matters - prefers 'trait'
          marker = c("marker", "markers", "id", "snp", "SNP"),
          chr = c("chr", "chrom", "chromosome"),
          pos = c("pos", "position", "bp", "location"),
          lod = c("lod", "LOD", "score") 
      )
      
      current_names <- colnames(peaks_data)
      for (std_name in names(col_mapping)) {
          # Find the first matching alternative name present in the data
          found_name <- NULL
          for (alt_name in col_mapping[[std_name]]) {
              # Use case-insensitive matching to find the column
              match_idx <- match(tolower(alt_name), tolower(current_names))
              if (!is.na(match_idx)) {
                  found_name <- current_names[match_idx]
                  break
              }
          }
          
          # If found and not already the standard name, rename it
          if (!is.null(found_name)) {
              if (found_name != std_name) {
                   message("Renaming peaks column ", sQuote(found_name), " to ", sQuote(std_name))
                   data.table::setnames(peaks_data, found_name, std_name)
                   current_names <- colnames(peaks_data) # Update current names after rename
              }
          } else {
              # Warn if essential columns (marker, chr, pos, lod) are missing
              if (std_name %in% c("marker", "chr", "pos", "lod")){
                  warning("Required peaks column ", sQuote(std_name), " or a known alternative was not found.")
              } 
              # Note: Missing 'trait' is handled later during filtering
          }
      }
      
      # Check for allele columns (A-H) and add NAs if missing
      allele_cols <- LETTERS[1:8]
      missing_alleles <- allele_cols[!allele_cols %in% colnames(peaks_data)]
      if (length(missing_alleles) > 0) {
          warning("Allele effect columns missing from peaks file: ", paste(missing_alleles, collapse=", "))
          for(col in missing_alleles) peaks_data[, (col) := NA_real_]
      }
      # --- End Standardization ---
      
      # Filter for the specific trait if provided (case-insensitive)
      if (!is.null(selected_trait)) {
           if ("trait" %in% colnames(peaks_data)) { # Check if trait column exists after standardization
              message("Filtering peaks for trait: ", selected_trait)
              original_rows <- nrow(peaks_data)
              # Use data.table filtering syntax
              peaks_data <- peaks_data[tolower(trait) == tolower(selected_trait)]
              message("Found ", nrow(peaks_data), " out of ", original_rows, " peaks for trait: ", sQuote(selected_trait))
          } else {
              warning("Cannot filter by trait (", sQuote(selected_trait), ") because 'trait' column is missing.")
          }
      }
      
      # Check if essential columns exist after standardization
      required_cols_final <- c("marker", "chr", "pos", "lod")
      missing_final <- required_cols_final[!required_cols_final %in% colnames(peaks_data)]
      if (length(missing_final) > 0) {
          # Stop if essential columns are missing, as downstream code likely depends on them
          stop("Essential columns missing after processing peaks file: ", paste(missing_final, collapse=", "))
      }
      
      # Optional: Convert column types if needed (fread usually does well)
      # peaks_data[, pos := as.numeric(pos)]
      # peaks_data[, lod := as.numeric(lod)]
      
      # Cache the result (as data.table)
      peaks_cache[[cache_key]] <- peaks_data
      message("Successfully processed and cached peaks data.")
      
  }, error = function(e) {
      warning("Error reading or processing peaks file ", csv_path, ": ", e$message)
      peaks_cache[[cache_key]] <- empty_peaks_dt # Cache empty table on error
  })
  
  # Return the cached data (might be empty data.table if errors occurred)
  return(peaks_cache[[cache_key]])
} 