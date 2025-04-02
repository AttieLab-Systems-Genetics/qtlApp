# Create a new environment for caching file paths to improve performance
file_path_cache <- new.env(parent = emptyenv())


trait_cache <- new.env(parent = emptyenv())


trait_scan <- function(file_dir, selected_dataset, selected_trait) {
   # Create a more unique cache key that includes the full dataset name
   cache_key <- paste(selected_dataset, tolower(selected_trait), sep = "_")
   
   if (!is.null(trait_cache[[cache_key]])) {
       message("Using cached data for trait: ", selected_trait, " in dataset: ", selected_dataset)
       return(trait_cache[[cache_key]])
   }
  
   message("Searching for trait: ", selected_trait, " in dataset: ", selected_dataset)
  
   # Filter the file directory for the selected dataset and scan type
   file_dir <- subset(file_dir, group == selected_dataset & file_type == "scans")
  
   # Check if we found any matching files
   if (nrow(file_dir) == 0) {
       stop("No matching files found for the selected dataset: ", selected_dataset)
   }
  
   # Initialize list to store data from each file
   all_data <- list()
  
   # Process each FST file (one per chromosome)
   for (i in 1:nrow(file_dir)) {
       # Extract chromosome number from ID_code
       chr_num <- file_dir$ID_code[i]
       fst_path <- file_dir$File_path[i]
      
       if (!file.exists(fst_path)) {
           warning("File does not exist: ", fst_path)
           next
       }
      
       # Ensure we are working with an FST file
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
                       setnames(data, possible_lod_cols[1], "LOD")
                   } else {
                       warning("LOD column not found in file: ", fst_path)
                       next
                   }
               }
              
               if (!"marker" %in% colnames(data)) {
                   possible_marker_cols <- grep("marker|id|snp", colnames(data), ignore.case = TRUE, value = TRUE)
                   if (length(possible_marker_cols) > 0) {
                       setnames(data, possible_marker_cols[1], "marker")
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
  
   # Check if we found any data
   if (length(all_data) == 0) {
       stop("Trait '", selected_trait, "' not found in any chromosome for dataset: ", selected_dataset)
   }
  
   # Combine all data
   combined_data <- data.table::rbindlist(all_data, fill = TRUE)
   message("Total rows in combined data: ", nrow(combined_data))
  
   # Cache the result
   trait_cache[[cache_key]] <- combined_data
   return(combined_data)
}

# Also add caching for peak_finder
peaks_cache <- new.env(parent = emptyenv())

peak_finder <- function(file_dir, selected_dataset, selected_trait = NULL) {
   # Create a unique cache key that includes the dataset
   cache_key <- if(is.null(selected_trait)) {
       selected_dataset  # Include full dataset name
   } else {
       paste(selected_dataset, tolower(selected_trait), sep = "_")  # Include full dataset name
   }
  
   # Check if we already have this data in cache
   if (is.null(peaks_cache[[cache_key]])) {
       message("Loading peaks data for dataset: ", selected_dataset)
       
       # Filter the file directory for the selected dataset and peaks files
       file_dir <- subset(file_dir, group == selected_dataset & file_type == "peaks")
      
       # Check if we found any matching files
       if (nrow(file_dir) == 0) {
           warning("No peaks files found for dataset: ", selected_dataset)
           # Return empty data frame with correct structure
           empty_peaks <- data.frame(
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
           peaks_cache[[cache_key]] <- empty_peaks
           return(empty_peaks)
       }
      
       # We now have a single consolidated peaks file
       message("Reading consolidated peaks file for dataset: ", selected_dataset)
       csv_path <- file_dir$File_path[1]
      
       if (!file.exists(csv_path)) {
           warning("Peaks file does not exist: ", csv_path)
           empty_peaks <- data.frame(
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
           if ("lodcolumn" %in% colnames(peaks_data)) {
               # We'll keep lodcolumn as is but also add a trait column for compatibility
               peaks_data$trait <- peaks_data$lodcolumn
               message("Added trait column based on lodcolumn")
           }
           if (any(grepl("phenotype", tolower(colnames(peaks_data))))) {
               phenotype_col <- grep("phenotype", tolower(colnames(peaks_data)), value = TRUE)[1]
               # We'll keep original column and add a trait column
               peaks_data$trait <- peaks_data[[phenotype_col]]
               message("Added trait column based on phenotype column")
           }
          
           # Filter for the specific trait if provided (case-insensitive)
           if (!is.null(selected_trait)) {
               message("Filtering for trait: ", selected_trait, " in dataset: ", selected_dataset)
               
               # Try various columns that might contain the trait
               rows_to_keep <- rep(FALSE, nrow(peaks_data))
               
               # Check trait column
               if ("trait" %in% colnames(peaks_data)) {
                   message("Checking trait column")
                   rows_to_keep <- rows_to_keep | (tolower(peaks_data$trait) == tolower(selected_trait))
               }
               
               # Check lodcolumn
               if ("lodcolumn" %in% colnames(peaks_data)) {
                   message("Checking lodcolumn column")
                   rows_to_keep <- rows_to_keep | (tolower(peaks_data$lodcolumn) == tolower(selected_trait))
               }
               
               # Check phenotype columns
               phenotype_cols <- grep("phenotype", tolower(colnames(peaks_data)), value = TRUE)
               for (col in phenotype_cols) {
                   message("Checking phenotype column: ", col)
                   rows_to_keep <- rows_to_keep | (tolower(peaks_data[[col]]) == tolower(selected_trait))
               }
               
               # Apply the filter
               peaks_data <- peaks_data[rows_to_keep]
               message("Found ", nrow(peaks_data), " peaks for trait: ", selected_trait, " in dataset: ", selected_dataset)
           }
          
           # Ensure we have all required columns
           required_cols <- c("marker", "chr", "pos", "lod")
          
           # Check and rename columns if needed
           col_mapping <- list(
               marker = c("marker", "markers", "id", "snp", "SNP"),
               chr = c("chr", "chrom", "chromosome"),
               pos = c("pos", "position", "bp", "location"),
               lod = c("lod", "LOD", "score", "pvalue")
           )
          
           # Try to map columns
           for (req_col in names(col_mapping)) {
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
              
               # Return empty data frame
               empty_peaks <- data.frame(
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
               peaks_cache[[cache_key]] <- empty_peaks
               return(empty_peaks)
           }
          
           # Print final column names for debug
           message("Final peaks data columns: ", paste(colnames(peaks_data), collapse=", "))
          
           # Convert to data.frame for compatibility
           peaks_data <- as.data.frame(peaks_data)
          
           # Cache the result
           peaks_cache[[cache_key]] <- peaks_data
          
       }, error = function(e) {
           warning("Error reading peaks file ", csv_path, ": ", e$message)
           # Return empty data frame
           empty_peaks <- data.frame(
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
           peaks_cache[[cache_key]] <- empty_peaks
       })
   } else {
       message("Using cached peaks data for ", if(is.null(selected_trait)) selected_dataset else paste(selected_trait, "in", selected_dataset))
   }
  
   return(peaks_cache[[cache_key]])
}


# set microfunctions==========================================================
QTL_plot_visualizer <- function(qtl.temp, phenotype, LOD_thr, mrkrs) {
  # Create a data.table for faster processing
  qtl.temp <- as.data.table(qtl.temp)
 
  # Select only the columns we need
  if ("marker" %in% colnames(qtl.temp) && "LOD" %in% colnames(qtl.temp)) {
      # If data already has the right columns, just select them
      qtl.temp <- qtl.temp[, .(markers = marker, LOD)]
  } else {
      # Otherwise, try to find the right columns
      qtl.temp <- qtl.temp[, .(markers = marker, LOD)]
  }
 
  # Convert markers to data.table for faster joins
  mrkrs2 <- as.data.table(mrkrs)[, .(markers = marker, chr, position = bp_grcm39/1e6)]
 
  # Use data.table join instead of merge (much faster)
  qtl.temp <- mrkrs2[qtl.temp, on = "markers", nomatch = NULL]
 
  # Convert chr to numeric efficiently
  qtl.temp[, chr := as.character(chr)]
  qtl.temp[chr == "X", chr := "20"]
  qtl.temp[chr == "Y", chr := "21"]
  qtl.temp[chr == "M", chr := "22"]
  qtl.temp[, chr := as.numeric(chr)]
 
  # Remove rows with NA in chr
  qtl.temp <- qtl.temp[!is.na(chr)]
 
  # Add order column (same as chr for numeric processing)
  qtl.temp[, order := chr]
 
  # Calculate chromosome lengths and cumulative positions
  chr_info <- qtl.temp[, .(chr_len = max(position)), by = chr]
  chr_info[, tot := cumsum(chr_len) - chr_len]
 
  # Join back to main data
  qtl.temp <- chr_info[qtl.temp, on = "chr", nomatch = NULL]
 
  # Calculate cumulative position
  qtl.temp[, BPcum := position + tot]
 
  # Sort by chromosome and position
  setorder(qtl.temp, chr, position)
 
  # Convert back to tibble for ggplot
  qtl_plot_obj <- as_tibble(qtl.temp)
 
  return(list(NULL, qtl_plot_obj))
}


# Function to convert CSV files to FST format for better performance
csv2fst <- function(csv_path, chunk_size = 50000) {
  # Check if the input file is a CSV
  if (stringr::str_detect(csv_path, "csv$")) {
      # Create FST filename by replacing .csv with .fst
      fst_path <- stringr::str_replace(csv_path, "csv$", "fst")
    
      # Only create FST file if it doesn't already exist
      if (!file.exists(fst_path)) {
          warning("Writing FST file in chunks: ", fst_path)
          # Initialize list to store chunks
          all_chunks <- list()
          skip <- 0
        
          # Read CSV in chunks to handle large files
          repeat {
              chunk <- data.table::fread(csv_path, skip = skip, nrows = chunk_size, showProgress = TRUE)
              if (nrow(chunk) == 0) break
              all_chunks[[length(all_chunks) + 1]] <- chunk
              skip <- skip + chunk_size
          }
        
          # Combine all chunks
          full_data <- data.table::rbindlist(all_chunks)
        
          # Sort by Phenotype if the column exists
          if ("Phenotype" %in% colnames(full_data)) {
              data.table::setorder(full_data, Phenotype)
          } else {
              warning("Column 'Phenotype' not found. Data will not be sorted.")
          }
        
          # Write to FST format with compression
          fst::write_fst(full_data, path = fst_path, compress = 50)
      }
  } else {
      # Check if the file is already an FST file
      if (!stringr::str_detect(csv_path, "fst$"))
          stop("No CSV or FST name provided: ", csv_path)
      fst_path <- csv_path
  }
  return(fst_path)
}

# Function to create row indices for FST files to speed up data access
fst_rows <- function(fst_path) {
  # Create row index filename
  row_path <- stringr::str_replace(fst_path, ".fst$", "_row.fst")
   # Only create row index if it doesn't exist
  if (!file.exists(row_path)) {
      # Read FST file and create row indices
      rows <- fst::read_fst(fst_path) |>
          # Select only Phenotype column
          dplyr::select(Phenotype) |>
          # Add row numbers
          dplyr::mutate(rown = dplyr::row_number()) |>
          # Group by Phenotype
          dplyr::group_by(Phenotype) |>
          # Get first and last row for each Phenotype
          dplyr::slice(c(1, dplyr::n())) |>
          # Create from/to columns
          dplyr::mutate(set = c("from", "to")) |>
          # Reshape to wide format
          tidyr::pivot_wider(names_from = "set", values_from = "rown")
      # Save row indices
      fst::write_fst(rows, row_path)
  }
  return(row_path)
}