# Helper functions for working with FST files, including indexing

# Function to convert CSV files to FST format for better performance
# Note: This function is less relevant if data is already pre-converted.
csv2fst <- function(csv_path, chunk_size = 50000) {
  # Check if the input file is a CSV
  if (!stringr::str_detect(csv_path, "\\.csv$")) {
    warning("Input path does not end with .csv: ", csv_path)
    # Check if it's already an FST file
    if (stringr::str_detect(csv_path, "\\.fst$")) {
        return(csv_path) # It's already FST, return path
    } else {
        stop("Input file is neither .csv nor .fst: ", csv_path)
    }
  }
  
  # Create FST filename by replacing .csv with .fst
  fst_path <- stringr::str_replace(csv_path, "\\.csv$", ".fst")
  
  # Only create FST file if it doesn't already exist
  if (!file.exists(fst_path)) {
      message("Converting CSV to FST (this might take a while): ", basename(csv_path))
      
      if (!requireNamespace("data.table", quietly = TRUE)) {
          stop("Package 'data.table' needed for csv2fst function. Please install it.", call. = FALSE)
      }
      if (!requireNamespace("fst", quietly = TRUE)) {
          stop("Package 'fst' needed for csv2fst function. Please install it.", call. = FALSE)
      }
      
      tryCatch({
        
          all_chunks <- list()
          skip <- 0
          
          
          repeat {
              message("Reading chunk starting at row ", skip + 1)
              chunk <- data.table::fread(csv_path, skip = skip, nrows = chunk_size, showProgress = FALSE) # Progress often noisy in logs
              if (nrow(chunk) == 0) break
              all_chunks[[length(all_chunks) + 1]] <- chunk
              skip <- skip + chunk_size
          }
          
          if (length(all_chunks) == 0) {
              stop("No data read from CSV file: ", csv_path)
          }
          
          # Combine all chunks
          full_data <- data.table::rbindlist(all_chunks, fill = TRUE)
          rm(all_chunks) # Free memory
          gc()
          
          
          if ("Phenotype" %in% colnames(full_data)) {
              message("Sorting data by Phenotype...")
              data.table::setorder(full_data, Phenotype)
          } else {
              warning("Column 'Phenotype' not found. Data will not be sorted. Row indexing might be incorrect.")
          }
          
          
          message("Writing FST file: ", basename(fst_path))
          fst::write_fst(full_data, path = fst_path, compress = 50)
          message("Successfully created FST file.")
          
      }, error = function(e) {
          stop("Error converting CSV to FST: ", e$message)
      })
      
  } else {
      message("FST file already exists: ", basename(fst_path))
  }
  
  return(fst_path)
}



fst_rows <- function(fst_path) {
  if (!stringr::str_detect(fst_path, "\\.fst$")) {
     stop("Input path is not an FST file: ", fst_path)
  }
  
  # Create row index filename
  row_path <- stringr::str_replace(fst_path, "\\.fst$", "_row.fst")
  
  # Only create row index if it doesn't exist
  if (!file.exists(row_path)) {
      message("Creating row index for: ", basename(fst_path))
      
      if (!requireNamespace("dplyr", quietly = TRUE)) {
          stop("Package 'dplyr' needed for fst_rows function. Please install it.", call. = FALSE)
      }
      if (!requireNamespace("fst", quietly = TRUE)) {
          stop("Package 'fst' needed for fst_rows function. Please install it.", call. = FALSE)
      }

      tryCatch({
          # Ensure Phenotype column exists before reading
          fst_metadata <- fst::metadata_fst(fst_path)
          if (!"Phenotype" %in% fst_metadata$columnNames) {
              stop("FST file must contain a 'Phenotype' column to create row index.")
          }
          
          # Read only the Phenotype column
          pheno_data <- fst::read_fst(fst_path, columns = "Phenotype")
          
          # Use dplyr for calculation (or data.table alternative)
          rows_data <- pheno_data %>%
              dplyr::mutate(rown = dplyr::row_number()) %>%
              dplyr::group_by(Phenotype) %>%
              dplyr::summarise(from = min(rown), to = max(rown), .groups = 'drop')
              
          # Save row indices
          fst::write_fst(rows_data, row_path)
          message("Successfully created row index: ", basename(row_path))
          
      }, error = function(e) {
          # Clean up potentially incomplete index file on error
          if (file.exists(row_path)) file.remove(row_path)
          stop("Error creating row index for ", basename(fst_path), ": ", e$message)
      })
  }
  return(row_path)
} 