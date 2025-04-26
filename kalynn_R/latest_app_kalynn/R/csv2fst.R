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

