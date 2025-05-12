#' Create row indices for FST files to speed up data access
#'
#' @param fst_path path to FST file
#'
#' @importFrom dplyr filter group_by mutate n row_number select slice
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_replace
#' @importFrom fst read_fst write_fst
#' @export
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

create_fst_rows_index <- function(fst_file_path) {
  # Create index for fst files for faster reading by trait
  # Stores min and max row for each trait in Phenotype column
  # Assumes fst file has columns marker, Phenotype, lod
  
  # These libraries are essential for the function's operation.
  # Consider loading them once globally if prepare_fst_data.R also uses them extensively,
  # or keep them here for encapsulation if this function is sourced by multiple scripts.
  library(fst)
  library(data.table)
  
  # Check if file exists
  if (!file.exists(fst_file_path)) {
    warning(paste("File not found:", fst_file_path))
    return(invisible(NULL)) # Return NULL invisibly to not print to console
  }
  
  message(paste("Generating row index for:", basename(fst_file_path)))
  
  tryCatch({
    # Read only necessary columns for efficiency
    DT <- fst::read_fst(fst_file_path, columns = c("marker", "Phenotype"), as.data.table = TRUE)
    
    # Ensure Phenotype column exists
    if (!"Phenotype" %in% colnames(DT)) {
      warning(paste("Phenotype column not found in:", basename(fst_file_path), "Skipping index generation."))
      return(invisible(NULL))
    }
    
    # Ensure there are rows to process
    if (nrow(DT) == 0) {
      warning(paste("No data rows found in:", basename(fst_file_path), "Skipping index generation."))
      return(invisible(NULL))
    }
    
    # Generate index: min and max row index (.I) for each unique Phenotype
    index_dt <- DT[, .(.row_min = min(.I), .row_max = max(.I)), by = Phenotype]
    
    # Define output path for the index file (e.g., path/to/file.fst -> path/to/file_rows.fst)
    out_path <- sub("\\.fst$", "_rows.fst", fst_file_path)
    
    # Safety check: if sub didn't change the path (e.g. no .fst extension), append _rows
    if (out_path == fst_file_path) { 
        out_path <- paste0(fst_file_path, "_rows.fst") 
    }

    fst::write_fst(index_dt, out_path)
    message(paste("Successfully created index file:", basename(out_path)))
    return(invisible(out_path)) # Return the path to the created index file invisibly
  }, error = function(e) {
    warning(paste("Error generating row index for", basename(fst_file_path), ":", e$message))
    return(invisible(NULL))
  })
}