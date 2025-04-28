#' Create row indices for FST files
#'
#' @param fst_path The path to the FST file.
#' @return The path to the created row index file, or NULL if creation fails.
#' @importFrom dplyr select mutate group_by slice row_number
#' @importFrom tidyr pivot_wider
#' @importFrom fst read_fst write_fst
#' @importFrom data.table as.data.table setnames
#' @export
fst_rows <- function(fst_path) {
  # Create row index filename
  row_path <- stringr::str_replace(fst_path, ".fst$", "_row.fst")
  
  # Only create row index if it doesn't exist
  if (!file.exists(row_path)) {
    message("Creating row index for: ", basename(fst_path))
    
    # Use tryCatch to handle potential errors during FST reading/processing
    result <- tryCatch({
      # Read FST file safely
      fst_data <- fst::read_fst(fst_path, as.data.table = TRUE)
      
      # Dynamically find the phenotype column (case-insensitive)
      phenotype_col <- grep("phenotype", colnames(fst_data), ignore.case = TRUE, value = TRUE)
      
      if (length(phenotype_col) == 0) {
        warning("No 'phenotype' column found in: ", basename(fst_path))
        return(NULL) # Return NULL if phenotype column is missing
      }
      
      # Use the first found phenotype column if multiple matches
      phenotype_col <- phenotype_col[1]
      
      # Rename the identified phenotype column to 'Phenotype' for consistency
      data.table::setnames(fst_data, old = phenotype_col, new = "Phenotype")
      
      # Create row indices
      rows <- fst_data |>
        dplyr::select(Phenotype) |> # Select the standardized 'Phenotype' column
        dplyr::mutate(rown = dplyr::row_number()) |>
        dplyr::group_by(Phenotype) |>
        dplyr::slice(c(1, dplyr::n())) |>
        dplyr::mutate(set = c("from", "to")) |>
        tidyr::pivot_wider(names_from = "set", values_from = "rown") |>
        dplyr::ungroup() # Ungroup after slicing
        
      # Save row indices
      fst::write_fst(rows, row_path)
      message("Successfully created row index: ", basename(row_path))
      return(row_path) # Return the path on success
      
    }, error = function(e) {
      warning("Error creating row index for ", basename(fst_path), ": ", e$message)
      return(NULL) # Return NULL on error
    })
    
    # If tryCatch returned NULL, return NULL from the function
    if (is.null(result)) {
      return(NULL)
    }
    
  } else {
      # If index already exists, just return the path
      return(row_path)
  }
}
