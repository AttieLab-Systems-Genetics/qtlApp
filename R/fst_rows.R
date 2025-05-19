#' Create row indices for FST files to speed up data access
#'
#' @param fst_path path to FST file
#'
#' @importFrom dplyr filter group_by mutate n row_number select slice
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_replace
#' @importFrom fst read_fst write_fst
#' @export
#' 
library(fst)
library(data.table)

fst_rows <- function(fst_path) {
  row_path <- stringr::str_replace(fst_path, "\\.fst$", "_row.fst")
  if (!file.exists(row_path)) {
    # Read FST file and create row indices
    rows <- fst::read_fst(fst_path) |>
      # Select only Phenotype column
      dplyr::select(Phenotype) |>
      # Add row numbers
      dplyr::mutate(rown = dplyr::row_number()) |>
      dplyr::group_by(Phenotype) |>
      dplyr::slice(c(1, dplyr::n())) |>
      dplyr::mutate(set = c("from", "to")) |>
      tidyr::pivot_wider(names_from = "set", values_from = "rown")
    
    fst::write_fst(rows, row_path)
  }
  return(row_path)
}

create_fst_rows_index <- function(fst_file_path) {
  # Create index for fst files for faster reading by trait
  if (!file.exists(fst_file_path)) {
    warning(paste("File not found:", fst_file_path))
    return(invisible(NULL)) 
  }
  message(paste("Generating row index for:", basename(fst_file_path)))
  
  tryCatch({
    DT <- fst::read_fst(fst_file_path, columns = c("marker", "Phenotype"), as.data.table = TRUE)
    if (!"Phenotype" %in% colnames(DT)) {
      warning(paste("Phenotype column not found in:", basename(fst_file_path), "Skipping index generation."))
      return(invisible(NULL))
    }
    
    if (nrow(DT) == 0) {
      warning(paste("No data rows found in:", basename(fst_file_path), "Skipping index generation."))
      return(invisible(NULL))
    }
    
    index_dt <- DT[, .(.row_min = min(.I), .row_max = max(.I)), by = Phenotype]
    out_path <- sub("\\.fst$", "_rows.fst", fst_file_path)
    
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