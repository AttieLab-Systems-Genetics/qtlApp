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