#' Convert CSV to FST File
#' 
#' @param csv_path path to CSV or FST file
#' 
#' @importFrom data.table fread
#' @importFrom dplyr arrange
#' @importFrom stringr str_detect str_replace
#' @importFrom fst write_fst
#' @export
csv2fst <- function(csv_path) {
  if(stringr::str_detect(csv_path, "csv$")) {
    fst_path <- stringr::str_replace(csv_path, "csv$", "fst")
    if(!file.exists(fst_path)) {
      # Write FST file.
      warning("Writing FST file: ", fst_path)
      fst::write_fst(
        data.table::fread(csv_path, drop = "Which_mice") |> dplyr::arrange(Phenotype),
        path = fst_path, compress = 100)
    }
  } else {
    if(!stringr::str_detect(csv_path, "fst$"))
      stop("No CSV or FST name provided: ", csv_path)
    fst_path <- csv_path
  }
  return(fst_path)
}
#' Convert CSV to FST File
#' 
#' @param fst_path path to FST file
#' 
#' @importFrom dplyr filter group_by mutate n row_number select slice
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_replace
#' @importFrom fst read_fst
#' @export
fst_rows <- function(fst_path) {
  row_path <- stringr::str_replace(fst_path, ".fst$", "_row.csv")
  if(!file.exists(row_path)) {
    # Database of first and last entries by phenotype
    rows <- fst::read_fst(fst_path) |>
      dplyr::select(Phenotype) |>
      dplyr::mutate(rown = dplyr::row_number()) |>
      dplyr::group_by(Phenotype) |>
      dplyr::slice(c(1, dplyr::n())) |>
      dplyr::mutate(set = c("from", "to")) |>
      tidyr::pivot_wider(names_from = "set", values_from = "rown")
    write.csv(rows, row_path, row.names = FALSE)
  }
  return(row_path)
}
