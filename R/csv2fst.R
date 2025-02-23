#' Convert CSV to FST File
#' 
#' @param file_path path to CSV or FST file
#' 
#' @importFrom data.table fread
#' @importFrom dplyr filter group_by mutate n row_number select slice
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_detect str_replace
#' @importFrom fst read_fst write_fst
#' @export
csv2fst <- function(file_path) {
  csv_name <- file_path
  if(stringr::str_detect(csv_name, "csv$")) {
    fst_name <- stringr::str_replace(csv_name, "csv$", "fst")
    if(!file.exists(fst_name)) {
      # Write FST file.
      warning("Writing FST file: ", fst_name)
      fst::write_fst(
        data.table::fread(csv_name, drop = "Which_mice") |> dplyr::arrange(Phenotype),
        path = fst_name, compress = 100)
    }
  } else {
    if(!stringr::str_detect(csv_name, "fst$"))
      stop("No CSV or FST name provided: ", csv_name)
    fst_name <- csv_name
  }
  row_name <- stringr::str_replace(fst_name, ".fst$", "_row.csv")
  if(!file.exists(row_name)) {
    # Database of first and last entries by phenotype
    rows <- fst::read_fst(fst_name) |>
      dplyr::select(Phenotype) |>
      dplyr::mutate(rown = dplyr::row_number()) |>
      dplyr::group_by(Phenotype) |>
      dplyr::slice(c(1, dplyr::n())) |>
      dplyr::mutate(set = c("from", "to")) |>
      tidyr::pivot_wider(names_from = "set", values_from = "rown")
    write.csv(rows, row_name, row.names = FALSE)
  }
  return(c(row_name = row_name, fst_name = fst_name))
}
