#' Find trait scan
#' 
#' @param file_dir data frame with file directory information
#' @param selected_dataset character string
#' @param selected_trait character string
#' 
#' @importFrom data.table fread
#' @importFrom dplyr filter group_by mutate n row_number select slice
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_detect str_replace
#' @importFrom fst read_fst write_fst
#' @export
trait_scan <- function(file_dir, selected_dataset, selected_trait) {
  # Subset the data
  file_dir <- subset(file_dir, group == selected_dataset)
  file_dir <- subset(file_dir, file_type == "scans")
  
  if (nrow(file_dir) == 0) {
    stop("No matching files found for the selected dataset")
  }
  
  if (!file.exists(file_dir$File_path[1])) {
    stop("File does not exist: ", file_dir$File_path[1])
  }
  
  scan_data <- tryCatch({
    # Read all data first
    csv_name <- file_dir$File_path[1]
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

    # Read FST data.
    # Filter for rows where Phenotype matches selected_trait
    rows <- dplyr::filter(read.csv(row_name), Phenotype == selected_trait)
    data <- fst::read_fst(fst_name, from = rows$from, to = rows$to)
    
    if (nrow(data) == 0) {
      stop("No data found for phenotype: ", selected_trait)
    }
    
    data
  }, error = function(e) {
    stop("Error reading file: ", e$message)
  })
  
  return(scan_data)
}
