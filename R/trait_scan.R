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
    # Convert CSV to FST (or retrieve FST name)
    browser()
    out <- csv2fst(file_dir$File_path[1])
    row_name = out["row_name"]
    fst_name = out["fst_name"]

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
