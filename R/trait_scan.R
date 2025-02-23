#' Find trait scan
#' 
#' @param file_dir data frame with file directory information
#' @param selected_dataset character string
#' @param selected_trait character string
#' 
#' @importFrom fst read_fst
#' @importFrom dplyr filter
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
    fst_path <- csv2fst(file_dir$File_path[1])
    row_path <- fst_rows(fst_path)

    # Read FST data.
    # Filter for rows where Phenotype matches selected_trait
    rows <- dplyr::filter(fst::read_fst(row_path), Phenotype == selected_trait)
    data <- fst::read_fst(fst_path, from = rows$from, to = rows$to)
    
    if (nrow(data) == 0) {
      stop("No data found for phenotype: ", selected_trait)
    }
    
    data
  }, error = function(e) {
    stop("Error reading file: ", e$message)
  })
  
  return(scan_data)
}
