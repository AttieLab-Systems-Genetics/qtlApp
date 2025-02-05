# find the trait
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
    data <- data.table::fread(file_dir$File_path[1], drop = "Which_mice")
    
    # Filter for rows where Phenotype matches selected_trait
    data <- data[data$Phenotype == selected_trait, ]
    
    if (nrow(data) == 0) {
      stop("No data found for phenotype: ", selected_trait)
    }
    
    data
  }, error = function(e) {
    stop("Error reading file: ", e$message)
  })
  
  return(scan_data)
}
