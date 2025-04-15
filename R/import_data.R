#' Import data from files
#'
#' This function imports data from files specified in the `import.csv` file
#' located in the package's `data` directory. The function reads files with
#' extensions `.csv`, `.rds`, or `.xlsx` and returns a list of data frames.
#'
#' @importFrom readxl read_excel
#' @importFrom tools file_ext
#' @export
import_data <- function() {
  import <- read.csv(system.file("data/import.csv", package = "qtlApp"))
  out <- list()
  for(i in 1:nrow(import)) {
    ext <- tools::file_ext(import$filename[i])
    switch(ext,
      csv =   out[[import$object[i]]] <- read.csv(import$filename[i]),
      rds =   out[[import$object[i]]] <- readRDS(import$filename[i]),
      xlsx =  out[[import$object[i]]] <- readxl::read_excel(import$filename[i])
    )
  }
  out$file_directory$group <- paste0(out$file_directory$diet, " ", 
    out$file_directory$trait_compartment, " ",
    out$file_directory$trait_type, ", ", 
    out$file_directory$scan_type)
  return(out)
}
