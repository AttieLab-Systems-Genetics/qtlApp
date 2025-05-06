#' Import data from files
#'
#' This function imports data from files specified in the `import.csv` file
#' located in the package's `data` directory. The function reads files with
#' extensions `.csv`, `.rds`, or `.xlsx` and returns a list of data frames.
#'
#' @importFrom readxl read_excel
#' @importFrom tools file_ext
#' @importFrom data.table fread
#' @export
import_data <- function() {
  # Define base data path within the container
  base_path <- "/data/dev/miniViewer_3.0"
  
  message("Importing data from: ", base_path)
  
  # Load file index
  file_index_path <- file.path(base_path, "file_index.csv")
  if (!file.exists(file_index_path)) {
    stop("Required file not found: ", file_index_path)
  }
  file_directory <- read.csv(file_index_path)
  message("Loaded file index: ", nrow(file_directory), " rows.")

  # Load gene symbols
  gene_symbols_path <- file.path(base_path, "gene_symbols.csv")
  gene_symbols <- c("Gnai3", "Cdc45", "Slc4a1", "Abca12", "Nadk", "Tfpi", "Scnn1b", "Cdc20", "Gpr89", "Cdc73") # Default
  if (file.exists(gene_symbols_path)) {
    message("Loading gene symbols from: ", gene_symbols_path)
    tryCatch({
      gene_symbols <- as.character(data.table::fread(gene_symbols_path)$gene_symbol)
    }, error = function(e) {
      warning("Error reading gene symbols file: ", e$message, ". Using default symbols.")
    })
  } else {
    warning("Gene symbols file not found at: ", gene_symbols_path, ". Using default symbols.")
  }
  gene_symbols <- sort(gene_symbols)
  message("Loaded ", length(gene_symbols), " gene symbols")

  # Load chromosome breaks
  chr_breaks_path <- file.path(base_path, "chromosomal_sep_mm11.csv")
  if (!file.exists(chr_breaks_path)) {
    stop("Required file not found: ", chr_breaks_path)
  }
  chr_breaks <- read.csv(chr_breaks_path)
  message("Loaded chromosome breaks.")
  
  # Load gene annotations
  annotation_list_path <- file.path(base_path, "annotation_list.rds")
  if (!file.exists(annotation_list_path)) {
    stop("Required file not found: ", annotation_list_path)
  }
  annotation_list <- readRDS(annotation_list_path)
  message("Loaded annotation list.")
  
  # Load markers
  markers_path <- file.path(base_path, "CHTC_dietDO_markers_RDSgrcm39.rds")
   if (!file.exists(markers_path)) {
    stop("Required file not found: ", markers_path)
  }
  markers <- readRDS(markers_path)
  message("Loaded markers.")
  
  # Create group identifier in file_directory
  file_directory$group <- paste0(file_directory$diet, " ", file_directory$trait_compartment, " ",
                           file_directory$trait_type, 
                           ifelse(file_directory$sexes == "Both", "", paste0(" (", file_directory$sexes, ")")),
                           ", ", file_directory$scan_type,
                           ifelse(file_directory$scan_type == "interactive",
                                  paste0(" (", file_directory$covars_interactive, ")"),
                                  ""))
  message("Created group identifier in file directory.")

  # --- Removed old logic relying on import.csv --- 
  # import <- read.csv(system.file("data/import.csv", package = "qtlApp"))
  # out <- list()
  # for (i in 1:nrow(import)) {
  #   ext <- tools::file_ext(import$filename[i])
  #   switch(ext,
  #     csv =   out[[import$object[i]]] <- read.csv(import$filename[i]),
  #     rds =   out[[import$object[i]]] <- readRDS(import$filename[i]),
  #     xlsx =  out[[import$object[i]]] <- readxl::read_excel(import$filename[i])
  #   )
  # }
  # out$file_directory$group <- paste0(out$file_directory$diet, " ", 
  #   out$file_directory$trait_compartment, " ",
  #   out$file_directory$trait_type, ", ", 
  #   out$file_directory$scan_type)
  # -----------------------------------------------
  
  # Return a list containing all loaded data objects
  out <- list(
    file_directory = file_directory,
    annotation_list = annotation_list,
    gene_symbols = gene_symbols,
    markers = markers,
    chr_breaks = chr_breaks
  )
  
  message("Data import complete.")
  return(out)
}
