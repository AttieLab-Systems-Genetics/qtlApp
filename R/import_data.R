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
  base_path <- "/data/dev/miniViewer_3.0"
  message("Importing data from: ", base_path)
  
  file_index_path <- file.path(base_path, "file_index.csv")
  if (!file.exists(file_index_path)) {
    stop("Required file not found: ", file_index_path)
  }
  file_directory <- read.csv(file_index_path)
  message("Loaded file index: ", nrow(file_directory), " rows.")

  file_directory$group <- paste0(file_directory$diet, " ", file_directory$trait_compartment, " ",
                           file_directory$trait_type, 
                           ifelse(file_directory$sexes == "Both", "", paste0(" (", file_directory$sexes, ")")),
                           ", ", file_directory$scan_type,
                           ifelse(file_directory$scan_type == "interactive",
                                  paste0(" (", file_directory$covars_interactive, ")"),
                                  ""))
  message("Created group identifier in file directory.")

  # Load gene symbols
  gene_symbols_path <- file.path(base_path, "gene_symbols.csv")
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
  
  if (!("lipids" %in% names(annotation_list))) {
    warning("annotation_list$lipids not found in the loaded annotation_list.rds. Lipid traits will not be available unless this is generated.")
    annotation_list$lipids <- data.frame(data_name = character(0), stringsAsFactors = FALSE)
    message("Created empty placeholder for annotation_list$lipids as it was not found in annotation_list.rds.")
  } else {
    if (is.data.frame(annotation_list$lipids) && "data_name" %in% colnames(annotation_list$lipids)){
        message(paste("Found annotation_list$lipids with", nrow(annotation_list$lipids), "entries."))
    } else {
        warning("annotation_list$lipids exists but is not a data.frame with a 'data_name' column. Re-initializing as empty.")
        annotation_list$lipids <- data.frame(data_name = character(0), stringsAsFactors = FALSE)
    }
  }

  markers_path <- file.path(base_path, "CHTC_dietDO_markers_RDSgrcm39.rds")
   if (!file.exists(markers_path)) {
    stop("Required file not found: ", markers_path)
  }
  markers <- readRDS(markers_path)
  message("Loaded markers.")
  
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
