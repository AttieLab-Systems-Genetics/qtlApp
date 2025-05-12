# R/load_initial_data.R

#' Load Initial Application Data
#'
#' Reads static data files required for the application at startup.
#'
#' @return A list containing:
#'   - file_directory: Data frame from file_index.csv.
#'   - gene_symbols: Character vector of gene symbols.
#'   - chr_breaks: Data frame with chromosome breaks.
#'   - annotation_list: List of annotation data frames.
#'   - markers: Data frame/object with marker information.
#'
#' @importFrom data.table fread
#' @importFrom stringr str_replace
#' @export
load_initial_data <- function() {
  message("--- Loading Initial Data ---")

  # Read the file index containing paths to all data files
  file_index_path <- "/data/dev/miniViewer_3.0/file_index.csv"
  if (!file.exists(file_index_path)) {
    stop("CRITICAL ERROR: File index not found at ", file_index_path)
  }
  file_directory <- read.csv(file_index_path)
  message("Loaded file directory.")

  # Try to load gene symbols file
  gene_symbols_path <- "/data/dev/miniViewer_3.0/gene_symbols.csv"
  # Default gene symbols in case loading fails
  default_gene_symbols <- c("Gnai3", "Cdc45", "Slc4a1", "Abca12", "Nadk", "Tfpi", "Scnn1b", "Cdc20", "Gpr89", "Cdc73")

  if (file.exists(gene_symbols_path)) {
    message("Loading gene symbols from: ", gene_symbols_path)
    gene_symbols <- tryCatch({
      as.character(data.table::fread(gene_symbols_path)$gene_symbol)
    }, error = function(e) {
      warning("Error reading gene symbols file: ", e$message, ". Using default symbols.")
      default_gene_symbols
    })
  } else {
    warning("Gene symbols file not found at: ", gene_symbols_path, ". Using default symbols.")
    gene_symbols <- default_gene_symbols
  }
  gene_symbols <- sort(unique(gene_symbols))
  message("Loaded ", length(gene_symbols), " gene symbols.")

  # Load chromosome break points for mm11 genome
  chr_breaks_path <- "/data/dev/miniViewer_3.0/chromosomal_sep_mm11.csv"
  if (!file.exists(chr_breaks_path)) {
    stop("CRITICAL ERROR: Chromosome breaks file not found at ", chr_breaks_path)
  }
  chr_breaks <- read.csv(chr_breaks_path)
  message("Loaded chromosome breaks.")

  # Load gene annotations
  annotation_list_path <- "/data/dev/miniViewer_3.0/annotation_list.rds"
  if (!file.exists(annotation_list_path)) {
    stop("CRITICAL ERROR: Annotation list file not found at ", annotation_list_path)
  }
  annotation_list <- readRDS(annotation_list_path)
  message("Loaded annotation list.")

  # Load marker information
  markers_path <- "/data/dev/miniViewer_3.0/CHTC_dietDO_markers_RDSgrcm39.rds"
   if (!file.exists(markers_path)) {
    stop("CRITICAL ERROR: Markers file not found at ", markers_path)
  }
  markers <- readRDS(markers_path)
  message("Loaded markers.")

  # Create a group identifier (needs to be done after loading file_directory)
  file_directory$group <- paste0(file_directory$diet, " ", file_directory$trait_compartment, " ",
                             file_directory$trait_type,
                             ifelse(file_directory$sexes == "Both", "", paste0(" (", file_directory$sexes, ")")),
                             ", ", file_directory$scan_type,
                             ifelse(file_directory$scan_type == "interactive",
                                    paste0(" (", file_directory$covars_interactive, ")"),
                                    ""))
  message("Created dataset group identifiers.")

  message("--- Initial Data Loading Complete ---")

  return(list(
    file_directory = file_directory,
    gene_symbols = gene_symbols,
    chr_breaks = chr_breaks,
    annotation_list = annotation_list,
    markers = markers
  ))
} 