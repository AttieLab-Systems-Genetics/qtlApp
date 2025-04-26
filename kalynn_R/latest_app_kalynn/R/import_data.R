#' Import data from files
#'
#' Reads necessary data files for the QTL viewer application.
#'
#' @return A list containing loaded data objects: 
#'         `file_directory`, `gene_symbols`, `chr_breaks`, 
#'         `annotation_list`, `markers`.
#'
#' @importFrom data.table fread
#' @importFrom utils head
#'
#' @export
import_data <- function() {
  # Create caches needed by other functions if they don't exist
  if (!exists("trait_cache", envir = .GlobalEnv)) {
    assign("trait_cache", new.env(parent = emptyenv()), envir = .GlobalEnv)
    message("Created trait_cache environment.")
  }
  if (!exists("peaks_cache", envir = .GlobalEnv)) {
    assign("peaks_cache", new.env(parent = emptyenv()), envir = .GlobalEnv)
    message("Created peaks_cache environment.")
  }
  
  # Define paths (adjust as necessary for your environment)
  base_path <- "/data/dev/miniViewer_3.0/"
  file_index_path <- file.path(base_path, "file_index.csv")
  gene_symbols_path <- file.path(base_path, "gene_symbols.csv")
  chr_breaks_path <- file.path(base_path, "chromosomal_sep_mm11.csv")
  annotation_list_path <- file.path(base_path, "annotation_list.rds")
  markers_path <- file.path(base_path, "CHTC_dietDO_markers_RDSgrcm39.rds")

  # Load file directory
  if (!file.exists(file_index_path)) stop("File index not found: ", file_index_path)
  file_directory <- read.csv(file_index_path)

  # Load gene symbols
  gene_symbols <- c("Gnai3", "Cdc45", "Slc4a1", "Abca12", "Nadk", "Tfpi", "Scnn1b", "Cdc20", "Gpr89", "Cdc73") # Default
  if (file.exists(gene_symbols_path)) {
    message("Loading gene symbols from: ", gene_symbols_path)
    gene_symbols <- tryCatch({
      as.character(data.table::fread(gene_symbols_path)$gene_symbol)
    }, error = function(e) {
      warning("Error reading gene symbols file: ", e$message, ". Using default symbols.")
      gene_symbols # Return default on error
    })
  } else {
    warning("Gene symbols file not found at: ", gene_symbols_path, ". Using default symbols.")
  }
  gene_symbols <- sort(unique(gene_symbols))
  message("Loaded ", length(gene_symbols), " unique gene symbols")

  # Load chromosome breaks
  if (!file.exists(chr_breaks_path)) stop("Chromosome breaks file not found: ", chr_breaks_path)
  chr_breaks <- read.csv(chr_breaks_path)

  # Load gene annotations
  if (!file.exists(annotation_list_path)) stop("Annotation list file not found: ", annotation_list_path)
  annotation_list <- readRDS(annotation_list_path)

  # Load markers
  if (!file.exists(markers_path)) stop("Markers file not found: ", markers_path)
  markers <- readRDS(markers_path)

  # Create group identifier in file_directory
  file_directory$group <- paste0(file_directory$diet, " ", 
                               file_directory$trait_compartment, " ",
                               file_directory$trait_type, 
                               ifelse(file_directory$sexes == "Both", "", paste0(" (", file_directory$sexes, ")")),
                               ", ", file_directory$scan_type,
                               ifelse(file_directory$scan_type == "interactive",
                                      paste0(" (", file_directory$covars_interactive, ")"),
                                      ""))

  # Return all loaded objects in a list
  return(list(
    file_directory = file_directory,
    gene_symbols = gene_symbols,
    chr_breaks = chr_breaks,
    annotation_list = annotation_list,
    markers = markers
  ))
} 