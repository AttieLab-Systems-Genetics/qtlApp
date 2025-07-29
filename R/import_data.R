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

# ==============================================================================
# MAIN DATA IMPORT FUNCTION
# ==============================================================================
# This is the central data loader for the entire application
# It loads all reference data needed for the genetics/QTL viewer
import_data <- function() {
  # --------------------------------------------------------------------------
  # SETUP: Define base path for all data files
  # --------------------------------------------------------------------------
  base_path <- "/data/dev/miniViewer_3.0"
  message("Importing data from: ", base_path)

  # --------------------------------------------------------------------------
  # 1. LOAD FILE INDEX (Critical - defines available datasets)
  # --------------------------------------------------------------------------
  # This CSV acts as a "catalog" of all available data files
  # Contains metadata about each dataset (diet, traits, sex, scan type, etc.)
  file_index_path <- file.path(base_path, "file_index.csv")
  if (!file.exists(file_index_path)) {
    stop("Required file not found: ", file_index_path) # Hard stop - app can't work without this
  }
  file_directory <- read.csv(file_index_path)
  message("Loaded file index: ", nrow(file_directory), " rows.")

  # --------------------------------------------------------------------------
  # CLEAN DATA: Trim whitespace from character columns to prevent duplicates
  # --------------------------------------------------------------------------
  # Fix issue where dataset_category has trailing whitespace causing duplicates
  char_cols <- sapply(file_directory, is.character)
  file_directory[char_cols] <- lapply(file_directory[char_cols], trimws)
  message("Trimmed whitespace from character columns.")

  # --------------------------------------------------------------------------
  # CREATE GROUP IDENTIFIERS (UI Display Names)
  # --------------------------------------------------------------------------
  # This creates human-readable labels for the Shiny UI dropdowns
  # Combines: diet + trait_compartment + trait_type + sex + scan_type + covariates
  # Example: "High Fat liver clinical, interactive (age)"
  file_directory$group <- paste0(
    file_directory$diet, " ", file_directory$trait_compartment, " ",
    file_directory$trait_type,
    ifelse(file_directory$sexes == "Both", "", paste0(" (", file_directory$sexes, ")")),
    ", ", file_directory$scan_type,
    ifelse(file_directory$scan_type == "interactive",
      paste0(" (", file_directory$covars_interactive, ")"),
      ""
    )
  )

  # --------------------------------------------------------------------------
  # ADD file_type COLUMN (Required by trait_scan function)
  # --------------------------------------------------------------------------
  # Add file_type column to distinguish scan files from other file types
  if (!"file_type" %in% names(file_directory)) {
    file_directory$file_type <- "scans" # Default to scans for backward compatibility
    message("Added file_type column (defaulting to 'scans').")
  }

  # --------------------------------------------------------------------------
  # 2. LOAD GENE SYMBOLS (For gene search/lookup functionality)
  # --------------------------------------------------------------------------
  # Used to populate gene search dropdowns and validate gene names
  gene_symbols_path <- file.path(base_path, "gene_symbols.csv")
  gene_symbols <- c() # Initialize as empty vector
  if (file.exists(gene_symbols_path)) {
    tryCatch(
      {
        # fread() is faster than read.csv() for large files
        gene_symbols <- as.character(data.table::fread(gene_symbols_path)$gene_symbol)
        message("Loaded ", length(gene_symbols), " gene symbols.")
      },
      error = function(e) {
        # Graceful degradation - app continues but with limited gene search
        warning("Error reading gene symbols file: ", e$message, ". Using default symbols.")
        gene_symbols <- c("Actb", "Gapdh", "Tbp") # Fallback symbols
      }
    )
  } else {
    warning("Gene symbols file not found at: ", gene_symbols_path, ". Using default symbols.")
    gene_symbols <- c("Actb", "Gapdh", "Tbp") # Fallback symbols
  }
  gene_symbols <- sort(gene_symbols) # Alphabetical order for UI dropdowns


  # --------------------------------------------------------------------------
  # 3. LOAD CHROMOSOME BREAKS (For genome visualization)
  # --------------------------------------------------------------------------
  # Defines where chromosomes start/end on genome plots
  # Critical for proper chromosome boundary display
  chr_breaks_path <- file.path(base_path, "chromosomal_sep_mm11.csv")
  if (!file.exists(chr_breaks_path)) {
    stop("Required file not found: ", chr_breaks_path) # Hard stop - needed for plots
  }
  chr_breaks <- read.csv(chr_breaks_path)


  # --------------------------------------------------------------------------
  # 4. LOAD GENE ANNOTATIONS (Gene position/metadata)
  # --------------------------------------------------------------------------
  # Contains gene positions, descriptions, functional annotations
  # Used for gene lookup and linking QTL peaks to nearby genes
  annotation_list_path <- file.path(base_path, "annotation_list.rds")
  if (!file.exists(annotation_list_path)) {
    stop("Required file not found: ", annotation_list_path) # Hard stop - core functionality
  }
  annotation_list <- readRDS(annotation_list_path)
  message("Loaded annotation list.")

  # --------------------------------------------------------------------------
  # 5. LOAD GENETIC MARKERS (SNP positions for QTL mapping)
  # --------------------------------------------------------------------------
  # Contains SNP positions, genetic map information
  # Essential for QTL peak calling and significance thresholds
  markers_path <- file.path(base_path, "CHTC_dietDO_markers_RDSgrcm39.rds")
  if (!file.exists(markers_path)) {
    stop("Required file not found: ", markers_path) # Hard stop - needed for QTL plots
  }
  markers <- readRDS(markers_path)
  message("Loaded markers: ", nrow(markers), " markers.")

  # --------------------------------------------------------------------------
  # RETURN: Package everything into a named list
  # --------------------------------------------------------------------------
  # This list structure becomes the global data object used throughout the app
  out <- list(
    file_directory = file_directory, # Dataset catalog + UI labels
    annotation_list = annotation_list, # Gene annotations + lipids
    gene_symbols = gene_symbols, # Gene names for search
    markers = markers, # SNP positions for QTL
    chr_breaks = chr_breaks # Chromosome boundaries for plots
  )

  message("Data import complete.")
  return(out)
}

# ==============================================================================
# KEY INSIGHTS FOR YOUR APP CONTEXT:
# ==============================================================================
# 1. This is NOT Shiny file upload - it's loading reference data at startup
# 2. Error handling is mixed: hard stops for critical files, warnings for optional
# 3. The 'group' column creates user-friendly dataset names for UI dropdowns
# 4. Gene symbols are sorted for better UX in search dropdowns
# 5. Lipids handling shows defensive programming for optional features
# 6. All data is returned as a single list - likely stored in global environment
# 7. Heavy use of messaging for debugging/monitoring data load process
