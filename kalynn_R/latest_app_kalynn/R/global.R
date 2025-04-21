# Global settings and data loading for the Shiny application

# --- 1. Load required packages --- 
suppressPackageStartupMessages({
    library(shiny)
    library(dplyr)
    library(stringr)
    library(ggplot2)
    library(plotly)
    library(data.table)
    library(fst)
    library(bslib)
    library(shinyjs)
    library(shinycssloaders)
    library(DT)
    library(reshape2)
    library(ggiraph)
    library(ggrepel)
    library(ggpubr)
    library(debounce)
})

# --- 2. Global Options --- 
options(shiny.maxRequestSize = 20000*1024^2) # Increased upload size

# --- 3. Load Static Data --- 
message("Loading application data...")
tryCatch({
    # Use data.table::fread for potentially faster reading, especially if large
    file_directory_path <- "/data/dev/miniViewer_3.0/file_index.csv"
    if (!file.exists(file_directory_path)) stop("File index not found: ", file_directory_path)
    file_directory <- data.table::fread(file_directory_path)
    message("Loaded file directory.")

    chr_breaks_path <- "/data/dev/miniViewer_3.0/chromosomal_sep_mm11.csv"
    if (!file.exists(chr_breaks_path)) warning("Chromosome breaks file not found: ", chr_breaks_path)
    
    chr_breaks <- data.table::fread(chr_breaks_path)
    message("Loaded chromosome breaks.")

    annotation_list_path <- "/data/dev/miniViewer_3.0/annotation_list.rds"
    if (!file.exists(annotation_list_path)) warning("Annotation list file not found: ", annotation_list_path)
    annotation_list <- readRDS(annotation_list_path)
    message("Loaded gene annotations.")

    markers_path <- "/data/dev/miniViewer_3.0/CHTC_dietDO_markers_RDSgrcm39.rds"
    if (!file.exists(markers_path)) stop("Markers file not found: ", markers_path)
    # Assuming markers might be large, convert to data.table if not already
    markers <- as.data.table(readRDS(markers_path))
    message("Loaded markers.")

}, error = function(e) {
    stop("Fatal Error: Could not load required application data. Details: ", e$message)
})

# Pre-process file_directory using data.table syntax for efficiency
message("Processing file directory group identifiers...")
tryCatch({
    file_directory[, group := paste0(
        diet, " ", 
        trait_compartment, " ",
        trait_type, 
        fifelse(sexes == "Both", "", paste0(" (", sexes, ")")),
        ", ", 
        scan_type,
        fifelse(scan_type == "interactive",
               paste0(" (", covars_interactive, ")"),
               "")
    )]
    message("Group identifiers created.")
}, error = function(e) {
    warning("Could not create group identifiers in file directory: ", e$message)
})


# --- 4. Caching Environments --- 
message("Initializing cache environments...")
trait_cache <- new.env(parent = emptyenv())
peaks_cache <- new.env(parent = emptyenv())
message("Cache environments ready.")

# --- 5. Source Helper Functions --- 
message("Sourcing helper functions...")
helper_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
if (length(helper_files) == 0) {
    warning("No R files found in the R/ subdirectory.")
} else {
    for (f in helper_files) {
        tryCatch({
            # Source into the current environment (global environment for global.R)
            sys.source(f, envir = .GlobalEnv)
            message("Successfully sourced: ", basename(f))
        }, error = function(e) {
            warning("Error sourcing file ", basename(f), ": ", e$message)
        })
    }
}

message("Global setup complete.")

# --- End of global.R ---

