# Global configuration and setup for the Shiny application

# --- Load Libraries --- 
suppressPackageStartupMessages({
    library(shiny)
    library(shinyjs)
    library(ggplot2)
    library(plotly)
    library(DT)
    library(fst)
    library(data.table) # Use data.table for efficient data handling
    library(dplyr)      # Keep for potential compatibility or user preference
    library(bslib)      # For theming
    library(reshape2)   # For melt/dcast if needed (server uses data.table::melt)
    library(shinycssloaders) # For loading spinners
    library(debounce)      # For debouncing inputs
})

# --- Shiny Options --- 
# Set maximum file upload size, if needed (adjust as necessary)
options(shiny.maxRequestSize = 300*1024^2) 

# --- Configuration Paths --- 
# Define base directory and file paths (adjust as needed)
# It's often better to use relative paths or environment variables
base_dir <- "/projects/compsci/qtlApp/data"
file_directory_path <- file.path(base_dir, "file_directory.csv")
chr_breaks_path <- file.path(base_dir, "chr_breaks.csv")

# --- Load Configuration Data --- 
# Use tryCatch for robust loading
file_directory <- tryCatch({
    dt <- data.table::fread(file_directory_path, 
                           stringsAsFactors = FALSE, 
                           data.table = TRUE) 
    message("Successfully loaded file directory from: ", file_directory_path)
    dt
}, error = function(e) {
    warning("Error loading file directory: ", e$message, "\nPath: ", file_directory_path)
    # Return an empty data.table with expected columns on error
    data.table::data.table(pheno_data = character(), 
                         peak_data = character(), 
                         fst_file = character(), 
                         group = character()) 
})

chr_breaks <- tryCatch({
    dt <- data.table::fread(chr_breaks_path, 
                           stringsAsFactors = FALSE, 
                           data.table = TRUE) 
    # Ensure correct column types if needed
    if("start" %in% colnames(dt)) dt[, start := as.numeric(start)]
    if("end" %in% colnames(dt)) dt[, end := as.numeric(end)]
    if("length" %in% colnames(dt)) dt[, length := as.numeric(length)]
    if("chr" %in% colnames(dt)) dt[, chr := as.character(chr)]
    message("Successfully loaded chromosome breaks from: ", chr_breaks_path)
    dt
}, error = function(e) {
    warning("Error loading chromosome breaks: ", e$message, "\nPath: ", chr_breaks_path)
    # Return an empty data.table with expected columns on error
    data.table::data.table(chr = character(), 
                         start = numeric(), 
                         end = numeric(), 
                         length = numeric()) 
})

# --- Prepare Marker Information --- 
# Create a combined marker object (ensure chr_breaks loaded successfully)
markers <- if (nrow(chr_breaks) > 0) {
    data.table::setorder(chr_breaks, chr) # Ensure correct order
    # Assuming columns are named correctly ('chr', 'length') in chr_breaks
    m <- chr_breaks[, .(marker = paste0(chr, "_0"), chr, pos = 0, cM = 0)]
    m2 <- chr_breaks[, .(marker = paste0(chr, "_", length), chr, pos = length, cM = 0)]
    data.table::rbindlist(list(m, m2), use.names = TRUE)
} else {
    warning("chr_breaks data is empty, cannot create markers object.")
    # Return an empty data.table with expected columns
    data.table::data.table(marker = character(), 
                         chr = character(), 
                         pos = numeric(), 
                         cM = numeric()) 
}


# --- Cache Environments --- 
# Use environments for caching results from trait_scan and peak_finder
trait_cache <- new.env(parent = emptyenv())
peaks_cache <- new.env(parent = emptyenv())

# --- Load Utility Functions --- 
# Source R scripts containing helper functions
r_scripts_path <- "R"
scripts_to_source <- c("data_access.R", "fst_helpers.R", "plotting.R")

for (script in scripts_to_source) {
    script_path <- file.path(r_scripts_path, script)
    if (file.exists(script_path)) {
        tryCatch({
            source(script_path, local = TRUE) # Source into the global environment
            message("Successfully sourced: ", script_path)
        }, error = function(e) {
            warning("Error sourcing ", script_path, ": ", e$message)
        })
    } else {
        warning("Script not found: ", script_path)
    }
}

message("Global setup complete.")
# --- End of global.R --- 