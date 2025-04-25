#!/usr/bin/env Rscript

# Script to compile data from gz files into separate FST files for each chromosome




# Load required libraries
library(dplyr)
library(data.table)
library(fst)
library(stringr)
library(parallel)

# Define input and output directories
INPUT_DIR <- "/mnt/rdrive/mkeller3/General/main_directory/scans/LIVER_DO1200/v5_simple_scan_genes_diet_interactive_female_mice/output"
OUTPUT_DIR <- "/data/dev/miniViewer_3.0"

# Set data.table specific temp directory
#options(datatable.tmpdir = my_temp_dir)

# Extract the basename of the input directory for use in output filenames
# This will extract "isoform_additive" from the path
input_dirname <- basename(dirname(INPUT_DIR))
if (grepl("^v5_simple_scan_", input_dirname)) {
  # Remove the v5_simple_scan_ prefix if it exists
  input_dirname <- sub("^v5_simple_scan_", "", input_dirname)
}
message(paste("Using directory basename for output files:", input_dirname))

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Define all chromosomes
#CHROMOSOMES <- c(1:19, "X", "Y", "M")
CHROMOSOMES <- c(1:8)
# Load marker information
message("Loading marker information...")
markers_file <- "/data/dev/miniViewer_3.0/CHTC_dietDO_markers_RDSgrcm39.rds"
if (file.exists(markers_file)) {
  markers <- readRDS(markers_file)
  # Create a list of markers for each chromosome
  chr_markers <- lapply(CHROMOSOMES, function(chr) {
    chr_markers <- markers[markers$chr == as.character(chr), "marker"]
    message(paste("Found", length(chr_markers), "markers on chromosome", chr))
    return(chr_markers)
  })
  names(chr_markers) <- CHROMOSOMES
} else {
  message("Markers file not found. Will extract chromosome information from marker names.")
  chr_markers <- NULL
}

# Function to process a single gz file for all chromosomes
process_file <- function(file_path) {
  tryCatch({
    # Explicitly use the custom temp dir for reading (via options)
    data <- fread(file_path, header = TRUE)
    
    # Create a list to store data for each chromosome
    chr_data <- vector("list", length(CHROMOSOMES))
    # Ensure names are characters for consistent access
    names(chr_data) <- as.character(CHROMOSOMES)
    
    if (is.null(chr_markers)) {
      # Extract chromosome from marker names if no marker file
      data$chr <- str_extract(data$marker, "(?<=chr)[0-9XYM]+")
      # Split data by chromosome using character names
      for (chr_char in as.character(CHROMOSOMES)) {
        chr_data[[chr_char]] <- data[data$chr == chr_char]
      }
    } else {
      # Use marker lists to filter data for each chromosome
      # Loop through numeric values but access lists using character names
      for (chr in CHROMOSOMES) { 
        chr_char <- as.character(chr) # Use character version for indexing
        # Check if the name exists in chr_markers before accessing
        if (chr_char %in% names(chr_markers)) {
          chr_data[[chr_char]] <- data[data$marker %in% chr_markers[[chr_char]]]
        } else {
          # Handle case where chromosome name might be missing (shouldn't happen with current setup)
          chr_data[[chr_char]] <- data.table() # Assign empty table
        }
      }
    }
    
    return(chr_data)
  }, error = function(e) {
    message(paste("Error processing file:", file_path, "-", e$message))
    return(NULL)
  })
}

# List all gz files in the input directory
message("Listing gz files...")
gz_files <- list.files(INPUT_DIR, pattern = "\\.gz$", full.names = TRUE)
message(paste("Found", length(gz_files), "gz files"))

# Process files in parallel
message("Processing files...")
num_cores <- 8 # Try reducing to 4 cores
message(paste("Using", num_cores, "cores for parallel processing"))

# Initialize data tables for each chromosome
all_chr_data <- vector("list", length(CHROMOSOMES))
names(all_chr_data) <- CHROMOSOMES
for (chr in CHROMOSOMES) {
  all_chr_data[[chr]] <- data.table()
}

# Process in chunks to avoid memory issues
chunk_size <- 100
num_chunks <- ceiling(length(gz_files) / chunk_size)

for (i in 1:num_chunks) {
  message(paste("Processing chunk", i, "of", num_chunks))
  start_idx <- (i-1) * chunk_size + 1
  end_idx <- min(i * chunk_size, length(gz_files))
  chunk_files <- gz_files[start_idx:end_idx]
  
  # Process this chunk of files
  chunk_results <- mclapply(chunk_files, process_file, mc.cores = num_cores)
  
  # Combine results from this chunk for each chromosome
  for (chr in CHROMOSOMES) {
    chr_char <- as.character(chr) # Use character version for indexing
    # Extract data for this chromosome from all files in the chunk
    # Use chr_char to access the named list elements from chunk_results
    chr_chunk_data <- rbindlist(
      lapply(chunk_results[!sapply(chunk_results, is.null)], 
             function(x) x[[chr_char]]), # Access using character name
      fill = TRUE
    )
    
    # Append to the main data table for this chromosome
    if (nrow(chr_chunk_data) > 0) {
      # Access all_chr_data using character name
      if (chr_char %in% names(all_chr_data)) { 
        all_chr_data[[chr_char]] <- rbindlist(list(all_chr_data[[chr_char]], chr_chunk_data), fill = TRUE)
        message(paste("  Accumulated", nrow(all_chr_data[[chr_char]]), "rows for chromosome", chr_char))
      } else {
        message(paste("Warning: Chromosome", chr_char, "not found in all_chr_data list during accumulation."))
      }
    }
  }
  
  # Clean up to free memory
  rm(chunk_results, chr_chunk_data)
  gc()
}

# Write separate FST files for each chromosome with the input directory name included
for (chr in CHROMOSOMES) {
  chr_char <- as.character(chr) # Use character version for indexing
  # Check if the name exists and if data exists before writing
  if (chr_char %in% names(all_chr_data) && nrow(all_chr_data[[chr_char]]) > 0) { 
    output_file <- file.path(OUTPUT_DIR, paste0("chromosome", chr_char, "_", input_dirname, "_data.fst"))
    message(paste("Writing chromosome", chr_char, "data to FST file..."))
    # Explicitly use the custom temp directory for writing
    write_fst(all_chr_data[[chr_char]], output_file, compress = 50)
    message(paste("Successfully wrote", nrow(all_chr_data[[chr_char]]), "rows to", output_file))
  } else {
    message(paste("No data found for chromosome", chr_char))
  }
}

message("Done!")
