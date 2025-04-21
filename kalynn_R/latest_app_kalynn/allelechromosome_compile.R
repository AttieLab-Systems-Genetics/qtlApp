#!/usr/bin/env Rscript

my_temp_dir <- "/data/dev/tmp_KW"
dir.create(my_temp_dir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
message(paste("Using custom temporary directory for operations:", my_temp_dir))

# Ensure custom temp directory is cleaned up on exit (even if script is interrupted)
cleanup <- function() {
  if (dir.exists(my_temp_dir)) {
    message(paste("Cleaning up custom temporary directory:", my_temp_dir))
    unlink(my_temp_dir, recursive = TRUE, force = TRUE)
  }
}
on.exit(cleanup(), add = TRUE)

# Load required libraries
library(dplyr)
library(data.table)
library(stringr)
library(parallel)

# Define input and output directories
INPUT_DIR <- "/mnt/rdrive/mkeller3/General/main_directory/scans/LIVER_DO1200/v5_simple_scan_genes_additive_female_mice/output"
OUTPUT_DIR <- "/data/dev/miniViewer_3.0/"

options(datatable.tmpdir = my_temp_dir)

# Extract the directory name from the input path
dir_name <- basename(dirname(INPUT_DIR))
message("Processing data from directory: ", dir_name)

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Define all chromosomes
CHROMOSOMES <- c(1:19, "X", "Y", "M")

# Function to process a single CSV file for all chromosomes
process_file <- function(file_path) {
  tryCatch({
    # Read the CSV file
    data <- fread(file_path, header = TRUE)
    
    # Create a list to store data for each chromosome
    chr_data <- vector("list", length(CHROMOSOMES))
    names(chr_data) <- CHROMOSOMES
    
    # Split data by chromosome
    for (chr in CHROMOSOMES) {
      chr_data[[chr]] <- data[chr == as.character(chr)]
    }
    
    return(chr_data)
  }, error = function(e) {
    message(paste("Error processing file:", file_path, "-", e$message))
    return(NULL)
  })
}

# List all CSV files in the input directory
message("Listing CSV files...")
csv_files <- list.files(INPUT_DIR, pattern = "\\.csv$", full.names = TRUE)
message(paste("Found", length(csv_files), "CSV files"))

# Process files in parallel
message("Processing files...")
num_cores <- min(detectCores() - 1, 8)  # Use at most 8 cores
message(paste("Using", num_cores, "cores for parallel processing"))

# Initialize data tables for each chromosome
all_chr_data <- vector("list", length(CHROMOSOMES))
names(all_chr_data) <- CHROMOSOMES
for (chr in CHROMOSOMES) {
  all_chr_data[[chr]] <- data.table()
}

# Process in chunks to avoid memory issues
chunk_size <- 100
num_chunks <- ceiling(length(csv_files) / chunk_size)

for (i in 1:num_chunks) {
  message(paste("Processing chunk", i, "of", num_chunks))
  start_idx <- (i-1) * chunk_size + 1
  end_idx <- min(i * chunk_size, length(csv_files))
  chunk_files <- csv_files[start_idx:end_idx]
  
  # Process this chunk of files
  chunk_results <- mclapply(chunk_files, process_file, mc.cores = num_cores)
  
  # Combine results from this chunk for each chromosome
  for (chr in CHROMOSOMES) {
    # Extract data for this chromosome from all files in the chunk
    chr_chunk_data <- rbindlist(
      lapply(chunk_results[!sapply(chunk_results, is.null)], 
             function(x) x[[chr]]), 
      fill = TRUE
    )
    
    # Append to the main data table for this chromosome
    if (nrow(chr_chunk_data) > 0) {
      all_chr_data[[chr]] <- rbindlist(list(all_chr_data[[chr]], chr_chunk_data), fill = TRUE)
      message(paste("  Accumulated", nrow(all_chr_data[[chr]]), "rows for chromosome", chr))
    }
  }
  
  # Clean up to free memory
  rm(chunk_results, chr_chunk_data)
  gc()
}

# Write separate CSV files for each chromosome and prepare for consolidation
message("Writing individual chromosome files and preparing for consolidation...")
consolidated_data <- data.table()

for (chr in CHROMOSOMES) {
  if (nrow(all_chr_data[[chr]]) > 0) {
    # Sort by position before writing
    setorder(all_chr_data[[chr]], pos)
    
    # Write individual chromosome file
    output_file <- file.path(OUTPUT_DIR, paste0("chromosome", chr, "_allele_effects.csv"))
    message(paste("Writing chromosome", chr, "data to CSV file..."))
    fwrite(all_chr_data[[chr]], output_file)
    message(paste("Successfully wrote", nrow(all_chr_data[[chr]]), "rows to", output_file))
    
    # Add to consolidated data
    consolidated_data <- rbindlist(list(consolidated_data, all_chr_data[[chr]]), fill = TRUE)
    message(paste("Added chromosome", chr, "data to consolidated file"))
  } else {
    message(paste("No data found for chromosome", chr))
  }
}

# Write consolidated file
if (nrow(consolidated_data) > 0) {
  # Sort consolidated data by chromosome and position
  setorder(consolidated_data, chr, pos)
  
  # Write consolidated file with directory name
  consolidated_file <- file.path(OUTPUT_DIR, paste0("consolidate_allele_effects_", dir_name, ".csv"))
  message("Writing consolidated data to CSV file...")
  fwrite(consolidated_data, consolidated_file)
  message(paste("Successfully wrote", nrow(consolidated_data), "rows to consolidated file"))
  
  # Delete individual chromosome files after successful consolidation
  message("Cleaning up individual chromosome files...")
  for (chr in CHROMOSOMES) {
    chr_file <- file.path(OUTPUT_DIR, paste0("chromosome", chr, "_allele_effects.csv"))
    if (file.exists(chr_file)) {
      file.remove(chr_file)
      message(paste("Removed", chr_file))
    }
  }
  message("Cleanup complete")
} else {
  message("No data to consolidate")
}

message("Done!")
