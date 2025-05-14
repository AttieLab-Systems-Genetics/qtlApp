#!/usr/bin/env Rscript

# Load required packages
library(data.table)
library(dplyr)
library(stringr)
library(fst)

# Source the fst_rows.R script to make create_fst_rows_index function available
# Assumes the script is run with the working directory at the project root (e.g., ~/qtlApp)
# and R/fst_rows.R is at the project root's R directory.
source("R/fst_rows.R")

# Function to read and process gene annotations
read_gene_annotations <- function(anno_file) {
    if (!file.exists(anno_file)) {
        warning("Gene annotation file not found: ", anno_file)
        return(list(symbol_to_id = data.table(), id_to_symbol = data.table(), full_annos = data.table()))
    }
    gene_annos <- fread(anno_file)
    symbol_to_id <- gene_annos[, .(gene_symbol, gene_id)] %>%
        filter(!is.na(gene_symbol)) %>%
        distinct()
    
    id_to_symbol <- gene_annos[, .(gene_id, gene_symbol)] %>%
        filter(!is.na(gene_symbol)) %>%
        distinct()
    
    return(list(
        symbol_to_id = symbol_to_id,
        id_to_symbol = id_to_symbol,
        full_annos = gene_annos
    ))
}

# Function to process FST files
process_fst_file <- function(fst_path, id_to_symbol, file_type) {
    message("Processing ", basename(fst_path), " as type: ", file_type)
    
    if (!file.exists(fst_path)) {
        warning("File does not exist: ", fst_path)
        return(NULL)
    }
    
    # Check if this is already a processed file with _with_symbols or _processed suffix
    if (grepl("_with_symbols\\.fst$", fst_path) || grepl("_processed\\.fst$", fst_path)) {
        message("File appears to be already processed, skipping: ", basename(fst_path))
        # Still, ensure index exists for these already processed files if script is re-run
        create_fst_rows_index(fst_path)
        return(fst_path)
    }
    
    tryCatch({
        data <- read_fst(fst_path, as.data.table = TRUE)
        
        if ("which_mice" %in% colnames(data)) {
            data[, which_mice := NULL]
            message("Removed which_mice column from ", basename(fst_path))
        }
        
        new_fst_path <- NULL

        if (file_type == "genes") {
            if ("Phenotype" %in% colnames(data) && nrow(id_to_symbol) > 0) {
                data[, clean_phenotype := gsub("^liver_", "", Phenotype)]
                data_orig_rows <- nrow(data)
                data <- merge(data, id_to_symbol, 
                             by.x = "clean_phenotype", 
                             by.y = "gene_id", 
                             all.x = TRUE)
                # If merge changed row count due to many-to-many, stop and warn
                if(nrow(data) != data_orig_rows && !isTRUE(all.equal(unique(data$clean_phenotype), data$clean_phenotype))) {
                    warning(paste("Merge resulted in changed row count for", basename(fst_path), " Check for many-to-many mapping in annotations for IDs in this file. Skipping."))
                    return(NULL)
                }

                data[!is.na(gene_symbol), Phenotype := gene_symbol]
                data[, clean_phenotype := NULL]
                data[, gene_symbol := NULL] # Remove the explicitly merged gene_symbol column
                message("Applied gene symbol mapping for ", basename(fst_path))
            } else {
                message("Skipping gene symbol mapping for ", basename(fst_path), " (Phenotype column missing or no gene annotations).")
            }
            new_fst_path <- paste0(tools::file_path_sans_ext(fst_path), "_with_symbols.fst")
        } else if (file_type == "clinical") {
            # No specific processing for clinical traits' Phenotype column, it's used as is.
            message("Clinical trait file, Phenotype column (", paste(head(unique(data$Phenotype)), collapse=", "), ",...) will be used as is for ", basename(fst_path))
            message("DEBUG: Original fst_path for clinical: ", fst_path) # DEBUG
            new_fst_path <- paste0(tools::file_path_sans_ext(fst_path), "_processed.fst")
            message("DEBUG: new_fst_path for clinical after paste0: ", new_fst_path) # DEBUG
        } else if (file_type == "liver_lipids") {
            message("Liver lipid trait file, Phenotype column will be used as is for ", basename(fst_path))
            new_fst_path <- paste0(tools::file_path_sans_ext(fst_path), "_processed.fst")
        } else {
            warning("Unknown file type: ", file_type, " for file: ", basename(fst_path))
            return(NULL)
        }
        
        write_fst(data, new_fst_path, compress = 50)
        message("Created ", basename(new_fst_path))
        
        # Generate _rows.fst index file
        create_fst_rows_index(new_fst_path)
        
        return(new_fst_path)
    }, error = function(e) {
        warning("Error processing file ", basename(fst_path), ": ", e$message)
        return(NULL)
    })
}

# Main execution
main <- function() {
    # Path to gene annotations, ensure this is correct for your environment
    gene_anno_file <- "/mnt/rdrive/mkeller3/General/main_directory/files_for_cross_object/gene_annotations.csv"
    gene_data <- read_gene_annotations(gene_anno_file)
    
    # Define input directory, ensure this is correct for your environment
    # This should be the path on your HOST machine where the raw FST files are located.
    input_dir <- "/data/dev/miniViewer_3.0"
    if (!dir.exists(input_dir)) {
        stop("Input directory not found: ", input_dir)
    }

    file_processing_configs <- list(
        list(type = "clinical", pattern = "chromosome[X]+_clinical_traits_all_mice_diet_interactive_data\\.fst$"),
        list(type = "liver_lipids", pattern = "chromosome[0-9XYM]+_liver_lipids_all_mice_additive_data\\.fst$")
        # Add new types and patterns here, e.g., for isoforms:
        # list(type = "isoforms", pattern = "isoform_pattern\\.fst$")
    )
    
    all_processed_paths <- character(0) # Initialize empty character vector

    for (config in file_processing_configs) {
        message(paste("\\nProcessing type:", config$type, "with pattern:", config$pattern))
        
        # Find files, excluding ones that already have the processed/symbol suffix
        # to prevent reading them as input if script is run multiple times.
        # The check inside process_fst_file also serves as a safeguard.
        fst_files <- list.files(input_dir, pattern = config$pattern, full.names = TRUE)
        fst_files <- fst_files[!grepl("_with_symbols\\.fst$", fst_files) & !grepl("_processed.fst$", fst_files)]

        if (length(fst_files) == 0) {
            message("No raw files found for type: ", config$type)
            next
        }
        
        message("Found ", length(fst_files), " raw ", config$type, " FST files to process.")
        
        # Use lapply instead of sapply if you expect NULLs or complex return objects frequently
        # and want to keep them as list elements.
        # Here, we expect paths or NULLs, sapply is fine and will simplify to a vector.
        processed_paths_for_type <- sapply(fst_files, process_fst_file, 
                                           id_to_symbol = gene_data$id_to_symbol, 
                                           file_type = config$type)
        
        all_processed_paths <- c(all_processed_paths, processed_paths_for_type[!sapply(processed_paths_for_type, is.null)])
    }
    
    message("\\nFST processing complete. Successfully processed files:")
    if (length(all_processed_paths) > 0) {
      for(p in all_processed_paths) message(p)
    } else {
      message("No files were processed in this run.")
    }
}

# Run the script
main() 