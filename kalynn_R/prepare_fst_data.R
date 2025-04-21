#!/usr/bin/env Rscript

# Load required packages
library(data.table)
library(dplyr)
library(stringr)
library(fst)

# Function to read and process gene annotations
read_gene_annotations <- function(anno_file) {
    gene_annos <- fread(anno_file)
    # Create mapping dataframes for both directions
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

# Function to process FST files and add gene symbols
process_fst_file <- function(fst_path, id_to_symbol) {
    message("Processing ", basename(fst_path))
    x
    # Check if file exists first
    if (!file.exists(fst_path)) {
        warning("File does not exist: ", fst_path)
        return(NULL)
    }
    
    # Check if this is already a processed file with _with_symbols suffix
    if (grepl("_with_symbols\\.fst$", fst_path)) {
        message("File appears to be already processed, skipping")
        return(fst_path)
    }
    
    tryCatch({
        # Read FST file
        data <- read_fst(fst_path, as.data.table = TRUE)
        
        # Remove which_mice column if it exists
        if ("which_mice" %in% colnames(data)) {
            data[, which_mice := NULL]
            message("Removed which_mice column")
        }
        
        # If Phenotype column exists, it likely contains ENSMUSG IDs
        if ("Phenotype" %in% colnames(data)) {
            # Remove 'liver_' prefix if it exists
            data[, clean_phenotype := gsub("^liver_", "", Phenotype)]
            
            # Merge with gene symbols
            data <- merge(data, id_to_symbol, 
                         by.x = "clean_phenotype", 
                         by.y = "gene_id", 
                         all.x = TRUE)
            
            # Create new phenotype column with gene symbol
            data[!is.na(gene_symbol), Phenotype := gene_symbol]
            
            # Clean up
            data[, clean_phenotype := NULL]
        }
        
        # Create new FST filename - maintain original pattern
        new_fst_path <- sub("\\.fst$", "_with_symbols.fst", fst_path)
        
        # Write new FST file
        write_fst(data, new_fst_path, compress = 50)
        
        message("Created ", basename(new_fst_path))
        return(new_fst_path)
    }, error = function(e) {
        warning("Error processing file ", basename(fst_path), ": ", e$message)
        return(NULL)
    })
}

# Main execution
main <- function() {
    # Read gene annotations
    gene_data <- read_gene_annotations("/mnt/rdrive/mkeller3/General/main_directory/files_for_cross_object_v5/gene_annotations.csv")
    
    # Define input directory
    input_dir <- "/data/dev/miniViewer_3.0"
    
    # Process chromosome FST files
    fst_pattern <- "chromosome[0-9XYM]+_genes_additive_female_mice_data\\.fst$"
    fst_files <- list.files(input_dir, pattern = fst_pattern, full.names = TRUE)
    
    message("Found ", length(fst_files), " chromosome FST files to process")
    new_fst_paths <- sapply(fst_files, process_fst_file, 
                           id_to_symbol = gene_data$id_to_symbol)
    
    # Remove NULL entries (failed files)
    new_fst_paths <- new_fst_paths[!sapply(new_fst_paths, is.null)]
    
    message("FST processing complete")
}

# Run the script
main() 