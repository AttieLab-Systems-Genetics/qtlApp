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
library(data.table)
library(dplyr)
library(stringr)

options(datatable.tmpdir = my_temp_dir)
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

# Function to process allele effects file
process_allele_effects <- function(allele_path, id_to_symbol) {
    message("Processing allele effects file: ", basename(allele_path))
    
    if (!file.exists(allele_path)) {
        warning("Allele effects file does not exist: ", allele_path)
        return(NULL)
    }
    
    tryCatch({
        # Read allele effects CSV
        data <- fread(allele_path)
        
        # Remove which_mice column if it exists
        if ("which_mice" %in% colnames(data)) {
            data[, which_mice := NULL]
            message("Removed which_mice column from allele effects file")
        }
        
        # If trait column exists and contains ENSMUSG IDs
        trait_col <- intersect(c("trait", "Phenotype", "lodcolumn"), colnames(data))[1]
        
        if (!is.null(trait_col)) {
            # Remove 'liver_' prefix if it exists
            data[, clean_trait := gsub("^liver_", "", get(trait_col))]
            
            # Merge with gene symbols
            data <- merge(data, id_to_symbol, 
                          by.x = "clean_trait", 
                          by.y = "gene_id", 
                          all.x = TRUE)
            
            # Update trait column with gene symbol where available
            data[!is.na(gene_symbol), (trait_col) := gene_symbol]
            
            # Clean up
            data[, clean_trait := NULL]
        }
        
        # Write processed file - maintain original pattern
        output_path <- sub("\\.csv$", "_with_symbols.csv", allele_path)
        fwrite(data, output_path)
        
        message("Created ", basename(output_path))
        return(output_path)
    }, error = function(e) {
        warning("Error processing allele effects file: ", e$message)
        return(NULL)
    })
}

# Main execution
main <- function() {
    # Read gene annotations
    gene_data <- read_gene_annotations("/mnt/rdrive/mkeller3/General/main_directory/files_for_cross_object_v5/gene_annotations.csv")
    
    # Define input directory
    input_dir <- "/data/dev/miniViewer_3.0"
    
    # Process allele effects file
    allele_file <- file.path(input_dir, "consolidate_allele_effects_v5_simple_scan_genes_additive_female_mice.csv")
    new_allele_path <- process_allele_effects(allele_file, gene_data$id_to_symbol)
    
    message("Allele effects processing complete")
}

# Run the script
main() 