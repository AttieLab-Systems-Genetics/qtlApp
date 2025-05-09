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
    
    # Determine if it's a clinical file early on
    is_clinical_file <- grepl("clinical_traits", basename(allele_path), ignore.case = TRUE)
    
    tryCatch({
        # Read allele effects CSV
        data <- fread(allele_path)
        
        # Remove which_mice column if it exists
        if ("which_mice" %in% colnames(data)) {
            data[, which_mice := NULL]
            message("Removed which_mice column from allele effects file: ", basename(allele_path))
        }
        
        # Identify the primary trait column
        trait_col <- intersect(c("trait", "Phenotype", "lodcolumn"), colnames(data))[1]
        
        if (!is.na(trait_col)) {
            if (is_clinical_file) {
                message("Clinical trait file. Preserving original trait identifiers in column: '", trait_col, "' for file: ", basename(allele_path))
                # For clinical traits, we use the trait_col as is. No gene symbol mapping.
                # Ensure the identified trait column is treated as character if it's going to hold diverse names
                if (!is.character(data[[trait_col]])) {
                    data[, (trait_col) := as.character(get(trait_col))]
                }
            } else {
                # Existing logic for non-clinical files (genes, isoforms)
                message("Non-clinical trait file (e.g., genes/isoforms). Attempting gene symbol mapping for column: '", trait_col, "' for file: ", basename(allele_path))
                
                # Store original data for potential revert and clean phenotype for merge
                original_data_for_revert <- copy(data) # Use copy for data.table
                data[, clean_trait_for_merge := gsub("^liver_", "", get(trait_col))]
                
                data_orig_rows <- nrow(data) # Store original row count before merge

                # Merge with gene symbols
                data <- merge(data, id_to_symbol, 
                              by.x = "clean_trait_for_merge", 
                              by.y = "gene_id", 
                              all.x = TRUE,
                              sort = FALSE) # Maintain original order as much as possible
                
                # Check for row changes due to merge
                if(nrow(data) != data_orig_rows) {
                    # Further check: did unique trait identifiers change count?
                    # If unique clean_trait_for_merge before merge != unique clean_trait_for_merge after merge (considering only matched ones)
                    # this indicates a many-to-many or one-to-many issue from the perspective of original data.
                    warning(paste("Merge with id_to_symbol resulted in changed row count for", basename(allele_path), 
                                  ". This could indicate complex mappings (e.g., one-to-many).",
                                  "Original rows:", data_orig_rows, "New rows:", nrow(data), 
                                  "Symbol replacement will proceed, but review mappings if issues arise."))
                    # Decide if symbol replacement should be skipped or data reverted.
                    # For now, let's allow it but warn heavily. If this causes downstream problems,
                    # one might choose to revert 'data' to 'original_data_for_revert' here.
                }

                # Update trait column with gene symbol where available
                data[!is.na(gene_symbol) & !is.na(get(trait_col)), (trait_col) := gene_symbol] # Ensure original col wasn't all NA
                message("Applied gene symbol mapping for ", basename(allele_path))
                
                # Clean up columns added for merging/mapping
                data[, clean_trait_for_merge := NULL]
                if ("gene_symbol" %in% colnames(data)) {
                    data[, gene_symbol := NULL] # Remove the merged gene_symbol column
                }
            }
        } else {
            warning("No standard trait column (trait, Phenotype, lodcolumn) found in: ", basename(allele_path))
        }
        
        # Determine output path based on input filename
        output_suffix <- if (is_clinical_file) "_processed.csv" else "_with_symbols.csv"
        message("For file '", basename(allele_path), "', using output suffix: '", output_suffix, "'")
        
        # Construct output path carefully to avoid issues if allele_path is in a subdirectory
        output_dir <- dirname(allele_path)
        output_filename <- sub("\\.csv$", output_suffix, basename(allele_path))
        output_path <- file.path(output_dir, output_filename)
        
        fwrite(data, output_path)
        
        message("Created ", basename(output_path), " in ", output_dir)
        return(output_path)
    }, error = function(e) {
        warning("Error processing allele effects file '", basename(allele_path), "': ", e$message)
        return(NULL)
    })
}

# Main execution
main <- function() {
    # Read gene annotations
    # Ensure the path is correct for your environment when running the script
    gene_anno_file <- "/mnt/rdrive/mkeller3/General/main_directory/files_for_cross_object/gene_annotations.csv"
    if (!file.exists(gene_anno_file)) {
        stop("Gene annotation file not found: ", gene_anno_file)
    }
    gene_data <- read_gene_annotations(gene_anno_file)
    
    # Define input directory
    # This should be the path on your HOST machine where the raw allele CSV files are.
    input_dir <- "/data/dev/miniViewer_3.0"
    if (!dir.exists(input_dir)) {
        stop("Input directory not found: ", input_dir)
    }
    
    # Example: Process a specific clinical allele effects file
    clinical_allele_file <- file.path(input_dir, "consolidate_allele_effects_clinical_traits_all_mice_additive.csv")
    if (file.exists(clinical_allele_file)) {
        process_allele_effects(clinical_allele_file, gene_data$id_to_symbol)
    } else {
        warning("Clinical allele file not found: ", clinical_allele_file)
    }

    # Example: Process a specific gene/isoform allele effects file (replace with actual filename if you have one)
    # gene_allele_file <- file.path(input_dir, "your_gene_or_isoform_allele_effects.csv")
    # if (file.exists(gene_allele_file)) {
    #     process_allele_effects(gene_allele_file, gene_data$id_to_symbol)
    # } else {
    #     warning("Gene/isoform allele file not found: ", gene_allele_file)
    # }
    
    message("\nAllele effects processing attempt complete.")
    message("Please check messages above for status of each file.")
}

# Run the script
main() 