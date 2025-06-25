#!/usr/bin/env Rscript

# Load required packages
library(data.table)
library(dplyr)
library(stringr)
library(fst)

# Source the fst_rows.R script to make create_fst_rows_index function available
# Try different paths to find fst_rows.R
if (file.exists("R/fst_rows.R")) {
    source("R/fst_rows.R")
} else if (file.exists("../../R/fst_rows.R")) {
    source("../../R/fst_rows.R")
} else {
    stop("Cannot find R/fst_rows.R. Please run this script from the project root.")
}

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

# Function to read and process transcript annotations
read_transcript_annotations <- function(anno_file) {
    if (!file.exists(anno_file)) {
        warning("Transcript annotation file not found: ", anno_file)
        return(list(id_to_symbol = data.table(), full_annos = data.table()))
    }
    transcript_annos <- fread(anno_file)
    # Ensure required columns exist
    if (!all(c("transcript_id", "transcript_symbol") %in% colnames(transcript_annos))) {
        warning("Transcript annotation file must contain 'transcript_id' and 'transcript_symbol' columns. Path: ", anno_file)
        return(list(id_to_symbol = data.table(), full_annos = data.table()))
    }

    id_to_symbol <- transcript_annos[, .(transcript_id, transcript_symbol)] %>%
        filter(!is.na(transcript_symbol) & !is.na(transcript_id)) %>%
        distinct(transcript_id, .keep_all = TRUE) # Ensure transcript_id is unique for mapping

    return(list(
        id_to_symbol = id_to_symbol, # This will be used for mapping
        full_annos = transcript_annos
    ))
}

# Function to process FST files
process_fst_file <- function(fst_path, gene_id_to_symbol_map, transcript_id_to_symbol_map, file_type) {
    message("Processing ", basename(fst_path), " as type: ", file_type)

    if (!file.exists(fst_path)) {
        warning("File does not exist: ", fst_path)
        return(NULL)
    }

    # Check if this is already a processed file with _with_symbols or _processed suffix
    if (grepl("_with_symbols\\.fst$", fst_path) || grepl("_processed\\.fst$", fst_path) || grepl("_with_transcript_symbols\\.fst$", fst_path)) {
        message("File appears to be already processed, skipping: ", basename(fst_path))
        # Still, ensure index exists for these already processed files if script is re-run
        create_fst_rows_index(fst_path) # Make sure this function is available
        return(fst_path)
    }

    tryCatch(
        {
            data <- read_fst(fst_path, as.data.table = TRUE)

            # Remove unwanted columns
            if ("which_mice" %in% colnames(data)) {
                data[, which_mice := NULL]
                message("Removed which_mice column from ", basename(fst_path))
            }

            # Remove chr_from_map column if it exists (leftover from chromosome_compile.R)
            if ("chr_from_map" %in% colnames(data)) {
                data[, chr_from_map := NULL]
                message("Removed chr_from_map column from ", basename(fst_path))
            }

            new_fst_path <- NULL

            if (file_type == "genes") {
                if ("Phenotype" %in% colnames(data) && nrow(gene_id_to_symbol_map) > 0) {
                    data[, clean_phenotype_id_genes := gsub("^liver_", "", Phenotype)] # Assuming 'liver_' prefix for gene IDs
                    data_orig_rows <- nrow(data)
                    data <- merge(data, gene_id_to_symbol_map,
                        by.x = "clean_phenotype_id_genes",
                        by.y = "gene_id",
                        all.x = TRUE
                    )
                    if (nrow(data) != data_orig_rows && !isTRUE(all.equal(unique(data$clean_phenotype_id_genes), data$clean_phenotype_id_genes))) {
                        warning(paste("Merge for gene symbols resulted in changed row count for", basename(fst_path), ". Check for many-to-many mapping. Skipping."))
                        return(NULL)
                    }
                    data[!is.na(gene_symbol), Phenotype := gene_symbol]
                    data[, clean_phenotype_id_genes := NULL]
                    data[, gene_symbol := NULL]
                    message("Applied gene symbol mapping for ", basename(fst_path))
                    new_fst_path <- paste0(tools::file_path_sans_ext(fst_path), "_with_symbols.fst")
                } else {
                    message("Skipping gene symbol mapping for ", basename(fst_path), " (Phenotype column missing or no gene annotations).")
                    new_fst_path <- paste0(tools::file_path_sans_ext(fst_path), "_processed.fst") # Still save as processed
                }
            } else if (file_type == "isoforms") {
                if ("Phenotype" %in% colnames(data) && nrow(transcript_id_to_symbol_map) > 0) {
                    # The Phenotype column for raw isoform FSTs should contain transcript_id
                    data_orig_rows <- nrow(data)
                    data <- merge(data, transcript_id_to_symbol_map,
                        by.x = "Phenotype", # Phenotype column is transcript_id here
                        by.y = "transcript_id",
                        all.x = TRUE
                    )
                    if (nrow(data) != data_orig_rows && !isTRUE(all.equal(unique(data$Phenotype[1:min(nrow(data), 10000)]), data$Phenotype[1:min(nrow(data), 10000)]))) { # Check on a subset for performance
                        warning(paste("Merge for transcript symbols resulted in changed row count for", basename(fst_path), ". Check for many-to-many mapping. Skipping."))
                        return(NULL)
                    }
                    data[!is.na(transcript_symbol), Phenotype := transcript_symbol]
                    data[, transcript_symbol := NULL] # Remove the explicitly merged transcript_symbol column
                    message("Applied transcript symbol mapping for ", basename(fst_path))
                    new_fst_path <- paste0(tools::file_path_sans_ext(fst_path), "_with_transcript_symbols.fst")
                } else {
                    message("Skipping transcript symbol mapping for ", basename(fst_path), " (Phenotype column missing or no transcript annotations).")
                    new_fst_path <- paste0(tools::file_path_sans_ext(fst_path), "_processed.fst") # Still save as processed
                }
            } else if (file_type == "clinical") {
                message("Clinical trait file, Phenotype column will be used as is for ", basename(fst_path))
                new_fst_path <- paste0(tools::file_path_sans_ext(fst_path), "_processed.fst")
            } else if (file_type == "liver_lipids") {
                message("Liver lipid trait file, Phenotype column will be used as is for ", basename(fst_path))
                new_fst_path <- paste0(tools::file_path_sans_ext(fst_path), "_processed.fst")
            } else if (file_type == "plasma_metabolites") {
                message("Plasma 2H Metabolites trait file, Phenotype column will be used as is for ", basename(fst_path))
                new_fst_path <- paste0(tools::file_path_sans_ext(fst_path), "_processed.fst")
            } else {
                warning("Unknown file type: ", file_type, " for file: ", basename(fst_path))
                return(NULL)
            }

            write_fst(data, new_fst_path, compress = 50)
            message("Created ", basename(new_fst_path))

            # Sort the data by Phenotype for efficient trait-based reading
            message("Sorting ", basename(new_fst_path), " by Phenotype...")
            setorder(data, Phenotype)
            write_fst(data, new_fst_path, compress = 50)
            message("Successfully sorted and saved: ", basename(new_fst_path))

            create_fst_rows_index(new_fst_path)

            return(new_fst_path)
        },
        error = function(e) {
            warning("Error processing file ", basename(fst_path), ": ", e$message)
            return(NULL)
        }
    )
}

# Main execution
main <- function() {
    gene_anno_file <- "/mnt/rdrive/mkeller3/General/main_directory/files_for_cross_object/gene_annotations.csv"
    gene_data <- read_gene_annotations(gene_anno_file)

    transcript_anno_file <- "/mnt/rdrive/mkeller3/General/main_directory/files_for_cross_object/transcript_annotations.csv"
    transcript_data <- read_transcript_annotations(transcript_anno_file)

    input_dir <- "/data/dev/miniViewer_3.0"
    if (!dir.exists(input_dir)) {
        stop("Input directory not found: ", input_dir)
    }

    file_processing_configs <- list(
        # list(type = "liver_lipids", pattern = "chromosome[0-9XYM]+_liver_lipids_HC_mice_additive_data\.fst$")
        # list(type = "genes", pattern = "chromosome[0-9XYM]+_DO_NOT_USE_data\.fst$"), # Example for gene files, adjust pattern
        # list(type = "isoforms", pattern = "chromosome[0-9XYM]+_liver_isoforms_.*_data\.fst$") # Adjust pattern as needed
        # list(type = "clinical", pattern = "chromosome[0-9XYM]+_clinical_traits_HF_mice_additive_data\.fst$"),
        list(type = "plasma_metabolites", pattern = "chromosome[0-9XYM]+_plasma_metabolites_all_mice_diet_interactive_data\\.fst$")
    )

    all_processed_paths <- character(0)

    for (config in file_processing_configs) {
        message(paste("\\nProcessing type:", config$type, "with pattern:", config$pattern))

        fst_files <- list.files(input_dir, pattern = config$pattern, full.names = TRUE)
        fst_files <- fst_files[!grepl("_with_symbols\\.fst$", fst_files) &
            !grepl("_with_transcript_symbols\\.fst$", fst_files) &
            !grepl("_processed\\.fst$", fst_files)]

        if (length(fst_files) == 0) {
            message("No raw files found for type: ", config$type)
            next
        }

        message("Found ", length(fst_files), " raw ", config$type, " FST files to process.")

        current_gene_map <- if (!is.null(gene_data$id_to_symbol)) gene_data$id_to_symbol else data.table()
        current_transcript_map <- if (!is.null(transcript_data$id_to_symbol)) transcript_data$id_to_symbol else data.table()

        processed_paths_for_type <- sapply(fst_files, process_fst_file,
            gene_id_to_symbol_map = current_gene_map,
            transcript_id_to_symbol_map = current_transcript_map,
            file_type = config$type
        )

        all_processed_paths <- c(all_processed_paths, processed_paths_for_type[!sapply(processed_paths_for_type, is.null)])
    }

    message("\\nFST processing complete. Successfully processed files:")
    if (length(all_processed_paths) > 0) {
        for (p in all_processed_paths) message(p)
    } else {
        message("No files were processed in this run.")
    }
}

# Run the script
main()
