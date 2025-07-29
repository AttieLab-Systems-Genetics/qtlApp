#' Pre-process Phenotype Data for Profile Plot by Category
#'
#' This script reads the "wide" phenotype data, splits it by dataset category using
#' the main annotation list, reshapes each subset into a "long" format, and saves
#' them as separate FST and RDS files for high-performance loading in the Shiny app.
#'
#' @details
#' It performs the following steps:
#' 1. Reads 'annotation_list.rds' which maps traits to trait types.
#' 2. Reads the main 'pheno_with_covar_for_plot.csv'.
#' 3. For each category defined in the annotation list (genes, clinical, etc.):
#'    a. Extracts the list of trait names from the annotation list.
#'    b. Subsets the wide phenotype data to include only those traits.
#'    c. Reshapes the subset from wide to long.
#'    d. Saves the data to 'data/pheno_data_long_<category>.fst' and 'data/trait_names_<category>.rds'.

# Load necessary libraries
if (!require("data.table")) install.packages("data.table")
if (!require("fst")) install.packages("fst")
if (!require("stringr")) install.packages("stringr")

library(data.table)
library(fst)
library(stringr)

# --- Configuration ---
pheno_csv_path <- "/data/dev/miniViewer_3.0/pheno_with_covar_for_plot.csv"
annotation_path <- "/data/dev/miniViewer_3.0/annotation_list.rds"
output_dir <- "/data/dev/miniViewer_3.0"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# --- Step 1: Load Annotation and Phenotype Data ---
message("Loading annotation list from: ", annotation_path)
if (!file.exists(annotation_path)) {
    stop("Fatal: annotation_list.rds not found at the specified path. This file is essential.")
}
annotation_list <- readRDS(annotation_path)

# --- NEW: Consolidate all plasma metabolite types into one category ---
message("Consolidating plasma metabolite annotation data...")
plasma_keys <- c("plasma_metabolites", "plasma_2H_metabolite", "plasma_metabolite")
existing_plasma_keys <- intersect(plasma_keys, names(annotation_list))

if (length(existing_plasma_keys) > 0) {
    # Use rbindlist to efficiently combine the data tables
    combined_plasma_df <- rbindlist(annotation_list[existing_plasma_keys], use.names = TRUE, fill = TRUE)

    # Remove the old, separate keys
    for (key in existing_plasma_keys) {
        annotation_list[[key]] <- NULL
    }

    # Add the new, combined data under a single key and remove any duplicates
    annotation_list[["plasma_metabolites"]] <- unique(combined_plasma_df)
    message("Plasma metabolite categories combined successfully.")
}


message("Reading wide phenotype data from: ", pheno_csv_path)
pheno_dt <- fread(pheno_csv_path)
message("Phenotype data loaded: ", nrow(pheno_dt), " rows, ", ncol(pheno_dt), " columns.")

# Define the mapping from internal annotation list names to user-facing category names
category_map <- c(
    "genes" = "Liver Genes",
    "isoforms" = "Liver Isoforms",
    "clinical" = "Clinical Traits",
    "lipids" = "Liver Lipids",
    "plasma_metabolites" = "Plasma Metabolites"
    # No longer need mappings for the other plasma types as they are now combined
)

# --- Step 2: Process Data for Each Category in the Annotation List ---
all_annotation_types <- names(annotation_list)
message("Found ", length(all_annotation_types), " annotation types to process: ", paste(all_annotation_types, collapse = ", "))

for (trait_type in all_annotation_types) {
    user_facing_category <- category_map[trait_type]

    if (is.null(user_facing_category) || is.na(user_facing_category)) {
        warning("No user-facing category name defined for trait type '", trait_type, "'. Skipping.")
        next
    }

    message("\n--- Processing category: ", user_facing_category, " (type: ", trait_type, ") ---")

    trait_df <- as.data.table(annotation_list[[trait_type]])
    id_vars <- c("Mouse", "GenLit", "Sex", "Diet")

    if (trait_type %in% c("genes", "isoforms")) {
        # Logic for Genes/Isoforms: Match by ID with "liver_" prefix, then map back to symbol
        id_col <- if (trait_type == "genes") "gene.id" else "transcript.id"
        symbol_col <- "symbol"
        prefix <- "liver_"

        # Check if required columns exist in the annotation data
        if (!all(c(id_col, symbol_col) %in% names(trait_df))) {
            warning("Annotation data for '", trait_type, "' is missing required columns ('", id_col, "', '", symbol_col, "'). Skipping.")
            next
        }

        # Create prefixed IDs to search for in the phenotype data
        trait_ids_from_annot <- na.omit(trait_df[[id_col]])
        prefixed_trait_ids <- paste0(prefix, trait_ids_from_annot)

        # Find which of these exist in the phenotype data
        available_prefixed_ids <- intersect(prefixed_trait_ids, colnames(pheno_dt))

        if (length(available_prefixed_ids) == 0) {
            warning("None of the ", length(prefixed_trait_ids), " prefixed trait IDs for category '", user_facing_category, "' were found as columns in the phenotype data file. Skipping.")
            next
        }
        message(length(available_prefixed_ids), " out of ", length(prefixed_trait_ids), " prefixed trait IDs are available in the phenotype data.")

        # Map available prefixed IDs back to their original IDs and then to symbols
        available_unprefixed_ids <- gsub(prefix, "", available_prefixed_ids)
        id_to_symbol_map <- trait_df[get(id_col) %in% available_unprefixed_ids, .SD, .SDcols = c(id_col, symbol_col)]
        id_to_symbol_map <- unique(id_to_symbol_map, by = symbol_col)

        # Subset the wide data using the available prefixed IDs
        subset_cols <- c(id_vars, available_prefixed_ids)
        category_dt_wide <- pheno_dt[, ..subset_cols]

        # Reshape data
        message("Reshaping data from wide to long format...")
        category_dt_long <- melt(category_dt_wide,
            id.vars = id_vars,
            measure.vars = available_prefixed_ids,
            variable.name = "ID_prefixed",
            value.name = "Value",
            na.rm = TRUE
        )

        # Map the prefixed IDs in the long data back to symbols for user-friendliness
        category_dt_long[, (id_col) := gsub(prefix, "", ID_prefixed)]
        setnames(id_to_symbol_map, c(id_col, symbol_col), c(id_col, "Trait_Name"))
        category_dt_long[id_to_symbol_map, on = id_col, Trait_Name := i.Trait_Name]
        category_dt_long[, c("ID_prefixed", id_col) := NULL] # Clean up intermediate columns

        # The final list of traits are the symbols that had corresponding IDs
        available_traits <- sort(unique(na.omit(category_dt_long$Trait_Name)))
    } else {
        # Logic for other data types (clinical, lipids, etc.)
        trait_id_col <- "data_name"
        if (!trait_id_col %in% names(trait_df)) {
            warning("Trait ID column '", trait_id_col, "' not found for '", trait_type, "'. Skipping.")
            next
        }

        category_traits <- na.omit(trait_df[[trait_id_col]])
        available_traits <- intersect(category_traits, colnames(pheno_dt))

        if (length(available_traits) == 0) {
            warning("None of the ", length(category_traits), " traits for category '", user_facing_category, "' were found. Skipping.")
            next
        }
        message(length(available_traits), " traits are available.")

        subset_cols <- c(id_vars, available_traits)
        category_dt_wide <- pheno_dt[, ..subset_cols]

        message("Reshaping data from wide to long format...")
        category_dt_long <- melt(category_dt_wide,
            id.vars = id_vars,
            measure.vars = available_traits,
            variable.name = "Trait_Name",
            value.name = "Value",
            na.rm = TRUE
        )
    }

    message("Reshaped data has ", nrow(category_dt_long), " rows.")

    # Generate sanitized file names
    category_sanitized <- tolower(user_facing_category)
    category_sanitized <- str_replace_all(category_sanitized, "[^a-z0-9]+", "_")

    # Define output paths
    output_fst_path <- file.path(output_dir, paste0("pheno_data_long_", category_sanitized, ".fst"))
    output_rds_path <- file.path(output_dir, paste0("trait_names_", category_sanitized, ".rds"))

    # Save the long-format data and the trait names
    message("Saving long-format data to: ", output_fst_path)
    # Set a key on Trait_Name for ultra-fast reading of specific traits
    setkey(category_dt_long, Trait_Name)
    write_fst(category_dt_long, output_fst_path)

    message("Saving unique trait names to: ", output_rds_path)
    saveRDS(available_traits, file = output_rds_path)
}

message("\nPreprocessing complete!")
