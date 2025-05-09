#!/usr/bin/env Rscript

library(data.table)
library(dplyr)

# Define base data path (consistent with other scripts)
base_path <- "/data/dev/miniViewer_3.0"

# Path to the existing annotation_list.rds
annotation_list_rds_path <- file.path(base_path, "annotation_list.rds")

# Path to the processed allele effects file for clinical traits
# Ensure this filename matches exactly what prepare_allele_data.R creates
allele_effects_csv_path <- file.path(base_path, "consolidate_allele_effects_clinical_traits_all_mice_additive_processed.csv")

message("Attempting to update: ", annotation_list_rds_path)
message("Using clinical trait names from: ", allele_effects_csv_path)

# 1. Read existing annotation_list.rds or initialize if not found
if (file.exists(annotation_list_rds_path)) {
  message("Loading existing annotation_list.rds...")
  annotation_list <- readRDS(annotation_list_rds_path)
} else {
  message("annotation_list.rds not found. Initializing a new list.")
  annotation_list <- list()
}

# 2. Read the processed allele effects CSV
if (!file.exists(allele_effects_csv_path)) {
  stop("Processed allele effects CSV not found: ", allele_effects_csv_path, 
       "\nPlease ensure prepare_allele_data.R has been run successfully.")
}
message("Reading allele effects CSV...")
allele_data <- data.table::fread(allele_effects_csv_path)

# 3. Extract unique, sorted trait names from the 'lodcolumn' column
# This CSV file uses 'lodcolumn' for trait names before symbol mapping (which doesn't happen for clinical).
if (!"lodcolumn" %in% colnames(allele_data)) {
  stop("'lodcolumn' column not found in ", allele_effects_csv_path)
}

clinical_trait_names <- sort(unique(allele_data$lodcolumn))

if (length(clinical_trait_names) == 0) {
  warning("No clinical trait names found in the 'lodcolumn' column of ", allele_effects_csv_path)
} else {
  message(paste("Found", length(clinical_trait_names), "unique clinical trait names. First few:", 
                paste(head(clinical_trait_names), collapse=", ")))
}

# 4. Ensure the structure annotation_list$clinical$data_name exists
if (!("clinical" %in% names(annotation_list))) {
  annotation_list$clinical <- list()
  message("Created $clinical element in annotation_list.")
} else {
  message("Existing annotation_list$clinical structure:")
  print(str(annotation_list$clinical))
}

# 5. Assign the trait names
# The key 'data_name' is what get_trait_id('clinical') and then get_trait_choices() will look for.
# We will overwrite annotation_list$clinical to ensure it has the correct structure and names.
annotation_list$clinical <- data.frame(
  data_name = clinical_trait_names,
  stringsAsFactors = FALSE
)
# If you need to preserve other columns that might have existed in annotation_list$clinical
# and can map them to the new clinical_trait_names, a more complex merge/join would be needed here.
# For now, this ensures the data_name column required by get_trait_choices is correct.

message("Replaced annotation_list$clinical with a new data frame containing updated clinical_trait_names.")

# 6. Save the updated annotation_list
tryCatch({
  saveRDS(annotation_list, file = annotation_list_rds_path)
  message("Successfully updated and saved annotation_list.rds to: ", annotation_list_rds_path)
}, error = function(e) {
  stop("Error saving updated annotation_list.rds: ", e$message)
})

message("Clinical trait names update complete.") 