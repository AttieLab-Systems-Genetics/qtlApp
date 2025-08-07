# Update annotation_list.rds for Clinical Traits and Liver Lipids using CSV headers
# - Reads existing annotation_list.rds
# - Replaces/updates 'clinical' and 'lipids' entries with full column lists from their CSVs
# - Preserves (and does not overwrite) the 'plasma_metabolite' key
# - Fills missing gene symbols for 'genes' by using the corresponding gene.id
# - Backs up the original RDS before writing
#
# Details:
# - Clinical/Lipids: Traits are taken directly from CSV column names (excluding
#   Mouse/GenLit/Sex/Diet), ensuring no traits are lost due to outdated annotations.
# - Genes: If 'symbol' is NA/empty for a row in annotation_list$genes, it will be
#   replaced with the 'gene.id' so users can search by gene id in the app.

if (!require("data.table")) install.packages("data.table")

suppressPackageStartupMessages(library(data.table))

base_dir <- "/data/dev/miniViewer_3.0"
annot_path <- file.path(base_dir, "annotation_list.rds")
backup_path <- file.path(base_dir, paste0("annotation_list_backup_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"))

# CSVs
clinical_csv <- file.path(base_dir, "clinical_traits_pheno_for_plot.csv")
lipids_csv <- file.path(base_dir, "liver_lipids_pheno_for_plot.csv")

if (!file.exists(annot_path)) stop("annotation_list.rds not found: ", annot_path)

message("Loading annotation list: ", annot_path)
annotation_list <- readRDS(annot_path)

# Helper: derive trait names from CSV header
get_trait_names_from_csv <- function(csv_path) {
    if (!file.exists(csv_path)) stop("CSV not found: ", csv_path)
    # Read only header
    header <- names(fread(csv_path, nrows = 0))
    # Remove id/covariate columns
    setdiff(header, c("Mouse", "GenLit", "Sex", "Diet"))
}

# Read names from CSVs
message("Extracting Clinical Traits from CSV header...")
clinical_traits <- get_trait_names_from_csv(clinical_csv)
message("Clinical trait count: ", length(clinical_traits))

message("Extracting Liver Lipids from CSV header...")
lipid_traits <- get_trait_names_from_csv(lipids_csv)
message("Liver lipid trait count: ", length(lipid_traits))

# Ensure list elements exist, then replace with data.frame having data_name column
if (!("clinical" %in% names(annotation_list))) {
    message("Adding missing 'clinical' entry to annotation list")
}
annotation_list[["clinical"]] <- data.frame(data_name = clinical_traits, stringsAsFactors = FALSE)

if (!("lipids" %in% names(annotation_list))) {
    message("Adding missing 'lipids' entry to annotation list")
}
annotation_list[["lipids"]] <- data.frame(data_name = lipid_traits, stringsAsFactors = FALSE)

# Preserve plasma metabolites: keep key named exactly 'plasma_metabolite' if present
if ("plasma_metabolite" %in% names(annotation_list)) {
    message("Preserving existing 'plasma_metabolite' entry (no changes applied)")
} else {
    warning("'plasma_metabolite' key not found in annotation list; leaving as-is")
}

# New: Fill missing gene symbols with gene.id so users can search by id
if ("genes" %in% names(annotation_list)) {
    message("Checking for missing gene symbols in 'genes'...")
    genes_dt <- as.data.table(annotation_list[["genes"]])
    if (all(c("gene.id", "symbol") %in% names(genes_dt))) {
        missing_mask <- is.na(genes_dt$symbol) | genes_dt$symbol == ""
        n_missing <- sum(missing_mask, na.rm = TRUE)
        message("Missing gene symbols before: ", n_missing)
        if (n_missing > 0) {
            genes_dt[missing_mask, symbol := gene.id]
            n_missing_after <- sum(is.na(genes_dt$symbol) | genes_dt$symbol == "", na.rm = TRUE)
            message("Missing gene symbols after fill: ", n_missing_after)
            annotation_list[["genes"]] <- as.data.frame(genes_dt)
        } else {
            message("No missing gene symbols detected.")
        }
    } else {
        warning("'genes' entry missing required columns 'gene.id' and/or 'symbol'; skipping symbol fill.")
    }
} else {
    message("No 'genes' entry found in annotation list; skipping symbol fill.")
}

# Write backup then updated list
message("Backing up original to: ", backup_path)
saveRDS(readRDS(annot_path), backup_path)

message("Writing updated annotation list: ", annot_path)
saveRDS(annotation_list, annot_path)

message("Done updating annotation list.")
