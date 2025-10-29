#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(fst)
})

# Quiet data.table NSE notes
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "data_name", "junction_id", "data_name_lower", "junction_id_lower", "target_phenotype"
    ))
}

# Try to load row index helper if available
try(
    {
        if (file.exists("R/fst_rows.R")) {
            source("R/fst_rows.R")
        } else if (file.exists("../../R/fst_rows.R")) {
            source("../../R/fst_rows.R")
        }
    },
    silent = TRUE
)

parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    base_dir <- "/data/dev/miniViewer_3.0"
    inplace <- TRUE
    allow_large <- TRUE

    for (a in args) {
        if (grepl("^--base_dir=", a)) {
            base_dir <- sub("^--base_dir=", "", a)
        } else if (grepl("^--inplace=", a)) {
            val <- tolower(sub("^--inplace=", "", a))
            inplace <- val %in% c("1", "true", "yes", "y")
        } else if (grepl("^--allow_large=", a)) {
            val <- tolower(sub("^--allow_large=", "", a))
            allow_large <- val %in% c("1", "true", "yes", "y")
        }
    }

    list(base_dir = base_dir, inplace = inplace, allow_large = allow_large)
}

load_splice_map <- function(annotation_list_path) {
    if (!file.exists(annotation_list_path)) {
        stop("annotation_list.rds not found at ", annotation_list_path)
    }
    annotation_list <- readRDS(annotation_list_path)
    splice <- annotation_list[["splice_junctions"]]
    if (is.null(splice) || !is.data.frame(splice)) {
        stop("annotation_list$splice_junctions is missing or not a data.frame")
    }
    req_cols <- c("data_name", "junction_id")
    if (!all(req_cols %in% colnames(splice))) {
        stop("splice_junctions must contain columns: ", paste(req_cols, collapse = ", "))
    }

    map_dt <- as.data.table(splice[, req_cols])
    map_dt[["data_name"]] <- as.character(map_dt[["data_name"]])
    map_dt[["junction_id"]] <- as.character(map_dt[["junction_id"]])
    map_dt[["data_name_lower"]] <- tolower(map_dt[["data_name"]])
    map_dt[["junction_id_lower"]] <- tolower(map_dt[["junction_id"]])
    map_dt[["target_phenotype"]] <- paste0("liver_", map_dt[["data_name"]])
    map_dt
}

find_splice_fsts <- function(base_dir) {
    pats <- c(
        "chromosome[0-9XYM]+_liver_splice_juncs_.*_data_processed_rows\\.fst$",
        "chromosome[0-9XYM]+_liver_splice_juncs_.*_data_processed\\.fst$",
        "chromosome[0-9XYM]+_liver_splice_juncs_.*_data\\.fst$"
    )
    files <- character(0)
    for (p in pats) {
        files <- c(files, list.files(base_dir, pattern = p, full.names = TRUE))
    }
    unique(files)
}

map_phenotypes <- function(phen_vec, map_dt) {
    phen_vec <- as.character(phen_vec)
    phen_lower <- tolower(phen_vec)

    idx1 <- match(phen_lower, map_dt$junction_id_lower)
    matched1 <- !is.na(idx1)

    idx2 <- match(phen_lower, map_dt$data_name_lower)
    matched2 <- !is.na(idx2) & !matched1

    stripped <- sub("^liver_", "", phen_lower)
    idx3 <- match(stripped, map_dt$data_name_lower)
    matched3 <- !is.na(idx3) & !(matched1 | matched2)

    new_phen <- phen_vec
    if (any(matched1)) new_phen[matched1] <- map_dt$target_phenotype[idx1[matched1]]
    if (any(matched2)) new_phen[matched2] <- map_dt$target_phenotype[idx2[matched2]]
    if (any(matched3)) new_phen[matched3] <- map_dt$target_phenotype[idx3[matched3]]

    list(
        new_phen = new_phen,
        num_mapped = sum(matched1 | matched2 | matched3),
        unmatched_examples = unique(utils::head(phen_vec[!(matched1 | matched2 | matched3)], 5))
    )
}

write_out <- function(dt, src_path, inplace) {
    out_path <- if (isTRUE(inplace)) src_path else paste0(tools::file_path_sans_ext(src_path), "_fixed_names.fst")
    write_fst(dt, out_path, compress = 50)
    # Rebuild row index if this is a processed_rows file and helper is available
    if (grepl("_processed_rows\\.fst$", out_path) && exists("create_fst_rows_index", mode = "function")) {
        try(create_fst_rows_index(out_path), silent = TRUE)
    }
    out_path
}

main <- function() {
    cfg <- parse_args()
    base_dir <- cfg$base_dir
    inplace <- cfg$inplace
    allow_large <- cfg$allow_large

    message("Base dir: ", base_dir)
    message("In-place: ", inplace)

    if (!dir.exists(base_dir)) stop("Base dir not found: ", base_dir)

    annot_path <- file.path(base_dir, "annotation_list.rds")
    map_dt <- load_splice_map(annot_path)
    message("Loaded splice junction map with ", nrow(map_dt), " entries.")

    files <- find_splice_fsts(base_dir)
    if (length(files) == 0) {
        message("No splice junction FST files found in ", base_dir)
        quit(status = 0)
    }
    message("Found ", length(files), " splice junction FST files to check/fix.")

    out_paths <- character(0)
    total_mapped <- 0L

    for (f in files) {
        message("\nProcessing ", basename(f), "...")
        # Quickly check available columns without loading entire file
        cols <- tryCatch(names(read_fst(f, from = 1, to = 1)), error = function(e) character(0))
        if (!("Phenotype" %in% cols)) {
            message("  Skipping (no Phenotype column)")
            next
        }

        is_rows <- grepl("_processed_rows\\.fst$", f)
        file_size <- tryCatch(file.info(f)$size, error = function(e) NA_real_)

        if (!is_rows && is.finite(file_size) && file_size > 5e8 && !isTRUE(allow_large)) {
            message("  Phenotype exists but file is large (", round(file_size / 1e9, 2), " GB). Pass --allow_large=yes to process.")
            next
        }

        dt <- tryCatch(read_fst(f, as.data.table = TRUE), error = function(e) NULL)
        if (is.null(dt)) {
            warning("  Failed to read FST: ", f)
            next
        }

        res <- map_phenotypes(dt[["Phenotype"]], map_dt)
        if (res$num_mapped == 0) {
            message("  No changes needed (already canonical or no matches)")
            next
        }

        dt[["Phenotype"]] <- res$new_phen
        out_path <- write_out(dt, f, inplace)
        message("  Mapped ", res$num_mapped, " Phenotype values -> ", basename(out_path))
        if (length(res$unmatched_examples) > 0) {
            message("  Unmatched examples (first 5): ", paste(res$unmatched_examples, collapse = ", "))
        }
        out_paths <- c(out_paths, out_path)
        total_mapped <- total_mapped + res$num_mapped
    }

    message("\nDone. Total phenotypes mapped: ", total_mapped)
    if (length(out_paths) > 0) {
        message("Updated files:")
        for (p in out_paths) message("  ", p)
    }
}

main()
