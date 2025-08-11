# Load required R packages for the application
require("rstudioapi")
require("dplyr")
require("stringr")
require("tidyverse")
require("BiocManager")
require("ggplot2")
require("qtl2")
require("grid")
require("ggrepel")
require("gridGraphics")
require("ggpubr")
require("shiny")
require("shinyFiles")
require("bslib")
require("spsComps")
require("DT")
require("shinyjs")
require("shinycssloaders")
require("data.table")
require("reshape2")
require("plotly")
require("ggiraph")
require("writexl")
require("fontawesome")
require("debounce")

# Set maximum file upload size for Shiny (20GB)
options(shiny.maxRequestSize = 20000 * 1024^2)

# Load required data files
# Read the file index containing paths to all data files
file_directory <- read.csv("/data/dev/miniViewer_3.0/file_index.csv")

# Try to load gene symbols file
gene_symbols_path <- "/data/dev/miniViewer_3.0/gene_symbols.csv"

# Default gene symbols in case loading fails
gene_symbols <- c("Gnai3", "Cdc45", "Slc4a1", "Abca12", "Nadk", "Tfpi", "Scnn1b", "Cdc20", "Gpr89", "Cdc73")

if (file.exists(gene_symbols_path)) {
    message("Loading gene symbols from: ", gene_symbols_path)
    gene_symbols <- tryCatch(
        {
            as.character(fread(gene_symbols_path)$gene_symbol)
        },
        error = function(e) {
            warning("Error reading gene symbols file: ", e$message)
            c("Gnai3", "Cdc45", "Slc4a1", "Abca12", "Nadk", "Tfpi", "Scnn1b", "Cdc20", "Gpr89", "Cdc73")
        }
    )
} else {
    warning("Gene symbols file not found at: ", gene_symbols_path, ". Using default symbols.")
}

# Sort gene symbols
gene_symbols <- sort(gene_symbols)
message("Loaded ", length(gene_symbols), " gene symbols")

# Load chromosome break points for mm11 genome
chr_breaks <- read.csv("/data/dev/miniViewer_3.0/chromosomal_sep_mm11.csv")


# Load gene annotations
annotation_list <- readRDS("/data/dev/miniViewer_3.0/annotation_list.rds")


# Load marker information for the diet DO study
markers <- readRDS(file.path("/data/dev/miniViewer_3.0/CHTC_dietDO_markers_RDSgrcm39.rds"))


# Create a group identifier by combining diet, trait compartment, trait type, and scan type
# Include 'sexes' to differentiate datasets like male-only vs both sexes, but only if not 'Both'
file_directory$group <- paste0(
    file_directory$diet, " ", file_directory$trait_compartment, " ",
    file_directory$trait_type,
    # Conditionally add sex info if it's not 'Both'
    ifelse(file_directory$sexes == "Both", "", paste0(" (", file_directory$sexes, ")")),
    ", ", file_directory$scan_type,
    ifelse(file_directory$scan_type == "interactive",
        paste0(" (", file_directory$covars_interactive, ")"),
        ""
    )
)


# Create a new environment for caching file paths to improve performance
file_path_cache <- new.env(parent = emptyenv())


trait_cache <- new.env(parent = emptyenv())


trait_scan <- function(file_dir, selected_dataset, selected_trait) {
    # Create a more unique cache key that includes the full dataset name
    cache_key <- paste(selected_dataset, tolower(selected_trait), sep = "_")

    if (!is.null(trait_cache[[cache_key]])) {
        message("Using cached data for trait: ", selected_trait, " in dataset: ", selected_dataset)
        return(trait_cache[[cache_key]])
    }

    message("Searching for trait: ", selected_trait, " in dataset: ", selected_dataset)

    # Filter the file directory for the selected dataset and scan type
    file_dir <- subset(file_dir, group == selected_dataset & file_type == "scans")

    # Check if we found any matching files
    if (nrow(file_dir) == 0) {
        stop("No matching files found for the selected dataset: ", selected_dataset)
    }

    # Initialize list to store data from each file
    all_data <- list()

    # Process each FST file (one per chromosome)
    for (i in 1:nrow(file_dir)) {
        # Extract chromosome number from ID_code
        chr_num <- file_dir$ID_code[i]
        fst_path <- file_dir$File_path[i]

        if (!file.exists(fst_path)) {
            warning("File does not exist: ", fst_path)
            next
        }

        # Ensure we are working with an FST file
        if (!stringr::str_detect(fst_path, "fst$")) {
            fst_path <- stringr::str_replace(fst_path, "csv$", "fst")
            if (!file.exists(fst_path)) {
                warning("FST file not found: ", fst_path)
                next
            }
        }

        message("Checking chromosome ", chr_num, " for trait: ", selected_trait, " in dataset: ", selected_dataset)

        # Create row index if it doesn't exist
        row_index_path <- fst_rows(fst_path)

        tryCatch(
            {
                # Read the row index to find the trait
                trait_index <- fst::read_fst(row_index_path, as.data.table = TRUE)

                # Convert Phenotype column to lowercase for case-insensitive matching
                trait_index[, Phenotype := tolower(Phenotype)]

                # Check if the trait is present in this chromosome (case-insensitive)
                trait_rows <- trait_index[Phenotype == tolower(selected_trait), ]

                if (nrow(trait_rows) > 0) {
                    message("Found trait in chromosome ", chr_num, " at rows ", trait_rows$from, "-", trait_rows$to)

                    # Read only the rows for this trait
                    data <- fst::read_fst(fst_path,
                        from = trait_rows$from,
                        to = trait_rows$to,
                        as.data.table = TRUE
                    )

                    # Ensure required columns are present
                    if (!"LOD" %in% colnames(data)) {
                        possible_lod_cols <- grep("lod|LOD|score", colnames(data), ignore.case = TRUE, value = TRUE)
                        if (length(possible_lod_cols) > 0) {
                            setnames(data, possible_lod_cols[1], "LOD")
                        } else {
                            warning("LOD column not found in file: ", fst_path)
                            next
                        }
                    }

                    if (!"marker" %in% colnames(data)) {
                        possible_marker_cols <- grep("marker|id|snp", colnames(data), ignore.case = TRUE, value = TRUE)
                        if (length(possible_marker_cols) > 0) {
                            setnames(data, possible_marker_cols[1], "marker")
                        } else {
                            warning("marker column not found in file: ", fst_path)
                            next
                        }
                    }

                    # Verify that we have the correct trait data (case-insensitive)
                    if ("Phenotype" %in% colnames(data)) {
                        # Double-check that all rows are for the requested trait
                        data <- data[tolower(Phenotype) == tolower(selected_trait)]
                        message("Verified ", nrow(data), " rows for trait: ", selected_trait)
                    }

                    if (nrow(data) > 0) {
                        message("Adding ", nrow(data), " rows from chromosome ", chr_num)
                        all_data[[length(all_data) + 1]] <- data
                    }
                } else {
                    message("Trait not found in chromosome ", chr_num)
                }
            },
            error = function(e) {
                warning("Error processing chromosome ", chr_num, ": ", e$message)
            }
        )
    }

    # Check if we found any data
    if (length(all_data) == 0) {
        stop("Trait '", selected_trait, "' not found in any chromosome for dataset: ", selected_dataset)
    }

    # Combine all data
    combined_data <- data.table::rbindlist(all_data, fill = TRUE)
    message("Total rows in combined data: ", nrow(combined_data))

    # Cache the result
    trait_cache[[cache_key]] <- combined_data
    return(combined_data)
}

# Also add caching for peak_finder
peaks_cache <- new.env(parent = emptyenv())

peak_finder <- function(file_dir, selected_dataset, selected_trait = NULL) {
    # Create a unique cache key that includes the dataset
    cache_key <- if (is.null(selected_trait)) {
        selected_dataset # Include full dataset name
    } else {
        paste(selected_dataset, tolower(selected_trait), sep = "_") # Include full dataset name
    }

    # Check if we already have this data in cache
    if (is.null(peaks_cache[[cache_key]])) {
        message("Loading peaks data for dataset: ", selected_dataset)

        # Filter the file directory for the selected dataset and peaks files
        file_dir <- subset(file_dir, group == selected_dataset & file_type == "peaks")

        # Check if we found any matching files
        if (nrow(file_dir) == 0) {
            warning("No peaks files found for dataset: ", selected_dataset)
            # Return empty data frame with correct structure
            empty_peaks <- data.frame(
                marker = character(0),
                trait = character(0),
                chr = character(0),
                pos = numeric(0),
                lod = numeric(0),
                A = numeric(0),
                B = numeric(0),
                C = numeric(0),
                D = numeric(0),
                E = numeric(0),
                F = numeric(0),
                G = numeric(0),
                H = numeric(0)
            )
            peaks_cache[[cache_key]] <- empty_peaks
            return(empty_peaks)
        }

        # We now have a single consolidated peaks file
        message("Reading consolidated peaks file for dataset: ", selected_dataset)
        csv_path <- file_dir$File_path[1]

        if (!file.exists(csv_path)) {
            warning("Peaks file does not exist: ", csv_path)
            empty_peaks <- data.frame(
                marker = character(0),
                trait = character(0),
                chr = character(0),
                pos = numeric(0),
                lod = numeric(0),
                A = numeric(0),
                B = numeric(0),
                C = numeric(0),
                D = numeric(0),
                E = numeric(0),
                F = numeric(0),
                G = numeric(0),
                H = numeric(0)
            )
            peaks_cache[[cache_key]] <- empty_peaks
            return(empty_peaks)
        }

        tryCatch(
            {
                # Read the consolidated CSV file - this might be large, so use data.table for efficiency
                message("Reading peaks file: ", basename(csv_path), " for dataset: ", selected_dataset)
                peaks_data <- data.table::fread(csv_path)

                # Print out column names for debug
                message("Original peaks file columns: ", paste(colnames(peaks_data), collapse = ", "))

                # Check for required columns and standardize names
                if ("lodcolumn" %in% colnames(peaks_data)) {
                    # We'll keep lodcolumn as is but also add a trait column for compatibility
                    peaks_data$trait <- peaks_data$lodcolumn
                    message("Added trait column based on lodcolumn")
                }
                if (any(grepl("phenotype", tolower(colnames(peaks_data))))) {
                    phenotype_col <- grep("phenotype", tolower(colnames(peaks_data)), value = TRUE)[1]
                    # We'll keep original column and add a trait column
                    peaks_data$trait <- peaks_data[[phenotype_col]]
                    message("Added trait column based on phenotype column")
                }

                # Filter for the specific trait if provided (case-insensitive)
                if (!is.null(selected_trait)) {
                    message("Filtering for trait: ", selected_trait, " in dataset: ", selected_dataset)

                    # Try various columns that might contain the trait
                    rows_to_keep <- rep(FALSE, nrow(peaks_data))

                    # Check trait column
                    if ("trait" %in% colnames(peaks_data)) {
                        message("Checking trait column")
                        rows_to_keep <- rows_to_keep | (tolower(peaks_data$trait) == tolower(selected_trait))
                    }

                    # Check lodcolumn
                    if ("lodcolumn" %in% colnames(peaks_data)) {
                        message("Checking lodcolumn column")
                        rows_to_keep <- rows_to_keep | (tolower(peaks_data$lodcolumn) == tolower(selected_trait))
                    }

                    # Check phenotype columns
                    phenotype_cols <- grep("phenotype", tolower(colnames(peaks_data)), value = TRUE)
                    for (col in phenotype_cols) {
                        message("Checking phenotype column: ", col)
                        rows_to_keep <- rows_to_keep | (tolower(peaks_data[[col]]) == tolower(selected_trait))
                    }

                    # Apply the filter
                    peaks_data <- peaks_data[rows_to_keep]
                    message("Found ", nrow(peaks_data), " peaks for trait: ", selected_trait, " in dataset: ", selected_dataset)
                }

                # Ensure we have all required columns
                required_cols <- c("marker", "chr", "pos", "lod")

                # Check and rename columns if needed
                col_mapping <- list(
                    marker = c("marker", "markers", "id", "snp", "SNP"),
                    chr = c("chr", "chrom", "chromosome"),
                    pos = c("pos", "position", "bp", "location"),
                    lod = c("lod", "LOD", "score", "pvalue")
                )

                # Try to map columns
                for (req_col in names(col_mapping)) {
                    if (!(req_col %in% colnames(peaks_data))) {
                        # Try to find a matching column
                        for (alt_name in col_mapping[[req_col]]) {
                            if (alt_name %in% colnames(peaks_data)) {
                                # Rename the column
                                message("Renaming column '", alt_name, "' to '", req_col, "'")
                                data.table::setnames(peaks_data, alt_name, req_col)
                                break
                            }
                        }
                    }
                }

                # Check if we have all required columns
                if (!all(required_cols %in% colnames(peaks_data))) {
                    missing_cols <- required_cols[!(required_cols %in% colnames(peaks_data))]
                    warning("Missing required columns in peaks file: ", paste(missing_cols, collapse = ", "))

                    # Return empty data frame
                    empty_peaks <- data.frame(
                        marker = character(0),
                        trait = character(0),
                        chr = character(0),
                        pos = numeric(0),
                        lod = numeric(0),
                        A = numeric(0),
                        B = numeric(0),
                        C = numeric(0),
                        D = numeric(0),
                        E = numeric(0),
                        F = numeric(0),
                        G = numeric(0),
                        H = numeric(0)
                    )
                    peaks_cache[[cache_key]] <- empty_peaks
                    return(empty_peaks)
                }

                # Print final column names for debug
                message("Final peaks data columns: ", paste(colnames(peaks_data), collapse = ", "))

                # Convert to data.frame for compatibility
                peaks_data <- as.data.frame(peaks_data)

                # Cache the result
                peaks_cache[[cache_key]] <- peaks_data
            },
            error = function(e) {
                warning("Error reading peaks file ", csv_path, ": ", e$message)
                # Return empty data frame
                empty_peaks <- data.frame(
                    marker = character(0),
                    trait = character(0),
                    chr = character(0),
                    pos = numeric(0),
                    lod = numeric(0),
                    A = numeric(0),
                    B = numeric(0),
                    C = numeric(0),
                    D = numeric(0),
                    E = numeric(0),
                    F = numeric(0),
                    G = numeric(0),
                    H = numeric(0)
                )
                peaks_cache[[cache_key]] <- empty_peaks
            }
        )
    } else {
        message("Using cached peaks data for ", if (is.null(selected_trait)) selected_dataset else paste(selected_trait, "in", selected_dataset))
    }

    return(peaks_cache[[cache_key]])
}


# set microfunctions==========================================================
QTL_plot_visualizer <- function(qtl.temp, phenotype, LOD_thr, mrkrs) {
    # Create a data.table for faster processing
    qtl.temp <- as.data.table(qtl.temp)

    # Select only the columns we need
    if ("marker" %in% colnames(qtl.temp) && "LOD" %in% colnames(qtl.temp)) {
        # If data already has the right columns, just select them
        qtl.temp <- qtl.temp[, .(markers = marker, LOD)]
    } else {
        # Otherwise, try to find the right columns
        qtl.temp <- qtl.temp[, .(markers = marker, LOD)]
    }

    # Convert markers to data.table for faster joins
    mrkrs2 <- as.data.table(mrkrs)[, .(markers = marker, chr, position = bp_grcm39 / 1e6)]

    # Use data.table join instead of merge (much faster)
    qtl.temp <- mrkrs2[qtl.temp, on = "markers", nomatch = NULL]

    # Convert chr to numeric
    qtl.temp[, chr := as.character(chr)]
    qtl.temp[chr == "X", chr := "20"]
    qtl.temp[chr == "Y", chr := "21"]
    qtl.temp[chr == "M", chr := "22"]
    qtl.temp[, chr := as.numeric(chr)]

    # Remove rows with NA in chr
    qtl.temp <- qtl.temp[!is.na(chr)]

    # Add order column (same as chr for numeric processing)
    qtl.temp[, order := chr]

    # Calculate chromosome lengths and cumulative positions
    chr_info <- qtl.temp[, .(chr_len = max(position)), by = chr]
    chr_info[, tot := cumsum(chr_len) - chr_len]

    # Join back to main data
    qtl.temp <- chr_info[qtl.temp, on = "chr", nomatch = NULL]

    # Calculate cumulative position
    qtl.temp[, BPcum := position + tot]

    # Sort by chromosome and position
    setorder(qtl.temp, chr, position)

    # Convert back to tibble for ggplot
    qtl_plot_obj <- as_tibble(qtl.temp)

    return(list(NULL, qtl_plot_obj))
}


# Function to convert CSV files to FST format for better performance
csv2fst <- function(csv_path, chunk_size = 50000) {
    # Check if the input file is a CSV
    if (stringr::str_detect(csv_path, "csv$")) {
        # Create FST filename by replacing .csv with .fst
        fst_path <- stringr::str_replace(csv_path, "csv$", "fst")

        # Only create FST file if it doesn't already exist
        if (!file.exists(fst_path)) {
            warning("Writing FST file in chunks: ", fst_path)
            # Initialize list to store chunks
            all_chunks <- list()
            skip <- 0

            # Read CSV in chunks to handle large files
            repeat {
                chunk <- data.table::fread(csv_path, skip = skip, nrows = chunk_size, showProgress = TRUE)
                if (nrow(chunk) == 0) break
                all_chunks[[length(all_chunks) + 1]] <- chunk
                skip <- skip + chunk_size
            }

            # Combine all chunks
            full_data <- data.table::rbindlist(all_chunks)

            # Sort by Phenotype if the column exists
            if ("Phenotype" %in% colnames(full_data)) {
                data.table::setorder(full_data, Phenotype)
            } else {
                warning("Column 'Phenotype' not found. Data will not be sorted.")
            }

            # Write to FST format with compression
            fst::write_fst(full_data, path = fst_path, compress = 50)
        }
    } else {
        # Check if the file is already an FST file
        if (!stringr::str_detect(csv_path, "fst$")) {
            stop("No CSV or FST name provided: ", csv_path)
        }
        fst_path <- csv_path
    }
    return(fst_path)
}

# Function to create row indices for FST files to speed up data access
fst_rows <- function(fst_path) {
    # Create row index filename
    row_path <- stringr::str_replace(fst_path, ".fst$", "_row.fst")

    # Only create row index if it doesn't exist
    if (!file.exists(row_path)) {
        message("Creating row index for: ", basename(fst_path))
        # Read FST file and create row indices
        rows <- fst::read_fst(fst_path) |>
            # Select only Phenotype column
            dplyr::select(Phenotype) |>
            # Add row numbers
            dplyr::mutate(rown = dplyr::row_number()) |>
            # Group by Phenotype
            dplyr::group_by(Phenotype) |>
            # Get first and last row for each Phenotype
            dplyr::slice(c(1, dplyr::n())) |>
            # Create from/to columns
            dplyr::mutate(set = c("from", "to")) |>
            # Reshape to wide format
            tidyr::pivot_wider(names_from = "set", values_from = "rown")
        # Save row indices
        fst::write_fst(rows, row_path)
        message("Created row index for: ", basename(fst_path))
    }
    return(row_path)
}

# set UI=======================================================================
ui <- fluidPage(
    useShinyjs(),
    # Add custom CSS
    tags$head(
        tags$style(HTML("
  .well {
    background-color: #ffffff;
    border: 1px solid #e3e3e3;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    margin-bottom: 20px;
  }
  /* Add lever switch styling */
  .lever-switch {
    position: relative;
    display: inline-block;
    width: 60px;
    height: 30px;
    background-color: #2c3e50;
    border-radius: 15px;
    cursor: pointer;
    transition: all 0.3s ease;
    overflow: hidden;
  }
  .lever-switch:before {
    content: '';
    position: absolute;
    width: 26px;
    height: 26px;
    left: 2px;
    bottom: 2px;
    background-color: white;
    border-radius: 50%;
    transition: all 0.3s ease;
  }
  .lever-switch.active {
    background-color: #3498db;
  }
  .lever-switch.active:before {
    transform: translateX(30px);
  }
  .lever-switch:hover {
    box-shadow: 0 0 5px rgba(0,0,0,0.2);
  }
  .control-label {
    color: #2c3e50;
    font-weight: 500;
    font-size: 14px;
  }
  .form-control {
    border-radius: 6px;
    border: 1px solid #dce4ec;
  }
  .btn {
    border-radius: 6px;
    text-transform: uppercase;
    font-size: 12px;
    font-weight: 600;
    padding: 8px 16px;
  }
  .btn-default {
    background-color: #3498db;
    color: white;
    border: none;
  }
  .btn-default:hover {
    background-color: #2980b9;
    color: white;
  }
  .selectize-input {
    border-radius: 6px;
    border: 1px solid #dce4ec;
    position: relative;
    z-index: 1;
  }
  .selectize-dropdown {
    z-index: 10000 !important;
  }
  /* Ensure the strain effects selectize appears above the spinner */
  #which_peak + .selectize-control .selectize-dropdown {
    z-index: 99999 !important;
  }
  /* Ensure all selectize dropdowns remain on top */
  .selectize-control {
    position: relative;
    z-index: 10;
  }
  .selectize-dropdown {
    z-index: 99999 !important;
  }
  /* Fix positioning of the spinner relative to the plot */
  .shiny-spinner-output-container {
    position: relative;
    z-index: 0;
  }
  .title-panel {
    background-color: #2c3e50;
    color: white;
    padding: 20px;
    margin-bottom: 30px;
    border-radius: 8px;
  }
  .plot-container {
    background-color: white;
    padding: 20px;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    margin-bottom: 20px;
  }
  .datatable-container {
    background-color: white;
    padding: 20px;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
  }
  #error_message {
    color: #e74c3c;
    padding: 10px;
    margin-top: 10px;
    font-size: 14px;
  }
   /* Add animation for hover effects */
  .well:hover {
    transform: translateY(-2px);
    transition: transform 0.2s ease;
  }
   .btn {
    transition: all 0.3s ease;
  }
   .btn:hover {
    transform: translateY(-1px);
    box-shadow: 0 4px 8px rgba(0,0,0,0.1);
  }
   /* Add pulsing effect to the spinner */
  .spinner-border {
    animation: spinner-pulse 1s ease infinite;
  }
   @keyframes spinner-pulse {
    0% { transform: scale(1); }
    50% { transform: scale(1.1); }
    100% { transform: scale(1); }
  }
   /* Add tooltip styling */
  .tooltip-inner {
    background-color: #2c3e50;
    color: white;
    border-radius: 4px;
    padding: 8px 12px;
  }
   /* Add success message styling */
  .success-message {
    color: #18bc9c;
    padding: 10px;
    margin-top: 10px;
    font-size: 14px;
    display: none;
  }
"))
    ),
    # Use bslib theme with custom settings
    theme = bs_theme(
        version = 4,
        bootswatch = "flatly",
        primary = "#3498db",
        secondary = "#2c3e50",
        success = "#18bc9c",
        info = "#3498db",
        warning = "#f39c12",
        danger = "#e74c3c"
    ),
    # Title Panel with modern styling
    div(
        class = "title-panel",
        h1("Pre-scanned QTL Visualizer for Diet DO Study",
            style = "font-size: 28px; margin: 0;"
        ),
        p("Interactive visualization tool for QTL analysis",
            style = "margin: 10px 0 0 0; opacity: 0.8;"
        )
    ),
    fluidRow(
        # Left column with inputs and allele effects
        column(
            3,
            wellPanel(
                style = "padding: 20px;",
                div(
                    style = "margin-bottom: 25px;",
                    selectizeInput("selected_dataset",
                        h4("Dataset Selection", style = "color: #2c3e50; margin-bottom: 15px;"),
                        choices = unique(file_directory$group),
                        multiple = FALSE,
                        options = list(
                            placeholder = "Select a dataset...",
                            onInitialize = I('function() { this.setValue(""); }')
                        )
                    )
                ),
                div(
                    style = "margin-bottom: 25px;",
                    selectizeInput("which_trait",
                        h4("Search Trait", style = "color: #2c3e50; margin-bottom: 15px;"),
                        choices = gene_symbols,
                        selected = NULL,
                        multiple = FALSE,
                        options = list(
                            placeholder = "Search gene symbol (e.g., Gnai3)",
                            maxOptions = 7,
                            create = FALSE,
                            maxItems = 1
                        )
                    )
                ),
                # Add error message display
                div(
                    id = "error_message_container", style = "margin-bottom: 15px; color: #e74c3c; font-weight: bold; display: none;",
                    textOutput("error_message")
                ),
                # Hide the search button since we'll auto-search
                div(
                    style = "margin-bottom: 25px; display: none;",
                    actionButton("search_trait", "Search",
                        class = "btn-primary",
                        style = "width: 100%;"
                    )
                )
            ),
            # Allele effects panel
            wellPanel(
                style = "padding: 20px; position: relative; overflow: visible;",
                h4("Strain Effects", style = "color: #2c3e50; margin-bottom: 15px;"),
                div(
                    style = "color: #7f8c8d; margin-bottom: 15px;",
                    "Select a peak to see strain effects (Only for additive scans)."
                ),
                div(
                    style = "margin-bottom: 20px; position: relative; z-index: 2;",
                    selectizeInput("which_peak", "Choose peak",
                        choices = NULL,
                        multiple = FALSE,
                        options = list(
                            placeholder = "Select a peak...",
                            onInitialize = I('function() { this.setValue(""); }')
                        )
                    )
                ),
                div(
                    style = "margin-bottom: 20px;",
                    downloadButton("download_effects_plot_png", "Download PNG"),
                    downloadButton("download_effects_plot_pdf", "Download PDF")
                ),
                div(
                    style = "margin-top: 20px; position: relative; z-index: 0;",
                    plotOutput("allele_effects", height = "400px") %>%
                        withSpinner(type = 8, color = "#3498db", proxy.height = "400px")
                )
            )
        ),
        # Right column with tabbed interface
        column(
            9,
            tabsetPanel(
                # LOD Plot Tab
                tabPanel(
                    "LOD Plot",
                    div(
                        style = "margin-bottom: 25px;",
                        sliderInput("LOD_thr",
                            h4("LOD Threshold", style = "color: #2c3e50; margin-bottom: 15px;"),
                            min = 4, max = 120, value = 7, step = 0.5,
                            ticks = TRUE
                        )
                    ),
                    # LOD plot section with clicked point info
                    div(
                        class = "plot-container",
                        div(
                            style = "display: flex; flex-direction: column; gap: 15px; margin-bottom: 20px;",
                            # Title and download buttons row
                            div(
                                style = "display: flex; justify-content: space-between; align-items: center;",
                                h3("LOD Score Plot",
                                    style = "margin: 0; color: #2c3e50; font-weight: 600; font-family: 'Montserrat', 'Helvetica Neue', sans-serif; font-size: 28px; letter-spacing: 0.5px;"
                                ),
                                div(
                                    style = "display: flex; gap: 10px;",
                                    downloadButton("download_qtl_plot_png", "Download PNG", class = "btn-sm"),
                                    downloadButton("download_qtl_plot_pdf", "Download PDF", class = "btn-sm")
                                )
                            ),
                            # Controls row
                            div(
                                style = "display: flex; gap: 20px; align-items: center;",
                                # Chromosome selector
                                div(
                                    style = "display: flex; align-items: center; gap: 10px;",
                                    selectInput("selected_chr", "Zoom to Chromosome:",
                                        choices = c(
                                            "All", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                            "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y", "M"
                                        ),
                                        selected = "All",
                                        width = "150px"
                                    )
                                ),
                                # Plot dimensions
                                div(
                                    style = "display: flex; align-items: center; gap: 15px;",
                                    div(
                                        style = "display: flex; align-items: center; gap: 10px;",
                                        numericInput("plot_width", "Width:",
                                            value = 1000, min = 400, max = 2000, step = 50,
                                            width = "100px"
                                        ),
                                        numericInput("plot_height", "Height:",
                                            value = 600, min = 300, max = 1200, step = 50,
                                            width = "100px"
                                        )
                                    ),
                                    div(
                                        style = "display: flex; gap: 5px;",
                                        actionButton("preset_1to1", "1:1", class = "btn-sm"),
                                        actionButton("preset_3to2", "3:2", class = "btn-sm"),
                                        actionButton("preset_16to9", "16:9", class = "btn-sm")
                                    ),
                                    # Add color toggle button
                                    div(
                                        style = "display: flex; align-items: center; gap: 10px;",
                                        div(
                                            class = "lever-switch", id = "color_toggle",
                                            style = "margin-left: 10px;",
                                            onclick = "Shiny.setInputValue('toggle_colors', Date.now())"
                                        )
                                    )
                                )
                            )
                        )
                    ),
                    plotlyOutput("scan_plot", width = "100%", height = "auto") %>%
                        withSpinner(type = 8, color = "#3498db", proxy.height = "600px"),
                    # Add clicked point info directly below plot
                    div(
                        style = "margin-top: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 4px; border: 1px solid #e9ecef;",
                        h5("Selected Point Information", style = "margin: 0 0 10px 0; color: #2c3e50; font-size: 14px;"),
                        DTOutput("clicked_point_info")
                    )
                ),
                # Mediation Analysis Tab (Placeholder)
                tabPanel(
                    "Mediation Analysis",
                    div(
                        style = "padding: 20px; text-align: center; color: #7f8c8d;",
                        h3("Mediation Analysis", style = "color: #2c3e50;"),
                        p("This feature will be implemented in a future update.",
                            style = "font-size: 16px; margin-top: 20px;"
                        )
                    )
                ),
                # Regression Analysis Tab (Placeholder)
                tabPanel(
                    "Regression Analysis",
                    div(
                        style = "padding: 20px; text-align: center; color: #7f8c8d;",
                        h3("Regression Analysis", style = "color: #2c3e50;"),
                        p("This feature will be implemented in a future update.",
                            style = "font-size: 16px; margin-top: 20px;"
                        )
                    )
                )
            )
        )
    )
)


# Server=======================================================================
server <- function(input, output, session) {
    # Initialize selectize with all gene symbols when the app starts
    # This needs to happen in an observe block
    observe({
        # Only initialize once on startup
        isolate({
            if (is.null(input$which_trait) || input$which_trait == "") {
                # Load initial gene symbols (limited for performance)
                top_genes <- head(gene_symbols, 100)
                updateSelectizeInput(session, "which_trait",
                    choices = top_genes,
                    server = TRUE,
                    selected = NULL
                )
            }
        })
    })

    # Server-side handler for the selectize search
    observeEvent(input$which_trait_search, {
        search_string <- input$which_trait_search
        if (!is.null(search_string) && nchar(search_string) >= 1) {
            # Find matching gene symbols (case-insensitive, match anywhere in string)
            matches <- gene_symbols[grep(search_string, gene_symbols, ignore.case = TRUE)]

            # Limit to top matches
            matches <- head(matches, 7)

            # Return the matches
            if (length(matches) > 0) {
                updateSelectizeInput(session, "which_trait",
                    choices = matches,
                    server = TRUE
                )
            }
        }
    })

    # Reactive value to store scan data
    scan_data <- reactiveVal(NULL)

    # Add reactive value to track click events
    clicked_data <- reactiveVal(NULL)

    # Add reactive value to track the current trait being searched
    current_trait <- reactiveVal("")

    # Create a debounced version of the trait input to prevent too many searches while typing
    # This will trigger the search 800ms after the user stops typing
    trait_debounced <- reactive({
        input$which_trait
    }) %>% debounce(800)

    # Observe the debounced trait input and trigger search automatically
    observeEvent(trait_debounced(), {
        req(input$selected_dataset)
        trait <- trait_debounced()

        # Only search if the trait is not empty and different from the current one
        if (nchar(trait) > 0 && trait != current_trait()) {
            # Update current trait
            current_trait(trait)

            # Clear error message and hide the container
            output$error_message <- renderText("")
            shinyjs::hide("error_message_container")

            # Clear clicked point data when a new trait is searched
            clicked_data(NULL)

            # Try to load the scan data for the specific trait
            tryCatch(
                {
                    message("Searching for trait: ", trait, " in dataset: ", input$selected_dataset)
                    data <- trait_scan(file_directory, input$selected_dataset, trait)
                    scan_data(data)
                },
                error = function(e) {
                    # Show error message
                    output$error_message <- renderText({
                        paste("Error:", e$message)
                    })
                    shinyjs::show("error_message_container")
                    scan_data(NULL)
                }
            )
        }
    })

    # Clear caches when dataset changes
    observeEvent(input$selected_dataset, {
        # Clear all caches to ensure we get fresh data for the new dataset
        rm(list = ls(envir = trait_cache), envir = trait_cache)
        rm(list = ls(envir = peaks_cache), envir = peaks_cache)

        # Remember the current trait
        current_trait_value <- input$which_trait

        # Clear scan data and clicked data
        scan_data(NULL)
        clicked_data(NULL)

        # Clear error message and hide the container
        output$error_message <- renderText("")
        shinyjs::hide("error_message_container")

        # Clear peak selection
        updateSelectizeInput(session, "which_peak", choices = NULL, selected = NULL)

        message("Cleared caches for dataset change to: ", input$selected_dataset)

        # If there was a trait already entered, automatically search for it in the new dataset
        if (nchar(current_trait_value) > 0) {
            message("Automatically searching for trait: ", current_trait_value, " in new dataset: ", input$selected_dataset)

            # Try to load the scan data for the same trait in the new dataset
            tryCatch(
                {
                    data <- trait_scan(file_directory, input$selected_dataset, current_trait_value)
                    current_trait(current_trait_value) # Update the current trait reactive
                    scan_data(data)
                    message("Successfully loaded trait data for new dataset")
                },
                error = function(e) {
                    # Show error message
                    output$error_message <- renderText({
                        paste("Error in new dataset:", e$message)
                    })
                    shinyjs::show("error_message_container")
                    # Don't clear the trait input even if it fails - let the user see what they searched for
                }
            )
        }
    })

    # Keep the search button functionality for backward compatibility
    observeEvent(input$search_trait, {
        req(input$selected_dataset, input$which_trait)

        # Clear error message
        output$error_message <- renderText("")

        # Clear clicked point data when a new trait is searched
        clicked_data(NULL)

        # Update current trait
        current_trait(input$which_trait)

        # Try to load the scan data for the specific trait
        tryCatch(
            {
                data <- trait_scan(file_directory, input$selected_dataset, input$which_trait)
                scan_data(data)
            },
            error = function(e) {
                output$error_message <- renderText({
                    paste("Error:", e$message)
                })
                scan_data(NULL)
            }
        )
    })

    # Also clear clicked data when trait input changes
    observeEvent(input$which_trait, {
        clicked_data(NULL)
    })

    # Main reactive expression for plot data
    plot_obj <- reactive({
        req(scan_data())
        QTL_plot_visualizer(scan_data(), current_trait(), input$LOD_thr, markers)
    })

    # Add reactive value to store the official gene symbol
    official_gene_symbol <- reactiveVal("")

    # Add reactive value to track color mode
    use_alternating_colors <- reactiveVal(TRUE)

    # Add observer for color toggle button
    observeEvent(input$toggle_colors, {
        use_alternating_colors(!use_alternating_colors())
        # Toggle the active class on the lever switch based on the new state
        runjs(sprintf("
        const toggle = document.getElementById('color_toggle');
        if (%s) {
          toggle.classList.add('active');
        } else {
          toggle.classList.remove('active');
        }
      ", tolower(use_alternating_colors())))
    })

    # Add observers for preset buttons
    observeEvent(input$preset_1to1, {
        updateNumericInput(session, "plot_width", value = 800)
        updateNumericInput(session, "plot_height", value = 800)
    })

    observeEvent(input$preset_3to2, {
        updateNumericInput(session, "plot_width", value = 900)
        updateNumericInput(session, "plot_height", value = 600)
    })

    observeEvent(input$preset_16to9, {
        updateNumericInput(session, "plot_width", value = 1600)
        updateNumericInput(session, "plot_height", value = 900)
    })

    # Update official gene symbol and automatically select highest peak when scan data changes
    observeEvent(scan_data(), {
        if (!is.null(scan_data()) && "Phenotype" %in% colnames(scan_data())) {
            # Set official gene symbol
            official_gene_symbol(unique(scan_data()$Phenotype)[1])

            # Automatically find and select the highest peak when a trait is loaded
            if (nrow(scan_data()) > 0) {
                message("Finding highest peak for trait: ", official_gene_symbol())

                # Get peaks data for the selected trait
                peaks_data <- peak_finder(file_directory, input$selected_dataset, input$which_trait)

                # Make sure we have data and the required columns
                if (nrow(peaks_data) > 0 && all(c("marker", "lod") %in% colnames(peaks_data))) {
                    # Filter by LOD threshold and sort to find highest
                    filtered_peaks <- peaks_data %>%
                        filter(lod >= input$LOD_thr) %>%
                        arrange(desc(lod))

                    if (nrow(filtered_peaks) > 0) {
                        # Get the highest peak
                        highest_peak <- filtered_peaks$marker[1]
                        message("Found highest peak: ", highest_peak, " with LOD: ", round(filtered_peaks$lod[1], 2))

                        # Update dropdown to select this peak
                        updateSelectizeInput(session, "which_peak", selected = highest_peak)

                        # Find this peak in the plot data
                        plot_data <- plot_obj()[[2]]
                        peak_point <- plot_data %>% filter(markers == highest_peak)

                        if (nrow(peak_point) > 0) {
                            # Create a fake click event at this peak's position
                            fake_click <- list(
                                x = if (input$selected_chr == "All") peak_point$BPcum[1] else peak_point$position[1],
                                y = peak_point$LOD[1],
                                curveNumber = 0,
                                pointNumber = which(plot_data$markers == highest_peak)[1]
                            )

                            # Update clicked data to automatically show point info for highest peak
                            clicked_data(fake_click)
                            message("Auto-selected highest peak with LOD: ", round(filtered_peaks$lod[1], 2))
                        } else {
                            message("Could not find highest peak marker in plot data")
                        }
                    } else {
                        message("No peaks found above threshold (", input$LOD_thr, ") for this trait")
                    }
                } else {
                    if (nrow(peaks_data) == 0) {
                        message("No peaks found for this trait")
                    } else {
                        message(
                            "Required columns missing from peaks data. Have: ",
                            paste(colnames(peaks_data), collapse = ", ")
                        )
                    }
                }
            }
        }
    })

    # Handle clicked points display
    observeEvent(event_data("plotly_click", source = "scan_plot"), {
        clicked_data(event_data("plotly_click", source = "scan_plot"))
    })

    # Also update clicked_data when a peak is selected in the dropdown
    observeEvent(input$which_peak, {
        req(input$which_peak, plot_obj())

        # Get the peak data
        peaks_data <- peak_finder(file_directory, input$selected_dataset, input$which_trait)
        selected_peak <- peaks_data %>% filter(marker == input$which_peak)

        if (nrow(selected_peak) > 0) {
            # Find this peak in the plot data
            plot_data <- plot_obj()[[2]]
            peak_point <- plot_data %>% filter(markers == selected_peak$marker[1])

            if (nrow(peak_point) > 0) {
                # Create a fake click event at this peak's position
                fake_click <- list(
                    x = if (input$selected_chr == "All") peak_point$BPcum[1] else peak_point$position[1],
                    y = peak_point$LOD[1],
                    curveNumber = 0,
                    pointNumber = which(plot_data$markers == selected_peak$marker[1])[1]
                )

                # Update the clicked data as if the plot was clicked
                clicked_data(fake_click)

                # Log that we're updating the clicked data
                message("Updated clicked data for peak: ", input$which_peak)
            }
        }
    })

    # Main reactive expression for plot data
    plot_base <- reactive({
        req(plot_obj())
        plot_data <- plot_obj()[[2]]

        # Filter data based on selected chromosome
        if (input$selected_chr != "All") {
            chr_num <- switch(input$selected_chr,
                "X" = 20,
                "Y" = 21,
                "M" = 22,
                as.numeric(input$selected_chr)
            )
            plot_data <- plot_data %>% filter(chr == chr_num)
        }

        # Calculate chromosome axis positions
        if (input$selected_chr == "All") {
            # For all chromosomes view
            axisdf <- plot_data %>%
                group_by(chr) %>%
                summarise(center = mean(BPcum))

            # Convert chromosome numbers to proper labels
            chr_labels <- as.character(axisdf$chr)
            chr_labels[chr_labels == "20"] <- "X"
            chr_labels[chr_labels == "21"] <- "Y"
            chr_labels[chr_labels == "22"] <- "M"

            # Create the base ggplot object for all chromosomes
            p <- ggplot(plot_data, aes(x = BPcum, y = LOD)) +
                geom_line(aes(color = if (use_alternating_colors()) as.factor(chr) else NULL), size = 0.75) +
                geom_hline(
                    yintercept = input$LOD_thr, color = "#e74c3c",
                    linetype = "dashed", size = 0.8
                ) +
                scale_color_manual(values = if (use_alternating_colors()) rep(c("#3498db", "#2c3e50"), 22) else "#3498db") +
                scale_x_continuous(
                    breaks = axisdf$center,
                    labels = chr_labels,
                    expand = expansion(mult = c(0.01, 0.01))
                )
        } else {
            # For single chromosome view
            # Create the base ggplot object for a single chromosome
            p <- ggplot(plot_data, aes(x = position, y = LOD)) +
                geom_line(color = "#3498db", size = 0.75) +
                geom_hline(
                    yintercept = input$LOD_thr, color = "#e74c3c",
                    linetype = "dashed", size = 0.8
                ) +
                scale_x_continuous(
                    expand = expansion(mult = c(0.02, 0.02))
                )
        }

        # Common theme and y-axis settings for both views
        p <- p +
            scale_y_continuous(
                expand = expansion(mult = c(0.02, 0.1))
            ) +
            theme_minimal() +
            theme(
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "#ecf0f1", size = 0.2),
                axis.line = element_line(color = "#2c3e50", size = 0.5),
                axis.text = element_text(size = 11, color = "#2c3e50"),
                axis.title = element_text(size = 12, face = "bold", color = "#2c3e50"),
                axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
                legend.position = "none",
                plot.title = element_text(size = 14, face = "bold", color = "#2c3e50"),
                plot.subtitle = element_text(size = 11, color = "#7f8c8d"),
                plot.background = element_rect(fill = "white", color = NA),
                panel.background = element_rect(fill = "white", color = NA),
                plot.margin = margin(t = 20, r = 20, b = 40, l = 40)
            ) +
            labs(
                x = if (input$selected_chr == "All") "Chromosome" else paste("Position on Chromosome", input$selected_chr, "(Mb)"),
                y = "LOD Score"
            )

        return(list(p = p, data = plot_data))
    })

    # Create the interactive plotly plot
    output$scan_plot <- renderPlotly({
        req(plot_base(), input$which_trait)

        # Get the base plot and data
        plot_result <- plot_base()
        p <- plot_result$p
        plot_data <- plot_result$data

        # Show the highest LOD peak
        peaks_info <- peak_finder(file_directory, input$selected_dataset, input$which_trait)
        if (nrow(peaks_info) > 0) {
            peaks_info <- peaks_info %>%
                arrange(desc(lod)) %>%
                slice(1)

            peak_point <- plot_data %>%
                filter(markers == peaks_info$marker)

            if (nrow(peak_point) > 0) {
                # Add the red diamond at the peak
                if (input$selected_chr == "All") {
                    p <- p + geom_point(
                        data = peak_point,
                        aes(x = BPcum, y = LOD),
                        color = "red",
                        size = 3,
                        shape = 20
                    )
                } else {
                    p <- p + geom_point(
                        data = peak_point,
                        aes(x = position, y = LOD),
                        color = "red",
                        size = 3,
                        shape = 20
                    )
                }
            }
        }

        # Create formatted trait text for plot title using the official gene symbol
        trait_text <- paste0("<b style='font-size: 24px;'>", official_gene_symbol(), "</b>")

        # Create subtitle with peak information
        subtitle <- if (nrow(peaks_info) > 0) {
            chr_label <- if (peaks_info$chr %in% c(20, 21, 22)) {
                c("X", "Y", "M")[peaks_info$chr - 19]
            } else {
                peaks_info$chr
            }
            paste0(
                "<span style='font-size: 16px;'>",
                "<b>Peak Marker:</b> ", peaks_info$marker,
                " (Chr", chr_label, ":", round(peaks_info$pos, 2), " Mb) | ",
                "<b>LOD:</b> ", round(peaks_info$lod, 2),
                "</span>"
            )
        } else {
            "<span style='font-size: 16px; color: #7f8c8d;'>No significant peaks</span>"
        }

        # Convert ggplot to plotly with custom dimensions and removed features
        plt <- ggplotly(p,
            source = "scan_plot", width = input$plot_width, height = input$plot_height,
            tooltip = c("x", "y", "chr")
        ) %>%
            layout(
                title = list(
                    text = paste0(trait_text, "<br>", subtitle),
                    font = list(family = "Arial"),
                    x = 0,
                    xanchor = "left",
                    y = 0.95,
                    yanchor = "top",
                    pad = list(b = 20)
                ),
                margin = list(t = 80),
                hoverlabel = list(
                    bgcolor = "white",
                    font = list(family = "Arial", size = 12, color = "#2c3e50"),
                    bordercolor = "#95a5a6"
                ),
                hovermode = "closest",
                # Add double click event to reset to full view
                doubleclick = if (input$selected_chr != "All") TRUE else FALSE
            )

        # Remove unwanted modebar buttons
        plt <- plt %>% config(
            displaylogo = FALSE,
            modeBarButtonsToRemove = c(
                "select2d", "lasso2d", "autoScale2d", "hoverClosestCartesian",
                "hoverCompareCartesian", "toggleSpikelines"
            )
        )

        return(plt)
    })

    # Add observer for plotly double click event
    observeEvent(event_data("plotly_doubleclick", source = "scan_plot"), {
        if (input$selected_chr != "All") {
            updateSelectInput(session, "selected_chr", selected = "All")
        }
    })

    # Update which_peak dropdown when trait is found
    observe({
        req(input$selected_dataset, input$which_trait, input$LOD_thr)

        # Get peaks data for the selected trait from the peaks CSV files
        peaks_data <- peak_finder(file_directory, input$selected_dataset, input$which_trait) %>%
            filter(lod >= input$LOD_thr) %>%
            arrange(desc(lod)) # Sort by LOD score descending

        if (nrow(peaks_data) > 0) {
            # Create named vector for dropdown
            marker_choices <- peaks_data$marker
            names(marker_choices) <- paste0(
                peaks_data$marker,
                " (Chr", peaks_data$chr, ": ", round(peaks_data$pos, 2),
                " Mb, LOD: ", round(peaks_data$lod, 2), ")"
            )

            # Select the maximum LOD peak by default
            updateSelectizeInput(session, "which_peak",
                choices = marker_choices,
                selected = marker_choices[1] # First one is max LOD due to arrange(desc(lod))
            )
        } else {
            updateSelectizeInput(session, "which_peak",
                choices = character(0),
                selected = NULL
            )
        }
    })

    output$clicked_point_info <- renderDT({
        event_data <- clicked_data()

        # If no click data, try to get information from the selected peak instead
        if (is.null(event_data) && !is.null(input$which_peak)) {
            message("No click data, using selected peak: ", input$which_peak)
            peaks_data <- peak_finder(file_directory, input$selected_dataset, input$which_trait)

            # Debug the structure of peaks_data
            message("Peaks data columns: ", paste(colnames(peaks_data), collapse = ", "))
            message("Number of rows in peaks_data: ", nrow(peaks_data))

            # Find peak by marker, considering possible column name differences
            peak_match <- NULL
            if ("marker" %in% colnames(peaks_data)) {
                peak_match <- peaks_data %>% filter(marker == input$which_peak)
            } else if ("markers" %in% colnames(peaks_data)) {
                peak_match <- peaks_data %>% filter(markers == input$which_peak)
            }

            if (!is.null(peak_match) && nrow(peak_match) > 0) {
                message("Found matching peak information row")

                # Debug the matched peak data
                message("Peak match columns: ", paste(colnames(peak_match), collapse = ", "))

                # Create peak info table with safer column access
                peak_info <- data.frame(
                    Marker = ifelse("marker" %in% colnames(peak_match),
                        peak_match$marker[1],
                        ifelse("markers" %in% colnames(peak_match),
                            peak_match$markers[1],
                            NA
                        )
                    ),
                    Chromosome = if ("chr" %in% colnames(peak_match)) {
                        chr_val <- peak_match$chr[1]
                        if (is.numeric(chr_val) && chr_val %in% c(20, 21, 22)) {
                            c("X", "Y", "M")[chr_val - 19]
                        } else {
                            chr_val
                        }
                    } else {
                        NA
                    },
                    Position = if ("pos" %in% colnames(peak_match)) round(peak_match$pos[1], 3) else NA,
                    LOD = if ("lod" %in% colnames(peak_match)) {
                        round(peak_match$lod[1], 3)
                    } else if ("LOD" %in% colnames(peak_match)) {
                        round(peak_match$LOD[1], 3)
                    } else {
                        NA
                    }
                )

                # Add trait column if available
                if ("trait" %in% colnames(peak_match)) {
                    peak_info$Trait <- peak_match$trait[1]
                } else if ("lodcolumn" %in% colnames(peak_match)) {
                    peak_info$Trait <- peak_match$lodcolumn[1]
                }

                # Add Cis column if available
                if ("cis" %in% colnames(peak_match)) {
                    peak_info$Cis <- peak_match$cis[1]
                }

                # Add confidence interval information if available
                if (all(c("ci_lo", "ci_hi") %in% colnames(peak_match))) {
                    peak_info$CI_Low <- round(peak_match$ci_lo[1], 3)
                    peak_info$CI_High <- round(peak_match$ci_hi[1], 3)
                }

                # Add strain effects if available
                if (all(c("A", "B", "C", "D", "E", "F", "G", "H") %in% colnames(peak_match))) {
                    # Add strain effects as columns, with safety checks
                    peak_info$AJ <- round(peak_match$A[1], 3)
                    peak_info$B6 <- round(peak_match$B[1], 3)
                    peak_info$`129` <- round(peak_match$C[1], 3)
                    peak_info$NOD <- round(peak_match$D[1], 3)
                    peak_info$NZO <- round(peak_match$E[1], 3)
                    peak_info$CAST <- round(peak_match$F[1], 3)
                    peak_info$PWK <- round(peak_match$G[1], 3)
                    peak_info$WSB <- round(peak_match$H[1], 3)
                }

                # Return the table
                return(datatable(
                    peak_info,
                    options = list(
                        dom = "t",
                        ordering = FALSE,
                        pageLength = 1
                    ),
                    rownames = FALSE,
                    class = "compact hover",
                    caption = htmltools::tags$caption(
                        style = "caption-side: top; text-align: left; color: #2c3e50; font-weight: bold; font-size: 14px;",
                        "Peak Information and Strain Effects"
                    )
                ))
            } else {
                message("No matching peak found in peaks data")
            }
        }

        if (is.null(event_data)) {
            return(NULL)
        }

        plot_data <- plot_obj()[[2]]
        if (is.null(plot_data)) {
            return(NULL)
        }

        # Find nearest point
        x_clicked <- event_data$x
        y_clicked <- event_data$y

        distances <- sqrt((plot_data$BPcum - x_clicked)^2 + (plot_data$LOD - y_clicked)^2)
        nearest_point <- plot_data[which.min(distances), ]

        if (nrow(nearest_point) == 0) {
            return(NULL)
        }

        # Check if this point is in the peaks file to get additional information
        peaks_data <- peak_finder(file_directory, input$selected_dataset, input$which_trait)

        # Debug
        message("Looking for marker: ", nearest_point$markers, " in peaks data")
        message("Peaks data columns: ", paste(colnames(peaks_data), collapse = ", "))

        # Find if the clicked marker is in the peaks data (check both marker and markers)
        peak_match <- NULL
        if ("marker" %in% colnames(peaks_data) && nrow(peaks_data) > 0) {
            peak_match <- peaks_data %>% filter(marker == nearest_point$markers)
        } else if ("markers" %in% colnames(peaks_data) && nrow(peaks_data) > 0) {
            peak_match <- peaks_data %>% filter(markers == nearest_point$markers)
        }

        if (!is.null(peak_match) && nrow(peak_match) > 0) {
            message("Found matching peak in peaks data")

            # Create base peak info with safer column access
            peak_info <- data.frame(
                Marker = ifelse("marker" %in% colnames(peak_match),
                    peak_match$marker[1],
                    ifelse("markers" %in% colnames(peak_match),
                        peak_match$markers[1],
                        nearest_point$markers
                    )
                ),
                Chromosome = if ("chr" %in% colnames(peak_match)) {
                    chr_val <- peak_match$chr[1]
                    if (is.numeric(chr_val) && chr_val %in% c(20, 21, 22)) {
                        c("X", "Y", "M")[chr_val - 19]
                    } else {
                        chr_val
                    }
                } else {
                    if (nearest_point$chr %in% c(20, 21, 22)) {
                        c("X", "Y", "M")[nearest_point$chr - 19]
                    } else {
                        nearest_point$chr
                    }
                },
                Position = if ("pos" %in% colnames(peak_match)) {
                    round(peak_match$pos[1], 3)
                } else {
                    round(nearest_point$position, 3)
                },
                LOD = if ("lod" %in% colnames(peak_match)) {
                    round(peak_match$lod[1], 3)
                } else if ("LOD" %in% colnames(peak_match)) {
                    round(peak_match$LOD[1], 3)
                } else {
                    round(nearest_point$LOD, 3)
                }
            )

            # Add trait column if available
            if ("trait" %in% colnames(peak_match)) {
                peak_info$Trait <- peak_match$trait[1]
            } else if ("lodcolumn" %in% colnames(peak_match)) {
                peak_info$Trait <- peak_match$lodcolumn[1]
            }

            # Add Cis column if available
            if ("cis" %in% colnames(peak_match)) {
                peak_info$Cis <- peak_match$cis[1]
            }

            # Add confidence interval information if available
            if (all(c("ci_lo", "ci_hi") %in% colnames(peak_match))) {
                peak_info$CI_Low <- round(peak_match$ci_lo[1], 3)
                peak_info$CI_High <- round(peak_match$ci_hi[1], 3)
            }

            # Add strain effects if available
            if (all(c("A", "B", "C", "D", "E", "F", "G", "H") %in% colnames(peak_match))) {
                # Add strain effects as columns
                peak_info$AJ <- round(peak_match$A[1], 3)
                peak_info$B6 <- round(peak_match$B[1], 3)
                peak_info$`129` <- round(peak_match$C[1], 3)
                peak_info$NOD <- round(peak_match$D[1], 3)
                peak_info$NZO <- round(peak_match$E[1], 3)
                peak_info$CAST <- round(peak_match$F[1], 3)
                peak_info$PWK <- round(peak_match$G[1], 3)
                peak_info$WSB <- round(peak_match$H[1], 3)
            }

            # Return the wide format table
            return(datatable(
                peak_info,
                options = list(
                    dom = "t",
                    ordering = FALSE,
                    pageLength = 1
                ),
                rownames = FALSE,
                class = "compact hover",
                caption = htmltools::tags$caption(
                    style = "caption-side: top; text-align: left; color: #2c3e50; font-weight: bold; font-size: 14px;",
                    "Peak Information and Strain Effects"
                )
            ))
        } else {
            message("No matching peak found in peaks data for clicked point")
        }

        # If no peak match found, return basic point info
        point_info <- data.frame(
            Marker = nearest_point$markers,
            Chromosome = if (nearest_point$chr %in% c(20, 21, 22)) {
                c("X", "Y", "M")[nearest_point$chr - 19]
            } else {
                nearest_point$chr
            },
            Position = round(nearest_point$position, 3),
            LOD = round(nearest_point$LOD, 3)
        )

        datatable(
            point_info,
            options = list(
                dom = "t",
                ordering = FALSE,
                pageLength = 1
            ),
            rownames = FALSE,
            class = "compact hover"
        )
    })

    # Automatically select first dataset on startup
    observe({
        req(input$selected_dataset == "")
        first_dataset <- unique(file_directory$group)[1]
        updateSelectizeInput(session, "selected_dataset", selected = first_dataset)
    })

    # Handle allele effects plot
    observeEvent(c(input$which_peak, input$which_trait), {
        req(input$selected_dataset, input$which_trait, input$which_peak)

        # Get peaks data from CSV files
        peaks <- peak_finder(file_directory, input$selected_dataset, input$which_trait) %>%
            filter(marker == input$which_peak)

        if (nrow(peaks) > 0 && all(c("A", "B", "C", "D", "E", "F", "G", "H") %in% colnames(peaks))) {
            peaks_plot <- peaks[1, c("marker", "A", "B", "C", "D", "E", "F", "G", "H")]
            colnames(peaks_plot)[2:9] <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
            peaks_plot <- reshape2::melt(peaks_plot, id.vars = "marker")

            newClrs <- c(
                "AJ" = "#000000", "B6" = "#96989A", "129" = "#E69F00",
                "NOD" = "#0072B2", "NZO" = "#619BFF", "CAST" = "#009E73",
                "PWK" = "#D55E00", "WSB" = "#CC79A7"
            )

            # Get the trait name - use official symbol if available, otherwise use input
            trait_name <- if (nchar(official_gene_symbol()) > 0) {
                official_gene_symbol()
            } else {
                input$which_trait
            }

            plot_alleles <- ggplot(data = peaks_plot, aes(x = marker, y = value, color = variable)) +
                geom_point(size = 4, alpha = 0.8) +
                scale_color_manual(values = newClrs) +
                theme_minimal() +
                theme(
                    legend.position = "right",
                    legend.title = element_text(size = 12, face = "bold"),
                    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5) # Make x-axis text horizontal
                ) +
                labs(
                    x = "Marker ID", y = "Founder Allele Effect", color = "Strain",
                    title = paste0("Strain Effects at ", input$which_peak),
                    subtitle = paste0("Trait: ", trait_name)
                )

            output$allele_effects <- renderPlot({
                plot_alleles
            })

            # Store the plot in a reactive value for download handlers
            strain_plot(plot_alleles)
        } else {
            # If no strain effects data available
            output$allele_effects <- renderPlot({
                plot.new()
                text(0.5, 0.5, "No strain effects data available for this peak", cex = 1.2)
            })
            strain_plot(NULL)
        }
    })


    # Add reactive value to store the strain effects plot
    strain_plot <- reactiveVal(NULL)


    # Download handler for strain effects PNG
    output$download_effects_plot_png <- downloadHandler(
        filename = function() {
            # Use official gene symbol instead of user input for consistency
            actual_trait_name <- official_gene_symbol()
            paste0("strain_effects_", actual_trait_name, "_", input$which_peak, "_", format(Sys.time(), "%Y%m%d"), ".png")
        },
        content = function(file) {
            req(strain_plot())
            ggsave(file, strain_plot(), width = 10, height = 7, dpi = 300)
        }
    )

    # Download handler for strain effects PDF
    output$download_effects_plot_pdf <- downloadHandler(
        filename = function() {
            # Use official gene symbol instead of user input for consistency
            actual_trait_name <- official_gene_symbol()
            paste0("strain_effects_", actual_trait_name, "_", input$which_peak, "_", format(Sys.time(), "%Y%m%d"), ".pdf")
        },
        content = function(file) {
            req(strain_plot())
            ggsave(file, strain_plot(), width = 10, height = 7, device = cairo_pdf)
        }
    )

    # Download handler for PNG
    output$download_qtl_plot_png <- downloadHandler(
        filename = function() {
            # Include chromosome info in filename if specific chromosome is selected
            chr_suffix <- if (input$selected_chr != "All") paste0("_chr", input$selected_chr) else ""
            paste0("lod_plot_", input$which_trait, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".png")
        },
        content = function(file) {
            # Get the base plot - this already has the chromosome filtering applied
            plot_result <- plot_base()
            p <- plot_result$p
            plot_data <- plot_result$data

            # Get peak information for the current trait
            peaks_info <- peak_finder(file_directory, input$selected_dataset, input$which_trait) %>%
                filter(lod >= input$LOD_thr) %>%
                arrange(desc(lod))

            # Filter peaks to only show those in the selected chromosome if applicable
            if (input$selected_chr != "All") {
                chr_num <- switch(input$selected_chr,
                    "X" = 20,
                    "Y" = 21,
                    "M" = 22,
                    as.numeric(input$selected_chr)
                )
                peaks_info <- peaks_info %>% filter(chr == chr_num)
            }

            # Get the highest peak (if any)
            if (nrow(peaks_info) > 0) {
                peaks_info <- peaks_info %>% slice(1)

                # Find the peak point in the filtered data
                peak_point <- plot_data %>%
                    filter(markers == peaks_info$marker)

                if (nrow(peak_point) > 0) {
                    # Add the red diamond at the peak
                    if (input$selected_chr == "All") {
                        p <- p + geom_point(
                            data = peak_point,
                            aes(x = BPcum, y = LOD),
                            color = "red",
                            size = 3,
                            shape = 18
                        )
                    } else {
                        p <- p + geom_point(
                            data = peak_point,
                            aes(x = position, y = LOD),
                            color = "red",
                            size = 3,
                            shape = 18
                        )
                    }
                }
            }

            # Create plot title and subtitle
            trait_title <- input$which_trait
            chr_info <- if (input$selected_chr != "All") paste0(" (Chromosome ", input$selected_chr, ")") else ""

            subtitle <- if (nrow(peaks_info) > 0) {
                chr_label <- if (peaks_info$chr %in% c(20, 21, 22)) {
                    c("X", "Y", "M")[peaks_info$chr - 19]
                } else {
                    peaks_info$chr
                }
                paste0(
                    "Peak Marker: ", peaks_info$marker,
                    " (Chr", chr_label, ":", round(peaks_info$pos, 2), " Mb) | ",
                    "LOD: ", round(peaks_info$lod, 2)
                )
            } else {
                if (input$selected_chr != "All") {
                    "No significant peaks in this chromosome"
                } else {
                    "No significant peaks"
                }
            }

            # Add title and subtitle to the plot
            p <- p + ggtitle(
                label = paste0(trait_title, chr_info),
                subtitle = subtitle
            )

            # Save the plot with high resolution
            ggsave(file, p, width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300)
        }
    )

    # Download handler for PDF
    output$download_qtl_plot_pdf <- downloadHandler(
        filename = function() {
            # Include chromosome info in filename if specific chromosome is selected
            chr_suffix <- if (input$selected_chr != "All") paste0("_chr", input$selected_chr) else ""
            paste0("lod_plot_", input$which_trait, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".pdf")
        },
        content = function(file) {
            # Get the base plot - this already has the chromosome filtering applied
            plot_result <- plot_base()
            p <- plot_result$p
            plot_data <- plot_result$data

            # Get peak information for the current trait
            peaks_info <- peak_finder(file_directory, input$selected_dataset, input$which_trait) %>%
                filter(lod >= input$LOD_thr) %>%
                arrange(desc(lod))

            # Filter peaks to only show those in the selected chromosome if applicable
            if (input$selected_chr != "All") {
                chr_num <- switch(input$selected_chr,
                    "X" = 20,
                    "Y" = 21,
                    "M" = 22,
                    as.numeric(input$selected_chr)
                )
                peaks_info <- peaks_info %>% filter(chr == chr_num)
            }

            # Get the highest peak (if any)
            if (nrow(peaks_info) > 0) {
                peaks_info <- peaks_info %>% slice(1)

                # Find the peak point in the filtered data
                peak_point <- plot_data %>%
                    filter(markers == peaks_info$marker)

                if (nrow(peak_point) > 0) {
                    # Add the red diamond at the peak
                    if (input$selected_chr == "All") {
                        p <- p + geom_point(
                            data = peak_point,
                            aes(x = BPcum, y = LOD),
                            color = "red",
                            size = 3,
                            shape = 18
                        )
                    } else {
                        p <- p + geom_point(
                            data = peak_point,
                            aes(x = position, y = LOD),
                            color = "red",
                            size = 3,
                            shape = 18
                        )
                    }
                }
            }

            # Create plot title and subtitle
            trait_title <- input$which_trait
            chr_info <- if (input$selected_chr != "All") paste0(" (Chromosome ", input$selected_chr, ")") else ""

            subtitle <- if (nrow(peaks_info) > 0) {
                chr_label <- if (peaks_info$chr %in% c(20, 21, 22)) {
                    c("X", "Y", "M")[peaks_info$chr - 19]
                } else {
                    peaks_info$chr
                }
                paste0(
                    "Peak Marker: ", peaks_info$marker,
                    " (Chr", chr_label, ":", round(peaks_info$pos, 2), " Mb) | ",
                    "LOD: ", round(peaks_info$lod, 2)
                )
            } else {
                if (input$selected_chr != "All") {
                    "No significant peaks in this chromosome"
                } else {
                    "No significant peaks"
                }
            }

            # Add title and subtitle to the plot
            p <- p + ggtitle(
                label = paste0(trait_title, chr_info),
                subtitle = subtitle
            )

            # Save the plot as PDF
            ggsave(file, p, width = input$plot_width / 72, height = input$plot_height / 72, device = cairo_pdf)
        }
    )
}

# Launch the Shiny application
shinyApp(ui = ui, server = server)
