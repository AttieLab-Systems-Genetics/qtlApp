# Suppress data.table NSE warnings
utils::globalVariables(c(
    "phenotype", "correlation_value", "p_value", "num_mice",
    "abs_correlation", "abs_corr_temp", ".N", ".SD", "target", "source",
    ".", ":="
))

#' Correlation Module Inputs
#'
#' @param id Module ID.
#' @return UI elements for correlation controls
#' @importFrom shiny NS selectInput tags div
#' @export
correlationInput <- function(id) {
    ns <- shiny::NS(id)
    shiny::div(
        shiny::div(
            style = "display: flex; align-items: center; gap: 12px;",
            shiny::tags$label(
                `for` = ns("correlation_dataset_selector"),
                style = "margin: 0; font-weight: 600;",
                "Correlate against dataset"
            ),
            shiny::selectInput(
                inputId = ns("correlation_dataset_selector"),
                label = NULL,
                choices = character(0),
                selected = NULL,
                width = "420px"
            ),
            shiny::textInput(
                inputId = ns("correlation_search"),
                label = NULL,
                placeholder = "Search trait...",
                width = "280px"
            ),
            shiny::div(
                style = "display: flex; align-items: center; gap: 6px;",
                shiny::checkboxInput(
                    inputId = ns("use_adjusted"),
                    label = "Covariate-Adjusted",
                    value = FALSE
                )
            )
        )
    )
}

#' Correlation Module UI Output
#'
#' @param id Module ID.
#' @return UI container for the correlation table
#' @importFrom bslib card card_header
#' @importFrom shinycssloaders withSpinner
#' @importFrom DT DTOutput
#' @export
correlationUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::div(
        bslib::card(
            bslib::card_header("Correlations"),
            shinycssloaders::withSpinner(
                DT::DTOutput(ns("correlation_table")),
                type = 4, color = "#3498db"
            )
        )
    )
}

#' Correlation Module Server
#'
#' Displays a table of correlations for the currently selected trait against
#' a chosen dataset. The table includes columns: trait, correlation_value, p_value, num_mice.
#'
#' Correlation CSV files are expected in the directory:
#'   data/correlations (relative to app working directory)
#' and named like: "liver_genes_vs_clinical_traits_corr.csv" for unadjusted correlations
#' or "liver_genes_vs_clinical_traits_genlitsexbydiet_adj_corr.csv" for covariate-adjusted.
#' Some correlation pairs (e.g., genes-vs-other datasets) only exist in adjusted form.
#'
#' @param id Module ID.
#' @param import_reactives Reactive that returns a list containing file_directory and annotation_list.
#' @param main_par Reactive that returns a list containing selected_dataset and which_trait.
#' @return Invisibly returns NULL. UI is rendered via outputs.
#' @importFrom data.table fread as.data.table melt setDT setnames setorder
#' @importFrom DT DTOutput renderDT datatable formatRound formatStyle
#' @importFrom shiny moduleServer NS reactive req validate need observeEvent updateSelectInput updateSelectizeInput debounce showNotification
#' @importFrom htmlwidgets JS
#' @export
correlationServer <- function(id, import_reactives, main_par) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Use the local data directory for correlations
        # Check multiple possible locations (Docker mount first!)
        possible_paths <- c(
            "/data/correlations", # Docker mount path (FIRST - this is where it's mounted!)
            "/home/khwillis@ad.wisc.edu/qtlApp/data/correlations", # Development local path
            file.path(getwd(), "data/correlations"), # Relative to working directory
            "/mnt/rdrive/mkeller3/General/main_directory/correlations" # Rdrive fallback
        )

        correlations_dir <- NULL
        for (path in possible_paths) {
            # Check if directory exists AND has files
            if (dir.exists(path)) {
                test_files <- list.files(path, pattern = "_corr\\.csv$")
                if (length(test_files) > 0) {
                    correlations_dir <- path
                    message("correlationServer: ✓ Found correlations directory with ", length(test_files), " files at: ", path)
                    break
                } else {
                    message("correlationServer: Directory exists but is empty: ", path)
                }
            }
        }

        if (is.null(correlations_dir)) {
            correlations_dir <- "/data/correlations" # Fallback to Docker mount path
            warning("correlationServer: No correlations directory with files found. Using fallback: ", correlations_dir)
        }

        # Helper function to convert file path to adjusted version
        # Handles both same-ordering and reversed-ordering adjusted files
        get_adjusted_file_path <- function(file_path, use_adjusted = FALSE) {
            if (!use_adjusted) {
                return(file_path)
            }

            # First try same-ordering adjusted file
            same_order <- sub("_corr\\.csv$", "_genlitsexbydiet_adj_corr.csv", file_path)
            if (file.exists(same_order)) {
                return(same_order)
            }

            # Try reversed-ordering adjusted file
            bn <- basename(file_path)
            m <- regexec("^([a-z0-9_]+)_vs_([a-z0-9_]+)_corr\\.csv$", tolower(bn))
            r <- regmatches(tolower(bn), m)[[1]]
            if (length(r) == 3) {
                source <- r[2]
                target <- r[3]
                reversed_adj <- file.path(
                    dirname(file_path),
                    paste0(target, "_vs_", source, "_genlitsexbydiet_adj_corr.csv")
                )
                if (file.exists(reversed_adj)) {
                    return(reversed_adj)
                }
            }

            # Fallback to same-ordering pattern even if it doesn't exist
            return(same_order)
        }

        # Map internal trait type to correlation file tokens and labels
        trait_type_to_token <- function(trait_type_char) {
            if (is.null(trait_type_char) || !nzchar(trait_type_char)) {
                return(NULL)
            }
            if (grepl("^gene", trait_type_char, ignore.case = TRUE)) {
                return("liver_genes")
            }
            if (grepl("isoform", trait_type_char, ignore.case = TRUE)) {
                return("liver_isoforms")
            }
            if (grepl("lipid", trait_type_char, ignore.case = TRUE)) {
                return("liver_lipids")
            }
            if (grepl("plasma.*metabol", trait_type_char, ignore.case = TRUE)) {
                return("plasma_metabolites")
            }
            if (grepl("clinical", trait_type_char, ignore.case = TRUE)) {
                return("clinical_traits")
            }
            if (grepl("splice.*junc", trait_type_char, ignore.case = TRUE)) {
                return("liver_splice_juncs")
            }
            return(NULL)
        }

        token_to_label <- function(token) {
            switch(token,
                liver_genes = "Liver Genes",
                liver_isoforms = "Liver Isoforms",
                liver_lipids = "Liver Lipids",
                plasma_metabolites = "Plasma Metabolites",
                clinical_traits = "Clinical Traits",
                liver_splice_juncs = "Liver Splice Junctions",
                token
            )
        }

        # Utility: List available correlation files and parse source/target tokens
        list_available_pairs <- shiny::reactive({
            # Expected/fallback files to show even if directory scan fails
            expected_files <- c(
                file.path(correlations_dir, "clinical_traits_vs_clinical_traits_corr.csv"),
                file.path(correlations_dir, "liver_genes_vs_clinical_traits_corr.csv"),
                file.path(correlations_dir, "liver_genes_vs_liver_genes_corr.csv"),
                file.path(correlations_dir, "liver_genes_vs_liver_lipids_corr.csv"),
                file.path(correlations_dir, "liver_genes_vs_plasma_metabolites_corr.csv"),
                file.path(correlations_dir, "liver_lipids_vs_liver_lipids_corr.csv"),
                file.path(correlations_dir, "plasma_metabolites_vs_plasma_metabolites_corr.csv"),
                file.path(correlations_dir, "clinical_traits_vs_liver_lipids_corr.csv"),
                file.path(correlations_dir, "clinical_traits_vs_plasma_metabolites_corr.csv"),
                file.path(correlations_dir, "liver_lipids_vs_plasma_metabolites_corr.csv"),
                # Add adjusted-only files that don't have base versions
                file.path(correlations_dir, "clinical_traits_vs_liver_genes_genlitsexbydiet_adj_corr.csv"),
                file.path(correlations_dir, "liver_lipids_vs_liver_genes_genlitsexbydiet_adj_corr.csv"),
                file.path(correlations_dir, "plasma_metabolites_vs_liver_genes_genlitsexbydiet_adj_corr.csv")
            )

            files <- character(0)
            message("correlationServer: Checking directory: ", correlations_dir)
            message("correlationServer: Directory exists: ", dir.exists(correlations_dir))
            if (dir.exists(correlations_dir)) {
                all_files <- list.files(correlations_dir, pattern = "_corr\\.csv$", full.names = TRUE)
                message("correlationServer: Found ", length(all_files), " total _corr.csv files")
                # Include base (unadjusted) files
                base_files <- all_files[!grepl("_adj_corr\\.csv$", all_files)]
                message("correlationServer: Found ", length(base_files), " base (unadjusted) files")
                # Include ALL adjusted files (both same-ordering and reversed-ordering)
                adj_files <- all_files[grepl("_adj_corr\\.csv$", all_files)]
                message("correlationServer: Found ", length(adj_files), " adjusted files")
                # Include all adjusted files (we'll determine adjusted-only status later)
                files <- c(base_files, adj_files)
                message("correlationServer: Total files to process: ", length(files))
            } else {
                warning("correlationServer: correlations_dir not found: ", correlations_dir)
            }

            # Include expected files regardless of current presence; parsing will still work
            files <- unique(c(files, expected_files))
            if (length(files) == 0) {
                return(data.frame(file = character(0), source = character(0), target = character(0), is_adj_only = logical(0)))
            }

            # Parse names like foo_vs_bar_corr.csv (including adjusted versions)
            parse_one <- function(fp) {
                bn <- basename(fp)
                # Normalize filename by removing adjustment suffix for parsing
                # Convert: foo_vs_bar_genlitsexbydiet_adj_corr.csv -> foo_vs_bar_corr.csv
                normalized_bn <- sub("_genlitsexbydiet_adj_corr\\.csv$", "_corr.csv", tolower(bn))

                # Parse the normalized filename
                m <- regexec("^([a-z0-9_]+)_vs_([a-z0-9_]+)_corr\\.csv$", normalized_bn)
                r <- regmatches(normalized_bn, m)[[1]]
                if (length(r) == 3) {
                    return(data.frame(file = fp, source = r[2], target = r[3], stringsAsFactors = FALSE))
                }
                return(NULL)
            }
            parsed <- lapply(files, parse_one)
            parsed <- parsed[!vapply(parsed, is.null, logical(1))]
            if (length(parsed) == 0) {
                warning("correlationServer: No correlation files could be parsed")
                return(data.frame(file = character(0), source = character(0), target = character(0), is_adj_only = logical(0)))
            }

            # Combine all parsed results
            result_df <- do.call(rbind, parsed)
            result_dt <- data.table::as.data.table(result_df)

            # Determine which files are adjusted-only by checking for corresponding base files
            # Need to check BOTH same-ordering and reversed-ordering base versions
            result_dt$is_adj_only <- vapply(result_dt$file, function(fp) {
                # If it's already a base file, it's not adjusted-only
                if (!grepl("_adj_corr\\.csv$", fp)) {
                    return(FALSE)
                }
                # If it's an adjusted file, check if base version exists (same ordering)
                base_version <- sub("_genlitsexbydiet_adj_corr\\.csv$", "_corr.csv", fp)
                if (file.exists(base_version)) {
                    return(FALSE) # Base version exists with same ordering
                }
                # Also check reversed ordering (e.g., clinical_traits_vs_liver_genes -> liver_genes_vs_clinical_traits)
                bn <- basename(fp)
                # Parse the adjusted filename
                m <- regexec("^([a-z0-9_]+)_vs_([a-z0-9_]+)_genlitsexbydiet_adj_corr\\.csv$", tolower(bn))
                r <- regmatches(tolower(bn), m)[[1]]
                if (length(r) == 3) {
                    source <- r[2]
                    target <- r[3]
                    # Check reversed ordering
                    reversed_base <- file.path(dirname(fp), paste0(target, "_vs_", source, "_corr.csv"))
                    if (file.exists(reversed_base)) {
                        return(FALSE) # Base version exists with reversed ordering
                    }
                }
                return(TRUE) # No base version found in either ordering
            }, logical(1))

            message("correlationServer: Found ", nrow(result_dt), " correlation file pairs")
            result_dt
        })

        # Determine current trait and its tokenized type
        current_trait_reactive <- shiny::reactive({
            mp <- main_par()
            shiny::req(mp, mp$which_trait)
            trait_val <- tryCatch(mp$which_trait(), error = function(e) NULL)
            if (is.null(trait_val) || !nzchar(trait_val)) {
                return(NULL)
            }
            trait_val
        })

        current_type_token <- shiny::reactive({
            mp <- main_par()
            shiny::req(mp, mp$selected_dataset)
            selected_group <- tryCatch(mp$selected_dataset(), error = function(e) NULL)
            shiny::req(selected_group)
            trait_type <- get_trait_type(import_reactives(), selected_group)
            trait_type_to_token(trait_type)
        })

        # Populate dataset choices based on available pairs that involve current trait type
        shiny::observeEvent(list(list_available_pairs(), current_type_token()),
            {
                shiny::req(current_type_token())
                pairs_dt <- list_available_pairs()
                if (nrow(pairs_dt) == 0) {
                    shiny::updateSelectInput(session, "correlation_dataset_selector",
                        choices = character(0), selected = character(0)
                    )
                    return()
                }
                token <- current_type_token()
                message("correlationServer: Current token = ", token)
                message("correlationServer: Available pairs with sources: ", paste(unique(pairs_dt$source), collapse = ", "))
                message("correlationServer: Available pairs with targets: ", paste(unique(pairs_dt$target), collapse = ", "))

                # Files where current type is on either side (using standard R subsetting)
                subset_dt <- pairs_dt[pairs_dt$source == token | pairs_dt$target == token, ]
                if (nrow(subset_dt) == 0) {
                    message("correlationServer: No matching pairs found for token: ", token)
                    shiny::updateSelectInput(session, "correlation_dataset_selector",
                        choices = character(0), selected = character(0)
                    )
                    return()
                }
                message("correlationServer: Found ", nrow(subset_dt), " matching pairs")
                # Build choices mapping to only show the target dataset label (the "other side")
                other_tokens <- ifelse(subset_dt$source == token, subset_dt$target, subset_dt$source)
                labels <- vapply(other_tokens, token_to_label, character(1))
                # Prefer files where the current token is on the source side; then deduplicate by label
                order_pref <- ifelse(subset_dt$source == token, 0L, 1L)
                ord <- order(order_pref, labels)
                subset_dt <- subset_dt[ord]
                labels <- labels[ord]
                keep_idx <- !duplicated(labels)
                choices <- stats::setNames(subset_dt$file[keep_idx], labels[keep_idx])
                shiny::updateSelectInput(session, "correlation_dataset_selector",
                    choices = choices,
                    selected = if (length(choices) > 0) choices[[1]] else NULL
                )
            },
            ignoreInit = FALSE
        )

        # Maximum number of correlations to return (for memory/performance optimization)
        MAX_CORRELATIONS <- 500

        # Resolve the proper column or row key for the current trait in a given file
        # For gene traits, correlation files typically use columns like "liver_<gene_id>"
        resolve_trait_keys <- function(trait_string, token_side, corr_file, import) {
            # token_side is one of c("liver_genes", "liver_isoforms", "liver_lipids", "plasma_metabolites", "clinical_traits")
            # Return list(mode = "column"|"row", key = <name>)
            header_only <- tryCatch(data.table::fread(corr_file, nrows = 0), error = function(e) NULL)
            if (is.null(header_only)) {
                return(NULL)
            }
            headers <- names(header_only)

            # Column candidates
            column_candidates <- character(0)

            if (token_side == "liver_genes") {
                # Try exact match
                if (trait_string %in% headers) column_candidates <- c(column_candidates, trait_string)
                # Try adding liver_ prefix if missing
                if (!grepl("^liver_", trait_string) && paste0("liver_", trait_string) %in% headers) {
                    column_candidates <- c(column_candidates, paste0("liver_", trait_string))
                }
                # Try mapping gene symbol -> gene.id via annotations available in import
                if (length(column_candidates) == 0) {
                    ann <- import$annotation_list
                    if (!is.null(ann) && !is.null(ann$genes)) {
                        genes_dt <- ann$genes
                        id_col <- if ("gene.id" %in% colnames(genes_dt)) "gene.id" else if ("gene_id" %in% colnames(genes_dt)) "gene_id" else NULL
                        sym_col <- if ("symbol" %in% colnames(genes_dt)) "symbol" else NULL
                        if (!is.null(id_col) && !is.null(sym_col)) {
                            match_row <- genes_dt[genes_dt[[sym_col]] == trait_string, , drop = FALSE]
                            if (nrow(match_row) > 0) {
                                gene_id <- match_row[[id_col]][1]
                                candidate <- paste0("liver_", gene_id)
                                if (candidate %in% headers) column_candidates <- c(column_candidates, candidate)
                            }
                        }
                    }
                }
            } else {
                # For non-gene types, first column should be phenotype; we'll usually match as a row
                if (trait_string %in% headers) column_candidates <- c(column_candidates, trait_string)
            }

            # If a direct column match exists, prefer column mode (efficient select)
            if (length(column_candidates) > 0) {
                return(list(mode = "column", key = column_candidates[1]))
            }

            # Otherwise, attempt row mode via Phenotype/phenotype column
            pheno_col <- if ("phenotype" %in% tolower(headers)) headers[which(tolower(headers) == "phenotype")[1]] else NULL
            if (!is.null(pheno_col)) {
                # For liver_genes in row mode, try gene symbol -> gene ID mapping
                row_key <- trait_string
                if (token_side == "liver_genes") {
                    ann <- import$annotation_list
                    if (!is.null(ann) && !is.null(ann$genes)) {
                        genes_dt <- ann$genes
                        id_col <- if ("gene.id" %in% colnames(genes_dt)) "gene.id" else if ("gene_id" %in% colnames(genes_dt)) "gene_id" else NULL
                        sym_col <- if ("symbol" %in% colnames(genes_dt)) "symbol" else NULL
                        if (!is.null(id_col) && !is.null(sym_col)) {
                            match_row <- genes_dt[genes_dt[[sym_col]] == trait_string, , drop = FALSE]
                            if (nrow(match_row) > 0) {
                                gene_id <- match_row[[id_col]][1]
                                # Try with liver_ prefix (adjusted files use this format)
                                row_key <- paste0("liver_", gene_id)
                            }
                        }
                    }
                }
                return(list(mode = "row", key = row_key, phenotype_col = pheno_col))
            }

            return(NULL)
        }

        # Build the correlation table for display
        correlation_table <- shiny::reactive({
            shiny::req(input$correlation_dataset_selector)
            trait <- current_trait_reactive()
            shiny::req(trait)

            # Get the appropriate file path based on adjusted toggle
            use_adjusted <- isTRUE(input$use_adjusted)

            # Check if the selected file is adjusted-only
            pairs_dt <- list_available_pairs()
            shiny::req(nrow(pairs_dt) > 0)
            base_file_path <- input$correlation_dataset_selector
            row <- pairs_dt[file == base_file_path]
            shiny::req(nrow(row) == 1)

            is_adj_only <- isTRUE(row$is_adj_only[1])

            # For adjusted-only files, always use the file as-is
            # For files with both versions, apply the toggle
            if (is_adj_only) {
                file_path <- base_file_path
            } else {
                file_path <- get_adjusted_file_path(base_file_path, use_adjusted)
                # If file doesn't exist, try adjusted version
                if (!file.exists(file_path) && !use_adjusted) {
                    adjusted_path <- get_adjusted_file_path(base_file_path, TRUE)
                    if (file.exists(adjusted_path)) {
                        file_path <- adjusted_path
                    }
                }
            }

            if (!file.exists(file_path)) {
                msg <- if (use_adjusted) {
                    paste("Covariate-adjusted correlation file not found:", basename(file_path))
                } else {
                    paste("Correlation file not found:", basename(file_path))
                }
                shiny::showNotification(msg, type = "error", duration = 4)
                return(data.frame(trait = character(0), correlation_value = numeric(0), p_value = numeric(0), num_mice = numeric(0)))
            }

            message("correlationServer: Loading correlation file: ", basename(file_path))
            message("correlationServer: Is adjusted: ", use_adjusted, ", Is adjusted-only: ", is_adj_only)

            # Determine which side (source/target) corresponds to the current trait type
            token <- current_type_token()

            side_token <- if (row$source[1] == token) row$source[1] else if (row$target[1] == token) row$target[1] else token

            message("correlationServer: Resolving trait keys for trait='", trait, "', side_token='", side_token, "'")
            key_info <- resolve_trait_keys(trait, side_token, file_path, import_reactives())
            message("correlationServer: key_info result: ", if (is.null(key_info)) "NULL" else paste("mode=", key_info$mode, ", key=", key_info$key))
            if (is.null(key_info)) {
                shiny::showNotification(paste("Trait not found in correlation file:", trait), type = "warning", duration = 4)
                return(data.frame(trait = character(0), correlation_value = numeric(0), p_value = numeric(0), num_mice = numeric(0)))
            }

            if (key_info$mode == "column") {
                # Check file size and show loading notification for large files
                file_size_mb <- file.info(file_path)$size / (1024^2)
                if (file_size_mb > 100) {
                    shiny::showNotification(
                        paste0("Loading correlation file (", round(file_size_mb, 0), " MB). This may take 30-60 seconds..."),
                        type = "message", duration = NULL, id = "loading_corr_col"
                    )
                }

                # Read only phenotype and the specific column for efficiency
                select_cols <- c("phenotype", "Phenotype", key_info$key)
                select_cols <- select_cols[select_cols %in% names(data.table::fread(file_path, nrows = 0))]
                dt <- data.table::fread(file_path, select = select_cols)

                if (file_size_mb > 100) {
                    shiny::removeNotification(id = "loading_corr_col")
                }
                # Normalize phenotype column name
                if ("Phenotype" %in% names(dt)) data.table::setnames(dt, "Phenotype", "phenotype")
                if (!"phenotype" %in% names(dt)) dt[, phenotype := seq_len(.N)]
                # Build output table
                out <- dt[, .(trait = phenotype, correlation_value = .SD[[1]]), .SDcols = key_info$key]

                # For large datasets (genes-vs-genes), filter to top N by absolute correlation
                total_rows <- nrow(out)
                truncated <- FALSE
                if (total_rows > MAX_CORRELATIONS) {
                    out[, abs_corr_temp := abs(correlation_value)]
                    data.table::setorder(out, -abs_corr_temp)
                    out <- out[1:MAX_CORRELATIONS]
                    out[, abs_corr_temp := NULL]
                    truncated <- TRUE
                }

                # Add p-values from companion file if available (only for filtered traits)
                pval_path <- sub("_genlitsexbydiet_adj_corr\\.csv$", "_genlitsexbydiet_adj_pval.csv", file_path)
                pval_path <- sub("_corr\\.csv$", "_pval.csv", pval_path)
                if (file.exists(pval_path)) {
                    pval_size_mb <- file.info(pval_path)$size / (1024^2)
                    if (pval_size_mb > 100) {
                        shiny::showNotification("Loading p-values...", type = "message", duration = NULL, id = "loading_pval_col")
                    }
                    psel <- c("phenotype", "Phenotype", key_info$key)
                    psel <- psel[psel %in% names(data.table::fread(pval_path, nrows = 0))]
                    pdt <- data.table::fread(pval_path, select = psel)
                    if (pval_size_mb > 100) {
                        shiny::removeNotification(id = "loading_pval_col")
                    }
                    if ("Phenotype" %in% names(pdt)) data.table::setnames(pdt, "Phenotype", "phenotype")
                    if (!"phenotype" %in% names(pdt)) pdt[, phenotype := seq_len(.N)]
                    pmerge <- pdt[, .(trait = phenotype, p_value = .SD[[1]]), .SDcols = key_info$key]
                    out <- merge(out, pmerge, by = "trait", all.x = TRUE, sort = FALSE)
                } else {
                    out[, p_value := as.numeric(NA)]
                }
                # Add num_mice from companion file if available (only for filtered traits)
                nmouse_path <- sub("_genlitsexbydiet_adj_corr\\.csv$", "_genlitsexbydiet_adj_num_mice.csv", file_path)
                nmouse_path <- sub("_corr\\.csv$", "_num_mice.csv", nmouse_path)
                if (file.exists(nmouse_path)) {
                    nmouse_size_mb <- file.info(nmouse_path)$size / (1024^2)
                    if (nmouse_size_mb > 100) {
                        shiny::showNotification("Loading sample sizes...", type = "message", duration = NULL, id = "loading_nmice_col")
                    }
                    nsel <- c("phenotype", "Phenotype", key_info$key)
                    nsel <- nsel[nsel %in% names(data.table::fread(nmouse_path, nrows = 0))]
                    ndt <- data.table::fread(nmouse_path, select = nsel)
                    if (nmouse_size_mb > 100) {
                        shiny::removeNotification(id = "loading_nmice_col")
                    }
                    if ("Phenotype" %in% names(ndt)) data.table::setnames(ndt, "Phenotype", "phenotype")
                    if (!"phenotype" %in% names(ndt)) ndt[, phenotype := seq_len(.N)]
                    nmerge <- ndt[, .(trait = phenotype, num_mice = .SD[[1]]), .SDcols = key_info$key]
                    out <- merge(out, nmerge, by = "trait", all.x = TRUE, sort = FALSE)
                } else {
                    out[, num_mice := as.numeric(NA)]
                }

                # Notify user if results were truncated
                if (truncated) {
                    shiny::showNotification(
                        paste0("Showing top ", MAX_CORRELATIONS, " of ", total_rows, " correlations (filtered by absolute correlation)"),
                        type = "message",
                        duration = 5
                    )
                }
            } else {
                # Row mode: read entire file, find the row, and transpose to long format
                # For large files (genes-vs-genes), show progress notification
                file_size_mb <- file.info(file_path)$size / (1024^2)
                if (file_size_mb > 100) {
                    shiny::showNotification(
                        paste0("Loading large correlation file (", round(file_size_mb, 0), " MB). This may take 30-60 seconds..."),
                        type = "message", duration = NULL, id = "loading_corr"
                    )
                }

                dt_full <- data.table::fread(file_path, showProgress = FALSE)

                if (file_size_mb > 100) {
                    shiny::removeNotification(id = "loading_corr")
                }

                # Normalize phenotype column name
                if ("Phenotype" %in% names(dt_full)) data.table::setnames(dt_full, "Phenotype", "phenotype")
                if (!"phenotype" %in% names(dt_full)) {
                    shiny::showNotification("Correlation file missing 'phenotype' column.", type = "error", duration = NULL)
                    return(data.frame(trait = character(0), correlation_value = numeric(0), p_value = numeric(0), num_mice = numeric(0)))
                }
                message("correlationServer: Row mode - looking for phenotype='", key_info$key, "' in file with ", nrow(dt_full), " rows")
                message("correlationServer: Row mode - first 5 phenotypes: ", paste(head(dt_full$phenotype, 5), collapse = ", "))
                row_dt <- dt_full[phenotype == key_info$key]
                message("correlationServer: Row mode - found ", nrow(row_dt), " matching rows")
                if (nrow(row_dt) == 0) {
                    shiny::showNotification(paste("Trait not found in 'phenotype' column:", trait), type = "warning", duration = 4)
                    return(data.frame(trait = character(0), correlation_value = numeric(0), p_value = numeric(0), num_mice = numeric(0)))
                }
                # Convert the single row (excluding phenotype) to long format
                numeric_cols <- setdiff(names(row_dt), "phenotype")
                long_dt <- data.table::melt(row_dt,
                    id.vars = "phenotype", measure.vars = numeric_cols,
                    variable.name = "trait", value.name = "correlation_value"
                )
                out <- long_dt[, .(trait, correlation_value)]

                # For large datasets (genes-vs-genes), filter to top N by absolute correlation
                total_rows <- nrow(out)
                truncated <- FALSE
                if (total_rows > MAX_CORRELATIONS) {
                    out[, abs_corr_temp := abs(correlation_value)]
                    data.table::setorder(out, -abs_corr_temp)
                    out <- out[1:MAX_CORRELATIONS]
                    out[, abs_corr_temp := NULL]
                    truncated <- TRUE
                }

                # Get the trait names we kept for efficient companion file reading
                kept_traits <- out$trait

                # Add p-values from companion pval file (same row, filtered columns)
                pval_path <- sub("_genlitsexbydiet_adj_corr\\.csv$", "_genlitsexbydiet_adj_pval.csv", file_path)
                pval_path <- sub("_corr\\.csv$", "_pval.csv", pval_path)
                if (file.exists(pval_path)) {
                    pval_size_mb <- file.info(pval_path)$size / (1024^2)
                    if (pval_size_mb > 100) {
                        shiny::showNotification("Loading p-values...", type = "message", duration = NULL, id = "loading_pval_row")
                    }
                    # Only read phenotype column + kept trait columns for efficiency
                    select_cols <- c("phenotype", as.character(kept_traits))
                    pdt_full <- data.table::fread(pval_path, select = select_cols)
                    if (pval_size_mb > 100) {
                        shiny::removeNotification(id = "loading_pval_row")
                    }
                    if ("Phenotype" %in% names(pdt_full)) data.table::setnames(pdt_full, "Phenotype", "phenotype")
                    if ("phenotype" %in% names(pdt_full)) {
                        prow <- pdt_full[phenotype == key_info$key]
                        if (nrow(prow) == 1) {
                            pnumeric <- setdiff(names(prow), "phenotype")
                            plong <- data.table::melt(prow,
                                id.vars = "phenotype", measure.vars = pnumeric,
                                variable.name = "trait", value.name = "p_value"
                            )
                            pmerge <- plong[, .(trait, p_value)]
                            out <- merge(out, pmerge, by = "trait", all.x = TRUE, sort = FALSE)
                        } else {
                            out[, p_value := as.numeric(NA)]
                        }
                    } else {
                        out[, p_value := as.numeric(NA)]
                    }
                } else {
                    out[, p_value := as.numeric(NA)]
                }

                # Add num_mice from companion num_mice file (same row, filtered columns)
                nmouse_path <- sub("_genlitsexbydiet_adj_corr\\.csv$", "_genlitsexbydiet_adj_num_mice.csv", file_path)
                nmouse_path <- sub("_corr\\.csv$", "_num_mice.csv", nmouse_path)
                if (file.exists(nmouse_path)) {
                    nmouse_size_mb <- file.info(nmouse_path)$size / (1024^2)
                    if (nmouse_size_mb > 100) {
                        shiny::showNotification("Loading sample sizes...", type = "message", duration = NULL, id = "loading_nmice_row")
                    }
                    # Only read phenotype column + kept trait columns for efficiency
                    select_cols <- c("phenotype", as.character(kept_traits))
                    nfull <- data.table::fread(nmouse_path, select = select_cols)
                    if (nmouse_size_mb > 100) {
                        shiny::removeNotification(id = "loading_nmice_row")
                    }
                    if ("Phenotype" %in% names(nfull)) data.table::setnames(nfull, "Phenotype", "phenotype")
                    if ("phenotype" %in% names(nfull)) {
                        nrow_dt <- nfull[phenotype == key_info$key]
                        if (nrow(nrow_dt) == 1) {
                            ncols <- setdiff(names(nrow_dt), "phenotype")
                            nlong <- data.table::melt(nrow_dt,
                                id.vars = "phenotype", measure.vars = ncols,
                                variable.name = "trait", value.name = "num_mice"
                            )
                            nmerge <- nlong[, .(trait, num_mice)]
                            out <- merge(out, nmerge, by = "trait", all.x = TRUE, sort = FALSE)
                        } else {
                            out[, num_mice := as.numeric(NA)]
                        }
                    } else {
                        out[, num_mice := as.numeric(NA)]
                    }
                } else {
                    out[, num_mice := as.numeric(NA)]
                }

                # Notify user if results were truncated
                if (truncated) {
                    shiny::showNotification(
                        paste0("Showing top ", MAX_CORRELATIONS, " of ", total_rows, " correlations (filtered by absolute correlation)"),
                        type = "message",
                        duration = 5
                    )
                }

                # Gene mapping happens after both column and row mode complete
            }

            # Map gene ids -> symbols for ALL genes-vs-genes correlations (applies to both column and row mode)
            has_liver_prefix <- any(grepl("^liver_ENSMUSG", out$trait, perl = TRUE))

            if (has_liver_prefix) {
                out$trait <- as.character(out$trait)
                gene_ids <- sub("^liver_", "", out$trait)

                # Use the same annotation source as resolve_trait_keys
                imp <- import_reactives()
                ann <- NULL
                if (!is.null(imp) && !is.null(imp$annotation_list)) {
                    ann <- imp$annotation_list
                }

                if (!is.null(ann) && !is.null(ann$genes)) {
                    genes_dt <- ann$genes
                    id_col <- if ("gene.id" %in% colnames(genes_dt)) "gene.id" else if ("gene_id" %in% colnames(genes_dt)) "gene_id" else NULL
                    sym_col <- if ("symbol" %in% colnames(genes_dt)) "symbol" else NULL

                    if (!is.null(id_col) && !is.null(sym_col)) {
                        # Create lookup: gene_id -> symbol
                        id_to_symbol <- stats::setNames(genes_dt[[sym_col]], genes_dt[[id_col]])

                        # Map gene IDs to symbols
                        mapped <- id_to_symbol[gene_ids]
                        display <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, gene_ids)
                        out$trait <- display
                    } else {
                        out$trait <- gene_ids
                    }
                } else {
                    out$trait <- gene_ids
                }
            }

            # Remove any rows with NA correlation values
            if ("correlation_value" %in% names(out)) {
                out <- out[!is.na(correlation_value)]
            }

            message("correlationServer: Correlation table prepared with ", nrow(out), " rows")
            message("correlationServer: Columns: ", paste(names(out), collapse = ", "))
            if (nrow(out) > 0) {
                message("correlationServer: First 3 traits: ", paste(head(out$trait, 3), collapse = ", "))
            }

            # Prepare display formatting and default sort by absolute correlation
            if ("p_value" %in% names(out)) {
                # For self-correlation (r ≈ 1.0), p-value may be missing (NA) from source data
                # This is expected because p-value calculation for perfect correlation is undefined
                out[, p_value := ifelse(is.na(p_value), "—", format(p_value, digits = 3, scientific = TRUE))]
            }
            out[, abs_correlation := abs(correlation_value)]
            out <- out[order(-abs_correlation, -correlation_value)]
            out <- data.table::as.data.table(out)

            message("correlationServer: Returning table with ", nrow(out), " rows")
            out
        }) |> shiny::debounce(150)

        # Debounced search query
        search_query <- shiny::reactive({
            val <- input$correlation_search
            if (is.null(val)) "" else trimws(as.character(val))
        }) |> shiny::debounce(150)

        # Filtered table by trait search; preserves ordering by absolute correlation
        filtered_table <- shiny::reactive({
            tbl <- correlation_table()
            message("correlationServer: filtered_table reactive - received ", nrow(tbl), " rows from correlation_table()")
            q <- search_query()
            if (!is.null(tbl) && nrow(tbl) > 0 && nzchar(q)) {
                keep <- tryCatch(grepl(q, tbl$trait, ignore.case = TRUE, perl = TRUE),
                    error = function(e) grepl(q, tbl$trait, ignore.case = TRUE, fixed = TRUE)
                )
                tbl <- tbl[keep]
                message("correlationServer: After search filter (query='", q, "'): ", nrow(tbl), " rows")
            }
            message("correlationServer: filtered_table returning ", nrow(tbl), " rows")
            tbl
        })

        output$correlation_table <- DT::renderDT({
            tbl <- filtered_table()
            if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0) {
                # Render an empty, stable table with correct columns but no rows
                empty_tbl <- data.table::as.data.table(data.frame(
                    abs_correlation = numeric(0),
                    Trait = character(0),
                    Correlation = numeric(0),
                    `Pvalue` = character(0),
                    `# Mice` = numeric(0),
                    check.names = FALSE
                ))
                return(DT::datatable(
                    empty_tbl,
                    rownames = FALSE,
                    width = "100%",
                    class = "compact",
                    escape = FALSE,
                    options = list(
                        pageLength = 25,
                        autoWidth = FALSE,
                        responsive = TRUE,
                        scrollX = FALSE,
                        order = list(list(0, "desc")),
                        columnDefs = list(
                            list(visible = FALSE, targets = 0),
                            list(width = "42%", targets = 1),
                            list(width = "16%", targets = 2),
                            list(width = "18%", targets = 3),
                            list(width = "12%", targets = 4)
                        ),
                        dom = "tip"
                    ),
                    callback = htmlwidgets::JS(
                        "table.settings()[0].aoColumns.forEach(function(col){",
                        "  if(col && col.nTh){col.nTh.style.whiteSpace='nowrap';col.nTh.style.wordBreak='normal';col.nTh.style.fontSize='90%';col.nTh.style.padding='6px 8px';}",
                        "});",
                        "table.on('draw.dt', function(){",
                        "  table.columns().every(function(){",
                        "    $(this.nodes()).css({'white-space':'normal','word-break':'break-word'});",
                        "  });",
                        "});"
                    )
                ))
            }
            tbl <- data.table::as.data.table(tbl)
            # Include a hidden abs_correlation column for default ordering by magnitude
            display_tbl <- tbl[, .(abs_correlation, trait, correlation_value, p_value, num_mice)]
            data.table::setnames(display_tbl,
                old = c("trait", "correlation_value", "p_value", "num_mice"),
                new = c("Trait", "Correlation", "Pvalue", "# Mice"),
                skip_absent = TRUE
            )
            dt <- DT::datatable(
                display_tbl,
                rownames = FALSE,
                width = "100%",
                class = "compact",
                escape = FALSE,
                options = list(
                    pageLength = 25,
                    autoWidth = FALSE,
                    responsive = TRUE,
                    scrollX = FALSE,
                    order = list(list(0, "desc")), # order by hidden abs_correlation
                    columnDefs = list(
                        list(visible = FALSE, targets = 0),
                        list(width = "42%", targets = 1), # Trait
                        list(width = "16%", targets = 2), # Correlation
                        list(width = "18%", targets = 3), # Pvalue
                        list(width = "12%", targets = 4) # # Mice
                    ),
                    dom = "tip"
                ),
                callback = htmlwidgets::JS(
                    "table.settings()[0].aoColumns.forEach(function(col){",
                    "  if(col && col.nTh){col.nTh.style.whiteSpace='nowrap';col.nTh.style.wordBreak='normal';col.nTh.style.fontSize='90%';col.nTh.style.padding='6px 8px';}",
                    "});",
                    "table.on('draw.dt', function(){",
                    "  table.columns().every(function(){",
                    "    $(this.nodes()).css({'white-space':'normal','word-break':'break-word'});",
                    "  });",
                    "});"
                )
            )
            dt <- DT::formatRound(dt, columns = "Correlation", digits = 4)
            dt <- DT::formatStyle(dt, columns = c("Trait", "Correlation", "Pvalue", "# Mice"), `white-space` = "normal", `word-wrap` = "break-word", fontSize = "90%")
            dt
        })

        return(invisible(NULL))
    })
}
