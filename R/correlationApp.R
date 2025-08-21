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
#' a chosen dataset. The table includes columns: trait, correlation_value, p_value.
#'
#' Correlation CSV files are expected in the directory:
#'   /mnt/rdrive/mkeller3/General/main_directory/correlations
#' and named like: "liver_genes_vs_clinical_traits_corr.csv".
#'
#' @param id Module ID.
#' @param import_reactives Reactive that returns a list containing file_directory and annotation_list.
#' @param main_par Reactive that returns a list containing selected_dataset and which_trait.
#' @return Invisibly returns NULL. UI is rendered via outputs.
#' @importFrom data.table fread as.data.table melt setDT setnames
#' @importFrom DT DTOutput renderDT datatable
#' @importFrom shiny moduleServer NS reactive req validate need observeEvent updateSelectInput updateSelectizeInput
#' @export
correlationServer <- function(id, import_reactives, main_par) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Prefer Docker-mounted correlation dir(s)
        correlations_dir <- "/data/dev/miniViewer_3.0"


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
            return(NULL)
        }

        token_to_label <- function(token) {
            switch(token,
                liver_genes = "Liver Genes",
                liver_isoforms = "Liver Isoforms",
                liver_lipids = "Liver Lipids",
                plasma_metabolites = "Plasma Metabolites",
                clinical_traits = "Clinical Traits",
                token
            )
        }

        # Load annotation list (for mapping gene ids -> symbols) once
        annotation_map_rv <- shiny::reactiveVal(NULL)
        get_gene_symbol_map <- function() {
            map <- annotation_map_rv()
            if (!is.null(map)) {
                return(map)
            }
            # Try Docker-mounted RDS first; if missing, fall back to import_reactives()$annotation_list
            ann <- NULL
            ann_candidates <- c(
                "/data/dev/miniViewer_3.0/annotationlist.RDS",
                "/data/dev/miniviewier_3.0/annotationlist.RDS"
            )
            for (p in ann_candidates) {
                if (file.exists(p)) {
                    try(
                        {
                            ann <- readRDS(p)
                        },
                        silent = TRUE
                    )
                    if (!is.null(ann)) break
                }
            }
            if (is.null(ann)) {
                imp <- import_reactives()
                if (!is.null(imp) && !is.null(imp$annotation_list)) ann <- imp$annotation_list
            }
            if (is.null(ann) || is.null(ann$genes)) {
                return(NULL)
            }
            genes_dt <- ann$genes
            # Robustly detect id and symbol columns
            id_candidates <- c("gene.id", "gene_id", "ensembl_id", "ensembl_gene_id", "ensemblid", "geneid")
            sym_candidates <- c("symbol", "gene_symbol", "mgi_symbol", "Gene.symbol", "gene.name")
            id_col <- id_candidates[id_candidates %in% colnames(genes_dt)][1]
            sym_col <- sym_candidates[sym_candidates %in% colnames(genes_dt)][1]
            if (is.na(id_col) || is.na(sym_col) || is.null(id_col) || is.null(sym_col)) {
                return(NULL)
            }
            # Build a named vector: names = gene ids, values = symbols
            valid_rows <- !is.na(genes_dt[[id_col]]) & nzchar(as.character(genes_dt[[id_col]]))
            ids <- as.character(genes_dt[[id_col]][valid_rows])
            syms <- as.character(genes_dt[[sym_col]][valid_rows])
            syms[is.na(syms)] <- ""
            out_map <- syms
            names(out_map) <- ids
            annotation_map_rv(out_map)
            out_map
        }

        # Utility: List available correlation files and parse source/target tokens
        list_available_pairs <- shiny::reactive({
            # Expected/fallback files to show even if directory scan fails
            expected_files <- c(
                file.path(correlations_dir, "clinical_traits_vs_clinical_traits_corr.csv"),
                file.path(correlations_dir, "liver_genes_vs_clinical_traits_corr.csv"),
                file.path(correlations_dir, "liver_genes_vs_liver_lipids_corr.csv"),
                file.path(correlations_dir, "liver_genes_vs_plasma_metabolites_corr.csv"),
                file.path(correlations_dir, "liver_lipids_vs_liver_lipids_corr.csv"),
                file.path(correlations_dir, "plasma_metabolites_vs_plasma_metabolites_corr.csv"),
                file.path(correlations_dir, "clinical_traits_vs_liver_lipids_corr.csv"),
                file.path(correlations_dir, "clinical_traits_vs_plasma_metabolites_corr.csv"),
                file.path(correlations_dir, "liver_lipids_vs_plasma_metabolites_corr.csv")
            )

            files <- character(0)
            if (dir.exists(correlations_dir)) {
                files <- list.files(correlations_dir, pattern = "_corr\\.csv$", full.names = TRUE)
            } else {
                warning("correlationServer: correlations_dir not found: ", correlations_dir)
            }

            # Include expected files regardless of current presence; parsing will still work
            files <- unique(c(files, expected_files))
            if (length(files) == 0) {
                return(data.frame(file = character(0), source = character(0), target = character(0)))
            }

            # Parse names like foo_vs_bar_corr.csv
            parse_one <- function(fp) {
                bn <- basename(fp)
                m <- regexec("^([a-z0-9_]+)_vs_([a-z0-9_]+)_corr\\.csv$", tolower(bn))
                r <- regmatches(tolower(bn), m)[[1]]
                if (length(r) == 3) {
                    return(data.frame(file = fp, source = r[2], target = r[3]))
                }
                return(NULL)
            }
            parsed <- lapply(files, parse_one)
            parsed <- parsed[!vapply(parsed, is.null, logical(1))]
            if (length(parsed) == 0) {
                return(data.frame(file = character(0), source = character(0), target = character(0)))
            }
            data.table::as.data.table(do.call(rbind, parsed))
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
                # Files where current type is on either side
                subset_dt <- pairs_dt[(source == token) | (target == token)]
                if (nrow(subset_dt) == 0) {
                    shiny::updateSelectInput(session, "correlation_dataset_selector",
                        choices = character(0), selected = character(0)
                    )
                    return()
                }
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
                return(list(mode = "row", key = trait_string, phenotype_col = pheno_col))
            }

            return(NULL)
        }

        # Build the correlation table for display
        correlation_table <- shiny::reactive({
            shiny::req(input$correlation_dataset_selector)
            trait <- current_trait_reactive()
            shiny::req(trait)

            file_path <- input$correlation_dataset_selector
            if (!file.exists(file_path)) {
                shiny::showNotification(paste("Correlation file not found:", basename(file_path)), type = "error", duration = 4)
                return(data.frame(trait = character(0), correlation_value = numeric(0), p_value = numeric(0), num_mice = numeric(0)))
            }
            # Determine which side (source/target) corresponds to the current trait type
            pairs_dt <- list_available_pairs()
            shiny::req(nrow(pairs_dt) > 0)
            token <- current_type_token()
            row <- pairs_dt[file == file_path]
            shiny::req(nrow(row) == 1)

            side_token <- if (row$source[1] == token) row$source[1] else if (row$target[1] == token) row$target[1] else token

            key_info <- resolve_trait_keys(trait, side_token, file_path, import_reactives())
            if (is.null(key_info)) {
                shiny::showNotification(paste("Trait not found in correlation file:", trait), type = "warning", duration = 4)
                return(data.frame(trait = character(0), correlation_value = numeric(0), p_value = numeric(0), num_mice = numeric(0)))
            }

            if (key_info$mode == "column") {
                # Read only phenotype and the specific column for efficiency
                select_cols <- c("phenotype", "Phenotype", key_info$key)
                select_cols <- select_cols[select_cols %in% names(data.table::fread(file_path, nrows = 0))]
                dt <- data.table::fread(file_path, select = select_cols)
                # Normalize phenotype column name
                if ("Phenotype" %in% names(dt)) data.table::setnames(dt, "Phenotype", "phenotype")
                if (!"phenotype" %in% names(dt)) dt[, phenotype := seq_len(.N)]
                # Build output table
                out <- dt[, .(trait = phenotype, correlation_value = .SD[[1]]), .SDcols = key_info$key]

                # Add p-values from companion file if available
                pval_path <- sub("_corr\\.csv$", "_pval.csv", file_path)
                if (file.exists(pval_path)) {
                    psel <- c("phenotype", "Phenotype", key_info$key)
                    psel <- psel[psel %in% names(data.table::fread(pval_path, nrows = 0))]
                    pdt <- data.table::fread(pval_path, select = psel)
                    if ("Phenotype" %in% names(pdt)) data.table::setnames(pdt, "Phenotype", "phenotype")
                    if (!"phenotype" %in% names(pdt)) pdt[, phenotype := seq_len(.N)]
                    pmerge <- pdt[, .(trait = phenotype, p_value = .SD[[1]]), .SDcols = key_info$key]
                    out <- merge(out, pmerge, by = "trait", all.x = TRUE, sort = FALSE)
                } else {
                    out[, p_value := as.numeric(NA)]
                }
                # Add num_mice from companion file if available
                nmouse_path <- sub("_corr\\.csv$", "_num_mice.csv", file_path)
                if (file.exists(nmouse_path)) {
                    nsel <- c("phenotype", "Phenotype", key_info$key)
                    nsel <- nsel[nsel %in% names(data.table::fread(nmouse_path, nrows = 0))]
                    ndt <- data.table::fread(nmouse_path, select = nsel)
                    if ("Phenotype" %in% names(ndt)) data.table::setnames(ndt, "Phenotype", "phenotype")
                    if (!"phenotype" %in% names(ndt)) ndt[, phenotype := seq_len(.N)]
                    nmerge <- ndt[, .(trait = phenotype, num_mice = .SD[[1]]), .SDcols = key_info$key]
                    out <- merge(out, nmerge, by = "trait", all.x = TRUE, sort = FALSE)
                } else {
                    out[, num_mice := as.numeric(NA)]
                }
            } else {
                # Row mode: read entire file, find the row, and transpose to long format
                dt_full <- data.table::fread(file_path)
                # Normalize phenotype column name
                if ("Phenotype" %in% names(dt_full)) data.table::setnames(dt_full, "Phenotype", "phenotype")
                if (!"phenotype" %in% names(dt_full)) {
                    shiny::showNotification("Correlation file missing 'phenotype' column.", type = "error", duration = NULL)
                    return(data.frame(trait = character(0), correlation_value = numeric(0), p_value = numeric(0), num_mice = numeric(0)))
                }
                row_dt <- dt_full[phenotype == key_info$key]
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

                # Add p-values from companion pval file (same row)
                pval_path <- sub("_corr\\.csv$", "_pval.csv", file_path)
                if (file.exists(pval_path)) {
                    pdt_full <- data.table::fread(pval_path)
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

                # Add num_mice from companion num_mice file (same row)
                nmouse_path <- sub("_corr\\.csv$", "_num_mice.csv", file_path)
                if (file.exists(nmouse_path)) {
                    nfull <- data.table::fread(nmouse_path)
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

                # Map gene ids -> symbols when the other side is genes (columns), i.e., when traits look like 'liver_<gene_id>'
                has_liver_prefix <- any(grepl("^liver_ENSMUSG", out$trait, perl = TRUE))
                if (has_liver_prefix) {
                    id_to_symbol <- get_gene_symbol_map()
                    gene_ids <- sub("^liver_", "", out$trait)
                    if (!is.null(id_to_symbol)) {
                        mapped <- id_to_symbol[gene_ids]
                    } else {
                        mapped <- rep(NA_character_, length(gene_ids))
                    }
                    # Display symbol when available; otherwise show bare gene id (no liver_ prefix)
                    display <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, gene_ids)
                    out$trait <- display
                }
            }

            # Prepare display formatting and default sort by absolute correlation
            if ("p_value" %in% names(out)) {
                out[, p_value := ifelse(is.na(p_value), NA_character_, format(p_value, digits = 3, scientific = TRUE))]
            }
            out[, abs_correlation := abs(correlation_value)]
            out <- out[order(-abs_correlation, -correlation_value)]
            out <- data.table::as.data.table(out)
            out
        }) %>% shiny::debounce(150)

        # Debounced search query
        search_query <- shiny::reactive({
            val <- input$correlation_search
            if (is.null(val)) "" else trimws(as.character(val))
        }) |> shiny::debounce(150)

        # Filtered table by trait search; preserves ordering by absolute correlation
        filtered_table <- shiny::reactive({
            tbl <- correlation_table()
            q <- search_query()
            if (!is.null(tbl) && nrow(tbl) > 0 && nzchar(q)) {
                keep <- tryCatch(grepl(q, tbl$trait, ignore.case = TRUE, perl = TRUE),
                    error = function(e) grepl(q, tbl$trait, ignore.case = TRUE, fixed = TRUE)
                )
                tbl <- tbl[keep]
            }
            tbl
        })

        output$correlation_table <- DT::renderDT({
            tbl <- filtered_table()
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
