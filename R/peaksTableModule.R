#' Peaks Table Module
#'
#' Lists peaks for the current trait and dataset context and allows selecting a peak
#' to drive allele effects. In additive mode, shows all additive peaks. In
#' interactive modes (sex/diet), shows split-by additive peaks when available.
#' Sex x Diet is attempted if files exist; otherwise, the table is empty.
#'
#' @param id Module ID
#' @param trait_reactive Reactive returning the current trait name (character)
#' @param dataset_group_reactive Reactive returning the selected dataset group (character)
#' @param interaction_type_reactive Reactive returning interaction type: "none", "sex", "diet", or "sex_diet"
#' @param import_reactives List-like object containing at least `file_directory`
#' @param set_selected_peak_fn A reactiveVal function from the scan plot module to set the selected peak (call with data.frame row)
#' @param clear_diff_peak_1_fn Optional reactiveVal function to clear diff peak 1 by calling with NULL
#' @param clear_diff_peak_2_fn Optional reactiveVal function to clear diff peak 2 by calling with NULL
#'
#' @return Renders a DT table UI and updates selection via `set_selected_peak_fn`
#' @importFrom DT DTOutput renderDT datatable
#' @importFrom data.table fread as.data.table
#' @export
peaksTableUI <- function(id) {
    ns <- shiny::NS(id)
    DT::DTOutput(ns("peaks_table"))
}

#' @rdname peaksTableUI
#' @export
peaksTableServer <- function(
    id,
    trait_reactive,
    dataset_group_reactive,
    interaction_type_reactive,
    import_reactives,
    set_selected_peak_fn,
    clear_diff_peak_1_fn = NULL,
    clear_diff_peak_2_fn = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Local cache for peaks table
        peaks_table_full_rv <- shiny::reactiveVal(NULL)

        # Helper to resolve import_reactives to a list
        resolve_import <- function() {
            ir <- import_reactives
            if (is.function(ir)) {
                ir <- ir()
            }
            return(ir)
        }

        # Helper to map interactive split-by peak files
        get_interactive_peak_filepaths <- function(dataset_group, interaction_type) {
            base_path <- "/data/dev/miniViewer_3.0/"
            base_name <- gsub(",[[:space:]]*(interactive|additive).*$", "", dataset_group, ignore.case = TRUE)
            base_name <- trimws(base_name)

            dataset_component <- NULL
            if (grepl("Liver Genes", base_name, ignore.case = TRUE)) {
                dataset_component <- "liver_genes"
            } else if (grepl("Liver Lipids", base_name, ignore.case = TRUE)) {
                dataset_component <- "liver_lipids"
            } else if (grepl("Clinical Traits", base_name, ignore.case = TRUE)) {
                dataset_component <- "clinical_traits"
            } else if (grepl("Plasma.*Metabol|plasma.*metabolite", base_name, ignore.case = TRUE)) {
                dataset_component <- "plasma_metabolites"
            } else if (grepl("Liver Isoforms", base_name, ignore.case = TRUE)) {
                dataset_component <- "liver_isoforms"
            } else {
                return(NULL)
            }

            if (interaction_type == "sex") {
                return(list(
                    file1 = file.path(base_path, paste0("DO1200_", dataset_component, "_qtlxsex_peaks_in_female_mice_additive.csv")),
                    file2 = file.path(base_path, paste0("DO1200_", dataset_component, "_qtlxsex_peaks_in_male_mice_additive.csv")),
                    labels = c("Female", "Male")
                ))
            } else if (interaction_type == "diet") {
                return(list(
                    file1 = file.path(base_path, paste0("DO1200_", dataset_component, "_qtlxdiet_peaks_in_HC_mice_additive.csv")),
                    file2 = file.path(base_path, paste0("DO1200_", dataset_component, "_qtlxdiet_peaks_in_HF_mice_additive.csv")),
                    labels = c("HC Diet", "HF Diet")
                ))
            } else if (interaction_type == "sex_diet") {
                return(NULL)
            }
            return(NULL)
        }

        # Helper to map QTLx summary peaks file (single file for interaction)
        get_qtlx_summary_filepath <- function(dataset_group, interaction_type) {
            base_path <- "/data/dev/miniViewer_3.0/"
            base_name <- gsub(",[[:space:]]*(interactive|additive).*$", "", dataset_group, ignore.case = TRUE)
            base_name <- trimws(base_name)

            dataset_component <- NULL
            if (grepl("Liver Genes", base_name, ignore.case = TRUE)) {
                dataset_component <- "liver_genes"
            } else if (grepl("Liver Lipids", base_name, ignore.case = TRUE)) {
                dataset_component <- "liver_lipids"
            } else if (grepl("Clinical Traits", base_name, ignore.case = TRUE)) {
                dataset_component <- "clinical_traits"
            } else if (grepl("Plasma.*Metabol|plasma.*metabolite", base_name, ignore.case = TRUE)) {
                dataset_component <- "plasma_metabolites"
            } else if (grepl("Liver Isoforms", base_name, ignore.case = TRUE)) {
                dataset_component <- "liver_isoforms"
            } else {
                return(NULL)
            }

            suffix <- if (interaction_type == "sex") {
                "qtlxsex_peaks.csv"
            } else if (interaction_type == "diet") {
                "qtlxdiet_peaks.csv"
            } else if (interaction_type == "sex_diet") {
                "qtlxsexbydiet_peaks.csv"
            } else {
                return(NULL)
            }

            file.path(base_path, paste0("DO1200_", dataset_component, "_all_mice_", suffix))
        }

        # Build peaks table for current state
        peaks_table_display <- shiny::reactive({
            trait <- trait_reactive()
            dataset_group <- dataset_group_reactive()
            shiny::req(trait, nzchar(trait), dataset_group)

            interaction_type <- interaction_type_reactive()
            ir <- resolve_import()

            # Additive mode
            if (is.null(interaction_type) || interaction_type == "none") {
                trait_type_val <- get_trait_type(ir, dataset_group)
                all_peaks <- peak_finder(
                    file_dir = ir$file_directory,
                    selected_dataset = dataset_group,
                    selected_trait = resolve_trait_aliases_for_peaks(ir, dataset_group, trait),
                    trait_type = trait_type_val,
                    cache_env = new.env(parent = emptyenv()),
                    use_cache = TRUE
                )
                if (is.null(all_peaks) || nrow(all_peaks) == 0) {
                    peaks_table_full_rv(NULL)
                    return(NULL)
                }
                peaks_table_full_rv(all_peaks)
                display <- all_peaks
                display$Source <- "Additive"
                display <- display[, c("Source", "marker", "qtl_chr", "qtl_pos", "qtl_lod", intersect(colnames(display), c("qtl_pval", "pval"))), drop = FALSE]
                # Explicitly rename core columns only; keep qtl_pval as-is
                names(display)[names(display) == "qtl_chr"] <- "chr"
                names(display)[names(display) == "qtl_pos"] <- "pos"
                names(display)[names(display) == "qtl_lod"] <- "lod"
                return(display)
            }

            # Interactive mode: combine any available sources (split-by and QTLx summary)
            if (interaction_type %in% c("sex", "diet", "sex_diet")) {
                pieces <- list()

                # Split-by rows removed from peaks table; defaults rendered separately

                # QTLx summary
                qtlx_file <- get_qtlx_summary_filepath(dataset_group, interaction_type)
                if (!is.null(qtlx_file) && file.exists(qtlx_file)) {
                    dfs <- data.table::fread(qtlx_file)
                    trait_col_s <- if ("gene_symbol" %in% colnames(dfs)) "gene_symbol" else if ("phenotype" %in% colnames(dfs)) "phenotype" else if ("metabolite" %in% colnames(dfs)) "metabolite" else if ("metabolite_name" %in% colnames(dfs)) "metabolite_name" else NULL
                    if (!is.null(trait_col_s)) {
                        dts <- data.table::as.data.table(dfs)[get(trait_col_s) == trait]
                        if (nrow(dts) > 0) {
                            label <- if (interaction_type == "sex") "QTLx Sex" else if (interaction_type == "diet") "QTLx Diet" else "QTLx Sex×Diet"
                            dts[, Source := label]
                            pieces[[length(pieces) + 1]] <- dts
                        }
                    }
                }

                if (length(pieces) == 0) {
                    peaks_table_full_rv(NULL)
                    return(NULL)
                }

                combined <- data.table::rbindlist(pieces, fill = TRUE)
                peaks_table_full_rv(as.data.frame(combined))
                display <- as.data.frame(combined)

                # Build display columns
                trait_cols <- intersect(colnames(display), c("gene_symbol", "phenotype"))
                lod_diff_candidates <- intersect(colnames(display), c("lod_diff", "qtl_lod_diff"))
                selected_cols <- c(
                    "Source",
                    intersect(colnames(display), c("marker")),
                    trait_cols,
                    intersect(colnames(display), c("qtl_chr", "qtl_pos", "qtl_lod", "qtl_pval", "pval_diff")),
                    lod_diff_candidates
                )
                display <- display[, unique(selected_cols), drop = FALSE]
                # Explicitly rename core columns only; keep qtl_pval as-is
                names(display)[names(display) == "qtl_chr"] <- "chr"
                names(display)[names(display) == "qtl_pos"] <- "pos"
                names(display)[names(display) == "qtl_lod"] <- "lod"
                if ("gene_symbol" %in% colnames(display)) {
                    colnames(display)[match("gene_symbol", colnames(display))] <- "trait"
                } else if ("phenotype" %in% colnames(display)) {
                    colnames(display)[match("phenotype", colnames(display))] <- "trait"
                }
                return(display)
            }

            peaks_table_full_rv(NULL)
            return(NULL)
        }) |> shiny::debounce(150)

        output$peaks_table <- DT::renderDT({
            tbl <- peaks_table_display()
            if (is.null(tbl) || nrow(tbl) == 0) {
                return(DT::datatable(data.frame(Info = "No peaks available for this selection"), options = list(dom = "t"), rownames = FALSE))
            }

            # Build prettier, structured view
            has_lod_diff <- "lod_diff" %in% colnames(tbl)

            # Normalize p-value column names across sources
            if (!("qtl_pval" %in% colnames(tbl)) && ("pval" %in% colnames(tbl))) {
                tbl$qtl_pval <- tbl$pval
            }
            if (!("pval_diff" %in% colnames(tbl)) && ("diff_pval" %in% colnames(tbl))) {
                tbl$pval_diff <- tbl$diff_pval
            }

            # Create single formatted position column and drop separate chr/pos
            if (all(c("chr", "pos") %in% colnames(tbl))) {
                # Safely format position with 2 decimals
                suppressWarnings({
                    numeric_pos <- as.numeric(tbl$pos)
                })
                formatted_pos <- ifelse(is.na(numeric_pos), as.character(tbl$pos), sprintf("%.2f", numeric_pos))
                tbl$position <- paste0(tbl$chr, ":", formatted_pos)
            }

            # Columns to show (no trait; use combined position)
            base_cols <- c("Source", "position", "lod", if (has_lod_diff) "lod_diff", "qtl_pval", "pval_diff", "marker")
            cols_present <- intersect(base_cols, colnames(tbl))
            tbl_view <- tbl[, cols_present, drop = FALSE]

            # Human-friendly column titles
            # Dynamic LOD title based on interactive state
            is_interactive <- {
                it <- interaction_type_reactive()
                !is.null(it) && it != "none"
            }
            lod_title <- if (is_interactive) "LOD Int" else "LOD Add"

            pretty_names <- c(
                Source = "Source",
                position = "Position",
                lod = lod_title,
                lod_diff = "LOD Diff",
                qtl_pval = "QtL Pval",
                pval_diff = "Diff Pval",
                marker = "Marker"
            )
            shown_titles <- unname(pretty_names[cols_present])

            # Determine numeric columns to format
            numeric_cols <- intersect(c("lod", "lod_diff"), colnames(tbl_view))
            pval_cols <- intersect(c("qtl_pval", "pval_diff"), colnames(tbl_view))

            # RowGroup by Source when present
            extensions <- c("RowGroup", "FixedHeader")
            row_group <- if ("Source" %in% colnames(tbl_view)) list(dataSrc = which(colnames(tbl_view) == "Source") - 1L) else NULL

            # Default order: by LOD desc when present
            lod_idx <- which(colnames(tbl_view) == "lod") - 1L
            order_opt <- if (length(lod_idx) == 1) list(list(lod_idx, "desc")) else list()

            dt <- DT::datatable(
                tbl_view,
                colnames = shown_titles,
                options = list(
                    pageLength = 8,
                    lengthChange = FALSE,
                    searching = FALSE,
                    ordering = TRUE,
                    dom = "tip",
                    rowGroup = row_group,
                    fixedHeader = TRUE,
                    order = order_opt,
                    columnDefs = list(
                        list(className = "dt-center", targets = which(colnames(tbl_view) %in% c("position", "lod", "lod_diff", "qtl_pval", "pval_diff")) - 1L),
                        list(className = "dt-left", targets = which(colnames(tbl_view) %in% c("marker")) - 1L)
                    )
                ),
                selection = "single",
                rownames = FALSE,
                escape = FALSE
            )

            # Numeric formatting
            if (length(numeric_cols) > 0) {
                dt <- DT::formatRound(dt, columns = numeric_cols, digits = 2)
            }
            if (length(pval_cols) > 0) {
                dt <- DT::formatSignif(dt, columns = pval_cols, digits = 3)
            }

            # Removed Cis/Trans styling and column per requirements

            dt
        })

        shiny::observeEvent(input$peaks_table_rows_selected,
            {
                sel <- input$peaks_table_rows_selected
                if (is.null(sel) || length(sel) == 0) {
                    return(NULL)
                }
                full <- peaks_table_full_rv()
                disp <- peaks_table_display()
                if (is.null(full) || is.null(disp)) {
                    return(NULL)
                }

                selected_row <- disp[sel, , drop = FALSE]
                marker_val <- if ("marker" %in% colnames(selected_row)) selected_row$marker else NA

                # Try to locate the corresponding full row (prefer marker match)
                idx <- integer(0)
                if (!is.na(marker_val) && "marker" %in% colnames(full)) {
                    idx <- which(full$marker == marker_val)
                }
                if (length(idx) == 0) {
                    chr_col <- if ("qtl_chr" %in% colnames(full)) "qtl_chr" else if ("chr" %in% colnames(full)) "chr" else NULL
                    pos_col <- if ("qtl_pos" %in% colnames(full)) "qtl_pos" else if ("pos" %in% colnames(full)) "pos" else NULL
                    lod_col <- if ("qtl_lod" %in% colnames(full)) "qtl_lod" else if ("lod" %in% colnames(full)) "lod" else NULL
                    if (!is.null(chr_col) && !is.null(pos_col) && !is.null(lod_col)) {
                        idx <- which(
                            as.character(full[[chr_col]]) == as.character(selected_row$chr) &
                                abs(full[[pos_col]] - as.numeric(selected_row$pos)) < 1e-6 &
                                abs(full[[lod_col]] - as.numeric(selected_row$lod)) < 1e-6
                        )
                    }
                }
                if (length(idx) == 0) {
                    return(NULL)
                }

                peak_row <- as.data.frame(full[idx[1], , drop = FALSE])

                # Clear any difference selections and set the additive/selected peak
                if (!is.null(clear_diff_peak_1_fn)) clear_diff_peak_1_fn(NULL)
                if (!is.null(clear_diff_peak_2_fn)) clear_diff_peak_2_fn(NULL)
                set_selected_peak_fn(peak_row)
            },
            ignoreInit = TRUE
        )
    })
}
