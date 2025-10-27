#' SNP Association Tab UI
#'
#' Creates UI for SNP association analysis at a selected QTL peak
#'
#' @param current_category Optional string of the current dataset category
#' @return Shiny UI for the SNP association tab
#' @importFrom shiny tags div h6 numericInput checkboxInput sliderInput
#' @export
snp_association_tab_ui <- function(current_category = NULL) {
    message("SNP Association Tab UI: initializing")

    # Check if current category is Liver Splice Junctions (no SNP association for these)
    if (!is.null(current_category) && grepl("Liver.*Splice.*Junction", current_category, ignore.case = TRUE)) {
        return(
            shiny::div(
                style = "padding: 20px; text-align: center; color: #7f8c8d; font-size: 1.1em;",
                shiny::tags$em("No SNP association for liver splice junctions")
            )
        )
    }

    shiny::div(
        style = "padding: 15px;",

        # Analysis Parameters Section
        shiny::tags$div(
            style = "margin-bottom: 20px; padding: 15px; background-color: #f8f9fa; border-radius: 5px;",
            shiny::h6("Analysis Parameters", style = "color: #2c3e50; font-weight: bold; margin-bottom: 10px;"),

            # Window size parameter
            shiny::div(
                style = "margin-bottom: 15px;",
                shiny::numericInput(
                    shiny::NS("app_controller", "snp_window"),
                    label = "Window size (Mb):",
                    value = 1.5,
                    min = 0.5,
                    max = 10,
                    step = 0.5,
                    width = "100%"
                ),
                shiny::tags$p(
                    "Window on either side of QTL peak for SNP association",
                    style = "font-size: 11px; color: #7f8c8d; margin: 5px 0 0 0;"
                )
            ),

            # LOD drop parameter
            shiny::div(
                style = "margin-bottom: 15px;",
                shiny::numericInput(
                    shiny::NS("app_controller", "snp_lod_drop"),
                    label = "LOD drop threshold:",
                    value = 1.5,
                    min = 0.5,
                    max = 5,
                    step = 0.5,
                    width = "100%"
                ),
                shiny::tags$p(
                    "LOD drop for highlighting top variants",
                    style = "font-size: 11px; color: #7f8c8d; margin: 5px 0 0 0;"
                )
            ),

            # Show SDP panel option
            shiny::div(
                style = "margin-bottom: 10px;",
                shiny::checkboxInput(
                    shiny::NS("app_controller", "snp_show_sdp"),
                    label = "Show strain distribution pattern panel",
                    value = TRUE,
                    width = "100%"
                )
            )
        ),

        # SNP Association Plot Container
        shiny::tags$div(
            style = "margin-top: 20px;",
            shiny::h6("SNP Association Plot", style = "color: #2c3e50; font-weight: bold; margin-bottom: 10px;"),
            shiny::uiOutput(shiny::NS("app_controller", "snp_plot_container"))
        ),

        # Top Variants Table
        shiny::tags$div(
            style = "margin-top: 20px;",
            shiny::h6("Top Variants", style = "color: #2c3e50; font-weight: bold; margin-bottom: 10px;"),
            shiny::uiOutput(shiny::NS("app_controller", "snp_table_container"))
        )
    )
}

#' Helper function to rank-z transform phenotype data
#'
#' @param x Numeric vector
#' @return Rank-z transformed vector
#' @export
rankz <- function(x) {
    x <- rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
    return(qnorm(x))
}

#' Load cross object for SNP association
#'
#' This function loads the cross object once and caches it
#'
#' @return Cross object or NULL if loading fails
#' @importFrom qtl2 create_variant_query_func create_gene_query_func
#' @export
load_cross_for_snp <- function() {
    data_dir <- "/data/dev/DO_mapping_files"
    cross_file <- file.path(data_dir, "cross_DO1200_grcm39.rds")

    message(sprintf("SNP Association: Checking for cross file at: %s", cross_file))
    message(sprintf(
        "SNP Association: File exists: %s, Readable: %s",
        file.exists(cross_file),
        if (file.exists(cross_file)) file.access(cross_file, 4) == 0 else "N/A"
    ))

    if (!file.exists(cross_file)) {
        warning(paste("SNP Association: Cross file not found:", cross_file))
        return(NULL)
    }

    tryCatch(
        {
            message("SNP Association: Loading cross object (this may take 10-20 seconds)...")
            cross <- readRDS(cross_file)
            message(sprintf(
                "SNP Association: Cross object loaded successfully. Components: %s",
                paste(names(cross)[1:5], collapse = ", ")
            ))
            message(sprintf("SNP Association: Cross has %d phenotypes", ncol(cross$pheno)))
            return(cross)
        },
        error = function(e) {
            warning(paste("SNP Association: Error loading cross object:", e$message))
            message(sprintf("SNP Association: Full error: %s", toString(e)))
            return(NULL)
        }
    )
}

#' Load genoprobs for SNP association
#'
#' @return Genoprobs object or NULL if loading fails
#' @export
load_genoprobs_for_snp <- function() {
    data_dir <- "/data/dev/DO_mapping_files"
    genoprobs_file <- file.path(data_dir, "genoprobs_DO1200_grcm39_fstindex.rds")

    message(sprintf("SNP Association: Checking for genoprobs file at: %s", genoprobs_file))
    message(sprintf(
        "SNP Association: File exists: %s, Readable: %s",
        file.exists(genoprobs_file),
        if (file.exists(genoprobs_file)) file.access(genoprobs_file, 4) == 0 else "N/A"
    ))

    if (!file.exists(genoprobs_file)) {
        warning(paste("SNP Association: Genoprobs file not found:", genoprobs_file))
        return(NULL)
    }

    tryCatch(
        {
            message("SNP Association: Loading genoprobs (this is quick, ~1 second)...")
            genoprobs <- readRDS(genoprobs_file)
            message(sprintf(
                "SNP Association: Genoprobs loaded successfully. Class: %s",
                paste(class(genoprobs), collapse = ", ")
            ))
            return(genoprobs)
        },
        error = function(e) {
            warning(paste("SNP Association: Error loading genoprobs:", e$message))
            message(sprintf("SNP Association: Full error: %s", toString(e)))
            return(NULL)
        }
    )
}

#' Perform SNP association analysis
#'
#' @param cross Cross object
#' @param genoprobs Genoprobs object
#' @param phenotype Phenotype name
#' @param chr Chromosome
#' @param pos Position (Mb)
#' @param window Window size (Mb)
#' @param interaction_type Type of interaction ("none", "sex", "diet")
#' @param ncores Number of cores (default 1)
#' @return List with lod and snpinfo, or NULL if analysis fails
#' @importFrom qtl2 scan1snps create_variant_query_func create_gene_query_func
#' @export
perform_snp_association <- function(cross, genoprobs, phenotype, chr, pos,
                                    window = 1.5, interaction_type = "none",
                                    ncores = 1) {
    if (is.null(cross) || is.null(genoprobs)) {
        warning("SNP Association: Cross or genoprobs object is NULL")
        return(list(error = "Required data files not available"))
    }

    # Validate inputs
    if (is.null(phenotype) || !nzchar(phenotype)) {
        warning("SNP Association: Invalid phenotype")
        return(list(error = "Invalid phenotype name"))
    }

    if (is.null(chr) || is.null(pos)) {
        warning("SNP Association: Invalid chromosome or position")
        return(list(error = "Invalid chromosome or position"))
    }

    # Check if phenotype exists in cross object, or try to map gene symbol to transcript
    pheno_to_use <- phenotype

    if (!phenotype %in% colnames(cross$pheno)) {
        message(sprintf("SNP Association: '%s' not found directly, attempting gene symbol -> transcript mapping...", phenotype))

        # Try to map gene symbol to transcript ID
        if (!is.null(cross$gene_annos) && !is.null(cross$transcript_annos)) {
            # Convert to data frame if it's a tibble to ensure proper indexing
            gene_annos_df <- as.data.frame(cross$gene_annos)
            transcript_annos_df <- as.data.frame(cross$transcript_annos)

            # Find gene by symbol using which() for safer indexing
            gene_idx <- which(gene_annos_df$gene_symbol == phenotype)

            if (length(gene_idx) > 0) {
                gene_id <- gene_annos_df$gene_id[gene_idx[1]]
                message(sprintf("SNP Association: Found gene ID %s for symbol %s", gene_id, phenotype))

                # Find transcripts for this gene
                trans_idx <- which(transcript_annos_df$gene_id == gene_id)

                if (length(trans_idx) > 0) {
                    # Use the first transcript
                    transcript_id <- transcript_annos_df$transcript_id[trans_idx[1]]
                    pheno_to_use <- paste0("liver_", transcript_id)

                    message(sprintf("SNP Association: Checking for transcript phenotype: %s", pheno_to_use))

                    if (pheno_to_use %in% colnames(cross$pheno)) {
                        message(sprintf("SNP Association: Successfully mapped to transcript: %s", pheno_to_use))
                    } else {
                        message(sprintf("SNP Association: Transcript %s not found in cross phenotypes", pheno_to_use))
                        return(list(error = paste0("Gene '", phenotype, "' found but its transcript data not available in cross object")))
                    }
                } else {
                    message(sprintf("SNP Association: No transcripts found for gene %s", gene_id))
                    return(list(error = paste0("Gene '", phenotype, "' found but has no transcripts in cross object")))
                }
            } else {
                message(sprintf(
                    "SNP Association: Gene symbol '%s' not found in annotations. Available phenotypes (first 10): %s",
                    phenotype, paste(head(colnames(cross$pheno), 10), collapse = ", ")
                ))
                return(list(error = paste0("Phenotype '", phenotype, "' not available for SNP association.")))
            }
        } else {
            message("SNP Association: No gene annotations available for mapping")
            return(list(error = paste0("Phenotype '", phenotype, "' not found and gene mapping not available")))
        }
    }

    tryCatch(
        {
            message(sprintf(
                "SNP Association: phenotype=%s, chr=%s, pos=%.2f, window=%.2f, interaction=%s",
                phenotype, chr, pos, window, interaction_type
            ))

            # Make a copy of cross to avoid modifying the original
            cross_copy <- cross

            # Use the mapped phenotype name (might be transcript ID)
            message(sprintf("SNP Association: Using phenotype column: %s", pheno_to_use))

            # Perform rankz transformation
            cross_copy$pheno[, pheno_to_use] <- rankz(cross_copy$pheno[, pheno_to_use])

            # Set up covariate matrices
            # Additive covariate matrix (always the same)
            addcovar <- model.matrix(as.formula("~GenLit+Sex*Diet"), cross_copy$covar)[, -1, drop = FALSE]

            # Interactive covariate matrix depends on interaction type
            intcovar <- NULL
            if (interaction_type == "sex") {
                intcovar <- model.matrix(as.formula("~Sex"), cross_copy$covar)[, -1, drop = FALSE]
            } else if (interaction_type == "diet") {
                intcovar <- model.matrix(as.formula("~Diet"), cross_copy$covar)[, -1, drop = FALSE]
            }

            # Set up variant query function
            data_dir <- "/data/dev/DO_mapping_files"
            dbfile <- file.path(data_dir, "variants_db_by_chr", paste0("founder_variants_chr", chr, ".sqlite"))

            if (!file.exists(dbfile)) {
                warning(paste("Variant database file not found:", dbfile))
                return(list(error = paste0("Variant database not found for chromosome ", chr)))
            }

            query_variants <- qtl2::create_variant_query_func(
                dbfile,
                table_name = "variants",
                chr_field = "chr",
                pos_field = "pos",
                id_field = "variant_id",
                sdp_field = "sdp_num"
            )

            # Query genes in the region
            genes <- tryCatch(
                {
                    query_genes <- qtl2::create_gene_query_func(
                        dbfile,
                        table_name = "genes",
                        chr_field = "gene_chr",
                        start_field = "gene_start",
                        stop_field = "gene_end",
                        name_field = "gene_symbol",
                        strand_field = "gene_strand"
                    )
                    query_genes(chr, pos - window, pos + window)
                },
                error = function(e) {
                    message(paste("Could not query genes:", e$message))
                    NULL
                }
            )

            # Ensure chromosome is a string for qtl2 functions
            chr_str <- as.character(chr)

            # For FST genoprobs, we need to manually subset to the chromosome
            # because the subset() method doesn't work properly
            message(sprintf("Subsetting genoprobs for chromosome %s...", chr_str))

            # Extract just the chromosome we need by reading the FST file directly
            tryCatch(
                {
                    fst_file <- paste0(genoprobs$fst, "_", chr_str, ".fst")
                    message(sprintf("Reading FST file: %s", fst_file))

                    if (!file.exists(fst_file)) {
                        stop(sprintf("FST file not found: %s", fst_file))
                    }

                    # Read the FST data for this chromosome
                    gp_data <- fst::read_fst(fst_file)
                    message(sprintf("FST data loaded: %d x %d", nrow(gp_data), ncol(gp_data)))

                    # Get dimensions from genoprobs metadata
                    chr_idx <- which(genoprobs$chr == chr_str)
                    n_ind <- genoprobs$dim[1, chr_idx]
                    n_alleles <- genoprobs$dim[2, chr_idx]
                    n_pos <- genoprobs$dim[3, chr_idx]

                    message(sprintf("Expected dimensions: %d ind x %d alleles x %d pos", n_ind, n_alleles, n_pos))
                    message(sprintf("FST dimensions: %d rows x %d cols", nrow(gp_data), ncol(gp_data)))

                    # FST stores as (n_ind * n_alleles) rows x n_pos columns
                    # We need to reshape to n_ind x n_alleles x n_pos
                    if (nrow(gp_data) != n_ind * n_alleles || ncol(gp_data) != n_pos) {
                        stop(sprintf(
                            "FST dimension mismatch: expected %d x %d, got %d x %d",
                            n_ind * n_alleles, n_pos, nrow(gp_data), ncol(gp_data)
                        ))
                    }

                    message(sprintf("Reshaping genoprobs: %d ind x %d alleles x %d pos", n_ind, n_alleles, n_pos))

                    # Convert to 3D array
                    gp_matrix <- as.matrix(gp_data)
                    gp_array <- array(dim = c(n_ind, n_alleles, n_pos))

                    # Reshape: FST has rows as (ind1_allele1, ind1_allele2, ..., ind2_allele1, ...)
                    for (i in 1:n_ind) {
                        for (j in 1:n_alleles) {
                            row_idx <- (i - 1) * n_alleles + j
                            gp_array[i, j, ] <- gp_matrix[row_idx, ]
                        }
                    }

                    # Set dimnames using the stored dimnames from genoprobs
                    # genoprobs$dimnames[[chr_idx]] is a list with 3 elements: [ind_names, allele_names, marker_names]
                    dimnames(gp_array) <- genoprobs$dimnames[[chr_idx]]

                    # Reconstruct a single-chromosome genoprobs object
                    gp_chr <- list()
                    gp_chr[[chr_str]] <- gp_array
                    class(gp_chr) <- c("calc_genoprob", "list")
                    attr(gp_chr, "crosstype") <- attr(genoprobs, "crosstype")
                    attr(gp_chr, "is_x_chr") <- stats::setNames(chr_str == "X", chr_str)
                    attr(gp_chr, "alleles") <- attr(genoprobs, "alleles")
                    attr(gp_chr, "alleleprobs") <- attr(genoprobs, "alleleprobs")

                    message(sprintf("Genoprobs subset successful for chr %s", chr_str))
                },
                error = function(e) {
                    stop(sprintf("Failed to subset genoprobs for chr %s: %s", chr_str, e$message))
                }
            )

            # Perform SNP association
            message(sprintf("Running scan1snps for chromosome %s...", chr_str))
            snp_assoc <- qtl2::scan1snps(
                genoprobs = gp_chr,
                map = cross_copy$pmap,
                pheno = cross_copy$pheno[, pheno_to_use, drop = FALSE],
                kinship = cross_copy$kinship_loco[[chr_str]],
                addcovar = addcovar,
                intcovar = intcovar,
                query_func = query_variants,
                chr = chr_str,
                start = pos - window,
                end = pos + window,
                keep_all_snps = TRUE,
                cores = ncores
            )

            if (is.null(snp_assoc) || is.null(snp_assoc$lod)) {
                warning("SNP Association: scan1snps returned NULL or missing LOD scores")
                return(list(error = "SNP association failed - no LOD scores returned. Check chromosome naming and data availability."))
            }

            message(sprintf("SNP association complete: %d SNPs analyzed", nrow(snp_assoc$snpinfo)))

            # Add genes to the result
            snp_assoc$genes <- genes

            return(snp_assoc)
        },
        error = function(e) {
            warning(paste("Error in SNP association analysis:", e$message))
            message(sprintf("Full error details: %s", toString(e)))
            return(list(error = paste0("Analysis error: ", e$message)))
        }
    )
}

#' Get gene information for SNP association plot
#'
#' @param chr Chromosome
#' @param pos Position (Mb)
#' @param window Window size (Mb)
#' @return Gene information data frame or NULL
#' @importFrom qtl2 create_gene_query_func
#' @export
get_genes_for_snp_plot <- function(chr, pos, window = 1.5) {
    data_dir <- "/data/dev/DO_mapping_files"
    dbfile <- file.path(data_dir, "variants_db_by_chr", paste0("founder_variants_chr", chr, ".sqlite"))

    if (!file.exists(dbfile)) {
        warning(paste("Variant database file not found:", dbfile))
        return(NULL)
    }

    tryCatch(
        {
            query_genes <- qtl2::create_gene_query_func(
                dbfile,
                table_name = "genes",
                chr_field = "gene_chr",
                start_field = "gene_start",
                stop_field = "gene_end",
                name_field = "gene_symbol",
                strand_field = "gene_strand"
            )

            genes <- query_genes(chr, pos - window, pos + window)
            return(genes)
        },
        error = function(e) {
            warning(paste("Error querying genes:", e$message))
            return(NULL)
        }
    )
}
