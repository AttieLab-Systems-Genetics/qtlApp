#!/usr/bin/env Rscript

# Helper utilities for on-the-fly liver isoform self-correlations.
# These functions are used by correlationApp.R when the user requests
# liver_isoforms vs liver_isoforms correlations. Instead of reading a
# huge precomputed CSV (~44 GB), we compute one-vs-all correlations
# on demand from the underlying phenotype matrix stored in the DO
# cross object.

isoform_corr_env <- new.env(parent = emptyenv())

# Load and cache the liver isoform phenotype matrix from the DO cross.
# This is equivalent to the data used by the batch correlation script,
# restricted to dataset "liver_isoforms".
load_isoform_pheno_matrix <- function() {
  if (!is.null(isoform_corr_env$pheno_mat) &&
    is.matrix(isoform_corr_env$pheno_mat)) {
    return(invisible(NULL))
  }

  cross_path <- "/data/dev/DO_mapping_files/cross_DO1200_grcm39.rds"
  if (!file.exists(cross_path)) {
    warning("load_isoform_pheno_matrix: cross file not found at ", cross_path)
    return(invisible(NULL))
  }

  message("load_isoform_pheno_matrix: Loading cross object from ", cross_path)
  cross <- readRDS(cross_path)

  if (is.null(cross$phenocovar) || is.null(cross$pheno)) {
    warning("load_isoform_pheno_matrix: cross object missing phenocovar or pheno.")
    return(invisible(NULL))
  }

  traits_iso <- subset(cross$phenocovar,
    dataset %in% "liver_isoforms",
    phenotype,
    drop = TRUE
  )

  if (length(traits_iso) == 0) {
    warning("load_isoform_pheno_matrix: No traits found for dataset 'liver_isoforms' in cross$phenocovar.")
    return(invisible(NULL))
  }

  pheno_mat <- cross$pheno[, traits_iso, drop = FALSE]
  if (!is.matrix(pheno_mat)) {
    pheno_mat <- as.matrix(pheno_mat)
  }
  mode(pheno_mat) <- "numeric"

  isoform_corr_env$pheno_mat <- pheno_mat
  isoform_corr_env$traits_iso <- traits_iso
  isoform_corr_env$covar <- cross$covar

  dims <- dim(pheno_mat)
  message(
    "load_isoform_pheno_matrix: Cached isoform phenotype matrix with ",
    dims[1], " samples x ", dims[2], " traits."
  )

  invisible(NULL)
}

# Ensure that a covariate-adjusted isoform phenotype matrix exists in the cache.
# This version hard-codes additive covariates ~GenLit+Sex, matching the app's
# desired behavior when the "Covariate-Adjusted" toggle is enabled.
ensure_isoform_adjusted_matrix <- function() {
  if (!is.null(isoform_corr_env$pheno_mat_adj) &&
    is.matrix(isoform_corr_env$pheno_mat_adj)) {
    return(invisible(NULL))
  }

  pheno_mat <- isoform_corr_env$pheno_mat
  covar <- isoform_corr_env$covar
  if (is.null(pheno_mat) || !is.matrix(pheno_mat) ||
    is.null(covar) || !"GenLit" %in% colnames(covar) || !"Sex" %in% colnames(covar)) {
    warning("ensure_isoform_adjusted_matrix: Missing phenotype or covariate data; cannot compute adjusted matrix.")
    return(invisible(NULL))
  }

  GenLit <- as.factor(covar[, "GenLit"])
  Sex <- as.factor(covar[, "Sex"])

  message("ensure_isoform_adjusted_matrix: Computing covariate-adjusted isoform matrix with ~GenLit+Sex.")
  pheno_adj <- pheno_mat

  for (j in seq_len(ncol(pheno_mat))) {
    y <- pheno_mat[, j]
    # Skip columns that are all NA or constant
    if (all(is.na(y)) || length(unique(y[!is.na(y)])) <= 1) {
      next
    }
    fit <- stats::lm(y ~ GenLit + Sex, na.action = stats::na.exclude)
    pheno_adj[, j] <- stats::resid(fit)
  }

  isoform_corr_env$pheno_mat_adj <- pheno_adj
  invisible(NULL)
}

# Resolve the matrix column name corresponding to a UI-selected isoform trait.
# The UI typically provides an isoform symbol (e.g. "0610005C13Rik-201").
# We map:
#   symbol -> transcript_id -> matrix column ("liver_<transcript_id>")
resolve_isoform_column_from_trait <- function(trait_string, import) {
  if (is.null(trait_string) || !nzchar(trait_string)) {
    return(NULL)
  }

  load_isoform_pheno_matrix()
  pheno_mat <- isoform_corr_env$pheno_mat
  if (is.null(pheno_mat) || !is.matrix(pheno_mat)) {
    return(NULL)
  }

  # If the trait string already matches a matrix column, use it directly.
  if (trait_string %in% colnames(pheno_mat)) {
    return(trait_string)
  }
  if (paste0("liver_", trait_string) %in% colnames(pheno_mat)) {
    return(paste0("liver_", trait_string))
  }

  # Use annotation_list$isoforms to map symbol -> transcript_id.
  ann <- if (!is.null(import) && !is.null(import$annotation_list)) {
    import$annotation_list
  } else {
    NULL
  }

  if (is.null(ann) || is.null(ann$isoforms)) {
    return(NULL)
  }
  iso_dt <- ann$isoforms

  get_col_chr <- function(df, nms) {
    for (nm in nms) {
      if (nm %in% colnames(df)) {
        return(as.character(df[[nm]]))
      }
    }
    return(NULL)
  }

  sym_vec <- get_col_chr(iso_dt, c("symbol", "gene.symbol"))
  tid_vec <- get_col_chr(iso_dt, c("transcript_id", "transcript.id"))

  trait_clean <- trimws(trait_string)
  trait_lower <- tolower(trait_clean)

  # Try matching by symbol first
  iso_tid <- NULL
  if (!is.null(sym_vec) && !is.null(tid_vec)) {
    idx_sym <- which(tolower(sym_vec) == trait_lower)
    if (length(idx_sym) >= 1) {
      cand <- tid_vec[idx_sym[1]]
      if (!is.na(cand) && nzchar(cand)) {
        iso_tid <- cand
      }
    }
  }

  # If that failed, try matching directly against transcript_id-like values
  if (is.null(iso_tid) && !is.null(tid_vec)) {
    tid_lower <- tolower(tid_vec)
    stripped <- sub("^liver_", "", trait_lower)
    idx_tid <- which(tid_lower == trait_lower | tid_lower == stripped)
    if (length(idx_tid) >= 1) {
      cand <- tid_vec[idx_tid[1]]
      if (!is.na(cand) && nzchar(cand)) {
        iso_tid <- cand
      }
    }
  }

  if (is.null(iso_tid) || !nzchar(iso_tid)) {
    return(NULL)
  }

  # Matrix columns are of the form "liver_<transcript_id>"
  col1 <- iso_tid
  col2 <- paste0("liver_", iso_tid)

  if (col1 %in% colnames(pheno_mat)) {
    return(col1)
  }
  if (col2 %in% colnames(pheno_mat)) {
    return(col2)
  }

  NULL
}

# Compute top-N absolute correlations for a selected isoform trait vs all
# isoforms, using the cached phenotype matrix. Returns a data.frame/table
# with columns:
#   trait, correlation_value, p_value, num_mice
compute_isoform_cor_top_n <- function(trait_string, import, top_n = 500, use_adjusted = FALSE) {
  if (is.null(trait_string) || !nzchar(trait_string)) {
    message("compute_isoform_cor_top_n: Empty trait_string; returning empty result.")
    return(data.frame(
      trait = character(0),
      correlation_value = numeric(0),
      p_value = numeric(0),
      num_mice = numeric(0)
    ))
  }

  load_isoform_pheno_matrix()

  # Choose unadjusted vs covariate-adjusted matrix
  if (isTRUE(use_adjusted)) {
    ensure_isoform_adjusted_matrix()
    pheno_mat <- isoform_corr_env$pheno_mat_adj
  } else {
    pheno_mat <- isoform_corr_env$pheno_mat
  }
  if (is.null(pheno_mat) || !is.matrix(pheno_mat)) {
    warning("compute_isoform_cor_top_n: isoform phenotype matrix not available.")
    return(data.frame(
      trait = character(0),
      correlation_value = numeric(0),
      p_value = numeric(0),
      num_mice = numeric(0)
    ))
  }

  col_name <- resolve_isoform_column_from_trait(trait_string, import)
  if (is.null(col_name)) {
    warning("compute_isoform_cor_top_n: Could not resolve matrix column for trait '", trait_string, "'.")
    return(data.frame(
      trait = character(0),
      correlation_value = numeric(0),
      p_value = numeric(0),
      num_mice = numeric(0)
    ))
  }

  if (!col_name %in% colnames(pheno_mat)) {
    warning("compute_isoform_cor_top_n: Column '", col_name, "' not found in phenotype matrix.")
    return(data.frame(
      trait = character(0),
      correlation_value = numeric(0),
      p_value = numeric(0),
      num_mice = numeric(0)
    ))
  }

  y <- pheno_mat[, col_name]

  # One-vs-all correlation
  message("compute_isoform_cor_top_n: Computing one-vs-all correlations for column '", col_name, "'.")
  r_vec <- stats::cor(pheno_mat, y, use = "pairwise.complete.obs")
  abs_r <- abs(r_vec)

  # Ensure correlation vector has names; fall back to phenotype column names
  if (is.null(names(r_vec))) {
    names(r_vec) <- colnames(pheno_mat)
  }

  # Compute per-partner sample sizes (number of mice with non-missing values)
  # for each correlation. This matches the logic used in the batch script,
  # but simplified for the one-vs-all case.
  valid_y <- !is.na(y)
  if (!any(valid_y)) {
    warning("compute_isoform_cor_top_n: No non-missing values for trait '", trait_string, "'.")
    return(data.frame(
      trait = character(0),
      correlation_value = numeric(0),
      p_value = numeric(0),
      num_mice = numeric(0)
    ))
  }
  # Restrict to rows where y is non-missing; then count non-missing per column
  valid_mat <- !is.na(pheno_mat[valid_y, , drop = FALSE])
  n_obs_vec <- colSums(valid_mat)

  ord <- order(abs_r, decreasing = TRUE)
  if (length(ord) > top_n) {
    ord <- ord[seq_len(top_n)]
  }

  if (length(ord) == 0L) {
    warning("compute_isoform_cor_top_n: No finite correlations found for trait '", trait_string, "'.")
    return(data.frame(
      trait = character(0),
      correlation_value = numeric(0),
      p_value = numeric(0),
      num_mice = numeric(0)
    ))
  }

  partners <- names(r_vec)[ord]
  cor_vals <- as.numeric(r_vec[ord])
  n_vals <- as.numeric(n_obs_vec[ord])

  # Vectorized p-value computation (same formula as in correlation script)
  p_vals <- rep(NA_real_, length(cor_vals))
  valid_p <- !is.na(cor_vals) & !is.na(n_vals) & n_vals >= 3 & abs(cor_vals) < 1
  if (any(valid_p)) {
    r_ok <- cor_vals[valid_p]
    n_ok <- n_vals[valid_p]
    t_stat <- r_ok * sqrt(n_ok - 2) / sqrt(1 - r_ok^2)
    p_vals[valid_p] <- 2 * stats::pt(-abs(t_stat), df = n_ok - 2)
  }

  # Map partner matrix column names back to isoform symbols when possible
  ann <- if (!is.null(import) && !is.null(import$annotation_list)) {
    import$annotation_list
  } else {
    NULL
  }
  display_traits <- partners

  if (!is.null(ann) && !is.null(ann$isoforms)) {
    iso_dt <- ann$isoforms

    get_col_chr <- function(df, nms) {
      for (nm in nms) {
        if (nm %in% colnames(df)) {
          return(as.character(df[[nm]]))
        }
      }
      return(NULL)
    }

    sym_vec <- get_col_chr(iso_dt, c("symbol", "gene.symbol"))
    tid_vec <- get_col_chr(iso_dt, c("transcript_id", "transcript.id"))

    if (!is.null(sym_vec) && !is.null(tid_vec)) {
      partner_tid <- sub("^liver_", "", partners)
      # Build a named lookup: transcript_id -> symbol
      sym_lookup <- stats::setNames(sym_vec, tid_vec)
      sym_match <- sym_lookup[partner_tid]
      # Prefer symbol when available; otherwise use transcript_id
      display_traits <- ifelse(!is.na(sym_match) & nzchar(sym_match), sym_match, partner_tid)
    }
  }

  # Defensive logging before constructing the data frame
  message(
    "compute_isoform_cor_top_n: Building result data.frame with lengths -> ",
    "trait: ", length(display_traits),
    ", correlation_value: ", length(cor_vals)
  )

  n_rows <- length(cor_vals)
  if (length(display_traits) != n_rows ||
    length(n_vals) != n_rows ||
    length(p_vals) != n_rows) {
    warning(
      "compute_isoform_cor_top_n: Length mismatch; trait len = ",
      length(display_traits), ", cor len = ", n_rows,
      ", n_vals len = ", length(n_vals),
      ", p_vals len = ", length(p_vals),
      ". Truncating to the smaller length for safety."
    )
    n_rows <- min(length(display_traits), n_rows, length(n_vals), length(p_vals))
    display_traits <- display_traits[seq_len(n_rows)]
    cor_vals <- cor_vals[seq_len(n_rows)]
    n_vals <- n_vals[seq_len(n_rows)]
    p_vals <- p_vals[seq_len(n_rows)]
  }

  data.frame(
    trait = display_traits,
    correlation_value = cor_vals,
    p_value = p_vals,
    num_mice = n_vals,
    stringsAsFactors = FALSE
  )
}
