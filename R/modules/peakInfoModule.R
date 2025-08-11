#' Peak Info Module
#'
#' Renders a detailed peak information panel including marker, position,
#' LOD, cis/trans status, confidence interval, and founder allele effects.
#'
#' @param id Module ID
#' @param selected_peak_reactive Reactive returning a data.frame with peak info
#'
#' @importFrom shiny moduleServer NS renderUI uiOutput req
#' @importFrom htmltools tags tagList
#' @export
peakInfoUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("peak_info_display"))
}

#' @rdname peakInfoUI
#' @export
peakInfoServer <- function(id, selected_peak_reactive) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$peak_info_display <- shiny::renderUI({
      peak_info <- selected_peak_reactive()
      # Debug: log structure
      if (is.null(peak_info)) {
        message("peakInfoModule: selected_peak is NULL")
      } else {
        message("peakInfoModule: selected_peak type=", paste(class(peak_info), collapse=","))
        if (is.data.frame(peak_info)) {
          message("peakInfoModule: selected_peak nrow=", nrow(peak_info), " cols=", paste(colnames(peak_info), collapse=","))
        }
      }
      if (is.null(peak_info) || nrow(peak_info) == 0) {
        return(htmltools::tags$div("No peak selected.", style = "color: #7f8c8d; text-align: center; padding-top: 20px;"))
      }

      # Ensure 1-row frame for simple extraction
      if (nrow(peak_info) > 1) {
        peak_info <- peak_info[1, , drop = FALSE]
      }
      # Ensure we are working with a 1-row data.frame
      if (inherits(peak_info, "data.table")) {
        peak_info <- as.data.frame(peak_info)
      }
      if (!is.data.frame(peak_info)) {
        # Handle named vector/list cases defensively
        if (is.atomic(peak_info) && !is.null(names(peak_info))) {
          peak_info <- as.data.frame(as.list(peak_info), stringsAsFactors = FALSE)
        } else if (is.list(peak_info)) {
          peak_info <- as.data.frame(peak_info, stringsAsFactors = FALSE)
        }
      }
      if (nrow(peak_info) == 0) {
        return(htmltools::tags$div("No peak selected.", style = "color: #7f8c8d; text-align: center; padding-top: 20px;"))
      }

      # Helper to get the first available non-NA value from candidate columns
      get_field <- function(cols) {
        for (c in cols) {
          if (c %in% colnames(peak_info)) {
            v <- peak_info[[c]][1]
            if (!is.null(v) && !is.na(v)) return(v)
          }
        }
        return(NA)
      }

      marker <- get_field(c("marker", "markers", "snp", "id", "marker_name"))
      chr_raw <- get_field(c("qtl_chr", "chr", "qtl_chr_char", "chromosome"))
      pos <- get_field(c("qtl_pos", "pos", "position", "mb", "bp", "bp_grcm39"))
      lod <- get_field(c("qtl_lod", "lod", "LOD", "lod_diff"))
      ci_lo <- get_field(c("qtl_ci_lo", "ci_lo"))
      ci_hi <- get_field(c("qtl_ci_hi", "ci_hi"))
      trait <- get_field(c("trait", "phenotype"))
      gene_symbol <- get_field(c("gene_symbol", "symbol", "gene"))
      gene_id <- get_field(c("gene_id", "ensembl_id", "geneid"))

      # Normalize chromosome label
      chr_label <- NULL
      if (!is.null(chr_raw) && !is.na(chr_raw)) {
        if (is.numeric(chr_raw)) {
          # Use helpers' numeric_to_chr if available
          if (exists("numeric_to_chr", mode = "function")) {
            chr_label <- numeric_to_chr(chr_raw)
          } else {
            chr_label <- as.character(chr_raw)
          }
        } else {
          chr_label <- as.character(chr_raw)
        }
      }

      info_elements <- list()

      # Gene/Trait line with preference for gene symbol
      if (!is.null(gene_symbol) && !is.na(gene_symbol)) {
        # Show Gene symbol and optional Ensembl ID
        gene_line <- gene_symbol
        if (!is.null(gene_id) && !is.na(gene_id)) {
          gene_line <- paste0(gene_symbol, " (", gene_id, ")")
        }
        info_elements <- c(info_elements, list(
          htmltools::tags$strong("Gene: "), gene_line, htmltools::tags$br()
        ))
      } else if (!is.null(trait) && !is.na(trait)) {
        # Fallback to trait/phenotype string
        info_elements <- c(info_elements, list(
          htmltools::tags$strong("Trait: "), trait, htmltools::tags$br()
        ))
      }

      # Basic info
      if (!is.null(marker) && !is.na(marker)) {
        info_elements <- c(info_elements, list(htmltools::tags$strong("Marker: "), marker, htmltools::tags$br()))
      }
      if (!is.null(chr_label) && !is.na(chr_label) && !is.null(pos) && !is.na(pos)) {
        info_elements <- c(info_elements, list(
          htmltools::tags$strong("Position: "), paste0("Chr", chr_label, ":", round(as.numeric(pos), 2), " Mb"), htmltools::tags$br()
        ))
      }
      if (!is.null(lod) && !is.na(lod)) {
        # If this came from lod_diff, show both signed and |lod|
        if ("lod_diff" %in% colnames(peak_info) && !is.na(peak_info$lod_diff[1])) {
          info_elements <- c(info_elements, list(
            htmltools::tags$strong("LOD Difference: "), paste0(round(as.numeric(peak_info$lod_diff[1]), 3),
            " (|LOD| = ", round(abs(as.numeric(peak_info$lod_diff[1])), 3), ")"), htmltools::tags$br()
          ))
        } else {
          info_elements <- c(info_elements, list(
            htmltools::tags$strong("LOD Score: "), round(as.numeric(lod), 3), htmltools::tags$br()
          ))
        }
      }

      # Cis/Trans status
      if ("cis" %in% colnames(peak_info)) {
        cis_val <- peak_info[["cis"]][1]
        cis_status <- if (is.logical(cis_val)) {
          ifelse(isTRUE(cis_val), "Cis", "Trans")
        } else if (is.character(cis_val)) {
          ifelse(toupper(cis_val) %in% c("TRUE", "1", "YES", "CIS"), "Cis", "Trans")
        } else if (is.numeric(cis_val)) {
          ifelse(cis_val == 1, "Cis", "Trans")
        } else {
          "Unknown"
        }
        cis_color <- if (cis_status == "Cis") "#27ae60" else "#e74c3c"
        info_elements <- c(info_elements, list(
          htmltools::tags$strong("Type: "),
          htmltools::tags$span(cis_status, style = paste0("color: ", cis_color, "; font-weight: bold;")),
          htmltools::tags$br()
        ))
      }

      # Confidence interval
      if (!is.null(ci_lo) && !is.na(ci_lo) && !is.null(ci_hi) && !is.na(ci_hi)) {
        info_elements <- c(info_elements, list(
          htmltools::tags$strong("95% CI: "),
          paste0(round(as.numeric(ci_lo), 2), " - ", round(as.numeric(ci_hi), 2), " Mb"),
          htmltools::tags$br()
        ))
      }

      # Founder allele effects display
      allele_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
      strain_names <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
      available_alleles <- allele_cols[allele_cols %in% colnames(peak_info)]
      if (length(available_alleles) > 0) {
        allele_effects <- list()
        for (i in seq_along(available_alleles)) {
          col <- available_alleles[i]
          value <- peak_info[[col]][1]
          if (!is.null(value) && !is.na(value)) {
            strain <- strain_names[i]
            allele_effects[[length(allele_effects) + 1]] <- paste0(strain, ": ", round(as.numeric(value), 3))
          }
        }
        if (length(allele_effects) > 0) {
          info_elements <- c(info_elements, list(
            htmltools::tags$strong("Founder Effects:"), htmltools::tags$br(),
            htmltools::tags$div(
              style = "margin-left: 10px; font-family: monospace; font-size: 11px;",
              lapply(allele_effects, function(effect) {
                htmltools::tags$div(effect, style = "margin: 2px 0;")
              })
            )
          ))
        }
      }

      # If nothing was captured (all fields missing), show a compact key-value dump as fallback
      if (length(info_elements) == 0) {
        message("peakInfoModule: No standard fields found. Available columns:", paste(colnames(peak_info), collapse = ", "))
        # Build a compact table of key fields if present
        key_cols <- c("trait", "phenotype", "marker", "markers", "qtl_chr", "chr", "qtl_pos", "pos", "qtl_lod", "lod", "lod_diff",
                      "A","B","C","D","E","F","G","H")
        present <- intersect(key_cols, colnames(peak_info))
        if (length(present) > 0) {
          rows <- lapply(present, function(nm) {
            htmltools::tags$tr(
              htmltools::tags$td(htmltools::tags$strong(nm), style = "padding:2px 6px;"),
              htmltools::tags$td(as.character(peak_info[[nm]][1]), style = "padding:2px 6px;")
            )
          })
          return(htmltools::tags$table(
            style = "font-size:11px; width:100%;",
            do.call(htmltools::tagList, rows)
          ))
        } else {
          return(htmltools::tags$div(
            "No metadata available for this peak.",
            style = "color: #7f8c8d; text-align: center; padding: 6px;"
          ))
        }
      }

      do.call(htmltools::tagList, info_elements)
    })
  })
}


