# R/cisTransPlotApp.R

#' Cis-Trans Plot Module UI components
#'
#' @param id Module ID
#' @export
library(ggplot2)
library(dplyr)
library(stringr)

# Source the helpers file to make get_qtlxcovar_file_paths available
if (file.exists("R/helpers.R")) {
  source("R/helpers.R")
}

#' @param id Module ID
#' @export
cisTransPlotInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList()
}

#' Cis-Trans Plot Module Output UI
#'
#' @param id Module ID
#' @export
cisTransPlotUI <- function(id) {
  ns <- shiny::NS(id)
  shinycssloaders::withSpinner(
    plotly::plotlyOutput(
      ns("cis_trans_plot_output"),
      height = "100%",
      width = "100%"
    )
  )
}

#' Helper function to get paths to qtlxcovar files for cis-trans plot
#' @param base_dataset The base name of the dataset (e.g., "HC_HF Liver Genes, additive")
#' @param interaction_type The type of interaction ("sex" or "diet")
#' @return A character vector of file paths, or NULL if not applicable.
get_qtlxcovar_file_paths <- function(base_dataset, interaction_type) {
  qtlxcovar_dir <- "/data/dev/miniViewer_3.0"

  file_prefix <- NULL
  if (grepl("Liver.*Genes", base_dataset, ignore.case = TRUE)) {
    file_prefix <- "DO1200_liver_genes_all_mice"
  } else if (grepl("Clinical", base_dataset, ignore.case = TRUE)) {
    file_prefix <- "DO1200_clinical_traits_all_mice"
  } else {
    return(NULL) # Not a dataset type supported by this plot's qtlxcovar files
  }

  file_suffix <- if (interaction_type == "sex") "qtlxsex_peaks.csv" else "qtlxdiet_peaks.csv"

  # Construct the full file path
  full_path <- file.path(qtlxcovar_dir, paste0(file_prefix, "_", file_suffix))

  # Return the path only if the file exists, otherwise return NULL to prevent errors
  if (file.exists(full_path)) {
    return(full_path)
  } else {
    warning(paste("qtlxcovar file not found at path:", full_path))
    return(NULL)
  }
}


#' Cis-Trans Plot Module Server
#'
#' @param id Module ID
#' @param import_reactives Reactive list containing file_directory and annotation_list
#' @param main_par A REACTIVE that returns a list, expected to contain selected_dataset (which is a reactive group name)
#' @param peaks_cache Environment for caching peak finder results.
#' @export
cisTransPlotServer <- function(id, import_reactives, main_par, peaks_cache, sidebar_interaction_type = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    selected_dataset <- shiny::reactive({
      main_par_val <- main_par()
      shiny::req(main_par_val, main_par_val$selected_dataset, main_par_val$selected_dataset())
      main_par_val$selected_dataset()
    })

    trait_type <- shiny::reactive({
      get_trait_type(import_reactives(), selected_dataset())
    })

    peaks_data <- shiny::reactive({
      shiny::req(import_reactives(), selected_dataset(), trait_type())

      interaction_type <- if (is.reactive(sidebar_interaction_type)) sidebar_interaction_type() else "none"
      base_dataset <- selected_dataset()
      use_qtlxcovar <- !is.null(interaction_type) && interaction_type != "none" && grepl("^HC_HF", base_dataset, ignore.case = TRUE)

      found_peaks <- NULL
      if (use_qtlxcovar) {
        qtlxcovar_files <- get_qtlxcovar_file_paths(base_dataset, interaction_type)
        if (!is.null(qtlxcovar_files) && length(qtlxcovar_files) > 0) {
          found_peaks <- tryCatch(
            {
              peak_data_list <- lapply(qtlxcovar_files, data.table::fread, stringsAsFactors = FALSE)
              data.table::rbindlist(peak_data_list, use.names = TRUE, fill = TRUE)
            },
            error = function(e) {
              warning(paste("Error reading qtlxcovar files:", e$message))
              NULL
            }
          )
        }
      } else {
        found_peaks <- peak_finder(
          file_dir = import_reactives()$file_directory,
          selected_dataset = selected_dataset(),
          trait_type = trait_type(),
          cache_env = peaks_cache,
          use_cache = TRUE
        )
      }

      if (use_qtlxcovar && is.null(found_peaks)) {
        return(NULL)
      }
      found_peaks
    })

    plot_data <- shiny::reactive({
      df <- peaks_data()
      shiny::req(main_par(), main_par()$LOD_thr)
      lod_threshold <- main_par()$LOD_thr()
      interaction_type <- if (is.reactive(sidebar_interaction_type)) sidebar_interaction_type() else "none"

      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }

      is_qtlxcovar_data <- "lod_diff" %in% colnames(df)
      use_lod_diff <- is_qtlxcovar_data && !is.null(interaction_type) && interaction_type != "none"
      lod_column <- if (use_lod_diff) "lod_diff" else "qtl_lod"

      req_cols <- c("qtl_chr", "qtl_pos", "cis", "gene_chr", "gene_start", lod_column)
      if (!all(req_cols %in% colnames(df))) {
        return(NULL)
      }

      # Using data.table for efficiency
      plot_data_dt <- data.table::as.data.table(df)

      # Filter NAs based on required columns
      plot_data_dt <- na.omit(plot_data_dt, cols = req_cols)

      if (use_lod_diff) {
        plot_data_dt <- plot_data_dt[abs(get(lod_column)) >= lod_threshold]
        plot_data_dt[, qtl_lod := abs(get(lod_column))]
      } else {
        plot_data_dt <- plot_data_dt[get(lod_column) >= lod_threshold]
      }

      if (nrow(plot_data_dt) == 0) {
        return(NULL)
      }

      markers_data <- shiny::req(import_reactives()$markers)
      interaction_type <- if (is.reactive(sidebar_interaction_type)) sidebar_interaction_type() else "none"

      chr_summary <- as.data.table(markers_data)[, .(chr_len = max(pos)), by = chr][
        , chr_num := chr_to_numeric(chr)
      ][
        order(chr_num)
      ][
        , tot := cumsum(as.numeric(chr_len)) - chr_len
      ][
        , center := tot + (chr_len / 2)
      ]

      plot_data_dt[, qtl_chr_char := chr_XYM(qtl_chr)]
      plot_data_dt[, gene_chr_char := chr_XYM(gene_chr)]

      plot_data_dt[chr_summary, on = .(qtl_chr_char = chr), qtl_BPcum := qtl_pos + i.tot]
      plot_data_dt[chr_summary, on = .(gene_chr_char = chr), gene_BPcum := gene_start + i.tot]

      list(plot_data = plot_data_dt, chr_summary = chr_summary, interaction_type = interaction_type)
    })

    clicked_phenotype <- shiny::reactiveVal(NULL)

    output$cis_trans_plot_output <- plotly::renderPlotly({
      plot_data_list <- plot_data()
      if (is.null(plot_data_list) || nrow(plot_data_list$plot_data) == 0) {
        return(plot_null("No significant peaks found above the LOD threshold."))
      }

      plot_data_dt <- plot_data_list$plot_data
      chr_summary <- plot_data_list$chr_summary
      interaction_type <- plot_data_list$interaction_type

      # Create a robust gene label for display/clicks: prefer symbol, else gene_id, else trait
      if ("gene_id" %in% names(plot_data_dt)) {
        plot_data_dt[, gene_label := ifelse(!is.na(gene_symbol) & gene_symbol != "" & gene_symbol != "N/A", gene_symbol, gene_id)]
      } else if ("trait" %in% names(plot_data_dt)) {
        plot_data_dt[, gene_label := ifelse(!is.na(gene_symbol) & gene_symbol != "" & gene_symbol != "N/A", gene_symbol, trait)]
      } else {
        plot_data_dt[, gene_label := gene_symbol]
      }

      # Vectorized conditional creation of the difference text part for the hover text
      diff_text <- if ("lod_diff" %in% names(plot_data_dt)) {
        ifelse(!is.na(plot_data_dt$lod_diff), paste0(" (Diff: ", round(plot_data_dt$lod_diff, 2), ")"), "")
      } else {
        ""
      }

      plot_data_dt[, hover_text := paste0(
        "Gene: ", gene_label, "<br>",
        "LOD: ", round(qtl_lod, 2),
        diff_text,
        "<br>Gene Pos: ", gene_chr_char, ":", round(gene_start, 2), " Mb",
        "<br>Marker Pos: ", qtl_chr_char, ":", round(qtl_pos, 2), " Mb"
      )]

      plot_title <- if (!is.null(interaction_type) && interaction_type != "none") {
        paste("Cis-Trans QTL Plot -", stringr::str_to_title(interaction_type), "Interaction Differences")
      } else {
        "Cis-Trans QTL Plot"
      }

      g <- ggplot2::ggplot(plot_data_dt, ggplot2::aes(
        x = qtl_BPcum, y = gene_BPcum,
        color = as.character(cis),
        text = hover_text,
        customdata = gene_label
      )) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::scale_color_manual(
          values = c("TRUE" = "blue", "FALSE" = "red"),
          name = "QTL Type",
          labels = c("TRUE" = "Cis", "FALSE" = "Trans")
        ) +
        ggplot2::labs(
          title = plot_title,
          x = "QTL Chromosome",
          y = "Gene Chromosome"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
        ggplot2::geom_vline(xintercept = chr_summary$tot, linetype = "dotted", color = "grey") +
        ggplot2::geom_hline(yintercept = chr_summary$tot, linetype = "dotted", color = "grey") +
        ggplot2::scale_x_continuous(labels = chr_summary$chr, breaks = chr_summary$center) +
        ggplot2::scale_y_continuous(labels = chr_summary$chr, breaks = chr_summary$center)

      plotly::ggplotly(g, tooltip = "text", source = ns("cistrans_plot")) %>%
        plotly::layout(dragmode = "pan") %>%
        plotly::config(scrollZoom = TRUE, displaylogo = FALSE, modeBarButtonsToRemove = c("select2d", "lasso2d")) %>%
        plotly::event_register("plotly_click")
    })

    shiny::observeEvent(plotly::event_data("plotly_click", source = ns("cistrans_plot")), {
      event_data <- plotly::event_data("plotly_click", source = ns("cistrans_plot"))
      if (!is.null(event_data$customdata)) {
        clicked_phenotype(event_data$customdata)
      }
    })

    return(
      list(
        clicked_phenotype_for_lod_scan = clicked_phenotype
      )
    )
  })
}
