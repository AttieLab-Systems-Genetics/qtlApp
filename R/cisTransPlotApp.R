# R/cisTransPlotApp.R

#' Cis-Trans Plot Module UI components
#'
#' @param id Module ID
#' @export
library(ggplot2)
library(dplyr)

#' @param id Module ID
#' @export
cisTransPlotInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    uiOutput(ns("dataset_selector_ui")) # Placeholder for dynamic UI
  )
}

#' Cis-Trans Plot Module Output UI
#'
#' @param id Module ID
#' @export
cisTransPlotUI <- function(id) {
  ns <- shiny::NS(id)
  shinycssloaders::withSpinner(plotly::plotlyOutput(ns("cis_trans_plot_output"), height = "700px"))
}

#' Cis-Trans Plot Module Server
#'
#' @param id Module ID
#' @param import_reactives Reactive list containing file_directory and annotation_list
#' @param peaks_cache Environment for caching peak finder results.
#' @export
cisTransPlotServer <- function(id, import_reactives, peaks_cache) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive for generating dataset choices for the selectInput
    cistrans_dataset_choices <- shiny::reactive({
      shiny::req(import_reactives(), import_reactives()$file_directory)
      file_dir <- import_reactives()$file_directory
      unique_groups <- unique(file_dir$group)
      
      valid_datasets <- character(0)
      for (grp in unique_groups) {
        group_subset <- file_dir[file_dir$group == grp, , drop = FALSE]
        if (nrow(group_subset) > 0) {
          # Create a temporary import_data structure for get_trait_type
          temp_import_data <- list(file_directory = group_subset)
          type <- get_trait_type(temp_import_data, grp) 
          if (!is.null(type) && type %in% c("genes", "isoforms")) {
            valid_datasets <- c(valid_datasets, grp)
          }
        }
      }
      sort(unique(valid_datasets))
    })
    
    output$dataset_selector_ui <- renderUI({
        ns <- session$ns
        choices <- cistrans_dataset_choices()
        if (length(choices) > 0) {
            create_select_input(ns("selected_cistrans_dataset"), 
                                label = "Select Dataset for Cis/Trans Plot:", 
                                choices = choices, 
                                selected = choices[1])
        } else {
            shiny::p("No suitable (gene/isoform) datasets found for Cis/Trans plot.")
        }
    })

    selected_dataset <- shiny::reactive({
      shiny::req(input$selected_cistrans_dataset)
      message("cisTransPlotServer: Selected Dataset: ", input$selected_cistrans_dataset) # DEBUG
      input$selected_cistrans_dataset
    })
    trait_type <- shiny::reactive({
      shiny::req(import_reactives(), selected_dataset())
      current_trait_type <- get_trait_type(import_reactives(), selected_dataset())
      message("cisTransPlotServer: Determined Trait Type: ", current_trait_type) # DEBUG
      current_trait_type
    })

    
    peaks_data <- shiny::reactive({
      shiny::req(import_reactives(), selected_dataset(), trait_type())
      message("cisTransPlotServer: Calling peak_finder with dataset: ", selected_dataset(), ", trait_type: ", trait_type()) # DEBUG
      
      
      current_peaks_cache <- if (is.environment(peaks_cache)) peaks_cache else NULL
      
      found_peaks <- peak_finder(import_reactives()$file_directory, 
                                 selected_dataset(), 
                                 trait_type = trait_type(), 
                                 cache_env = current_peaks_cache, 
                                 use_cache = TRUE)
      message("cisTransPlotServer: Data returned by peak_finder:") # DEBUG
      if(is.data.frame(found_peaks)){
        message("peak_finder returned a data.frame with ", nrow(found_peaks), " rows and columns: ", paste(colnames(found_peaks), collapse=", "))
        if(nrow(found_peaks) > 0){
            message("Head of peak_finder output:")
            print(utils::head(found_peaks))
        } else {
            message("peak_finder returned an empty data.frame.")
        }
      } else {
        message("peak_finder did not return a data.frame. Object type: ", class(found_peaks))
        print(found_peaks)
      }
      found_peaks
    })

 
    plot_data <- shiny::reactive({
      df <- peaks_data()
      message("cisTransPlotServer: Preparing plot_data. Initial df from peaks_data() has ", if(is.null(df)) "NULL" else nrow(df), " rows.") # DEBUG
      if (is.null(df) || nrow(df) == 0) {
        message("cisTransPlotServer: plot_data - peaks_data() is NULL or empty, returning NULL.") # DEBUG
        return(NULL)
      }
     
      req_cols <- c("qtl_chr", "qtl_pos", "cis", "gene_chr", "gene_start")
      message("cisTransPlotServer: Required columns for plot: ", paste(req_cols, collapse=", ")) # DEBUG
      message("cisTransPlotServer: Columns in df from peaks_data(): ", paste(colnames(df), collapse=", ")) # DEBUG
      
      if (!all(req_cols %in% colnames(df))) {
        missing_cols <- req_cols[!req_cols %in% colnames(df)]
        message("cisTransPlotServer: plot_data - Missing required columns: ", paste(missing_cols, collapse=", "), ". Returning NULL.") # DEBUG
        return(NULL)
      }
    
      df_filtered <- dplyr::filter(df, !is.na(qtl_chr) & !is.na(qtl_pos) & !is.na(cis) & !is.na(gene_chr) & !is.na(gene_start))
      message("cisTransPlotServer: plot_data - After filtering NAs, df_filtered has ", nrow(df_filtered), " rows.") # DEBUG
      
      if(nrow(df_filtered) == 0){
          message("cisTransPlotServer: plot_data - After NA filtering, data frame is empty. Returning NULL.") # DEBUG
          return(NULL)
      }
      
     
      chrom_levels <- c(as.character(1:19), "X")
      df_filtered$qtl_chr <- factor(chr_XYM(df_filtered$qtl_chr), levels = chrom_levels, ordered = TRUE)
      df_filtered$gene_chr <- factor(chr_XYM(df_filtered$gene_chr), levels = chrom_levels, ordered = TRUE)
      # Ensure cis is logical or character
      if (is.logical(df_filtered$cis)) {
        df_filtered$cis <- as.character(df_filtered$cis)
      } else if (!is.character(df_filtered$cis)) {
        df_filtered$cis <- as.character(as.logical(df_filtered$cis))
      }
      df_filtered
    })

    output$cis_trans_plot_output <- plotly::renderPlotly({
      df <- plot_data()
      if (is.null(df) || nrow(df) == 0) {
        return(plotly::ggplotly(plot_null("No cis/trans data available for this dataset.")))
      }
      
      
      df$hover_text <- paste(
        "Gene ID:", df$gene_id,
        "<br>Symbol:", df$gene_symbol,
        "<br>QTL Pos:", round(df$qtl_pos, 2),
        "<br>Marker:", df$marker,
        "<br>A:", round(df$A, 3), " B:", round(df$B, 3),
        "<br>C:", round(df$C, 3), " D:", round(df$D, 3),
        "<br>E:", round(df$E, 3), " F:", round(df$F, 3),
        "<br>G:", round(df$G, 3), " H:", round(df$H, 3)
      )
      
      cis.colors <- c("FALSE" = "#E41A1C", "TRUE" = "blue")
      p <- ggplot2::ggplot(df, ggplot2::aes(x = qtl_pos, y = gene_start, color = cis, text = hover_text)) +
        ggplot2::geom_point(alpha = 0.5, size = 1.2) +
        ggplot2::scale_color_manual(values = cis.colors, labels = c("FALSE" = "Trans", "TRUE" = "Cis"), drop = FALSE) +
        ggplot2::labs(
          x = NULL,
          y = NULL,
          color = "Cis QTL"
        ) +
        theme_minimal() +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = "grey70", fill = NA),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          panel.spacing = ggplot2::unit(0.1, "lines")
        )
      plotly::ggplotly(p, tooltip = "text") %>%
        plotly::layout(dragmode = "pan", hovermode = "closest") %>%
        plotly::config(displaylogo = FALSE, modeBarButtonsToRemove = c("select2d", "lasso2d", "hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines", "sendDataToCloud"))
    })
  })
}
