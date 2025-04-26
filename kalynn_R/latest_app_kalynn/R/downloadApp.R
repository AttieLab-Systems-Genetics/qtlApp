#' Download Module UI for Allele Effects Plot
#'
#' Creates download buttons (PNG, PDF) for the allele effects plot.
#'
#' @param id Module ID.
#'
#' @importFrom shiny NS downloadButton tagList
#' @export
downloadAlleleUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::downloadButton(ns("download_effects_plot_png"), "Download PNG"),
    shiny::downloadButton(ns("download_effects_plot_pdf"), "Download PDF")
  )
}

#' Download Module UI for QTL Plot
#'
#' Creates download buttons (PNG, PDF) for the main QTL plot.
#'
#' @param id Module ID.
#'
#' @importFrom shiny NS downloadButton tagList
#' @export
downloadQtlUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
      shiny::downloadButton(ns("download_qtl_plot_png"), "Download PNG", class = "btn-sm"),
      shiny::downloadButton(ns("download_qtl_plot_pdf"), "Download PDF", class = "btn-sm")
  )
}

#' Download Module Server
#'
#' Handles the server-side logic for downloading plots.
#'
#' @param id Module ID.
#' @param main_par Reactive containing main parameters (`selected_chr`, `which_trait`).
#' @param plot_base_obj Reactive containing the base ggplot object for the QTL plot (`p`) and its data (`data`).
#' @param peak_mod_reactives Reactive containing outputs from the peak module (`allele_plot`, `selected_peak`).
#' @param official_gene_symbol Reactive containing the official trait symbol string.
#' @param plot_dims Reactive containing `width` and `height` for the QTL plot download.
#'
#' @importFrom shiny moduleServer NS downloadHandler req
#' @importFrom ggplot2 ggsave ggtitle geom_point aes
#' @importFrom dplyr filter arrange slice
#' @importFrom rlang .data
#' @export
downloadServer <- function(id, main_par, plot_base_obj, peak_mod_reactives, official_gene_symbol, plot_dims) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # == Allele Effects Plot Downloads ==

    output$download_effects_plot_png <- shiny::downloadHandler(
        filename = function() {
            shiny::req(peak_mod_reactives()$selected_peak, official_gene_symbol())
            actual_trait_name <- official_gene_symbol()
            peak_marker <- peak_mod_reactives()$selected_peak
            paste0("strain_effects_", actual_trait_name, "_", peak_marker, "_", format(Sys.time(), "%Y%m%d"), ".png")
        },
        content = function(file) {
            plot_obj <- peak_mod_reactives()$allele_plot
            shiny::req(plot_obj)
            ggplot2::ggsave(file, plot_obj, width = 10, height = 7, dpi = 300)
        }
    )

    output$download_effects_plot_pdf <- shiny::downloadHandler(
        filename = function() {
            shiny::req(peak_mod_reactives()$selected_peak, official_gene_symbol())
            actual_trait_name <- official_gene_symbol()
            peak_marker <- peak_mod_reactives()$selected_peak
            paste0("strain_effects_", actual_trait_name, "_", peak_marker, "_", format(Sys.time(), "%Y%m%d"), ".pdf")
        },
        content = function(file) {
            plot_obj <- peak_mod_reactives()$allele_plot
            shiny::req(plot_obj)
            # Ensure cairo_pdf device is available or handle error/alternative
             tryCatch({
                 ggplot2::ggsave(file, plot_obj, width = 10, height = 7, device = cairo_pdf)
             }, error = function(e){
                 warning("cairo_pdf device not available for PDF generation. Error: ", e$message)
                 # Fallback or stop?
                 # Fallback to standard pdf device if needed, though quality might differ
                 ggplot2::ggsave(file, plot_obj, width = 10, height = 7, device = "pdf")
             })
        }
    )

    # == QTL Plot Downloads ==
    
    # Helper function to generate the downloadable QTL plot (ggplot object)
    generate_qtl_download_plot <- reactive({ 
        shiny::req(plot_base_obj(), main_par(), peak_mod_reactives())
        
        plot_result <- plot_base_obj()
        p <- plot_result$p # The base ggplot object
        plot_data <- plot_result$data # Data used for the plot
            
        # Get peak information (already filtered/ordered in peak module if logic is sound)
        # Re-filter based on current LOD threshold just to be sure?
        peaks_full <- peak_mod_reactives()$peak_table
        peaks_info <- highest_peaks(peaks_full, main_par()$LOD_thr)
            
        # Filter peaks to only show those in the selected chromosome if applicable
        if(main_par()$selected_chr != "All") {
            chr_num <- switch(main_par()$selected_chr,
                "X" = 20, "Y" = 21, "M" = 22,
                as.numeric(main_par()$selected_chr)
            )
            if(!is.null(peaks_info)) {
               peaks_info <- peaks_info %>% dplyr::filter(.data$chr == chr_num)
            }
        }
            
        # Get the highest peak for annotation
        highest_peak_for_plot <- NULL
        if(!is.null(peaks_info) && nrow(peaks_info) > 0) {
            highest_peak_for_plot <- peaks_info %>% dplyr::slice(1)
                
            # Find the peak point in the plot data
            peak_point <- plot_data %>%
                dplyr::filter(.data$markers == highest_peak_for_plot$marker)
                
            # Add marker point to plot
            if (nrow(peak_point) > 0) {
                 xvar <- if (main_par()$selected_chr == "All") "BPcum" else "position"
                 p <- p + ggplot2::geom_point(data = peak_point,
                                           ggplot2::aes(x = .data[[xvar]], y = .data$LOD),
                                           color = "red", size = 3, shape = 18) # Use diamond shape
            }
        }
            
        # Create plot title and subtitle
        trait_title <- main_par()$which_trait # Use the value from the input
        chr_info_text <- if(main_par()$selected_chr != "All") paste0(" (Chromosome ", main_par()$selected_chr, ")") else ""
            
        subtitle <- if(!is.null(highest_peak_for_plot) && nrow(highest_peak_for_plot) > 0) {
            chr_label <- chr_XYM(highest_peak_for_plot$chr)
            paste0(
                "Peak Marker: ", highest_peak_for_plot$marker,
                " (Chr", chr_label, ":", round(highest_peak_for_plot$pos, 2), " Mb) | ",
                "LOD: ", round(highest_peak_for_plot$lod, 2)
            )
        } else {
            if(main_par()$selected_chr != "All") {
                "No significant peaks in this chromosome"
            } else {
                "No significant peaks"
            }
        }
            
        # Add title and subtitle to the plot
        p <- p + ggplot2::ggtitle(
                label = paste0(trait_title, chr_info_text),
                subtitle = subtitle
            )
            
        return(p)
    })

    output$download_qtl_plot_png <- shiny::downloadHandler(
        filename = function() {
            shiny::req(main_par()$selected_chr, main_par()$which_trait)
            chr_suffix <- if(main_par()$selected_chr != "All") paste0("_chr", main_par()$selected_chr) else ""
            paste0("lod_plot_", main_par()$which_trait, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".png")
        },
        content = function(file) {
            plot_obj <- generate_qtl_download_plot()
            dims <- plot_dims() # Get reactive dimensions
            shiny::req(plot_obj, dims$width, dims$height)
            # Save the plot with high resolution, using provided dimensions
            ggplot2::ggsave(file, plot_obj, width = dims$width/72, height = dims$height/72, dpi = 300, units = "in")
        }
    )

    output$download_qtl_plot_pdf <- shiny::downloadHandler(
        filename = function() {
            shiny::req(main_par()$selected_chr, main_par()$which_trait)
            chr_suffix <- if(main_par()$selected_chr != "All") paste0("_chr", main_par()$selected_chr) else ""
            paste0("lod_plot_", main_par()$which_trait, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".pdf")
        },
        content = function(file) {
            plot_obj <- generate_qtl_download_plot()
            dims <- plot_dims() # Get reactive dimensions
            shiny::req(plot_obj, dims$width, dims$height)
            # Save the plot as PDF, using provided dimensions
            tryCatch({
                ggplot2::ggsave(file, plot_obj, width = dims$width/72, height = dims$height/72, device = cairo_pdf, units = "in")
            }, error = function(e){
                 warning("device not available for PDF generation. Error: ", e$message)
                 ggplot2::ggsave(file, plot_obj, width = dims$width/72, height = dims$height/72, device = "pdf", units = "in")
            })
        }
    )

  })
} 