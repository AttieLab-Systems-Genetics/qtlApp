#' Split Allele Effects Module (Single or Side-by-Side)
#'
#' Renders allele effects for a selected additive peak or, when available,
#' renders two side-by-side allele effects cards for interactive splits
#' (e.g., Sex: Female vs Male, Diet: HC vs HF).
#'
#' @param id Module ID
#' @param selected_peak_reactive Reactive that returns a data.frame (single additive peak) or NULL
#' @param diff_peak_1_reactive Reactive that returns a data.frame for first split peak (or NULL)
#' @param diff_peak_2_reactive Reactive that returns a data.frame for second split peak (or NULL)
#'
#' @return A list of reactives (invisibly) that may be useful for testing
#'
#' @importFrom shiny moduleServer NS reactive renderUI uiOutput plotOutput renderPlot req tagList
#' @importFrom shinycssloaders withSpinner
#' @importFrom bslib card card_header card_body layout_columns
#' @importFrom htmltools tags
#' @export
splitAlleleEffectsUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("allele_effects_section"))
}

#' @rdname splitAlleleEffectsUI
#' @export
splitAlleleEffectsServer <- function(id,
                                     selected_peak_reactive,
                                     diff_peak_1_reactive,
                                     diff_peak_2_reactive) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Prepare data for single additive allele plot
    allele_effects_data <- shiny::reactive({
      peak_info <- selected_peak_reactive()
      if (is.null(peak_info)) {
        return(NULL)
      }
      reshaped_data <- pivot_peaks(peak_info, peak_info$marker)
      if (is.null(reshaped_data) || nrow(reshaped_data) == 0) {
        return(NULL)
      }
      reshaped_data$trait <- peak_info$trait
      reshaped_data
    })

    # Data for side-by-side split plots
    diff_allele_data_1 <- shiny::reactive({
      peak_info <- diff_peak_1_reactive()
      shiny::req(peak_info)
      reshaped <- pivot_peaks(peak_info, peak_info$marker)
      reshaped$trait <- peak_info$trait
      reshaped$plot_label <- peak_info$plot_label
      reshaped
    })

    diff_allele_data_2 <- shiny::reactive({
      peak_info <- diff_peak_2_reactive()
      shiny::req(peak_info)
      reshaped <- pivot_peaks(peak_info, peak_info$marker)
      reshaped$trait <- peak_info$trait
      reshaped$plot_label <- peak_info$plot_label
      reshaped
    })

    # Detailed peak info display (for single additive view)
    output$peak_info_display <- shiny::renderUI({
      peak_info <- selected_peak_reactive()
      if (is.null(peak_info) || nrow(peak_info) == 0) {
        return(htmltools::tags$div("No peak selected.", style = "color: #7f8c8d; text-align: center; padding-top: 20px;"))
      }

      info_elements <- list()
      info_elements <- c(info_elements, list(
        htmltools::tags$strong("Marker: "), peak_info$marker, htmltools::tags$br(),
        htmltools::tags$strong("Position: "), paste0("Chr", peak_info$qtl_chr, ":", round(peak_info$qtl_pos, 2), " Mb"), htmltools::tags$br(),
        htmltools::tags$strong("LOD Score: "), round(peak_info$qtl_lod, 3), htmltools::tags$br()
      ))

      # Cis/Trans status
      if ("cis" %in% colnames(peak_info)) {
        cis_status <- if (is.logical(peak_info$cis)) {
          ifelse(peak_info$cis, "Cis", "Trans")
        } else if (is.character(peak_info$cis)) {
          ifelse(toupper(peak_info$cis) %in% c("TRUE", "1", "YES"), "Cis", "Trans")
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
      if ("qtl_ci_lo" %in% colnames(peak_info) && "qtl_ci_hi" %in% colnames(peak_info)) {
        if (!is.na(peak_info$qtl_ci_lo) && !is.na(peak_info$qtl_ci_hi)) {
          info_elements <- c(info_elements, list(
            htmltools::tags$strong("95% CI: "),
            paste0(round(peak_info$qtl_ci_lo, 2), " - ", round(peak_info$qtl_ci_hi, 2), " Mb"),
            htmltools::tags$br()
          ))
        }
      }

      # Founder allele effects display (A-H mapped to strains)
      allele_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
      strain_names <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
      available_alleles <- allele_cols[allele_cols %in% colnames(peak_info)]
      if (length(available_alleles) > 0) {
        allele_effects <- list()
        for (i in seq_along(available_alleles)) {
          col <- available_alleles[i]
          value <- peak_info[[col]]
          if (!is.na(value) && !is.null(value)) {
            strain <- strain_names[i]
            allele_effects[[length(allele_effects) + 1]] <- paste0(strain, ": ", round(value, 3))
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

      do.call(htmltools::tagList, info_elements)
    })

    # Single additive allele effects plot
    output$allele_effects_plot_output <- shiny::renderPlot({
      effects_data <- allele_effects_data()
      if (is.null(effects_data)) {
        return(ggplot2::ggplot() +
                 ggplot2::theme_void() +
                 ggplot2::labs(title = "No strain effects data available for this peak") +
                 ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, color = "#7f8c8d")))
      }
      ggplot_alleles(effects_data)
    })

    # Titles for split plots
    output$diff_plot_title_1 <- shiny::renderText({
      data <- diff_allele_data_1(); shiny::req(data)
      unique(data$plot_label)
    })
    output$diff_plot_title_2 <- shiny::renderText({
      data <- diff_allele_data_2(); shiny::req(data)
      unique(data$plot_label)
    })

    # Split plots
    output$diff_allele_plot_1 <- shiny::renderPlot({
      data <- diff_allele_data_1(); shiny::req(data)
      ggplot_alleles(data)
    })
    output$diff_allele_plot_2 <- shiny::renderPlot({
      data <- diff_allele_data_2(); shiny::req(data)
      ggplot_alleles(data)
    })

    # High level UI switch (single vs split)
    output$allele_effects_section <- shiny::renderUI({
      additive_peak <- selected_peak_reactive()
      diff_peak_1 <- diff_peak_1_reactive()
      diff_peak_2 <- diff_peak_2_reactive()

      # Single additive
      if (!is.null(additive_peak)) {
        return(shiny::tagList(
          htmltools::tags$hr(style = "margin: 20px 0; border-top: 2px solid #3498db;"),
          htmltools::tags$div(
            style = "margin-bottom: 15px;",
            htmltools::tags$h5("Strain Effects", style = "color: #2c3e50; font-weight: bold;"),
            htmltools::tags$p("Showing strain effects for the selected peak.", style = "font-size: 12px;")
          ),
          bslib::layout_columns(
            col_widths = c(5, 7),
            htmltools::tags$div(
              style = "background: #f8f9fa; padding: 10px; border-radius: 5px; height: 450px; overflow-y: auto;",
              shiny::uiOutput(ns("peak_info_display"))
            ),
            shiny::plotOutput(ns("allele_effects_plot_output"), height = "450px", width = "450px") |>
              shinycssloaders::withSpinner(type = 8, color = "#3498db")
          )
        ))
      }

      # Side-by-side split when difference peaks exist
      if (!is.null(diff_peak_1) || !is.null(diff_peak_2)) {
        ui_elements <- list()
        if (!is.null(diff_peak_1)) {
          ui_elements <- c(ui_elements, list(
            bslib::card(
              bslib::card_header(shiny::textOutput(ns("diff_plot_title_1"))),
              bslib::card_body(
                shiny::plotOutput(ns("diff_allele_plot_1"), height = "400px") |>
                  shinycssloaders::withSpinner(type = 8, color = "#e74c3c")
              )
            )
          ))
        }
        if (!is.null(diff_peak_2)) {
          ui_elements <- c(ui_elements, list(
            bslib::card(
              bslib::card_header(shiny::textOutput(ns("diff_plot_title_2"))),
              bslib::card_body(
                shiny::plotOutput(ns("diff_allele_plot_2"), height = "400px") |>
                  shinycssloaders::withSpinner(type = 8, color = "#e74c3c")
              )
            )
          ))
        }

        info_message <- NULL
        if (is.null(diff_peak_1) || is.null(diff_peak_2)) {
          info_message <- htmltools::tags$p(
            "Note: A corresponding peak was found in only one of the comparison datasets within the search window.",
            style = "font-style: italic; font-size: 12px; color: #7f8c8d; margin-top: 0;"
          )
        }

        return(shiny::tagList(
          htmltools::tags$hr(style = "margin: 20px 0; border-top: 2px solid #e74c3c;"),
          htmltools::tags$h5("Comparative Strain Effects", style = "color: #2c3e50; font-weight: bold; margin-bottom: 15px;"),
          htmltools::tags$p("Showing strain effects for the peak found in each respective dataset based on your click on the difference plot.", style = "font-size: 12px; margin-bottom: 2px;"),
          info_message,
          do.call(bslib::layout_columns, c(list(col_widths = 6), ui_elements))
        ))
      }

      # Default: nothing to render
      return(NULL)
    })

    invisible(list(
      additive_data = allele_effects_data,
      diff_data_1 = diff_allele_data_1,
      diff_data_2 = diff_allele_data_2
    ))
  })
}


