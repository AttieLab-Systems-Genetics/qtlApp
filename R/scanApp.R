#' Refactored Scan App Module
#'
#' This is a streamlined version that uses existing modules:
#' - mainParApp for dataset selection and trait search
#' - chromosomeControlsApp for chromosome zoom functionality
#' - downloadApp for plot downloads
#' - existing plot modules (manhattanPlotApp, cisTransPlotApp, scanlyApp)
#'
#' @importFrom shiny moduleServer reactive renderText req selectInput shinyApp
#'             uiOutput observeEvent tags
#' @importFrom bslib page_sidebar sidebar layout_columns card card_header card_body
#' @importFrom shinyjs useShinyjs
#' @importFrom htmltools tagList div h4 hr h5
#' @export
scanApp <- function() {
  # Ensure ui_styles.R is sourced
  if (!exists("custom_css", mode = "character")) {
    source("R/ui_styles.R")
  }

  # Source required modules
  if (!exists("mainParServer", mode = "function")) {
    source("R/mainParApp.R")
  }
  if (!exists("chromosomeControlsServer", mode = "function")) {
    source("R/chromosomeControlsApp.R")
  }
  if (!exists("downloadServer", mode = "function")) {
    source("R/downloadApp.R")
  }
  if (!exists("plot_null", mode = "function")) {
    source("R/plot_null.R")
  }

  ui <- bslib::page_sidebar(
    shinyjs::useShinyjs(),
    tags$head(tags$style(custom_css)),
    create_title_panel(
      "QTL Scan Visualizer",
      "Interactive visualization tool using modular architecture"
    ),
    sidebar = bslib::sidebar(
      # Use existing mainPar module for dataset selection and trait search
      h4("📊 Dataset & Trait Selection", style = "color: #2c3e50;"),
      mainParInput("main_par"),
      mainParUI("main_par"),
      hr(style = "border-top: 1px solid #bdc3c7; margin: 20px 0;"),

      # Add chromosome controls module
      h4("🔍 Chromosome Controls", style = "color: #2c3e50;"),
      chromosomeControlsInput("chr_controls"),
      hr(style = "border-top: 1px solid #bdc3c7; margin: 20px 0;"),

      # Use existing download module
      h4("💾 Download Plots", style = "color: #2c3e50;"),
      downloadInput("download"),
      downloadOutput("download")
    ),

    # Main content area with plot outputs
    bslib::layout_columns(
      col_widths = bslib::breakpoints(
        sm = c(12, 12),
        md = c(6, 6)
      ),

      # Primary plot (Manhattan or Cis/Trans based on dataset category)
      bslib::card(
        id = "primary_plot_card",
        bslib::card_header(shiny::textOutput("plot_title")),
        bslib::card_body(
          shiny::uiOutput("conditional_plot_ui")
        )
      ),

      # LOD scan plot (when trait is selected)
      bslib::card(
        id = "lod_scan_card",
        bslib::card_header("LOD Scan"),
        bslib::card_body(
          scanlyOutput("scanly"),
          scanlyUI("scanly")
        )
      )
    ),

    # Download preview section
    bslib::card(
      id = "download_preview_card",
      bslib::card_header("📋 Download Preview"),
      bslib::card_body(
        downloadUI("download")
      )
    )
  )

  server <- function(input, output, session) {
    # Initialize data import
    import_reactives <- importServer("import")

    # Initialize main parameter controls (dataset selection, trait search, LOD threshold)
    main_par <- mainParServer("main_par", import_reactives)

    # Initialize chromosome controls
    chr_controls <- chromosomeControlsServer("chr_controls")

    # Determine plot type based on dataset category
    plot_type <- shiny::reactive({
      shiny::req(main_par$dataset_category())
      category <- main_par$dataset_category()
      if (category %in% c("Liver Lipids", "Clinical Traits")) {
        "manhattan"
      } else if (category %in% c("Liver Genes", "Liver Isoforms", "Plasma 2H Metabolites")) {
        "cistrans"
      } else {
        "unknown"
      }
    })

    # Dynamic plot title
    output$plot_title <- shiny::renderText({
      shiny::req(main_par$dataset_category(), main_par$selected_dataset())
      category <- main_par$dataset_category()
      dataset <- main_par$selected_dataset()

      plot_type_text <- switch(plot_type(),
        "manhattan" = "Manhattan Plot",
        "cistrans" = "Cis/Trans Plot",
        "Plot"
      )

      paste0(plot_type_text, " for: ", dataset, " (", category, ")")
    })

    # Dynamic plot UI based on dataset category
    output$conditional_plot_ui <- shiny::renderUI({
      switch(plot_type(),
        "manhattan" = manhattanPlotUI("manhattan_plot"),
        "cistrans" = tagList(
          cisTransPlotInput("cistrans_plot"),
          cisTransPlotUI("cistrans_plot")
        ),
        shiny::p("No plot type configured for this dataset category.")
      )
    })

    # Initialize plot servers with shared cache
    peaks_cache <- new.env(parent = emptyenv())
    trait_cache <- new.env(parent = emptyenv())

    # Combined main_par that includes chromosome selection
    enhanced_main_par <- shiny::reactive({
      list(
        dataset_category = main_par$dataset_category,
        selected_dataset = main_par$selected_dataset,
        which_trait = main_par$which_trait,
        LOD_thr = main_par$LOD_thr,
        selected_chr = chr_controls$selected_chr # Add chromosome selection
      )
    })

    # Initialize plot servers
    manhattan_outputs <- manhattanPlotServer("manhattan_plot",
      import_reactives = import_reactives,
      main_par = enhanced_main_par
    )

    cistrans_outputs <- cisTransPlotServer("cistrans_plot",
      import_reactives = import_reactives,
      main_par = enhanced_main_par,
      peaks_cache = peaks_cache
    )

    # Initialize scan server for LOD scans
    scan_outputs <- scanServer("scan_list", enhanced_main_par, import_reactives, trait_cache)
    peak_outputs <- peakServer("peak_list", enhanced_main_par, import_reactives, peaks_cache)

    # Initialize scanly server for interactive LOD plots
    scanlyServer("scanly", enhanced_main_par, scan_outputs, peak_outputs)

    # Create download list structure
    download_list <- shiny::reactiveValues(
      filename = shiny::reactive({
        dataset_name <- main_par$selected_dataset()
        if (is.null(dataset_name)) dataset_name <- "qtl_analysis"
        trait_name <- main_par$which_trait()
        base_name <- if (!is.null(trait_name) && nzchar(trait_name)) {
          paste(dataset_name, trait_name, sep = "_")
        } else {
          dataset_name
        }
        gsub("[^A-Za-z0-9_-]", "_", base_name)
      }),
      plots = shiny::reactiveValues(
        manhattan = shiny::reactive({
          if (plot_type() == "manhattan" && !is.null(manhattan_outputs)) {
            # We need to get the actual ggplot object from manhattanPlotServer
            # This would require modifying manhattanPlotServer to return plot objects
            plot_null("Manhattan plot available but needs ggplot export")
          } else {
            plot_null("Manhattan plot not available")
          }
        }),
        cistrans = shiny::reactive({
          if (plot_type() == "cistrans" && !is.null(cistrans_outputs)) {
            plot_null("Cis/Trans plot available but needs ggplot export")
          } else {
            plot_null("Cis/Trans plot not available")
          }
        }),
        lod_scan = shiny::reactive({
          if (!is.null(scan_outputs$plots$scan)) {
            scan_outputs$plots$scan()
          } else {
            plot_null("LOD scan plot not available")
          }
        })
      ),
      tables = shiny::reactiveValues(
        scan_data = shiny::reactive({
          if (!is.null(scan_outputs$tables$scan)) {
            scan_outputs$tables$scan()
          } else {
            data.frame(Info = "No scan data available")
          }
        })
      )
    )

    # Initialize download server
    downloadServer("download", download_list)

    # Clear caches when dataset changes
    shiny::observeEvent(main_par$selected_dataset(), {
      rm(list = ls(envir = trait_cache), envir = trait_cache)
      rm(list = ls(envir = peaks_cache), envir = peaks_cache)
    })
  }

  shiny::shinyApp(ui = ui, server = server)
}
