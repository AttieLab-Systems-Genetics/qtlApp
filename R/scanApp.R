#' Modular Scan App
#'
#' @param id shiny identifier
#'
#' @importFrom DT DTOutput renderDT datatable
#' @importFrom shiny actionButton h4 moduleServer NS plotOutput reactive reactiveVal
#'             reactiveValues renderPlot renderUI req setProgress shinyApp selectInput
#'             uiOutput withProgress div downloadButton numericInput tagList
#'             observeEvent updateNumericInput downloadHandler tabsetPanel tabPanel
#'             observe
#' @importFrom plotly plotlyOutput renderPlotly ggplotly event_data layout config
#' @importFrom bslib card card_header page_sidebar sidebar layout_columns navset_tab nav_panel card_body
#' @importFrom shinycssloaders withSpinner
#' @importFrom dplyr across mutate where filter select
#' @importFrom stringr str_split str_remove
#' @importFrom ggplot2 ggsave
#' @importFrom shinyjs useShinyjs
#' @importFrom htmltools tagList tags
#' @importFrom stats setNames
#'
#' @export
scanApp <- function() {
  # Source required module files and helpers first
  source("R/helpers.R")
  source("R/data_handling.R")
  source("R/scanApp_monolithic_backup.R")
  source("R/importApp.R")
  source("R/datasetSelectionModule.R")
  source("R/traitSearchModule.R")
  source("R/interactiveAnalysisModule.R")
  source("R/scanPlotModule.R")
  source("R/alleleEffectsModule.R")
  source("R/manhattanPlotApp.R")
  source("R/cisTransPlotApp.R")
  source("R/profilePlotApp.R")
  source("R/correlationApp.R")

  # Ensure ui_styles.R is sourced
  if (!exists("custom_css", mode = "character")) {
    source("R/ui_styles.R")
  }

  ui <- bslib::page_sidebar(
    shinyjs::useShinyjs(),
    tags$head(
      # Viewport meta tag for responsive design
      tags$meta(name = "viewport", content = "width=device-width, initial-scale=1.0, user-scalable=yes"),

      # Autosizing CSS and JavaScript
      tags$link(rel = "stylesheet", type = "text/css", href = "autosize.css"),
      tags$script(src = "autosize.js"),

      # Custom app styles
      tags$style(custom_css)
    ),

    # Dataset selection UI (includes top navigation bar)
    datasetSelectionUI("dataset_selection"),
    sidebar = bslib::sidebar(
      width = 600, # Increased sidebar width for better screen coverage
      # Tabbed sidebar content
      bslib::navset_pill(
        id = "sidebar_tabs",

        # Tab 1: Dataset Selection and Controls
        bslib::nav_panel(
          "Data Search",
          # Trait search module
          traitSearchUI("trait_search"),
        ),

        # Tab 2: Overview Plot (Manhattan/Cis-Trans)
        bslib::nav_panel(
          "LOD peaks",
          div(
            style = "padding: 10px;",

            # LOD Threshold Control (moved here from Data Search tab)
            h5("Peak Filtering", style = "color: #2c3e50; margin-bottom: 10px; font-weight: bold;"),
            # Dynamic LOD threshold slider that updates based on scan type
            uiOutput("lod_threshold_slider"),
            p("Filters peaks shown in the plot below",
              style = "font-size: 11px; color: #7f8c8d; margin: 5px 0 15px 0;"
            ),

            # Peak Selection Dropdown for future peak differences analysis
            hr(style = "border-top: 1px solid #bdc3c7; margin: 15px 0;"),
            h5("🎯 Peak Analysis", style = "color: #2c3e50; margin-bottom: 10px; font-weight: bold;"),
            p("Select peaks for detailed analysis and future comparison features",
              style = "font-size: 11px; color: #7f8c8d; margin: 5px 0 15px 0;"
            ),
            hr(style = "border-top: 1px solid #bdc3c7; margin: 15px 0;"),
            div(
              id = "overview-plot-container",
              class = "overview-plot-container",
              style = "height: 65vh; min-height: 400px; max-height: 800px; border: 1px solid #bdc3c7; border-radius: 5px; overflow: hidden;",
              shiny::uiOutput("conditional_plot_ui")
            ),
            p("Click on points to view detailed LOD scans. Plot titles show dataset and analysis type.",
              style = "font-size: 11px; color: #7f8c8d; margin: 10px 0 0 0; text-align: center;"
            )
          )
        )
      ),
    ),

    # Simplified main area - just the LOD scan plot
    bslib::card(
      id = "lod_scan_card",
      bslib::card_header(uiOutput("main_plot_title")),
      bslib::card_body(
        shiny::uiOutput("lod_scan_plot_ui_placeholder")
      )
    )
  )

  server <- function(input, output, session) {
    # Source helper for allele plots if not already available
    if (!exists("ggplot_alleles", mode = "function")) {
      source("R/ggplot_alleles.R")
    }

    # Initialize data import
    import_reactives <- importServer("import")

    # Initialize dataset selection module
    dataset_selection <- datasetSelectionServer("dataset_selection", import_reactives)

    # Initialize trait search module
    trait_search <- traitSearchServer("trait_search", import_reactives, dataset_selection$selected_dataset)

    # Initialize interactive analysis module
    interactive_analysis <- interactiveAnalysisServer("interactive_analysis", dataset_selection$selected_dataset)

    # Create cache environments
    trait_cache <- new.env(parent = emptyenv())
    peaks_cache <- new.env(parent = emptyenv())

    # Clear caches when dataset changes
    shiny::observeEvent(dataset_selection$selected_dataset(),
      {
        selected_ds_val <- dataset_selection$selected_dataset()

        if (!is.null(selected_ds_val)) {
          message(paste("Clearing caches for dataset:", selected_ds_val))
          rm(list = ls(envir = trait_cache), envir = trait_cache)
          rm(list = ls(envir = peaks_cache), envir = peaks_cache)
        } else {
          message("No dataset selected or selected_dataset is NULL, not clearing caches.")
        }
      },
      ignoreNULL = FALSE,
      ignoreInit = TRUE
    )

    # Detect if current dataset is additive or interactive based on dataset name and interaction type
    scan_type <- shiny::reactive({
      dataset_name <- interactive_analysis$mapped_dataset()

      if (is.null(dataset_name) || dataset_name == "") {
        return("additive") # Default to additive
      }

      # Check if dataset name contains "interactive" or if interaction type is selected
      if (grepl("interactive", dataset_name, ignore.case = TRUE)) {
        return("interactive")
      } else {
        return("additive")
      }
    })

    # Dynamic LOD threshold slider based on scan type
    output$lod_threshold_slider <- shiny::renderUI({
      current_scan_type <- scan_type()

      # Set different minimums based on scan type
      min_val <- if (current_scan_type == "interactive") 10.5 else 7.5
      default_val <- min_val # Start at minimum value

      sliderInput("LOD_thr",
        label = paste("LOD Threshold (", current_scan_type, "scan):"),
        min = min_val, max = 20, value = default_val, step = 0.5,
        width = "100%"
      )
    })

    # Our own LOD threshold reactive
    lod_threshold_rv <- shiny::reactive({
      current_scan_type <- scan_type()
      default_threshold <- if (current_scan_type == "interactive") 10.5 else 7.5
      input$LOD_thr %||% default_threshold
    }) %>% shiny::debounce(300) # Debounce LOD threshold to prevent rapid re-firing

    # Reactive for selected chromosome (for zooming into specific chromosomes)
    selected_chromosome_rv <- shiny::reactive({
      input$selected_chr %||% "All" # Default to "All" chromosomes
    })

    # Create our own main_par structure (compatible with existing modules)
    active_main_par <- shiny::reactive({
      # Create a compatible structure with the expected reactives
      list(
        selected_dataset = interactive_analysis$mapped_dataset, # Mapped dataset for interaction
        LOD_thr = lod_threshold_rv, # Our LOD threshold
        selected_chr = selected_chromosome_rv, # Selected chromosome
        which_trait = trait_search$trait_for_lod_scan, # Currently searched trait
        dataset_category = dataset_selection$dataset_category # Our dataset category
      )
    })

    # Instantiate plot modules and capture their outputs
    manhattan_plot_outputs <- manhattanPlotServer("manhattan_plot_module",
      import_reactives = import_reactives,
      main_par = active_main_par,
      sidebar_interaction_type = interactive_analysis$interaction_type
    )

    cistrans_plot_outputs <- cisTransPlotServer("cistrans_plot_module",
      import_reactives = import_reactives,
      main_par = active_main_par,
      peaks_cache = peaks_cache,
      sidebar_interaction_type = interactive_analysis$interaction_type
    )

    # Observe clicks from Manhattan plot
    shiny::observeEvent(manhattan_plot_outputs$clicked_phenotype_for_lod_scan(),
      {
        clicked_trait <- manhattan_plot_outputs$clicked_phenotype_for_lod_scan()
        if (!is.null(clicked_trait)) {
          trait_search$trait_for_lod_scan(clicked_trait)
          message(paste("scanApp: Manhattan plot click detected. Trait for LOD scan:", clicked_trait))
        }
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE
    )

    # Observe clicks from Cis/Trans plot
    shiny::observeEvent(cistrans_plot_outputs$clicked_phenotype_for_lod_scan(),
      {
        clicked_trait <- cistrans_plot_outputs$clicked_phenotype_for_lod_scan()
        if (!is.null(clicked_trait)) {
          trait_search$trait_for_lod_scan(clicked_trait)
          message(paste("scanApp: Cis/Trans plot click detected. Trait for LOD scan:", clicked_trait))
        }
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE
    )



    output$conditional_plot_ui <- shiny::renderUI({
      category <- dataset_selection$dataset_category()
      shiny::req(category)

      if (category %in% c("Liver Lipids", "Clinical Traits", "Plasma Metabolites")) {
        tagList(
          manhattanPlotUI("manhattan_plot_module")
        )
      } else if (category %in% c("Liver Genes", "Liver Isoforms")) {
        tagList(
          cisTransPlotInput("cistrans_plot_module"),
          cisTransPlotUI("cistrans_plot_module")
        )
      } else {
        shiny::p(paste("No specific plot type configured for category:", category))
      }
    })

    # Initialize scan plot module
    scan_plot_outputs <- scanServer(
      id = "scan_plot_module",
      trait_to_scan = trait_search$trait_for_lod_scan,
      selected_dataset_group = interactive_analysis$mapped_dataset,
      import_reactives = import_reactives,
      main_par_inputs = active_main_par,
      interaction_type_reactive = interactive_analysis$interaction_type
    )

    # Initialize allele effects module, now driven by the scan plot's selected peak
    allele_effects_outputs <- alleleEffectsServer(
      id = "allele_effects",
      selected_peak_reactive = scan_plot_outputs$selected_peak
    )

    # Dynamic title for the main LOD scan card
    output$main_plot_title <- shiny::renderUI({
      trait <- trait_search$trait_for_lod_scan()
      title_text <- if (!is.null(trait) && nzchar(trait)) {
        paste("LOD Scan for Trait:", trait)
      } else {
        "LOD Scan - Detailed View"
      }
      h4(title_text, style = "font-weight: bold; margin-bottom: 0;")
    })

    # UI for LOD Scan plot - only appears when a trait is selected
    output$lod_scan_plot_ui_placeholder <- shiny::renderUI({
      if (!is.null(trait_search$trait_for_lod_scan()) && nzchar(trait_search$trait_for_lod_scan())) {
        tagList(
          # Combined controls row
          div(
            style = "margin-bottom: 15px; background: #f8f9fa; padding: 10px 15px; border-radius: 4px; border: 1px solid #bdc3c7;",
            div(
              style = "display: flex; align-items: flex-end; gap: 15px; flex-wrap: wrap;",

              # Interaction analysis dropdown (now on the left)
              div(
                style = "flex: 1 1 180px; min-width: 180px;",
                interactiveAnalysisUI("interactive_analysis")
              ),

              # Chromosome selector
              div(
                style = "flex: 1 1 120px; min-width: 120px;",
                shiny::selectInput(
                  "selected_chr",
                  label = "Chromosome:",
                  choices = c(
                    "All" = "All",
                    setNames(as.character(1:19), paste("Chr", 1:19)),
                    "X" = "X", "Y" = "Y", "M" = "M"
                  ),
                  selected = "All",
                  width = "100%"
                )
              ),

              # Zoom/Reset buttons
              div(
                style = "flex: 0 0 auto;",
                shiny::actionButton(
                  "reset_chr_view", "🌐 Reset Zoom",
                  class = "btn btn-sm btn-secondary",
                  style = "background: #7f8c8d; border: none; color: white; font-size: 11px; padding: 4px 8px;"
                )
              )
            )
          ),

          # Conditional panel for interaction info
          shiny::conditionalPanel(
            condition = "input['interactive_analysis-interaction_type'] != 'none'",
            div(
              style = "margin-bottom: 15px; padding: 10px; background-color: #e8f4fd; border-radius: 5px; border-left: 4px solid #3498db;",
              p("ℹ️ Interactive analysis will show stacked plots: Interactive LOD scan (top) and Difference plot (Interactive - Additive, bottom).",
                style = "font-size: 12px; color: #2c3e50; margin: 0;"
              )
            )
          ),

          # LOD scan plot module
          scanPlotUI("scan_plot_module"),
          # Clicked point details table
          div(
            style = "margin-top: 15px;",
            DT::DTOutput("lod_scan_click_table")
          ),
          # Allele effects module
          alleleEffectsUI("allele_effects")
        )
      } else {
        # Show a placeholder message when no trait is selected
        div(
          style = "text-align: center; padding-top: 50px; color: #7f8c8d;",
          h5("No trait selected for LOD scan"),
          p("Select a trait from the 'Data Search' tab or click on a point in an overview plot to begin.")
        )
      }
    })

    # Render the LOD scan click details table
    output$lod_scan_click_table <- DT::renderDT({
      # Check if we have scan module outputs and clicked point details
      if (!is.null(scan_plot_outputs) && !is.null(scan_plot_outputs$clicked_point_details)) {
        clicked_details <- scan_plot_outputs$clicked_point_details()

        if (!is.null(clicked_details) && nrow(clicked_details) > 0) {
          # Format the clicked point data for display
          display_data <- data.frame(
            Chromosome = if ("chr" %in% colnames(clicked_details)) clicked_details$chr else "N/A",
            Marker = if ("markers" %in% colnames(clicked_details)) clicked_details$markers else "N/A",
            Position = if ("position" %in% colnames(clicked_details)) paste0(clicked_details$position, " Mb") else "N/A",
            LOD = if ("LOD" %in% colnames(clicked_details)) clicked_details$LOD else "N/A"
          )

          message(paste("scanApp: Displaying clicked point details for marker:", display_data$Marker))

          return(DT::datatable(
            display_data,
            options = list(
              dom = "t",
              paging = FALSE,
              searching = FALSE,
              columnDefs = list(list(targets = "_all", className = "dt-center"))
            ),
            rownames = FALSE,
            selection = "none",
            class = "compact hover"
          ))
        }
      }

      # Default message when no point is clicked
      return(DT::datatable(
        data.frame(Info = "Click on a point in the LOD scan plot to see details"),
        options = list(dom = "t", paging = FALSE, searching = FALSE),
        rownames = FALSE,
        selection = "none"
      ))
    })

    # CHROMOSOME ZOOM FUNCTIONALITY
    observeEvent(input$selected_chr,
      {
        selected_chr <- input$selected_chr
        if (!is.null(selected_chr) && selected_chr != "All") {
          shiny::showNotification(
            paste("Zoomed to chromosome", selected_chr),
            type = "message",
            duration = 2
          )
        }
      },
      ignoreInit = TRUE
    )

    # Observer for "Show All Chromosomes" button
    observeEvent(input$reset_chr_view, {
      message("scanApp: Resetting to show all chromosomes")
      shiny::updateSelectInput(session, "selected_chr",
        selected = "All"
      )
      shiny::showNotification(
        "Showing all chromosomes",
        type = "message",
        duration = 2
      )
    })

    # Debug output for interaction type
    output$debug_interaction_type <- shiny::renderText({
      interaction_type <- interactive_analysis$interaction_type()
      mapped_dataset <- interactive_analysis$mapped_dataset()

      paste(
        "Interaction Type:", interaction_type %||% "NULL",
        "| Mapped Dataset:", mapped_dataset %||% "NULL",
        "| Scan Type:", interactive_analysis$scan_type() %||% "NULL"
      )
    })

    # Define the %||% operator for null coalescing
    `%||%` <- function(a, b) if (!is.null(a)) a else b
  }

  shiny::shinyApp(ui = ui, server = server)
}
