#' Scan App Module
#'
#' @param id shiny identifier
#' @param main_par reactive list with selected_dataset, LOD_thr and which_trait
#' @param import reactive list with file_directory and markers
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
#' @importFrom stats setNames # Added stats for setNames
#'
#' @export
scanApp <- function() {
  # Ensure ui_styles.R is sourced
  if (!exists("custom_css", mode = "character")) {
    source("R/ui_styles.R")
  }

  ui <- bslib::page_sidebar(
    shinyjs::useShinyjs(),
    tags$head(tags$style(custom_css)),

    # Top navigation bar for dataset category selection
    div(
      style = "background: linear-gradient(135deg, #2c3e50, #3498db); padding: 15px; margin-bottom: 20px; border-radius: 8px;",
      div(
        style = "display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap;",
        h3("QTL Scan Visualizer",
          style = "color: white; margin: 0; font-weight: bold;"
        ),
        div(
          style = "display: flex; align-items: center; gap: 15px;",
          h5("Dataset Category:",
            style = "color: white; margin: 0; font-weight: bold;"
          ),
          shiny::selectInput(shiny::NS("app_controller", "dataset_category_selector"),
            NULL,
            choices = c("Loading..." = ""),
            width = "200px"
          )
        )
      )
    ),
    sidebar = bslib::sidebar(
      width = 600, # Increased sidebar width a bit more for better screen coverage
      # Tabbed sidebar content
      bslib::navset_pill(
        id = "sidebar_tabs",

        # Tab 1: Dataset Selection and Controls
        bslib::nav_panel(
          "Data Search",
          div(
            style = "padding: 10px;",

            # Dataset selection
            h5("Dataset Selection", style = "color: #2c3e50; margin-bottom: 15px; font-weight: bold;"),
            shiny::selectInput(shiny::NS("app_controller", "specific_dataset_selector"),
              "Select Specific Dataset:",
              choices = c("Loading..." = ""),
              width = "100%"
            ),

            # Trait search section
            hr(style = "border-top: 2px solid #3498db; margin: 20px 0;"),
            h5("🔍 Trait Search", style = "color: #2c3e50; margin-bottom: 15px; font-weight: bold;"),
            selectizeInput(shiny::NS("app_controller", "trait_search_input"),
              "Search for traits/genes:",
              choices = NULL,
              selected = NULL,
              multiple = FALSE,
              options = list(
                placeholder = "Type to search (e.g., Gapdh, Insulin, PI_38_3)",
                maxItems = 1,
                maxOptions = 10,
                create = FALSE
              ),
              width = "100%"
            ),
            div(
              style = "text-align: center; margin-top: 10px;",
              actionButton(shiny::NS("app_controller", "trait_search_button"),
                "🚀 Search & Plot LOD Scan",
                icon = icon("search"),
                class = "btn-primary",
                style = "background: #3498db; border: none; font-weight: bold; width: 100%;"
              )
            ),

            # Back button
            hr(style = "border-top: 1px solid #bdc3c7; margin: 20px 0;"),
            shiny::actionButton(shiny::NS("app_controller", "clear_lod_scan_btn"),
              "← Back to Overview Plot",
              icon = shiny::icon("arrow-left"),
              class = "btn-secondary",
              style = "width: 100%;"
            )
          )
        ),

        # Tab 2: Overview Plot (Manhattan/Cis-Trans)
        bslib::nav_panel(
          "LOD peaks",
          div(
            style = "padding: 10px;",

            # LOD Threshold Control (moved here from Data Search tab)
            h5("Peak Filtering", style = "color: #2c3e50; margin-bottom: 10px; font-weight: bold;"),
            # Dynamic LOD threshold slider that updates based on scan type
            uiOutput(shiny::NS("app_controller", "lod_threshold_slider")),
            p("Filters peaks shown in the plot below",
              style = "font-size: 11px; color: #7f8c8d; margin: 5px 0 15px 0;"
            ),
            hr(style = "border-top: 1px solid #bdc3c7; margin: 15px 0;"),
            h5(shiny::textOutput(shiny::NS("app_controller", "plot_title")),
              style = "color: #2c3e50; margin-bottom: 15px; font-weight: bold; text-align: center;"
            ),
            div(
              style = "height: 65vh; border: 1px solid #bdc3c7; border-radius: 5px; overflow: hidden;", # Increased to 65% viewport height
              shiny::uiOutput(shiny::NS("app_controller", "conditional_plot_ui"))
            ),
            p("Click on points to view detailed LOD scans",
              style = "font-size: 11px; color: #7f8c8d; margin: 10px 0 0 0; text-align: center;"
            )
          )
        )
      ),

      # Horizontal separator
      hr(style = "border-top: 2px solid #3498db; margin: 20px 0;"),

      # Additional Analyses section below the main tabs
      h5("📈 Additional Analyses", style = "color: #2c3e50; margin-bottom: 15px; font-weight: bold;"),
      bslib::navset_pill(
        id = "additional_analyses_tabs",

        # Profile Plot tab
        bslib::nav_panel(
          "Profile Plot",
          div(
            style = "padding: 10px;",
            div(
              style = "height: 50vh; border: 1px solid #bdc3c7; border-radius: 5px; overflow: hidden;",
              shinycssloaders::withSpinner(
                plotly::plotlyOutput(shiny::NS("app_controller", "profile_plot_output"),
                  height = "calc(50vh - 40px)"
                )
              )
            ),
            p("Profile plot visualization coming soon",
              style = "font-size: 11px; color: #7f8c8d; margin: 10px 0 0 0; text-align: center;"
            )
          )
        ),

        # Correlation tab
        bslib::nav_panel(
          "Correlation",
          div(
            style = "padding: 10px;",
            div(
              style = "height: 50vh; border: 1px solid #bdc3c7; border-radius: 5px; overflow: hidden;",
              shinycssloaders::withSpinner(
                plotly::plotlyOutput(shiny::NS("app_controller", "correlation_plot_output"),
                  height = "calc(50vh - 40px)"
                )
              )
            ),
            p("Correlation analysis visualization coming soon",
              style = "font-size: 11px; color: #7f8c8d; margin: 10px 0 0 0; text-align: center;"
            )
          )
        )
      )
    ),

    # Simplified main area - just the LOD scan plot
    bslib::card(
      id = "lod_scan_card",
      bslib::card_header("LOD Scan - Detailed View"),
      bslib::card_body(
        shiny::uiOutput(shiny::NS("app_controller", "lod_scan_plot_ui_placeholder"))
      )
    )
  )

  server <- function(input, output, session) {
    ns_app_controller <- shiny::NS("app_controller")

    # Source helper for allele plots if not already available
    if (!exists("ggplot_alleles", mode = "function")) {
      source("R/ggplot_alleles.R")
    }

    # Define the %||% operator for null coalescing
    `%||%` <- function(a, b) if (!is.null(a)) a else b

    trait_cache <- new.env(parent = emptyenv())
    peaks_cache <- new.env(parent = emptyenv())

    import_reactives <- importServer("import")

    # Detect if current dataset is additive or interactive based on dataset name
    scan_type <- shiny::reactive({
      shiny::req(input[[ns_app_controller("specific_dataset_selector")]])
      dataset_name <- input[[ns_app_controller("specific_dataset_selector")]]

      if (is.null(dataset_name) || dataset_name == "") {
        return("additive") # Default to additive
      }

      # Check if dataset name contains "interactive" or "diet_interactive"
      if (grepl("interactive|diet_interactive", dataset_name, ignore.case = TRUE)) {
        return("interactive")
      } else {
        return("additive")
      }
    })

    # Dynamic LOD threshold slider based on scan type
    output[[ns_app_controller("lod_threshold_slider")]] <- shiny::renderUI({
      current_scan_type <- scan_type()

      # Set different minimums based on scan type
      min_val <- if (current_scan_type == "interactive") 10.5 else 7.5
      default_val <- min_val # Start at minimum value

      sliderInput(shiny::NS("app_controller", "LOD_thr"),
        label = paste("LOD Threshold (", current_scan_type, "scan):"),
        min = min_val, max = 20, value = default_val, step = 0.5,
        width = "100%"
      )
    })

    # Our own LOD threshold reactive (no longer from mainParServer)
    lod_threshold_rv <- shiny::reactive({
      current_scan_type <- scan_type()
      default_threshold <- if (current_scan_type == "interactive") 10.5 else 7.5
      input[[ns_app_controller("LOD_thr")]] %||% default_threshold
    })

    # Reactive for selected chromosome (for zooming into specific chromosomes)
    selected_chromosome_rv <- shiny::reactive({
      input[[ns_app_controller("selected_chr")]] %||% "All" # Default to "All" chromosomes
    })

    # Reactive to find peak data for the selected trait
    peaks_data_for_trait <- shiny::reactive({
      trait_val <- trait_for_lod_scan_rv()
      if (is.null(trait_val)) {
        return(NULL)
      }

      shiny::req(main_selected_dataset_group(), import_reactives())

      dataset_val <- main_selected_dataset_group()

      message(paste("scanApp: Finding peaks for trait:", trait_val, "in dataset:", dataset_val))

      # Determine trait type for this dataset to pass to peak_finder
      trait_type_val <- get_trait_type(import_reactives(), dataset_val)

      # Use peak_finder to get peaks for this specific trait
      tryCatch(
        {
          peaks <- peak_finder(
            file_dir = import_reactives()$file_directory,
            selected_dataset = dataset_val,
            selected_trait = trait_val,
            trait_type = trait_type_val,
            cache_env = peaks_cache,
            use_cache = TRUE
          )
          message(paste("scanApp: Found", if (is.null(peaks)) 0 else nrow(peaks), "peaks for trait:", trait_val))
          peaks
        },
        error = function(e) {
          message(paste("scanApp: Error finding peaks for trait", trait_val, ":", e$message))
          NULL
        }
      )
    })

    # Reactive to get ALL peaks above threshold for the selected trait
    available_peaks_for_trait <- shiny::reactive({
      peaks <- peaks_data_for_trait()
      if (is.null(peaks) || nrow(peaks) == 0) {
        return(NULL)
      }

      # Filter by LOD threshold
      lod_thr <- lod_threshold_rv()
      shiny::req(lod_thr)

      # Check what columns are available
      message(paste("scanApp: Peak data columns:", paste(colnames(peaks), collapse = ", ")))

      # Use qtl_lod column (from peak_finder) instead of lod
      if (!"qtl_lod" %in% colnames(peaks)) {
        message("scanApp: No qtl_lod column found in peaks data")
        return(NULL)
      }

      # Filter by LOD threshold and sort by LOD descending
      filtered_peaks <- peaks[peaks$qtl_lod >= lod_thr, ]
      if (nrow(filtered_peaks) == 0) {
        message(paste("scanApp: No peaks above LOD threshold", lod_thr, "for trait:", trait_for_lod_scan_rv()))
        return(NULL)
      }

      # Sort by qtl_lod descending
      filtered_peaks <- filtered_peaks[order(-filtered_peaks$qtl_lod), ]

      message(paste("scanApp: Found", nrow(filtered_peaks), "peaks above threshold for trait:", trait_for_lod_scan_rv()))
      filtered_peaks
    })

    # Reactive value to store the selected peak marker from dropdown
    selected_peak_from_dropdown <- shiny::reactiveVal(NULL)

    # Reactive to get the currently selected peak marker (either from dropdown or default to highest)
    selected_peak_marker <- shiny::reactive({
      available_peaks <- available_peaks_for_trait()
      if (is.null(available_peaks) || nrow(available_peaks) == 0) {
        return(NULL)
      }

      # Use dropdown selection if available, otherwise default to highest peak
      dropdown_selection <- selected_peak_from_dropdown()
      if (!is.null(dropdown_selection) && dropdown_selection %in% available_peaks$marker) {
        selected_marker <- dropdown_selection
        selected_peak_info <- available_peaks[available_peaks$marker == selected_marker, ]
        message(paste("scanApp: Using dropdown-selected peak:", selected_marker, "with LOD:", selected_peak_info$qtl_lod[1]))
      } else {
        # Default to highest peak
        selected_marker <- available_peaks$marker[1]
        message(paste("scanApp: Using highest peak (default):", selected_marker, "with LOD:", available_peaks$qtl_lod[1]))
      }

      selected_marker
    })

    # Reactive to prepare allele effects data
    allele_effects_data <- shiny::reactive({
      peaks <- peaks_data_for_trait()
      marker <- selected_peak_marker()

      if (is.null(peaks) || is.null(marker)) {
        return(NULL)
      }

      # Use pivot_peaks helper function to reshape data for plotting
      reshaped_data <- pivot_peaks(peaks, marker)

      if (is.null(reshaped_data) || nrow(reshaped_data) == 0) {
        message(paste("scanApp: No allele effects data available for marker:", marker))
        return(NULL)
      }

      # Add trait name to the data for plot labeling
      reshaped_data$trait <- trait_for_lod_scan_rv()
      message(paste("scanApp: Prepared allele effects data for marker:", marker, "with", nrow(reshaped_data), "strain effects"))
      reshaped_data
    })

    file_index_dt <- shiny::reactive({
      shiny::req(import_reactives()$file_directory)
      dt <- data.table::as.data.table(import_reactives()$file_directory)
      shiny::validate(
        shiny::need("dataset_category" %in% names(dt), "Error: 'dataset_category' column missing in file_index.csv."),
        shiny::need("group" %in% names(dt), "Error: 'group' column missing in file_index.csv.")
      )
      return(dt)
    })

    shiny::observe({
      shiny::req(file_index_dt())
      categories <- unique(file_index_dt()$dataset_category)
      if (length(categories) > 0) {
        shiny::updateSelectInput(session, ns_app_controller("dataset_category_selector"),
          choices = stats::setNames(categories, categories),
          selected = categories[1]
        )
      } else {
        shiny::updateSelectInput(session, ns_app_controller("dataset_category_selector"),
          choices = c("No categories found" = ""), selected = ""
        )
      }
    })

    shiny::observe({
      shiny::req(file_index_dt(), input[[ns_app_controller("dataset_category_selector")]])
      selected_cat <- input[[ns_app_controller("dataset_category_selector")]]

      if (!is.null(selected_cat) && nzchar(selected_cat) && selected_cat != "No categories found") {
        datasets_in_category <- file_index_dt()[dataset_category == selected_cat, ]
        specific_datasets_choices <- unique(datasets_in_category$group)

        if (length(specific_datasets_choices) > 0) {
          shiny::updateSelectInput(session, ns_app_controller("specific_dataset_selector"),
            choices = stats::setNames(specific_datasets_choices, specific_datasets_choices),
            selected = specific_datasets_choices[1]
          )
        } else {
          shiny::updateSelectInput(session, ns_app_controller("specific_dataset_selector"),
            choices = c("No datasets in category" = ""), selected = ""
          )
        }
      } else {
        shiny::updateSelectInput(session, ns_app_controller("specific_dataset_selector"),
          choices = c("Select category first" = ""), selected = ""
        )
      }
    })

    main_selected_dataset_group <- shiny::reactive({
      shiny::req(input[[ns_app_controller("specific_dataset_selector")]])
      selected_group <- input[[ns_app_controller("specific_dataset_selector")]]
      if (is.null(selected_group) || !nzchar(selected_group) ||
        selected_group %in% c("Select category first", "No datasets in category")) {
        return(NULL)
      }
      message(paste("Main selected dataset group:", selected_group))
      return(selected_group)
    })

    selected_dataset_category_reactive <- shiny::reactive({
      shiny::req(main_selected_dataset_group(), file_index_dt())
      info <- file_index_dt()[group == main_selected_dataset_group()]
      if (nrow(info) > 0) {
        return(unique(info$dataset_category)[1])
      }
      return(NULL)
    })

    output[[ns_app_controller("plot_title")]] <- shiny::renderText({
      category <- selected_dataset_category_reactive()
      group <- main_selected_dataset_group()
      if (is.null(category) || is.null(group)) {
        return("Select Dataset Category and Specific Dataset")
      }

      plot_type_text <- "Plot"
      if (category %in% c("Liver Lipids", "Clinical Traits")) {
        plot_type_text <- "Manhattan Plot"
      } else if (category %in% c("Liver Genes", "Liver Isoforms")) {
        plot_type_text <- "Cis/Trans Plot"
      }
      return(paste0(plot_type_text, " for: ", group, " (", category, ")"))
    })

    output[[ns_app_controller("conditional_plot_ui")]] <- shiny::renderUI({
      category <- selected_dataset_category_reactive()
      shiny::req(category)

      # Ensure the module IDs are unique if using the same ns_app_controller
      if (category %in% c("Liver Lipids", "Clinical Traits")) {
        tagList(
          manhattanPlotUI(ns_app_controller("manhattan_plot_module"))
        )
      } else if (category %in% c("Liver Genes", "Liver Isoforms")) {
        tagList(
          cisTransPlotInput(ns_app_controller("cistrans_plot_module")),
          cisTransPlotUI(ns_app_controller("cistrans_plot_module"))
        )
      } else {
        shiny::p(paste("No specific plot type configured for category:", category))
      }
    })

    # Create our own main_par structure (no longer using mainParServer)
    active_main_par <- shiny::reactive({
      # Create a compatible structure with the expected reactives
      list(
        selected_dataset = main_selected_dataset_group, # Our dataset selection
        LOD_thr = lod_threshold_rv, # Our LOD threshold
        selected_chr = selected_chromosome_rv, # Selected chromosome
        which_trait = shiny::reactive(trait_for_lod_scan_rv()), # Currently searched trait
        dataset_category = selected_dataset_category_reactive # Our dataset category
      )
    })

    # Instantiate plot modules and capture their outputs
    manhattan_plot_outputs <- manhattanPlotServer(ns_app_controller("manhattan_plot_module"),
      import_reactives = import_reactives,
      main_par = active_main_par
    )

    cistrans_plot_outputs <- cisTransPlotServer(ns_app_controller("cistrans_plot_module"),
      import_reactives = import_reactives,
      main_par = active_main_par,
      peaks_cache = peaks_cache
    )

    # Reactive value to store the trait selected from either plot for LOD scanning
    trait_for_lod_scan_rv <- shiny::reactiveVal(NULL)

    # Observe clicks from Manhattan plot
    shiny::observeEvent(manhattan_plot_outputs$clicked_phenotype_for_lod_scan(),
      {
        clicked_trait <- manhattan_plot_outputs$clicked_phenotype_for_lod_scan()
        if (!is.null(clicked_trait)) {
          trait_for_lod_scan_rv(clicked_trait)
          message(paste("scanApp: Manhattan plot click detected. Trait for LOD scan:", clicked_trait))
          # Potentially switch to LOD scan tab/view here in the future
        }
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE
    ) # ignoreNULL=FALSE if we want to clear on background click

    # Observe clicks from Cis/Trans plot
    shiny::observeEvent(cistrans_plot_outputs$clicked_phenotype_for_lod_scan(),
      {
        clicked_trait <- cistrans_plot_outputs$clicked_phenotype_for_lod_scan()
        if (!is.null(clicked_trait)) {
          trait_for_lod_scan_rv(clicked_trait)
          message(paste("scanApp: Cis/Trans plot click detected. Trait for LOD scan:", clicked_trait))
          # Potentially switch to LOD scan tab/view here in the future
        }
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE
    )

    # When the main selected dataset changes, clear the trait_for_lod_scan_rv
    # to prevent a scan from an old selection on a new plot type.
    shiny::observeEvent(main_selected_dataset_group(),
      {
        message("scanApp: Main dataset group changed. Clearing trait_for_lod_scan_rv.")
        trait_for_lod_scan_rv(NULL) # Clear any active LOD scan
        # Clear the search input selection (choices will be updated by the other observer)
        updateSelectizeInput(session, ns_app_controller("trait_search_input"),
          selected = character(0)
        )
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE
    )

    # For debugging: observe the final trait selected for LOD scan
    shiny::observeEvent(trait_for_lod_scan_rv(),
      {
        message(paste("scanApp: trait_for_lod_scan_rv is now:", trait_for_lod_scan_rv()))
      },
      ignoreNULL = FALSE,
      ignoreInit = TRUE
    )

    # Clear caches when dataset changes (updated for new structure)
    shiny::observeEvent(main_selected_dataset_group(),
      {
        selected_ds_val <- main_selected_dataset_group()

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

    # Observer for the clear LOD scan button
    shiny::observeEvent(input[[ns_app_controller("clear_lod_scan_btn")]], {
      message("scanApp: Clear LOD scan button clicked.")
      trait_for_lod_scan_rv(NULL)
      # Also clear the search input field
      updateSelectizeInput(session, ns_app_controller("trait_search_input"),
        selected = character(0)
      )
    })

    # Call scanServer, passing the necessary reactives
    # scan_module_outputs will be a list of reactives/values returned by scanServer
    scan_module_outputs <- scanServer(
      id = ns_app_controller("scan_plot_module"),
      trait_to_scan = trait_for_lod_scan_rv, # Pass the reactive directly
      selected_dataset_group = main_selected_dataset_group, # Pass the reactive for the group
      import_reactives = import_reactives, # Pass the whole list of import reactives
      main_par_inputs = active_main_par # Pass the combined main parameters (for LOD_thr, etc.)
    )

    # Render the LOD scan click details table
    output[[ns_app_controller("lod_scan_click_table")]] <- DT::renderDT({
      # Check if we have scan module outputs and clicked point details
      if (!is.null(scan_module_outputs) && !is.null(scan_module_outputs$clicked_point_details)) {
        clicked_details <- scan_module_outputs$clicked_point_details()

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

    # UI for LOD Scan plot - only appears when a trait is selected
    output[[ns_app_controller("lod_scan_plot_ui_placeholder")]] <- shiny::renderUI({
      if (!is.null(trait_for_lod_scan_rv())) {
        bslib::card(
          id = "lod_scan_plot_card",
          bslib::card_header(paste("LOD Scan for Trait:", trait_for_lod_scan_rv())),
          bslib::card_body(
            # Chromosome selector
            div(
              style = "margin-bottom: 15px; background: #f8f9fa; padding: 10px; border-radius: 6px; border: 1px solid #bdc3c7;",
              div(
                style = "display: flex; align-items: center; gap: 15px; flex-wrap: wrap;",
                div(
                  style = "flex: 1; min-width: 200px;",
                  h6("🔍 Chromosome View", style = "color: #2c3e50; margin: 0 0 8px 0; font-weight: bold;"),
                  shiny::selectInput(
                    ns_app_controller("selected_chr"),
                    label = NULL,
                    choices = c(
                      "All" = "All",
                      setNames(as.character(1:19), paste("Chr", 1:19)),
                      "X" = "X", "Y" = "Y", "M" = "M"
                    ),
                    selected = "All",
                    width = "100%"
                  )
                ),
                div(
                  style = "display: flex; align-items: end; gap: 10px;",
                  shiny::actionButton(
                    ns_app_controller("zoom_to_chr"),
                    "🔍 Zoom to Chromosome",
                    class = "btn btn-primary",
                    style = "background: #3498db; border: none; color: white; font-weight: bold;"
                  ),
                  shiny::actionButton(
                    ns_app_controller("reset_chr_view"),
                    "🌐 Show All",
                    class = "btn btn-secondary",
                    style = "background: #7f8c8d; border: none; color: white;"
                  )
                )
              )
            ),
            # LOD scan plot
            scanOutput(ns_app_controller("scan_plot_module")),
            # Clicked point details table
            div(
              style = "margin-top: 15px;",
              h6("Click on plot to see point details:",
                style = "color: #2c3e50; margin-bottom: 10px; font-weight: bold;"
              ),
              DT::DTOutput(ns_app_controller("lod_scan_click_table"))
            ),
            # Conditional allele effects plot
            shiny::uiOutput(ns_app_controller("allele_effects_section"))
          )
        )
      } else {
        NULL # Don't show the card if no trait is selected for scanning
      }
    })

    # Render allele effects section conditionally
    output[[ns_app_controller("allele_effects_section")]] <- shiny::renderUI({
      # Check if we have allele effects data
      available_peaks <- available_peaks_for_trait()

      message(paste("scanApp: Checking allele effects section. Available peaks:", !is.null(available_peaks)))

      if (!is.null(available_peaks) && nrow(available_peaks) > 0) {
        # Create choices for dropdown - just show marker names
        peak_choices <- setNames(
          available_peaks$marker,
          available_peaks$marker
        )

        # Set default selection to highest peak if not already set
        current_selection <- selected_peak_from_dropdown()
        if (is.null(current_selection) || !current_selection %in% available_peaks$marker) {
          default_selection <- available_peaks$marker[1] # Highest peak
        } else {
          default_selection <- current_selection
        }

        # Get current peak info for display
        current_peak_info <- NULL
        if (!is.null(default_selection)) {
          current_peak_info <- available_peaks[available_peaks$marker == default_selection, ][1, ]
        }

        message(paste("scanApp: Rendering allele effects section with", length(peak_choices), "peak choices"))

        tagList(
          hr(style = "margin: 20px 0; border-top: 2px solid #3498db;"),
          div(
            style = "margin-bottom: 15px;",
            h5("Strain Effects",
              style = "color: #2c3e50; margin-bottom: 10px; font-weight: bold;"
            ),
            p(paste("Select a peak to view strain effects (", length(peak_choices), "peaks found):"),
              style = "color: #7f8c8d; margin-bottom: 10px; font-size: 12px;"
            ),
            shiny::selectInput(
              ns_app_controller("peak_selection_dropdown"),
              label = NULL,
              choices = peak_choices,
              selected = default_selection,
              width = "100%"
            ),
            # Add peak summary info
            if (!is.null(current_peak_info)) {
              div(
                id = ns_app_controller("peak_summary_info"),
                style = "background: #f8f9fa; padding: 10px; border-radius: 5px; margin: 10px 0; border-left: 4px solid #3498db;",
                shiny::uiOutput(ns_app_controller("peak_info_display"))
              )
            } else {
              NULL
            }
          ),
          shiny::plotOutput(ns_app_controller("allele_effects_plot_output"), height = "350px") %>%
            shinycssloaders::withSpinner(type = 8, color = "#3498db")
        )
      } else {
        message("scanApp: No allele effects section - no peaks above threshold")
        NULL
      }
    })

    # Dynamic peak info display
    output[[ns_app_controller("peak_info_display")]] <- shiny::renderUI({
      current_peaks <- available_peaks_for_trait()
      selected_marker <- selected_peak_marker()

      if (is.null(current_peaks) || is.null(selected_marker)) {
        return(NULL)
      }

      peak_info <- current_peaks[current_peaks$marker == selected_marker, ][1, ]
      if (nrow(peak_info) == 0) {
        return(NULL)
      }

      # Build summary info
      info_elements <- list()

      # Basic info
      info_elements <- c(info_elements, list(
        tags$strong("Marker: "), peak_info$marker, tags$br(),
        tags$strong("Position: "), paste0("Chr", peak_info$qtl_chr, ":", round(peak_info$qtl_pos, 2), " Mb"), tags$br(),
        tags$strong("LOD Score: "), round(peak_info$qtl_lod, 3), tags$br()
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
          tags$strong("Type: "),
          tags$span(cis_status, style = paste0("color: ", cis_color, "; font-weight: bold;")),
          tags$br()
        ))
      }

      # Confidence interval
      if ("qtl_ci_lo" %in% colnames(peak_info) && "qtl_ci_hi" %in% colnames(peak_info)) {
        if (!is.na(peak_info$qtl_ci_lo) && !is.na(peak_info$qtl_ci_hi)) {
          info_elements <- c(info_elements, list(
            tags$strong("95% CI: "),
            paste0(round(peak_info$qtl_ci_lo, 2), " - ", round(peak_info$qtl_ci_hi, 2), " Mb"),
            tags$br()
          ))
        }
      }

      # Add founder allele effects summary
      allele_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
      available_alleles <- allele_cols[allele_cols %in% colnames(peak_info)]

      if (length(available_alleles) > 0) {
        non_na_alleles <- available_alleles[!is.na(peak_info[available_alleles])]
        if (length(non_na_alleles) > 0) {
          info_elements <- c(info_elements, list(
            tags$strong("Founder Effects: "),
            paste0(length(non_na_alleles), " available (", paste(non_na_alleles, collapse = ", "), ")"),
            tags$br()
          ))
        }
      }

      do.call(tagList, info_elements)
    })

    # Observer to update selected peak when dropdown changes
    shiny::observeEvent(input[[ns_app_controller("peak_selection_dropdown")]],
      {
        new_selection <- input[[ns_app_controller("peak_selection_dropdown")]]
        if (!is.null(new_selection)) {
          selected_peak_from_dropdown(new_selection)
          message(paste("scanApp: Peak dropdown selection changed to:", new_selection))
        }
      },
      ignoreNULL = FALSE
    )

    # Observer to reset peak selection when trait changes
    shiny::observeEvent(trait_for_lod_scan_rv(),
      {
        selected_peak_from_dropdown(NULL) # Reset selection when trait changes
        message("scanApp: Reset peak selection due to trait change")
      },
      ignoreNULL = FALSE
    )

    # Render the allele effects plot
    output[[ns_app_controller("allele_effects_plot_output")]] <- shiny::renderPlot({
      effects_data <- allele_effects_data()

      message(paste("scanApp: Rendering allele effects plot. Data available:", !is.null(effects_data)))

      if (is.null(effects_data)) {
        message("scanApp: No allele effects data - showing placeholder")
        # Create a placeholder plot when no data
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::labs(title = "No strain effects data available for this peak") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, color = "#7f8c8d"))
      } else {
        message(paste("scanApp: Creating allele effects plot with", nrow(effects_data), "data points"))
        # Use the ggplot_alleles function to create the plot
        ggplot_alleles(effects_data)
      }
    })

    # --- FIXED Trait Search Logic ---
    observeEvent(input[[ns_app_controller("trait_search_button")]], {
      # FIX: Use main_selected_dataset_group() directly (it's a string), not main_selected_dataset_group()$group
      shiny::req(main_selected_dataset_group()) # Ensure a dataset is selected

      searched_trait <- input[[ns_app_controller("trait_search_input")]] # Get selected value from selectizeInput

      if (is.null(searched_trait) || !nzchar(searched_trait)) {
        shiny::showNotification("Please select a trait/gene to search.", type = "warning", duration = 3)
        return()
      }

      # Check if the searched trait is different from the current one to avoid re-triggering for no reason
      # Or if current is NULL, then definitely update.
      if (is.null(trait_for_lod_scan_rv()) || !identical(trait_for_lod_scan_rv(), searched_trait)) {
        message(paste(
          "scanApp: Trait search triggered. Trait for LOD scan set to:", searched_trait,
          "for dataset:", main_selected_dataset_group()
        ))
        trait_for_lod_scan_rv(searched_trait)

        # Show success notification
        shiny::showNotification(
          paste("Searching for trait:", searched_trait),
          type = "message",
          duration = 2
        )
      } else {
        message(paste("scanApp: Trait search for already selected trait:", searched_trait, "- no change."))
      }
    })

    # Store a flag to prevent auto-search immediately after dataset changes
    dataset_just_changed <- shiny::reactiveVal(FALSE)

    # Reset the flag after a short delay
    observe({
      if (dataset_just_changed()) {
        invalidateLater(1000) # Wait 1 second
        dataset_just_changed(FALSE)
      }
    })

    # --- AUTO-TRIGGER SEARCH WHEN TRAIT IS SELECTED FROM DROPDOWN ---
    # This allows users to skip clicking the button - just select from dropdown and it searches automatically
    # But NOT immediately after dataset changes to prevent slow loading
    observeEvent(input[[ns_app_controller("trait_search_input")]],
      {
        searched_trait <- input[[ns_app_controller("trait_search_input")]]

        # Only auto-search if:
        # 1. A valid trait is selected and dataset is available
        # 2. Dataset didn't just change (to prevent slow loading)
        # 3. This is a new selection by the user (not just preservation from dataset change)
        if (!is.null(searched_trait) && nzchar(searched_trait) &&
          !is.null(main_selected_dataset_group()) &&
          !dataset_just_changed()) {
          if (is.null(trait_for_lod_scan_rv()) || !identical(trait_for_lod_scan_rv(), searched_trait)) {
            message(paste("scanApp: Auto-search triggered for trait:", searched_trait))
            trait_for_lod_scan_rv(searched_trait)

            # Show notification for auto-search
            shiny::showNotification(
              paste("Auto-searching for trait:", searched_trait),
              type = "message",
              duration = 2
            )
          }
        } else if (dataset_just_changed()) {
          message(paste("scanApp: Trait selection preserved but auto-search skipped due to recent dataset change. Use search button to trigger LOD scan."))
        }
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE
    )

    # ====== TRAIT SEARCH DROPDOWN UPDATE LOGIC ======
    # Updates trait search choices when dataset changes (server-side selectize for performance)
    shiny::observeEvent(shiny::req(main_selected_dataset_group()), {
      shiny::req(import_reactives()) # Ensure import data is available
      current_ds <- main_selected_dataset_group()

      message(paste("scanApp: Updating trait search choices for dataset:", current_ds))

      # SET FLAG to prevent immediate auto-search after dataset change
      dataset_just_changed(TRUE)

      # PRESERVE CURRENT SELECTION: Store the current trait selection before updating
      current_trait_selection <- input[[ns_app_controller("trait_search_input")]]
      message(paste("scanApp: Current trait selection before dataset change:", current_trait_selection %||% "None"))

      # RESET TRAIT SEARCH INPUT: Prevent conflicts during update
      shiny::freezeReactiveValue(input, ns_app_controller("trait_search_input")) # Temporarily freeze reactive
      shiny::updateSelectizeInput(session, ns_app_controller("trait_search_input"),
        choices = character(0),
        selected = character(0),
        options = list(placeholder = "Loading traits...")
      )

      # GET NEW TRAIT CHOICES: Based on selected dataset
      choices <- get_trait_choices(import_reactives(), current_ds) # External helper function

      if (!is.null(choices) && length(choices) > 0) {
        message(paste("scanApp: Found", length(choices), "traits for dataset:", current_ds))

        # CHECK IF PREVIOUS SELECTION EXISTS IN NEW DATASET
        trait_to_select <- NULL
        if (!is.null(current_trait_selection) && nzchar(current_trait_selection)) {
          if (current_trait_selection %in% choices) {
            trait_to_select <- current_trait_selection
            message(paste("scanApp: Preserving trait selection:", current_trait_selection, "- found in new dataset"))
          } else {
            message(paste("scanApp: Previous trait selection:", current_trait_selection, "- not found in new dataset, clearing selection"))
          }
        }

        # UPDATE TRAIT SEARCH DROPDOWN: With server-side processing for large lists
        shiny::updateSelectizeInput(session, ns_app_controller("trait_search_input"),
          choices = choices,
          selected = trait_to_select, # Preserve selection if trait exists in new dataset
          options = list(
            placeholder = "Type to search traits/genes...",
            maxItems = 1, # Only allow single selection
            maxOptions = 10 # Limit displayed options for performance
          ),
          server = TRUE # Enable server-side processing
        )

        # Notify user if trait was preserved but auto-search was disabled
        if (!is.null(trait_to_select)) {
          shiny::showNotification(
            paste("Trait '", trait_to_select, "' preserved. Click 'Search & Plot LOD Scan' to analyze in new dataset."),
            type = "message",
            duration = 4
          )
        }
      } else {
        message(paste("scanApp: No traits found for dataset:", current_ds))
        # No traits available for this dataset
        shiny::updateSelectizeInput(session, ns_app_controller("trait_search_input"),
          choices = character(0),
          selected = character(0),
          options = list(placeholder = "No traits available for this dataset")
        )
      }
    })

    # Observer for the "Back to Overview Plot" button
    observeEvent(input[[ns_app_controller("clear_lod_scan_btn")]], {
      message("scanApp: clear_lod_scan_btn clicked. Clearing trait_for_lod_scan_rv.")
      trait_for_lod_scan_rv(NULL)
      # Also clear the search input field
      updateSelectizeInput(session, ns_app_controller("trait_search_input"),
        selected = character(0)
      )
    })

    # ====== CHROMOSOME ZOOM FUNCTIONALITY ======
    # Observer for "Zoom to Chromosome" button
    observeEvent(input[[ns_app_controller("zoom_to_chr")]], {
      selected_chr <- input[[ns_app_controller("selected_chr")]]
      if (!is.null(selected_chr) && selected_chr != "All") {
        message(paste("scanApp: Zooming to chromosome:", selected_chr))
        # The reactive selected_chromosome_rv will automatically pick up this change
        # and the plot will update via scan_table_chr reactive
        shiny::showNotification(
          paste("Zoomed to chromosome", selected_chr),
          type = "message",
          duration = 2
        )
      } else {
        shiny::showNotification(
          "Please select a specific chromosome to zoom to",
          type = "warning",
          duration = 3
        )
      }
    })

    # Observer for "Show All Chromosomes" button
    observeEvent(input[[ns_app_controller("reset_chr_view")]], {
      message("scanApp: Resetting to show all chromosomes")
      shiny::updateSelectInput(session, ns_app_controller("selected_chr"),
        selected = "All"
      )
      shiny::showNotification(
        "Showing all chromosomes",
        type = "message",
        duration = 2
      )
    })

    # ====== ADDITIONAL ANALYSES PLOT OUTPUTS ======
    # Profile Plot output
    output[[ns_app_controller("profile_plot_output")]] <- plotly::renderPlotly({
      plotly::plot_ly() %>%
        plotly::add_annotations(
          text = "Profile Plot Coming Soon",
          x = 0.5,
          y = 0.5,
          xref = "paper",
          yref = "paper",
          showarrow = FALSE,
          font = list(size = 20, color = "#2c3e50")
        ) %>%
        plotly::layout(
          title = "Profile Plot",
          showlegend = FALSE,
          xaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
          yaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
        )
    })

    # Correlation Plot output
    output[[ns_app_controller("correlation_plot_output")]] <- plotly::renderPlotly({
      plotly::plot_ly() %>%
        plotly::add_annotations(
          text = "Correlation Analysis Coming Soon",
          x = 0.5,
          y = 0.5,
          xref = "paper",
          yref = "paper",
          showarrow = FALSE,
          font = list(size = 20, color = "#2c3e50")
        ) %>%
        plotly::layout(
          title = "Correlation Analysis",
          showlegend = FALSE,
          xaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
          yaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
        )
    })
  }
  shiny::shinyApp(ui = ui, server = server)
}

# Server logic for scan related inputs and plot
scanServer <- function(id, trait_to_scan, selected_dataset_group, import_reactives, main_par_inputs) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Create local cache for peaks data
    local_peaks_cache <- new.env(parent = emptyenv())

    plot_width_rv <- shiny::reactiveVal(1200)
    plot_height_rv <- shiny::reactiveVal(600)
    use_alternating_colors_rv <- shiny::reactiveVal(TRUE)
    clicked_plotly_point_details_lod_scan_rv <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$plot_width,
      {
        plot_width_rv(input$plot_width)
      },
      ignoreNULL = TRUE
    )
    shiny::observeEvent(input$plot_height,
      {
        plot_height_rv(input$plot_height)
      },
      ignoreNULL = TRUE
    )

    shiny::observeEvent(input$preset_1to1, {
      shiny::updateNumericInput(session, "plot_width", value = 800)
      shiny::updateNumericInput(session, "plot_height", value = 800)
    })
    shiny::observeEvent(input$preset_3to2, {
      shiny::updateNumericInput(session, "plot_width", value = 900)
      shiny::updateNumericInput(session, "plot_height", value = 600)
    })
    shiny::observeEvent(input$preset_16to9, {
      shiny::updateNumericInput(session, "plot_width", value = 1280)
      shiny::updateNumericInput(session, "plot_height", value = 720)
    })

    shiny::observeEvent(input$color_toggle, {
      if (is.logical(input$color_toggle)) {
        use_alternating_colors_rv(input$color_toggle)
      } else {
        use_alternating_colors_rv(!use_alternating_colors_rv())
      }
    })

    current_trait_for_scan <- shiny::reactive({
      shiny::req(trait_to_scan())
      trait_val <- trait_to_scan()
      message(paste("scanServer: Received trait_to_scan value:", trait_val, "(this will be passed to trait_scan)."))
      return(trait_val)
    })

    scans <- shiny::reactive({
      req(current_trait_for_scan(), selected_dataset_group())
      trait_val <- current_trait_for_scan()
      dataset_group_val <- selected_dataset_group()

      message(paste0("scanServer (within scans reactive): ABOUT TO CALL trait_scan. Trait: '", trait_val, "', Dataset Group: '", dataset_group_val, "'"))

      file_dir_val <- import_reactives()$file_directory
      req(file_dir_val)
      if (!is.data.frame(file_dir_val) || !("group" %in% names(file_dir_val))) {
        stop("scanServer: import_reactives()$file_directory is not a valid data frame or missing 'group' column before calling trait_scan.")
      }

      result <- tryCatch(
        {
          trait_scan(
            file_dir = file_dir_val,
            selected_dataset = dataset_group_val,
            selected_trait = trait_val,
            cache_env = NULL
          )
        },
        error = function(e) {
          message(paste0("scanServer (within scans reactive): ERROR DURING trait_scan CALL. Trait: '", trait_val, "', Dataset Group: '", dataset_group_val, "'. Error: ", e$message))
          return(NULL)
        }
      )

      message(paste0("scanServer (within scans reactive): RETURNED FROM trait_scan. Trait: '", trait_val, "', Dataset Group: '", dataset_group_val, "'. Result class: ", class(result), ", Result nrows: ", if (!is.null(result) && (is.data.frame(result) || is.data.table(result))) nrow(result) else "N/A"))

      # Additional check: if result is NULL due to error, or if it's an empty data frame, handle appropriately.
      if (is.null(result) || ((is.data.frame(result) || is.data.table(result)) && nrow(result) == 0)) {
        message(paste0("scanServer: trait_scan returned NULL or empty for Trait: '", trait_val, "', Dataset: '", dataset_group_val, "'. Propagating as empty result."))
        # Return an empty data.table or data.frame as expected by downstream reactives to prevent crashes
        # Make sure it has the columns expected by QTL_plot_visualizer if possible, or handle this there.
        return(data.table::data.table())
      }

      result
    })

    scan_table <- shiny::reactive({
      # Access the list of reactives first
      main_par_list <- main_par_inputs()
      shiny::req(
        scans(),
        current_trait_for_scan(),
        main_par_list, # Ensure the list itself is available
        main_par_list$LOD_thr, # Ensure the reactive for LOD_thr is in the list
        main_par_list$LOD_thr(), # Ensure the value of LOD_thr is available
        import_reactives()$markers
      )
      QTL_plot_visualizer(scans(), current_trait_for_scan(), main_par_list$LOD_thr(), import_reactives()$markers)
    })

    scan_table_chr <- shiny::reactive({
      main_par_list <- main_par_inputs()
      shiny::req(
        scan_table(),
        main_par_list,
        main_par_list$selected_chr,
        main_par_list$selected_chr()
      )
      current_scan_table <- scan_table()
      selected_chromosome <- main_par_list$selected_chr()
      if (selected_chromosome == "All") {
        current_scan_table
      } else {
        sel_chr_num <- selected_chromosome
        if (selected_chromosome == "X") sel_chr_num <- 20
        if (selected_chromosome == "Y") sel_chr_num <- 21
        if (selected_chromosome == "M") sel_chr_num <- 22
        sel_chr_num <- as.numeric(sel_chr_num)

        dplyr::filter(current_scan_table, chr == sel_chr_num)
      }
    })

    current_scan_plot_gg <- shiny::reactive({
      main_par_list <- main_par_inputs()
      shiny::req(
        scan_table_chr(),
        main_par_list,
        main_par_list$LOD_thr,
        main_par_list$LOD_thr(),
        main_par_list$selected_chr,
        main_par_list$selected_chr()
      )
      ggplot_qtl_scan(scan_table_chr(), main_par_list$LOD_thr(), main_par_list$selected_chr())
    })

    output$scan_plot_ui_render <- shiny::renderUI({
      shiny::req(current_scan_plot_gg()) # This will now wait until ggplot_qtl_scan is successful
      plotly::plotlyOutput(ns("render_plotly_plot"),
        width = paste0(plot_width_rv(), "px"),
        height = paste0(plot_height_rv(), "px")
      ) |>
        shinycssloaders::withSpinner(type = 8, color = "#3498db")
    })

    output$render_plotly_plot <- plotly::renderPlotly({
      shiny::req(current_scan_plot_gg()) # Ensure ggplot object is ready
      plt <- plotly::ggplotly(current_scan_plot_gg(),
        source = ns("qtl_scan_plotly"),
        tooltip = c("x", "y", "chr")
      )

      plt <- plt %>% # Use the new pipe operator
        plotly::layout(
          dragmode = "zoom",
          hovermode = "closest",
          title = list(
            text = NULL
          ),
          xaxis = list(
            title = current_scan_plot_gg()$labels$x,
            fixedrange = FALSE # Allow x-axis zooming
          ),
          yaxis = list(
            title = current_scan_plot_gg()$labels$y,
            fixedrange = TRUE # Prevent y-axis zooming - always show full height
          )
        ) %>% # Use the new pipe operator
        plotly::config(
          displaylogo = FALSE,
          modeBarButtonsToRemove = c(
            "select2d", "lasso2d", "hoverClosestCartesian",
            "hoverCompareCartesian", "toggleSpikelines", "sendDataToCloud"
          )
        )

      # Register the plotly_click event here for the LOD scan plot
      plt <- plotly::event_register(plt, "plotly_click")
      plt
    })

    shiny::observeEvent(plotly::event_data("plotly_click", source = ns("qtl_scan_plotly")), {
      ev_data <- plotly::event_data("plotly_click", source = ns("qtl_scan_plotly"))

      if (!is.null(ev_data) && !is.null(ev_data$customdata) && length(ev_data$customdata) > 0) {
        clicked_marker_id <- ev_data$customdata[1]
        current_scan_data <- scan_table_chr()

        if (!"markers" %in% colnames(current_scan_data)) {
          warning("scanServer plotly_click: 'markers' column not found in scan_table_chr(). Cannot lookup clicked marker.")
          clicked_plotly_point_details_lod_scan_rv(data.frame(Info = "Marker column missing in scan data."))
          return()
        }

        # Get basic scan data
        selected_point_df_raw <- dplyr::filter(current_scan_data, markers == clicked_marker_id)

        if (nrow(selected_point_df_raw) > 0) {
          selected_point_df_raw <- selected_point_df_raw[1, , drop = FALSE]

          # Start with basic scan information
          basic_cols <- c("markers", "chr", "position", "LOD")
          basic_cols_present <- basic_cols[basic_cols %in% names(selected_point_df_raw)]

          if (length(basic_cols_present) > 0) {
            selected_point_df_processed <- dplyr::select(selected_point_df_raw, dplyr::all_of(basic_cols_present))
            if ("chr" %in% names(selected_point_df_processed)) {
              selected_point_df_processed$chr <- chr_XYM(selected_point_df_processed$chr)
            }
            selected_point_df_processed <- dplyr::mutate(selected_point_df_processed, dplyr::across(dplyr::where(is.numeric), \(x) signif(x, 4)))

            # Try to get additional peak information from peaks file
            tryCatch(
              {
                current_trait <- current_trait_for_scan()
                current_dataset <- selected_dataset_group()

                if (!is.null(current_trait) && !is.null(current_dataset)) {
                  # Get trait type for this dataset to pass to peak_finder
                  trait_type_val <- get_trait_type(import_reactives(), current_dataset)

                  # Get all peaks for this trait
                  trait_peaks <- peak_finder(
                    file_dir = import_reactives()$file_directory,
                    selected_dataset = current_dataset,
                    selected_trait = current_trait,
                    trait_type = trait_type_val,
                    cache_env = local_peaks_cache,
                    use_cache = TRUE
                  )

                  if (!is.null(trait_peaks) && nrow(trait_peaks) > 0) {
                    # Find the peak matching this marker
                    marker_peak <- dplyr::filter(trait_peaks, marker == clicked_marker_id)

                    if (nrow(marker_peak) > 0) {
                      marker_peak <- marker_peak[1, , drop = FALSE]

                      # Add cis/trans information if available
                      if ("cis" %in% colnames(marker_peak)) {
                        cis_status <- if (is.logical(marker_peak$cis)) {
                          ifelse(marker_peak$cis, "Cis", "Trans")
                        } else if (is.character(marker_peak$cis)) {
                          ifelse(toupper(marker_peak$cis) %in% c("TRUE", "1", "YES"), "Cis", "Trans")
                        } else {
                          "Unknown"
                        }
                        selected_point_df_processed$CisOrTrans <- cis_status
                      }

                      # Add confidence interval if available
                      if ("qtl_ci_lo" %in% colnames(marker_peak) && "qtl_ci_hi" %in% colnames(marker_peak)) {
                        ci_lo <- signif(marker_peak$qtl_ci_lo, 4)
                        ci_hi <- signif(marker_peak$qtl_ci_hi, 4)
                        selected_point_df_processed$CI_Range <- paste0("[", ci_lo, " - ", ci_hi, "]")
                      }

                      # Add founder allele effects A-H if available
                      allele_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
                      available_alleles <- allele_cols[allele_cols %in% colnames(marker_peak)]

                      # Debug: print all column names and A-H values
                      message("scanServer: Peak data columns: ", paste(colnames(marker_peak), collapse = ", "))
                      message("scanServer: Available A-H columns: ", paste(available_alleles, collapse = ", "))

                      if (length(available_alleles) > 0) {
                        for (allele in available_alleles) {
                          allele_value <- marker_peak[[allele]]
                          message(paste("scanServer: Column", allele, "value:", allele_value, "is.na:", is.na(allele_value)))
                          if (!is.na(allele_value)) {
                            selected_point_df_processed[[paste0("Founder_", allele)]] <- signif(allele_value, 4)
                          }
                        }
                      } else {
                        # Check if allele columns might have different names
                        potential_allele_cols <- grep("^[A-H]$|founder.*[A-H]|allele.*[A-H]", colnames(marker_peak), ignore.case = TRUE, value = TRUE)
                        message("scanServer: No standard A-H columns found. Potential allele columns: ", paste(potential_allele_cols, collapse = ", "))
                      }

                      message(paste("scanServer: Enhanced click details with peak info for marker:", clicked_marker_id))
                    } else {
                      message(paste("scanServer: No peak found matching marker:", clicked_marker_id))
                    }
                  } else {
                    message("scanServer: No peaks data available for trait:", current_trait)
                  }
                }
              },
              error = function(e) {
                message(paste("scanServer: Error getting peak info for clicked marker:", e$message))
              }
            )

            clicked_plotly_point_details_lod_scan_rv(selected_point_df_processed)
          } else {
            clicked_plotly_point_details_lod_scan_rv(data.frame(Info = "Selected point data columns not found."))
          }
        } else {
          clicked_plotly_point_details_lod_scan_rv(data.frame(Info = paste("Details for marker", clicked_marker_id, "not found.")))
        }
      } else {
        clicked_plotly_point_details_lod_scan_rv(NULL)
      }
    })

    output$plot_click_dt <- DT::renderDT({
      details <- clicked_plotly_point_details_lod_scan_rv()
      if (is.null(details) || nrow(details) == 0) {
        return(DT::datatable(data.frame(Info = "Click on the plot to see point details."),
          options = list(dom = "t"), rownames = FALSE
        ))
      }

      # Transpose the data for better display
      if (ncol(details) > 1) {
        # Create a two-column table with Property and Value
        transposed_details <- data.frame(
          Property = names(details),
          Value = as.character(unlist(details[1, ])),
          stringsAsFactors = FALSE
        )

        # Clean up column names for better display
        transposed_details$Property <- gsub("_", " ", transposed_details$Property)
        transposed_details$Property <- gsub("Founder ", "Founder ", transposed_details$Property)
        transposed_details$Property <- ifelse(transposed_details$Property == "markers", "Marker",
          ifelse(transposed_details$Property == "chr", "Chromosome",
            ifelse(transposed_details$Property == "position", "Position (Mb)",
              ifelse(transposed_details$Property == "CisOrTrans", "Cis/Trans",
                ifelse(transposed_details$Property == "CI Range", "Confidence Interval",
                  transposed_details$Property
                )
              )
            )
          )
        )

        return(DT::datatable(
          transposed_details,
          options = list(
            dom = "t",
            paging = FALSE,
            searching = FALSE,
            columnDefs = list(
              list(targets = 0, className = "dt-left", width = "40%"),
              list(targets = 1, className = "dt-left", width = "60%")
            )
          ),
          rownames = FALSE,
          selection = "none",
          colnames = c("Property", "Value"),
          class = "compact hover"
        ) %>%
          DT::formatStyle(
            "Property",
            fontWeight = "bold",
            backgroundColor = "#f8f9fa"
          ))
      } else {
        return(DT::datatable(details, options = list(dom = "t", paging = FALSE), rownames = FALSE, selection = "none"))
      }
    })

    output$download_qtl_plot_png <- shiny::downloadHandler(
      filename = function() {
        main_par_list <- main_par_inputs()
        trait_name <- current_trait_for_scan() %||% "plot" # Use the current scanned trait
        chr_suffix <- if (main_par_list$selected_chr() != "All") paste0("_chr", main_par_list$selected_chr()) else ""
        paste0("lod_plot_", trait_name, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".png")
      },
      content = function(file) {
        shiny::req(current_scan_plot_gg())
        ggplot2::ggsave(file,
          plot = current_scan_plot_gg(),
          width = plot_width_rv() / 96,
          height = plot_height_rv() / 96,
          dpi = 300, units = "in"
        )
      }
    )

    output$download_qtl_plot_pdf <- shiny::downloadHandler(
      filename = function() {
        main_par_list <- main_par_inputs()
        trait_name <- current_trait_for_scan() %||% "plot"
        chr_suffix <- if (main_par_list$selected_chr() != "All") paste0("_chr", main_par_list$selected_chr()) else ""
        paste0("lod_plot_", trait_name, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".pdf")
      },
      content = function(file) {
        shiny::req(current_scan_plot_gg())
        ggplot2::ggsave(file,
          plot = current_scan_plot_gg(),
          width = plot_width_rv() / 96,
          height = plot_height_rv() / 96,
          device = cairo_pdf, units = "in"
        )
      }
    )

    file_name_reactive <- shiny::reactive({
      main_par_list <- main_par_inputs()
      trait_name <- current_trait_for_scan() %||% "scan"
      instanceID <- trait_name
      selected_chromosome <- main_par_list$selected_chr()
      if (!is.null(selected_chromosome) && selected_chromosome != "All") {
        instanceID <- paste0(instanceID, "_chr", selected_chromosome)
      }
      paste("scan", instanceID, sep = "_")
    })

    # Return reactive values that might be useful for other modules (e.g., download peak info)
    # This is currently not used by the main app, but good practice for a module.
    return(shiny::reactiveValues(
      filename = file_name_reactive,
      tables = shiny::reactiveValues(scan = scan_table_chr),
      plots = shiny::reactiveValues(scan = current_scan_plot_gg),
      clicked_point_details = clicked_plotly_point_details_lod_scan_rv
    ))
  })
}

scanUI <- function(id) {
  ns <- shiny::NS(id)
  DT::DTOutput(ns("plot_click_dt"))
}

scanOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("scan_plot_ui_render"))
}
