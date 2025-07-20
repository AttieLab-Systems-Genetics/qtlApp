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
#' @importFrom stringr str_split str_remove str_to_title
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
    tags$head(
      # Viewport meta tag for responsive design
      tags$meta(name = "viewport", content = "width=device-width, initial-scale=1.0, user-scalable=yes"),

      # Autosizing CSS and JavaScript
      tags$link(rel = "stylesheet", type = "text/css", href = "autosize.css"),
      tags$script(src = "autosize.js"),

      # Custom app styles
      tags$style(custom_css)
    ),

    # Top navigation bar
    div(
      style = "background: linear-gradient(135deg, #2c3e50, #3498db); padding: 15px; margin-bottom: 20px; border-radius: 8px; text-align: center;",
      h3("QTL Scan Visualizer",
        style = "color: white; margin: 0; font-weight: bold;"
      )
    ),
    sidebar = bslib::sidebar(
      width = 600, # Increased sidebar width a bit more for better screen coverage

      # Dataset Category Selection - above all tabs
      div(
        style = "padding: 15px; margin-bottom: 20px; background: linear-gradient(135deg, #f8f9fa, #e9ecef); border-radius: 8px; border: 2px solid #3498db;",
        h4("Dataset Category", style = "color: #2c3e50; margin-bottom: 10px; font-weight: bold; text-align: center;"),
        shiny::selectInput(shiny::NS("app_controller", "dataset_category_selector"),
          NULL,
          choices = c("Loading..." = ""),
          width = "100%"
        ),
        p("Select the type of biological data to analyze",
          style = "font-size: 12px; color: #6c757d; margin: 5px 0 0 0; text-align: center;"
        )
      ),

      # Tabbed sidebar content
      bslib::navset_pill(
        id = "sidebar_tabs",
        selected = "LOD peaks", # Set default tab to LOD peaks

        # Tab 1: Overview Plot (Manhattan/Cis-Trans) - now the default tab
        bslib::nav_panel(
          "LOD peaks",
          div(
            style = "padding: 10px;",

            # LOD Threshold Control
            h5("Peak Filtering", style = "color: #2c3e50; margin-bottom: 10px; font-weight: bold;"),
            # Dynamic LOD threshold slider that updates based on scan type
            uiOutput(shiny::NS("app_controller", "lod_threshold_slider")),
            p("Filters peaks shown in the plot below",
              style = "font-size: 11px; color: #7f8c8d; margin: 5px 0 15px 0;"
            ),

            # Peak Selection Dropdown for future peak differences analysis
            hr(style = "border-top: 1px solid #bdc3c7; margin: 15px 0;"),
            h5("🎯 Peak Analysis", style = "color: #2c3e50; margin-bottom: 10px; font-weight: bold;"),
            shiny::uiOutput(shiny::NS("app_controller", "peak_selection_sidebar")),
            p("Control interaction analysis for the overview plots below",
              style = "font-size: 11px; color: #7f8c8d; margin: 5px 0 15px 0;"
            ),
            hr(style = "border-top: 1px solid #bdc3c7; margin: 15px 0;"),
            div(
              id = "overview-plot-container",
              class = "overview-plot-container",
              style = "height: 65vh; min-height: 400px; max-height: 800px; border: 1px solid #bdc3c7; border-radius: 5px; overflow: hidden;",
              shiny::uiOutput(shiny::NS("app_controller", "conditional_plot_ui"))
            ),
            p("Click on points to view detailed LOD scans. Plot titles show dataset and analysis type.",
              style = "font-size: 11px; color: #7f8c8d; margin: 10px 0 0 0; text-align: center;"
            )
          )
        ),

        # Tab 2: Data Search and Controls - now secondary
        bslib::nav_panel(
          "Data Search",
          div(
            style = "padding: 10px;",

            # Dataset selection
            h5("Dataset Selection", style = "color: #2c3e50; margin-bottom: 15px; font-weight: bold;"),

            # Conditional dataset selector - hidden for Genes category, shown for others
            shiny::uiOutput(shiny::NS("app_controller", "dataset_selection_ui")),

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
              id = "profile-plot-container",
              class = "sidebar-plot-container",
              style = "height: 50vh; min-height: 300px; max-height: 500px; border: 1px solid #bdc3c7; border-radius: 5px; overflow: hidden;",
              shinycssloaders::withSpinner(
                plotly::plotlyOutput(shiny::NS("app_controller", "profile_plot_output"),
                  height = "100%", width = "100%"
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
              id = "correlation-plot-container",
              class = "sidebar-plot-container",
              style = "height: 50vh; min-height: 300px; max-height: 500px; border: 1px solid #bdc3c7; border-radius: 5px; overflow: hidden;",
              shinycssloaders::withSpinner(
                plotly::plotlyOutput(shiny::NS("app_controller", "correlation_plot_output"),
                  height = "100%", width = "100%"
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

    # Simplified main area - now with a dynamic header
    bslib::card(
      id = "lod_scan_card",
      bslib::card_header(shiny::uiOutput(shiny::NS("app_controller", "main_plot_title"))),
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

    # Helper function to map HC_HF dataset names to their interactive versions
    # Only maps to datasets that actually exist in the file system
    get_interactive_dataset_name <- function(base_dataset, interaction_type) {
      if (is.null(interaction_type) || interaction_type == "none") {
        return(base_dataset)
      }

      # HC_HF Liver Genes (supports both Sex and Diet interactions)
      if (grepl("HC_HF Liver Genes", base_dataset, ignore.case = TRUE)) {
        if (interaction_type == "sex") {
          return("HC_HF Liver Genes, interactive (Sex)")
        } else if (interaction_type == "diet") {
          return("HC_HF Liver Genes, interactive (Diet)")
        }
      }
      # HC_HF Liver Lipids (only supports Diet interaction)
      else if (grepl("HC_HF.*Liver.*Lipid", base_dataset, ignore.case = TRUE)) {
        if (interaction_type == "diet") {
          return("HC_HF Liver Lipids, interactive (Diet)")
        }
        # No Sex interaction available for Liver Lipids - return original
      }
      # HC_HF Clinical Traits (supports both Sex and Diet interactions)
      else if (grepl("HC_HF.*Clinical", base_dataset, ignore.case = TRUE)) {
        if (interaction_type == "sex") {
          return("HC_HF Systemic Clinical Traits, interactive (Sex)")
        } else if (interaction_type == "diet") {
          return("HC_HF Systemic Clinical Traits, interactive (Diet)")
        }
      } else if (grepl("HC_HF.*Plasma.*Metabol", base_dataset, ignore.case = TRUE)) {
        if (interaction_type == "diet") {
          return("HC_HF Plasma plasma_metabolite, interactive (Diet)")
        }
        # No Sex interaction available for Plasma Metabolites - return original
      }

      # Fallback to original dataset if no mapping found
      return(base_dataset)
    }

    trait_cache <- new.env(parent = emptyenv())
    peaks_cache <- new.env(parent = emptyenv())

    import_reactives <- importServer("import")

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
    }) %>% shiny::debounce(300) # Debounce LOD threshold to prevent rapid re-firing

    # Reactive for selected chromosome (for zooming into specific chromosomes)
    selected_chromosome_rv <- shiny::reactive({
      input[[ns_app_controller("selected_chr")]] %||% "All" # Default to "All" chromosomes
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

    # Conditional dataset selection UI - show info for HC_HF categories, selector for others
    output[[ns_app_controller("dataset_selection_ui")]] <- shiny::renderUI({
      selected_cat <- input[[ns_app_controller("dataset_category_selector")]]

      if (is.null(selected_cat) || !nzchar(selected_cat)) {
        return(div(
          style = "padding: 15px; text-align: center; color: #7f8c8d; background: #f8f9fa; border-radius: 5px; border: 1px solid #e9ecef;",
          p("Select a dataset category above", style = "margin: 0; font-style: italic;")
        ))
      }

      # Show information panel for categories that have HC_HF auto-selection
      if (selected_cat %in% c("Liver Genes", "Liver Lipids", "Clinical Traits", "Plasma Metabolites", "Liver Isoforms")) {
        # Determine the dataset name and interaction info based on category
        dataset_info <- switch(selected_cat,
          "Liver Genes" = list(
            name = "HC_HF Liver Genes (Additive)",
            interaction_note = "Use interaction controls for Sex/Diet effects."
          ),
          "Liver Lipids" = list(
            name = "HC_HF Liver Lipids (Additive)",
            interaction_note = "Use interaction controls for Diet effects."
          ),
          "Clinical Traits" = list(
            name = "HC_HF Clinical Traits (Additive)",
            interaction_note = "Use interaction controls for Sex/Diet effects."
          ),
          "Plasma Metabolites" = list(
            name = "HC_HF Plasma Metabolites (Additive)",
            interaction_note = "Use interaction controls for Diet effects."
          ),
          "Liver Isoforms" = list(
            name = "HC_HF Liver Isoforms",
            interaction_note = "No interactive analysis available for this dataset."
          )
        )

        return(div(
          style = "padding: 15px; background: #e8f5e8; border-radius: 5px; border-left: 4px solid #28a745;",
          div(
            style = "display: flex; align-items: center; gap: 10px;",
            span("✓", style = "color: #28a745; font-weight: bold; font-size: 16px;"),
            div(
              h6(dataset_info$name, style = "color: #155724; margin: 0; font-weight: bold;"),
              p(paste("Auto-selected for streamlined analysis.", dataset_info$interaction_note),
                style = "color: #155724; margin: 5px 0 0 0; font-size: 12px;"
              )
            )
          )
        ))
      } else {
        # For other categories, show normal selector
        return(shiny::selectInput(
          ns_app_controller("specific_dataset_selector"),
          "Select Specific Dataset:",
          choices = c("Loading..." = ""),
          width = "100%"
        ))
      }
    })

    shiny::observe({
      shiny::req(file_index_dt(), input[[ns_app_controller("dataset_category_selector")]])
      selected_cat <- input[[ns_app_controller("dataset_category_selector")]]

      if (!is.null(selected_cat) && nzchar(selected_cat) && selected_cat != "No categories found") {
        datasets_in_category <- file_index_dt()[dataset_category == selected_cat, ]
        specific_datasets_choices <- unique(datasets_in_category$group)

        # Auto-select HC_HF datasets for all categories that have them
        hc_hf_dataset <- NULL

        if (selected_cat == "Liver Genes") {
          # Look for HC_HF Liver Genes dataset (additive version)
          hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Liver.*Genes", specific_datasets_choices, ignore.case = TRUE) &
            !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
        } else if (selected_cat == "Liver Lipids") {
          # Look for HC_HF Liver Lipids dataset (additive version)
          hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Liver.*Lipid", specific_datasets_choices, ignore.case = TRUE) &
            !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
        } else if (selected_cat == "Clinical Traits") {
          # Look for HC_HF Clinical Traits dataset (additive version)
          hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Clinical", specific_datasets_choices, ignore.case = TRUE) &
            !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
        } else if (selected_cat == "Plasma Metabolites") {
          # Look for HC_HF Plasma Metabolites dataset
          hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Plasma.*Metabol", specific_datasets_choices, ignore.case = TRUE) &
            !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
        } else if (selected_cat == "Liver Isoforms") {
          # Look for HC_HF Liver Isoforms dataset
          hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Liver.*Isoform", specific_datasets_choices, ignore.case = TRUE) &
            !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
        }

        if (!is.null(hc_hf_dataset) && length(hc_hf_dataset) > 0) {
          message(paste("Auto-selecting HC_HF dataset for", selected_cat, ":", hc_hf_dataset[1]))
          # Set choices to just the HC_HF dataset (hide the dropdown essentially)
          shiny::updateSelectInput(session, ns_app_controller("specific_dataset_selector"),
            choices = stats::setNames(hc_hf_dataset[1], hc_hf_dataset[1]),
            selected = hc_hf_dataset[1]
          )
        } else {
          # Fallback to normal behavior if HC_HF dataset not found
          message(paste("Warning: No HC_HF dataset found for", selected_cat, ", showing all options"))
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
        }
      } else {
        shiny::updateSelectInput(session, ns_app_controller("specific_dataset_selector"),
          choices = c("Select category first" = ""), selected = ""
        )
      }
    })

    main_selected_dataset_group <- shiny::reactive({
      selected_cat <- input[[ns_app_controller("dataset_category_selector")]]

      # Auto-select HC_HF datasets for categories that have them
      if (!is.null(selected_cat) && selected_cat %in% c("Liver Genes", "Liver Lipids", "Clinical Traits", "Plasma Metabolites", "Liver Isoforms")) {
        shiny::req(file_index_dt())
        datasets_in_category <- file_index_dt()[dataset_category == selected_cat, ]
        specific_datasets_choices <- unique(datasets_in_category$group)

        # Find the appropriate HC_HF dataset based on category
        hc_hf_dataset <- NULL

        if (selected_cat == "Liver Genes") {
          hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Liver.*Genes", specific_datasets_choices, ignore.case = TRUE) &
            !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
        } else if (selected_cat == "Liver Lipids") {
          hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Liver.*Lipid", specific_datasets_choices, ignore.case = TRUE) &
            !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
        } else if (selected_cat == "Clinical Traits") {
          hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Clinical", specific_datasets_choices, ignore.case = TRUE) &
            !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
        } else if (selected_cat == "Plasma Metabolites") {
          hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Plasma.*Metabol", specific_datasets_choices, ignore.case = TRUE) &
            !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
        } else if (selected_cat == "Liver Isoforms") {
          hc_hf_dataset <- specific_datasets_choices[grepl("^HC_HF.*Liver.*Isoform", specific_datasets_choices, ignore.case = TRUE) &
            !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
        }

        if (!is.null(hc_hf_dataset) && length(hc_hf_dataset) > 0) {
          message(paste("Auto-selected dataset for", selected_cat, "category:", hc_hf_dataset[1]))
          return(hc_hf_dataset[1])
        } else {
          message(paste("Warning: No HC_HF dataset found in", selected_cat, "category"))
          return(NULL)
        }
      }

      # Normal handling for other categories
      shiny::req(input[[ns_app_controller("specific_dataset_selector")]])
      selected_group <- input[[ns_app_controller("specific_dataset_selector")]]

      # Remove the forced dependency on interaction_type_rv to prevent circular reactive chain
      # Instead, we'll handle interaction type mapping in a separate reactive

      if (is.null(selected_group) || !nzchar(selected_group) ||
        selected_group %in% c("Select category first", "No datasets in category")) {
        return(NULL)
      }

      message(paste("Main selected dataset group:", selected_group))
      return(selected_group)
    })

    # New reactive that handles the dataset name mapping for interactive analysis
    mapped_dataset_for_interaction <- shiny::reactive({
      base_dataset <- main_selected_dataset_group()
      interaction_type <- interaction_type_rv()

      if (is.null(base_dataset)) {
        return(NULL)
      }

      # Check if this is an HC_HF dataset that supports interactive analysis
      is_hc_hf_dataset <- grepl("^HC_HF", base_dataset, ignore.case = TRUE)

      if (is_hc_hf_dataset && !is.null(interaction_type) && interaction_type != "none") {
        message(paste("mapped_dataset_for_interaction: HC_HF dataset detected:", base_dataset, "interaction_type is:", interaction_type))

        # Use helper function to get the appropriate dataset name
        interactive_dataset <- get_interactive_dataset_name(base_dataset, interaction_type)

        if (interactive_dataset != base_dataset) {
          message(paste("Interactive analysis mode: Using dataset", interactive_dataset, "for interaction type:", interaction_type))
          return(interactive_dataset)
        } else {
          message("mapped_dataset_for_interaction: Interaction type is none or null, using additive dataset")
        }
      }

      return(base_dataset)
    })

    # Detect if current dataset is additive or interactive based on dataset name and interaction type
    scan_type <- shiny::reactive({
      dataset_name <- main_selected_dataset_group()

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

    # Reactive to find peak data for the selected trait
    peaks_data_for_trait <- shiny::reactive({
      trait_val <- trait_for_lod_scan_rv()
      if (is.null(trait_val)) {
        return(NULL)
      }

      shiny::req(mapped_dataset_for_interaction(), import_reactives())

      dataset_val <- mapped_dataset_for_interaction()

      message(paste("scanApp: Finding peaks for trait:", trait_val, "in dataset:", dataset_val))

      # Determine trait type for this dataset to pass to peak_finder
      # Use the base dataset name for trait type determination
      base_dataset <- main_selected_dataset_group()
      trait_type_val <- get_trait_type(import_reactives(), base_dataset)

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

    # Reactive to prepare allele effects data, now driven by the selected peak
    allele_effects_data <- shiny::reactive({
      # This now correctly uses the output from the scanServer module
      peak_info <- scan_module_outputs$selected_peak()
      if (is.null(peak_info)) {
        return(NULL)
      }

      # Use pivot_peaks helper function to reshape data for plotting
      reshaped_data <- pivot_peaks(peak_info, peak_info$marker)

      if (is.null(reshaped_data) || nrow(reshaped_data) == 0) {
        message(paste("scanApp: No allele effects data available for marker:", peak_info$marker))
        return(NULL)
      }

      # Add trait name to the data for plot labeling
      reshaped_data$trait <- peak_info$trait
      message(paste("scanApp: Prepared allele effects data for marker:", peak_info$marker))
      reshaped_data
    })

    selected_dataset_category_reactive <- shiny::reactive({
      shiny::req(main_selected_dataset_group(), file_index_dt())
      info <- file_index_dt()[group == main_selected_dataset_group()]
      if (nrow(info) > 0) {
        return(unique(info$dataset_category)[1])
      }
      return(NULL)
    })





    output[[ns_app_controller("conditional_plot_ui")]] <- shiny::renderUI({
      category <- selected_dataset_category_reactive()
      shiny::req(category)

      # Ensure the module IDs are unique if using the same ns_app_controller
      if (category %in% c("Liver Lipids", "Clinical Traits", "Plasma Metabolites")) {
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
      main_par = active_main_par,
      sidebar_interaction_type = sidebar_interaction_type_rv
    )

    cistrans_plot_outputs <- cisTransPlotServer(ns_app_controller("cistrans_plot_module"),
      import_reactives = import_reactives,
      main_par = active_main_par,
      peaks_cache = peaks_cache,
      sidebar_interaction_type = sidebar_interaction_type_rv
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
    # BUT preserve trait selection when switching between HC_HF variants (additive <-> interactive)
    shiny::observeEvent(main_selected_dataset_group(),
      {
        current_dataset <- main_selected_dataset_group()
        message("scanApp: Main dataset group changed. Clearing trait_for_lod_scan_rv.")

        # Check if this is just switching between HC_HF dataset variants (any HC_HF dataset type)
        is_hc_hf_dataset_switch <- !is.null(current_dataset) &&
          grepl("^HC_HF", current_dataset, ignore.case = TRUE)

        if (!is_hc_hf_dataset_switch) {
          # Only clear trait selection for real dataset changes, not interaction type changes
          trait_for_lod_scan_rv(NULL) # Clear any active LOD scan
          # Clear the search input selection (choices will be updated by the other observer)
          updateSelectizeInput(session, ns_app_controller("trait_search_input"),
            selected = character(0)
          )
        } else {
          message("scanApp: HC_HF dataset switch detected - preserving trait selection")
        }
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

    # Interactive Analysis section - show for all HC_HF datasets
    output[[ns_app_controller("interactive_analysis_section")]] <- shiny::renderUI({
      dataset_group <- main_selected_dataset_group()

      # Show interactive analysis controls for all HC_HF datasets (Genes, Lipids, Clinical Traits, Metabolites)
      if (!is.null(dataset_group) && grepl("^HC_HF", dataset_group, ignore.case = TRUE)) {
        # Preserve the current selection when re-rendering
        current_selection <- current_interaction_type_rv()

        # Determine what interaction types are available for this dataset
        available_interactions <- c("None (Additive only)" = "none")

        # Check what interactions are actually available based on dataset type
        if (grepl("HC_HF Liver Genes", dataset_group, ignore.case = TRUE)) {
          available_interactions <- c(available_interactions,
            "Sex interaction" = "sex",
            "Diet interaction" = "diet"
          )
        } else if (grepl("HC_HF.*Liver.*Lipid", dataset_group, ignore.case = TRUE)) {
          available_interactions <- c(available_interactions,
            "Diet interaction" = "diet"
          )
          # No Sex interaction for Liver Lipids
        } else if (grepl("HC_HF.*Clinical", dataset_group, ignore.case = TRUE)) {
          available_interactions <- c(available_interactions,
            "Sex interaction" = "sex",
            "Diet interaction" = "diet"
          )
        } else if (grepl("HC_HF.*Plasma.*Metabol", dataset_group, ignore.case = TRUE)) {
          available_interactions <- c(available_interactions,
            "Diet interaction" = "diet"
          )
          # No Sex interaction for Plasma Metabolites
        }

        tagList(
          hr(style = "border-top: 2px solid #e74c3c; margin: 15px 0;"),
          h5("🧬 Interactive Analysis", style = "color: #2c3e50; margin-bottom: 15px; font-weight: bold;"),
          shiny::selectInput(
            ns_app_controller("interaction_type"),
            label = "Select interaction analysis:",
            choices = available_interactions,
            selected = if (current_selection %in% available_interactions) current_selection else "none",
            width = "100%"
          ),
          shiny::conditionalPanel(
            condition = paste0("input['", ns_app_controller("interaction_type"), "'] != 'none'"),
            div(
              style = "margin-top: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 5px; border-left: 4px solid #e74c3c;",
              p("ℹ️ Interactive analysis will show stacked plots: Interactive LOD scan (top) and Difference plot (Interactive - Additive, bottom).",
                style = "font-size: 12px; color: #6c757d; margin: 0;"
              )
            )
          )
        )
      } else {
        NULL # Don't show for other datasets
      }
    })

    # Reactive to store the selected interaction type
    interaction_type_rv <- shiny::reactive({
      interaction_value <- input[[ns_app_controller("interaction_type")]] %||% "none"
      message(paste("interaction_type_rv: Current value is:", interaction_value))
      return(interaction_value)
    })

    # Observer to update stored interaction type when input changes
    shiny::observeEvent(input[[ns_app_controller("interaction_type")]],
      {
        new_value <- input[[ns_app_controller("interaction_type")]]
        if (!is.null(new_value)) {
          current_interaction_type_rv(new_value)
          message(paste("Updated main UI interaction type to:", new_value))
        }
      },
      ignoreNULL = TRUE
    )

    # Observer to update sidebar interaction type when input changes (independent of main UI)
    shiny::observeEvent(input[[ns_app_controller("sidebar_interaction_type")]],
      {
        new_value <- input[[ns_app_controller("sidebar_interaction_type")]]
        if (!is.null(new_value)) {
          sidebar_interaction_type_rv(new_value)
          message(paste("Updated sidebar interaction type to:", new_value, "(for sidebar plots only)"))
        }
      },
      ignoreNULL = TRUE
    )

    # Observer to preserve sidebar interaction type during UI updates
    shiny::observeEvent(trait_for_lod_scan_rv(),
      {
        # When a new trait is selected, update the sidebar selectInput to preserve the current selection
        current_sidebar_type <- sidebar_interaction_type_rv()
        if (!is.null(current_sidebar_type) && current_sidebar_type != "none") {
          shiny::updateSelectInput(session, ns_app_controller("sidebar_interaction_type"),
            selected = current_sidebar_type
          )
          message(paste("Preserved sidebar interaction type:", current_sidebar_type, "during trait change"))
        }
      },
      ignoreNULL = FALSE,
      ignoreInit = TRUE
    )

    # NOTE: Removed the observer that re-triggered LOD scan when interaction type changes
    # The mapped_dataset_for_interaction reactive now handles dataset switching automatically

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
      selected_dataset_group = mapped_dataset_for_interaction, # Use mapped dataset for interactive analysis
      import_reactives = import_reactives, # Pass the whole list of import reactives
      main_par_inputs = active_main_par, # Pass the combined main parameters (for LOD_thr, etc.)
      interaction_type_reactive = interaction_type_rv # Pass the interaction type reactive
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

    # UI for LOD Scan plot - refactored for clarity and dynamic content
    output[[ns_app_controller("lod_scan_plot_ui_placeholder")]] <- shiny::renderUI({
      if (!is.null(trait_for_lod_scan_rv())) {
        # Logic for interactive analysis dropdown (moved from its own renderUI)
        dataset_group <- main_selected_dataset_group()
        interaction_analysis_ui <- NULL

        if (!is.null(dataset_group) && grepl("^HC_HF", dataset_group, ignore.case = TRUE)) {
          current_selection <- current_interaction_type_rv()
          available_interactions <- c("None (Additive only)" = "none")

          if (grepl("HC_HF Liver Genes", dataset_group, ignore.case = TRUE)) {
            available_interactions <- c(available_interactions, "Sex interaction" = "sex", "Diet interaction" = "diet")
          } else if (grepl("HC_HF.*Liver.*Lipid", dataset_group, ignore.case = TRUE)) {
            available_interactions <- c(available_interactions, "Diet interaction" = "diet")
          } else if (grepl("HC_HF.*Clinical", dataset_group, ignore.case = TRUE)) {
            available_interactions <- c(available_interactions, "Sex interaction" = "sex", "Diet interaction" = "diet")
          } else if (grepl("HC_HF.*Plasma.*Metabol", dataset_group, ignore.case = TRUE)) {
            available_interactions <- c(available_interactions, "Diet interaction" = "diet")
          }

          interaction_analysis_ui <- div(
            style = "flex: 1 1 180px; min-width: 180px;",
            shiny::selectInput(
              ns_app_controller("interaction_type"),
              label = "Interaction Analysis:",
              choices = available_interactions,
              selected = if (current_selection %in% available_interactions) current_selection else "none",
              width = "100%"
            )
          )
        }

        tagList(
          # Combined controls row
          div(
            style = "margin-bottom: 15px; background: #f8f9fa; padding: 10px 15px; border-radius: 4px; border: 1px solid #bdc3c7;",
            div(
              style = "display: flex; align-items: flex-end; gap: 15px; flex-wrap: wrap;",

              # Interaction analysis dropdown (now on the left)
              interaction_analysis_ui,

              # Chromosome selector
              div(
                style = "flex: 1 1 120px; min-width: 120px;",
                shiny::selectInput(
                  ns_app_controller("selected_chr"),
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
                style = "flex: 0 0 auto; display: flex; gap: 5px;",
                shiny::actionButton(
                  ns_app_controller("reset_chr_view"),
                  "🌐 Reset Zoom",
                  class = "btn btn-sm btn-secondary",
                  style = "background: #7f8c8d; border: none; color: white; font-size: 11px; padding: 4px 8px;"
                )
              )
            )
          ),

          # Conditional panel for interaction info
          if (!is.null(interaction_analysis_ui)) {
            shiny::conditionalPanel(
              condition = paste0("input['", ns_app_controller("interaction_type"), "'] != 'none'"),
              div(
                style = "margin-bottom: 15px; padding: 10px; background-color: #e8f4fd; border-radius: 5px; border-left: 4px solid #3498db;",
                p("ℹ️ Interactive analysis will show stacked plots: Interactive LOD scan (top) and Difference plot (Interactive - Additive, bottom).",
                  style = "font-size: 12px; color: #2c3e50; margin: 0;"
                )
              )
            )
          },

          # LOD scan plot
          scanOutput(ns_app_controller("scan_plot_module")),
          # Clicked point details table
          div(
            style = "margin-top: 15px;",
            DT::DTOutput(ns_app_controller("lod_scan_click_table"))
          ),
          # Conditional allele effects plot
          shiny::uiOutput(ns_app_controller("allele_effects_section"))
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

    # Dynamic title for the main LOD scan card
    output[[ns_app_controller("main_plot_title")]] <- shiny::renderUI({
      trait <- trait_for_lod_scan_rv()
      title_text <- if (!is.null(trait)) {
        paste("LOD Scan for Trait:", trait)
      } else {
        "LOD Scan - Detailed View"
      }
      h4(title_text, style = "font-weight: bold; margin-bottom: 0;")
    })

    # Render allele effects section conditionally
    output[[ns_app_controller("allele_effects_section")]] <- shiny::renderUI({
      # Check if we have allele effects data
      # Use the selected_peak from the module output as the source of truth
      current_peak_info <- scan_module_outputs$selected_peak()

      message(paste("scanApp: Checking allele effects section. Selected peak available:", !is.null(current_peak_info)))

      if (!is.null(current_peak_info)) {
        message(paste("scanApp: Rendering allele effects section for peak:", current_peak_info$marker))

        tagList(
          hr(style = "margin: 20px 0; border-top: 2px solid #3498db;"),
          div(
            style = "margin-bottom: 15px;",
            h5("Strain Effects",
              style = "color: #2c3e50; margin-bottom: 10px; font-weight: bold;"
            ),
            p("Showing strain effects for the selected peak. Click another peak on the LOD plot to update.",
              style = "color: #7f8c8d; margin-bottom: 10px; font-size: 12px;"
            )
          ),
          # Use layout_columns to place peak info and plot side-by-side
          bslib::layout_columns(
            col_widths = c(5, 7), # Allocate space for info and plot

            # Peak summary info
            div(
              id = ns_app_controller("peak_summary_info"),
              style = "background: #f8f9fa; padding: 10px; border-radius: 5px; border-left: 4px solid #3498db; height: 450px; overflow-y: auto;",
              shiny::uiOutput(ns_app_controller("peak_info_display"))
            ),

            # Allele effects plot (now square)
            shiny::plotOutput(ns_app_controller("allele_effects_plot_output"), height = "450px", width = "450px") %>%
              shinycssloaders::withSpinner(type = 8, color = "#3498db")
          )
        )
      } else {
        message("scanApp: No allele effects section - no peak selected or available")
        NULL
      }
    })

    # Dynamic peak info display
    output[[ns_app_controller("peak_info_display")]] <- shiny::renderUI({
      # This now correctly uses the output from the scanServer module
      peak_info <- scan_module_outputs$selected_peak()

      if (is.null(peak_info) || nrow(peak_info) == 0) {
        return(tags$div("No peak selected.", style = "color: #7f8c8d; text-align: center; padding-top: 20px;"))
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

      # Add founder allele effects with actual values
      allele_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
      strain_names <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
      available_alleles <- allele_cols[allele_cols %in% colnames(peak_info)]

      if (length(available_alleles) > 0) {
        # Get non-NA allele effects with their values
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
          # Create a more detailed display of founder effects
          info_elements <- c(info_elements, list(
            tags$strong("Founder Effects:"), tags$br(),
            # Display effects in a more readable format
            tags$div(
              style = "margin-left: 10px; font-family: monospace; font-size: 11px;",
              lapply(allele_effects, function(effect) {
                tags$div(effect, style = "margin: 2px 0;")
              })
            )
          ))
        } else {
          # Fallback to old display if no valid effects found
          info_elements <- c(info_elements, list(
            tags$strong("Founder Effects: "),
            paste0(length(available_alleles), " available (", paste(available_alleles, collapse = ", "), ")"),
            tags$br()
          ))
        }
      }

      do.call(tagList, info_elements)
    })

    # Observer to update selected peak when dropdown changes
    # shiny::observeEvent(input[[ns_app_controller("peak_selection_dropdown")]],
    #   {
    #     new_selection <- input[[ns_app_controller("peak_selection_dropdown")]]
    #     if (!is.null(new_selection)) {
    #       selected_peak_from_dropdown(new_selection)
    #       message(paste("scanApp: Peak dropdown selection changed to:", new_selection))
    #     }
    #   },
    #   ignoreNULL = FALSE
    # )

    # Observer to reset peak selection when trait changes
    # shiny::observeEvent(trait_for_lod_scan_rv(),
    #   {
    #     selected_peak_from_dropdown(NULL) # Reset selection when trait changes
    #     message("scanApp: Reset peak selection due to trait change")
    #   },
    #   ignoreNULL = FALSE
    # )

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
    # Observer for chromosome selection dropdown
    observeEvent(input[[ns_app_controller("selected_chr")]],
      {
        selected_chr <- input[[ns_app_controller("selected_chr")]]
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

    # Store the current interaction type to preserve across UI re-renders
    current_interaction_type_rv <- shiny::reactiveVal("none")

    # Store the sidebar interaction type separately (for independent sidebar plot control)
    sidebar_interaction_type_rv <- shiny::reactiveVal("none")

    # Observer to reset interaction types to default when dataset category changes
    # shiny::observeEvent(input[[ns_app_controller("dataset_category_selector")]],
    #   {
    #     message("scanApp: Dataset category changed. Resetting interaction analysis to default (additive).")
    #     current_interaction_type_rv("none")
    #     sidebar_interaction_type_rv("none")
    #   },
    #   ignoreInit = TRUE
    # )

    # Peak Analysis dropdown for sidebar - separate from main UI interaction controls
    output[[ns_app_controller("peak_selection_sidebar")]] <- shiny::renderUI({
      dataset_group <- main_selected_dataset_group()

      # Check if there's a current trait selected before attempting to get peaks
      current_trait <- trait_for_lod_scan_rv()
      available_peaks <- NULL

      # Only get peaks if trait is selected and use isolate to prevent reactive invalidation during tab switches
      if (!is.null(current_trait)) {
        available_peaks <- shiny::isolate(available_peaks_for_trait())
      }

      # Show interaction analysis controls for all HC_HF datasets (independent of main UI)
      if (!is.null(dataset_group) && grepl("^HC_HF", dataset_group, ignore.case = TRUE)) {
        # Determine what interaction types are available for this dataset
        available_interactions <- c("None (Additive only)" = "none")

        # Check what interactions are actually available based on dataset type
        if (grepl("HC_HF Liver Genes", dataset_group, ignore.case = TRUE)) {
          available_interactions <- c(available_interactions,
            "Sex interaction" = "sex",
            "Diet interaction" = "diet"
          )
        } else if (grepl("HC_HF.*Liver.*Lipid", dataset_group, ignore.case = TRUE)) {
          available_interactions <- c(available_interactions,
            "Diet interaction" = "diet"
          )
          # No Sex interaction for Liver Lipids
        } else if (grepl("HC_HF.*Clinical", dataset_group, ignore.case = TRUE)) {
          available_interactions <- c(available_interactions,
            "Sex interaction" = "sex",
            "Diet interaction" = "diet"
          )
        } else if (grepl("HC_HF.*Plasma.*Metabol", dataset_group, ignore.case = TRUE)) {
          available_interactions <- c(available_interactions,
            "Diet interaction" = "diet"
          )
          # No Sex interaction for Plasma Metabolites
        }

        # Remove peak selection dropdown - it's handled in the main LOD scan area

        tagList(
          h6("Sidebar Plot Analysis:", style = "color: #2c3e50; margin-bottom: 8px; font-weight: bold; font-size: 12px;"),
          shiny::selectInput(
            ns_app_controller("sidebar_interaction_type"),
            label = NULL,
            choices = available_interactions,
            selected = if (sidebar_interaction_type_rv() %in% available_interactions) sidebar_interaction_type_rv() else "none",
            width = "100%"
          ),
          div(
            style = "margin-top: 10px; padding: 8px; background-color: #e8f4fd; border-radius: 3px; border-left: 3px solid #3498db;",
            p("🔬 Independent control for sidebar Manhattan/Cis-Trans plots",
              style = "font-size: 10px; color: #2c3e50; margin: 0; font-style: italic;"
            )
          )
        )
      } else {
        div(
          style = "padding: 10px; text-align: center; color: #7f8c8d;",
          p("No peak analysis available", style = "margin: 0; font-size: 12px;"),
          p("Select an HC_HF dataset to enable analysis", style = "margin: 5px 0 0 0; font-size: 10px; font-style: italic;")
        )
      }
    })

    # Interactive Analysis section - show for all HC_HF datasets
  }
  shiny::shinyApp(ui = ui, server = server)
}

# Server logic for scan related inputs and plot
scanServer <- function(id, trait_to_scan, selected_dataset_group, import_reactives, main_par_inputs, interaction_type_reactive = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Create local cache for peaks data
    local_peaks_cache <- new.env(parent = emptyenv())

    plot_width_rv <- shiny::reactiveVal(1200)
    plot_height_rv <- shiny::reactiveVal(600)
    use_alternating_colors_rv <- shiny::reactiveVal(TRUE)
    clicked_plotly_point_details_lod_scan_rv <- shiny::reactiveVal(NULL)
    selected_peak_rv <- shiny::reactiveVal(NULL) # NEW: To store the full selected peak data

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

      # DEBUG: Show what files will be used for this dataset
      debug_files <- subset(file_dir_val, group == dataset_group_val & file_type == "scans")
      message(paste0("scanServer: Files that will be used for dataset '", dataset_group_val, "':"))
      for (i in 1:min(3, nrow(debug_files))) {
        message(paste0("  - Chr ", debug_files$ID_code[i], ": ", basename(debug_files$File_path[i])))
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
    }) %>% shiny::debounce(150) # Add debouncing to prevent rapid re-computation

    scan_table <- shiny::reactive({
      # Access the list of reactives first
      main_par_list <- main_par_inputs()
      shiny::req(
        scans(),
        current_trait_for_scan(),
        main_par_list,
        main_par_list$LOD_thr,
        main_par_list$LOD_thr(),
        import_reactives()$markers
      )

      # Streamlined processing - reduce debug messages
      scan_data <- scans()
      trait_val <- current_trait_for_scan()
      lod_val <- main_par_list$LOD_thr()
      markers_data <- import_reactives()$markers

      result <- QTL_plot_visualizer(scan_data, trait_val, lod_val, markers_data)
      result
    }) %>% shiny::debounce(150) # Add debouncing to prevent rapid re-computation

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
    }) %>% shiny::debounce(200) # Debounce scan table chr filtering

    # Observer to find the highest peak and set it as the default selected peak
    # This triggers whenever a new scan is successfully completed and loaded
    shiny::observeEvent(scan_table_chr(),
      {
        # Ensure we have a trait and scan data before proceeding
        shiny::req(current_trait_for_scan(), nrow(scan_table_chr()) > 0)

        trait_val <- current_trait_for_scan()
        dataset_group_val <- selected_dataset_group()

        message(paste("scanServer: New scan data available for trait", trait_val, ". Finding highest peak to set as default."))

        # It's safer to wrap the function call in a tryCatch in case helpers aren't loaded
        trait_type_val <- tryCatch(
          {
            get_trait_type(import_reactives(), dataset_group_val)
          },
          error = function(e) {
            warning(paste("Could not get trait type for default peak finding:", e$message))
            "genes" # Fallback to a sensible default
          }
        )

        # Get all peaks for the current trait
        all_peaks <- peak_finder(
          file_dir = import_reactives()$file_directory,
          selected_dataset = dataset_group_val,
          selected_trait = trait_val,
          trait_type = trait_type_val,
          cache_env = local_peaks_cache
        )

        if (!is.null(all_peaks) && nrow(all_peaks) > 0) {
          # Find the highest peak and set it as the selected one
          highest_peak <- all_peaks[which.max(all_peaks$qtl_lod), ]
          selected_peak_rv(highest_peak)
          message(paste("scanServer: Default peak set to highest LOD peak:", highest_peak$marker))
        } else {
          selected_peak_rv(NULL) # Clear selection if no peaks are found
          message("scanServer: No peaks found for this trait, clearing selected peak.")
        }
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE
    )


    # Reactive to store additive plot data for difference calculations
    additive_scan_data_rv <- shiny::reactiveVal(NULL)

    # Observer to clear the additive data cache whenever the selected trait changes
    shiny::observeEvent(current_trait_for_scan(),
      {
        # This ensures we don't use stale additive data from a previous trait
        message("scanServer: New trait selected, clearing additive data cache.")
        additive_scan_data_rv(NULL)
      },
      ignoreNULL = FALSE,
      ignoreInit = TRUE
    )


    # Reactive to automatically load additive data in background for any trait
    # This replaces the old background loader and observer pattern
    additive_scan_data_loader <- shiny::reactive({
      trait_val <- current_trait_for_scan()
      dataset_group_val <- selected_dataset_group()
      interaction_type <- if (!is.null(interaction_type_reactive)) interaction_type_reactive() else "none"

      # Only load if we're in interactive mode
      if (is.null(trait_val) || is.null(dataset_group_val) ||
        !grepl("interactive", dataset_group_val, ignore.case = TRUE) ||
        interaction_type == "none") {
        return(NULL)
      }

      # Derive the corresponding additive dataset name from the interactive one
      base_name <- gsub(",\\s*interactive\\s*\\([^)]+\\)", "", dataset_group_val)
      additive_dataset_name <- paste0(trimws(base_name), ", additive")

      message("scanServer: Automatically loading additive data for '", trait_val, "' from dataset '", additive_dataset_name, "'")

      tryCatch(
        {
          file_dir_val <- import_reactives()$file_directory
          req(file_dir_val)

          scan_data <- trait_scan(
            file_dir = file_dir_val,
            selected_dataset = additive_dataset_name,
            selected_trait = trait_val,
            cache_env = NULL
          )

          if (is.null(scan_data) || nrow(scan_data) == 0) {
            return(NULL)
          }

          main_par_list <- main_par_inputs()
          req(main_par_list$LOD_thr, import_reactives()$markers)

          processed_data <- QTL_plot_visualizer(scan_data, trait_val, main_par_list$LOD_thr(), import_reactives()$markers)

          if (!is.null(processed_data) && nrow(processed_data) > 0) {
            message("scanServer: Background loaded additive data with ", nrow(processed_data), " rows.")
            return(processed_data)
          }
          return(NULL)
        },
        error = function(e) {
          message("scanServer: Error automatically loading additive data: ", e$message)
          return(NULL)
        }
      )
    }) %>% shiny::debounce(250)

    # New observer to trigger the automatic loader only when necessary
    shiny::observe({
      interaction_type <- if (!is.null(interaction_type_reactive)) interaction_type_reactive() else "none"

      # If we are in interactive mode AND the cache is empty, trigger the loader to run and populate the cache.
      if (interaction_type != "none" && is.null(additive_scan_data_rv())) {
        loaded_data <- additive_scan_data_loader()
        if (!is.null(loaded_data)) {
          additive_scan_data_rv(loaded_data)
          message("scanServer: Additive data cache populated by automatic loader.")
        }
      }
    })

    # Observer to store additive data when interaction type is "none" (for manual workflow)
    shiny::observeEvent(scan_table_chr(),
      {
        interaction_type <- if (!is.null(interaction_type_reactive)) interaction_type_reactive() else "none"

        # Only store if this is a genuine additive scan, not an interactive one
        # This prevents overwriting the additive data when switching to an interactive view
        if (interaction_type == "none" && !grepl("interactive", selected_dataset_group(), ignore.case = TRUE)) {
          plot_data <- scan_table_chr()
          if (!is.null(plot_data) && nrow(plot_data) > 0) {
            additive_scan_data_rv(plot_data)
            message("scanServer: Stored additive scan data with ", nrow(plot_data), " rows (from 'none' mode).")
          }
        }
      },
      ignoreInit = TRUE,
      ignoreNULL = TRUE
    )


    current_scan_plot_gg <- shiny::reactive({
      # Get inputs with early validation
      scan_data <- scan_table_chr()
      main_par_list <- main_par_inputs()

      # Early returns for missing data
      if (is.null(scan_data) || nrow(scan_data) == 0 ||
        is.null(main_par_list) || is.null(main_par_list$LOD_thr) || is.null(main_par_list$selected_chr)) {
        return(NULL)
      }

      lod_threshold <- main_par_list$LOD_thr()
      selected_chr <- main_par_list$selected_chr()

      # Streamlined plot creation
      tryCatch(
        ggplot_qtl_scan(scan_data, lod_threshold, selected_chr),
        error = function(e) {
          message("scanServer: Error creating plot: ", e$message)
          return(NULL)
        }
      )
    }) %>% shiny::debounce(150) # Optimized debounce timing

    # Reactive to create difference plot (interactive - additive)
    difference_plot_gg <- shiny::reactive({
      if (is.null(interaction_type_reactive)) {
        return(NULL)
      }

      interaction_type <- interaction_type_reactive()

      # Only create difference plot for interactive datasets
      if (is.null(interaction_type) || interaction_type == "none") {
        return(NULL)
      }

      interactive_data <- scan_table_chr()

      # This now exclusively reads from the reactiveVal, which is populated by either the manual or automatic workflow.
      # This decouples it from the interactive scan's reactive chain.
      additive_data <- additive_scan_data_rv()

      # Early exit if data is missing
      if (is.null(interactive_data) || is.null(additive_data) ||
        nrow(interactive_data) == 0 || nrow(additive_data) == 0) {
        return(NULL)
      }

      # Optimized difference calculation
      tryCatch(
        {
          # Handle datasets with different numbers of rows by joining on markers
          if (nrow(interactive_data) != nrow(additive_data)) {
            # Check if both datasets have markers column
            if (!("markers" %in% colnames(interactive_data)) || !("markers" %in% colnames(additive_data))) {
              return(NULL)
            }

            # Efficient join using data.table if available
            if (requireNamespace("data.table", quietly = TRUE)) {
              # Use data.table for faster joins
              interactive_dt <- data.table::as.data.table(interactive_data)
              additive_dt <- data.table::as.data.table(additive_data)

              aligned_data <- merge(interactive_dt, additive_dt, by = "markers", suffixes = c("_int", "_add"))

              if (nrow(aligned_data) == 0) {
                return(NULL)
              }

              # Create difference efficiently
              diff_plot_data <- aligned_data[, .(
                markers = markers,
                chr = chr_int,
                position = position_int,
                BPcum = BPcum_int,
                LOD = LOD_int - LOD_add
              )]
            } else {
              # Fallback to dplyr
              aligned_data <- interactive_data %>%
                dplyr::inner_join(additive_data, by = "markers", suffix = c("_int", "_add"))

              if (nrow(aligned_data) == 0) {
                return(NULL)
              }

              diff_plot_data <- aligned_data %>%
                dplyr::mutate(LOD = LOD_int - LOD_add) %>%
                dplyr::select(markers,
                  chr = chr_int, position = position_int,
                  BPcum = BPcum_int, LOD
                )
            }
          } else {
            # Simple case - same number of rows
            diff_plot_data <- interactive_data
            diff_plot_data$LOD <- interactive_data$LOD - additive_data$LOD
          }

          # Handle NA values efficiently
          diff_plot_data$LOD[is.na(diff_plot_data$LOD)] <- 0

          # Create plot
          main_par_list <- main_par_inputs()
          selected_chr <- main_par_list$selected_chr()

          diff_plot <- ggplot_qtl_scan(diff_plot_data, -Inf, selected_chr)

          if (!is.null(diff_plot)) {
            diff_plot <- diff_plot + ggplot2::labs(title = paste("LOD Difference:", stringr::str_to_title(interaction_type), "- Additive"))
          }

          return(diff_plot)
        },
        error = function(e) {
          message("scanServer: Error creating difference plot: ", e$message)
          return(NULL)
        }
      )
    }) %>% shiny::debounce(300) # Optimized debounce timing

    # Simple UI state check
    show_stacked_plots <- shiny::reactive({
      if (is.null(interaction_type_reactive)) {
        return(FALSE)
      }
      interaction_type <- interaction_type_reactive()
      !is.null(interaction_type) && interaction_type != "none"
    })

    output$scan_plot_ui_render <- shiny::renderUI({
      # Minimal dependencies for faster rendering
      plot_gg <- current_scan_plot_gg()
      shiny::req(plot_gg)

      use_stacked <- show_stacked_plots()

      if (use_stacked) {
        # Calculate height for each plot (split the total height)
        individual_plot_height <- plot_height_rv() / 2

        shiny::tagList(
          # Interactive LOD plot (top)
          shiny::div(
            style = "margin-bottom: 10px;",
            shiny::h5("Interactive LOD Scan", style = "text-align: center; margin-bottom: 5px;"),
            plotly::plotlyOutput(ns("render_plotly_plot"),
              width = paste0(plot_width_rv(), "px"),
              height = paste0(individual_plot_height, "px")
            ) |>
              shinycssloaders::withSpinner(type = 8, color = "#3498db")
          ),
          # Difference plot (bottom)
          shiny::div(
            shiny::h5("LOD Difference (Interactive - Additive)", style = "text-align: center; margin-bottom: 5px;"),
            plotly::plotlyOutput(ns("render_difference_plot"),
              width = paste0(plot_width_rv(), "px"),
              height = paste0(individual_plot_height, "px")
            ) |>
              shinycssloaders::withSpinner(type = 8, color = "#e74c3c")
          )
        )
      } else {
        # Show only the main plot (additive or regular datasets)
        plotly::plotlyOutput(ns("render_plotly_plot"),
          width = paste0(plot_width_rv(), "px"),
          height = paste0(plot_height_rv(), "px")
        ) |>
          shinycssloaders::withSpinner(type = 8, color = "#3498db")
      }
    })

    output$render_plotly_plot <- plotly::renderPlotly({
      plot_gg <- current_scan_plot_gg()
      shiny::req(plot_gg)

      # Streamlined plotly creation
      plt <- plotly::ggplotly(plot_gg,
        source = ns("qtl_scan_plotly"),
        tooltip = c("x", "y", "chr")
      ) %>%
        plotly::layout(
          dragmode = "zoom",
          hovermode = "closest",
          title = list(text = NULL),
          xaxis = list(
            title = plot_gg$labels$x,
            fixedrange = FALSE
          ),
          yaxis = list(
            title = plot_gg$labels$y,
            fixedrange = TRUE
          )
        ) %>%
        plotly::config(
          displaylogo = FALSE,
          modeBarButtonsToRemove = c(
            "select2d", "lasso2d", "hoverClosestCartesian",
            "hoverCompareCartesian", "toggleSpikelines", "sendDataToCloud"
          )
        ) %>%
        plotly::event_register("plotly_click")

      plt
    })

    # Render function for the difference plot
    output$render_difference_plot <- plotly::renderPlotly({
      # Only render for interactive datasets
      shiny::req(show_stacked_plots())

      diff_plot_gg <- difference_plot_gg()

      if (is.null(diff_plot_gg)) {
        # Show a placeholder when no difference plot is available
        placeholder_plot <- ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::labs(title = "No difference plot available") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, color = "#7f8c8d"))

        return(plotly::ggplotly(placeholder_plot))
      }

      # Streamlined difference plot creation
      plotly::ggplotly(diff_plot_gg,
        source = ns("difference_plotly"),
        tooltip = c("x", "y", "chr")
      ) %>%
        plotly::layout(
          dragmode = "zoom",
          hovermode = "closest",
          title = list(text = NULL),
          xaxis = list(
            title = diff_plot_gg$labels$x,
            fixedrange = FALSE
          ),
          yaxis = list(
            title = diff_plot_gg$labels$y,
            fixedrange = TRUE
          )
        ) %>%
        plotly::config(
          displaylogo = FALSE,
          modeBarButtonsToRemove = c(
            "select2d", "lasso2d", "hoverClosestCartesian",
            "hoverCompareCartesian", "toggleSpikelines", "sendDataToCloud"
          )
        )
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

          # --- Start of new logic: Find closest peak ---
          # Get all peaks for this trait
          current_trait <- current_trait_for_scan()
          current_dataset <- selected_dataset_group()
          all_peaks <- NULL
          if (!is.null(current_trait) && !is.null(current_dataset)) {
            trait_type_val <- get_trait_type(import_reactives(), current_dataset)
            all_peaks <- peak_finder(
              file_dir = import_reactives()$file_directory,
              selected_dataset = current_dataset,
              selected_trait = current_trait,
              trait_type = trait_type_val,
              cache_env = local_peaks_cache,
              use_cache = TRUE
            )
          }

          if (!is.null(all_peaks) && nrow(all_peaks) > 0) {
            # Find the peak on the same chromosome that is closest in position
            clicked_chr <- selected_point_df_raw$chr
            clicked_pos <- selected_point_df_raw$position

            closest_peak <- all_peaks %>%
              dplyr::filter(qtl_chr == chr_XYM(clicked_chr)) %>%
              dplyr::mutate(pos_diff = abs(qtl_pos - clicked_pos)) %>%
              dplyr::filter(pos_diff == min(pos_diff))

            if (nrow(closest_peak) > 0) {
              selected_peak_rv(closest_peak[1, ])
              message(paste("scanServer: Click detected. Closest peak is:", closest_peak[1, ]$marker))
              # Use this closest peak for enhancing the details table
              marker_peak <- closest_peak[1, , drop = FALSE]
            } else {
              # Fallback if no peak is found on the same chromosome
              marker_peak <- NULL
            }
          } else {
            marker_peak <- NULL
          }
          # --- End of new logic ---

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
                # The logic to find the peak is now handled above,
                # we just need to use the 'marker_peak' data frame if it exists.
                if (!is.null(marker_peak) && nrow(marker_peak) > 0) {
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

                  message(paste("scanServer: Enhanced click details with peak info for marker:", marker_peak$marker))
                } else {
                  message("scanServer: No peak found for this click location.")
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
      clicked_point_details = clicked_plotly_point_details_lod_scan_rv,
      selected_peak = selected_peak_rv # NEW: Pass selected peak out
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
