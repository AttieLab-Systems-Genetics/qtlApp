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
#' @importFrom tags tagList
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
    
    create_title_panel(
      "QTL Scan Visualizer", 
      "Interactive visualization tool for QTL analysis"
    ),
    
    sidebar = bslib::sidebar("control_panel", 
      # Phase 1: New selection controls
      shiny::selectInput(shiny::NS("app_controller", "dataset_category_selector"),
                         "Select Dataset Category:",
                         choices = c("Loading..." = "")),
      shiny::selectInput(shiny::NS("app_controller", "specific_dataset_selector"),
                         "Select Specific Dataset:",
                         choices = c("Loading..." = "")), 
      hr(), # Separator
      # Existing controls (may be conditionally shown or modified later)
      mainParInput("main_par"),
      mainParUI("main_par"),
      peakInput("peak"), # Likely tied to a specific scan, for Phase 3
      peakUI("peak"),    # Likely tied to a specific scan, for Phase 3
      downloadInput("download"),
      downloadOutput("download"),
      # Action button to clear/go back from LOD scan view
      shiny::actionButton(shiny::NS("app_controller", "clear_lod_scan_btn"), "Back to Overview Plot", icon = shiny::icon("arrow-left")),
      # Trait search input
      hr(), # Add a horizontal rule for visual separation
      h4("Direct Trait/Gene LOD Scan"),
      textInput(shiny::NS("app_controller", "trait_search_input"), "Enter Trait/Gene ID:", placeholder = "e.g., Gapdh or PI_38_3"),
      actionButton(shiny::NS("app_controller", "trait_search_button"), "Plot LOD for Trait/Gene", icon = icon("search")),
      hr(), # Add another horizontal rule
    ),
    
    bslib::layout_columns(
      col_widths = bslib::breakpoints(
        sm = c(12, 12), # On small screens, stack them
        md = c(6, 6)    # On medium and larger, side-by-side if LOD scan is active
      ),
      bslib::card(
        id = "primary_plot_card",
        bslib::card_header(shiny::textOutput(shiny::NS("app_controller", "plot_title"))), 
        bslib::card_body(
          shiny::uiOutput(shiny::NS("app_controller", "conditional_plot_ui"))
        )
      ),
      # Conditional panel for the LOD scan plot
      shiny::uiOutput(shiny::NS("app_controller", "lod_scan_plot_ui_placeholder"))
    )
  )
  
  server <- function(input, output, session) {
      ns_app_controller <- shiny::NS("app_controller")

      trait_cache <- new.env(parent = emptyenv())
      peaks_cache <- new.env(parent = emptyenv())

      import_reactives <- importServer("import")
      
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
                                   selected = categories[1])
        } else {
          shiny::updateSelectInput(session, ns_app_controller("dataset_category_selector"),
                                   choices = c("No categories found" = ""), selected = "")
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
                                       selected = specific_datasets_choices[1])
            } else {
              shiny::updateSelectInput(session, ns_app_controller("specific_dataset_selector"), 
                                       choices = c("No datasets in category" = ""), selected = "")
            }
        } else {
            shiny::updateSelectInput(session, ns_app_controller("specific_dataset_selector"), 
                                     choices = c("Select category first" = ""), selected = "")
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
        if(nrow(info) > 0){
            return(unique(info$dataset_category)[1]) 
        }
        return(NULL)
     })

      output[[ns_app_controller("plot_title")]] <- shiny::renderText({
        category <- selected_dataset_category_reactive()
        group <- main_selected_dataset_group()
        if (is.null(category) || is.null(group)) return("Select Dataset Category and Specific Dataset")
        
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
      
      main_par_outputs <- mainParServer("main_par", import_reactives) 

      active_main_par <- shiny::reactive({
        # main_par_outputs is ALREADY the list of reactives returned by mainParServer
        
        params_list <- list()
        if(is.list(main_par_outputs)){ # Check if it's a list (it should be)
            for(name in names(main_par_outputs)){
                # Each element main_par_outputs[[name]] is itself a reactive expression
                params_list[[name]] <- main_par_outputs[[name]] 
            }
        }
        # Override or set selected_dataset with our main_selected_dataset_group reactive
        params_list$selected_dataset <- main_selected_dataset_group 
        
        return(params_list)
      })

      # Instantiate plot modules and capture their outputs
      manhattan_plot_outputs <- manhattanPlotServer(ns_app_controller("manhattan_plot_module"), 
                                                 import_reactives = import_reactives, 
                                                 main_par = active_main_par)

      cistrans_plot_outputs <- cisTransPlotServer(ns_app_controller("cistrans_plot_module"), 
                                                 import_reactives = import_reactives, 
                                                 main_par = active_main_par, 
                                                 peaks_cache = peaks_cache)

      # Reactive value to store the trait selected from either plot for LOD scanning
      trait_for_lod_scan_rv <- shiny::reactiveVal(NULL)

      # Observe clicks from Manhattan plot
      shiny::observeEvent(manhattan_plot_outputs$clicked_phenotype_for_lod_scan(), {
        clicked_trait <- manhattan_plot_outputs$clicked_phenotype_for_lod_scan()
        if (!is.null(clicked_trait)) {
          trait_for_lod_scan_rv(clicked_trait)
          message(paste("scanApp: Manhattan plot click detected. Trait for LOD scan:", clicked_trait))
          # Potentially switch to LOD scan tab/view here in the future
        }
      }, ignoreNULL = TRUE, ignoreInit = TRUE) # ignoreNULL=FALSE if we want to clear on background click

      # Observe clicks from Cis/Trans plot
      shiny::observeEvent(cistrans_plot_outputs$clicked_phenotype_for_lod_scan(), {
        clicked_trait <- cistrans_plot_outputs$clicked_phenotype_for_lod_scan()
        if (!is.null(clicked_trait)) {
          trait_for_lod_scan_rv(clicked_trait)
          message(paste("scanApp: Cis/Trans plot click detected. Trait for LOD scan:", clicked_trait))
          # Potentially switch to LOD scan tab/view here in the future
        }
      }, ignoreNULL = TRUE, ignoreInit = TRUE)

      # When the main selected dataset changes, clear the trait_for_lod_scan_rv
      # to prevent a scan from an old selection on a new plot type.
      shiny::observeEvent(main_selected_dataset_group(), {
        message(paste("scanApp: Main dataset changed to", main_selected_dataset_group(),"Clearing trait_for_lod_scan_rv."))
        trait_for_lod_scan_rv(NULL) 
      }, ignoreNULL = FALSE, ignoreInit = TRUE)

      # For debugging: observe the final trait selected for LOD scan
      shiny::observeEvent(trait_for_lod_scan_rv(), {
        message(paste("scanApp: trait_for_lod_scan_rv is now:", trait_for_lod_scan_rv()))
      }, ignoreNULL = FALSE, ignoreInit = TRUE)

      shiny::observeEvent(active_main_par()$selected_dataset(), {
        selected_ds_val <- active_main_par()$selected_dataset()
        
        if(!is.null(selected_ds_val)) { 
            message(paste("Clearing caches for dataset:", selected_ds_val))
            rm(list = ls(envir = trait_cache), envir = trait_cache)
            rm(list = ls(envir = peaks_cache), envir = peaks_cache)
        } else {
            message("No dataset selected or selected_dataset is NULL, not clearing caches.")
        }
      }, ignoreNULL = FALSE, ignoreInit = TRUE)
      
      # Observer for the clear LOD scan button
      shiny::observeEvent(input[[ns_app_controller("clear_lod_scan_btn")]], {
        message("scanApp: Clear LOD scan button clicked.")
        trait_for_lod_scan_rv(NULL)
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

      # UI for LOD Scan plot - only appears when a trait is selected
      output[[ns_app_controller("lod_scan_plot_ui_placeholder")]] <- shiny::renderUI({
        if(!is.null(trait_for_lod_scan_rv())){
          bslib::card(
            id = "lod_scan_plot_card",
            bslib::card_header(paste("LOD Scan for Trait:", trait_for_lod_scan_rv())),
            bslib::card_body(
              scanOutput(ns_app_controller("scan_plot_module"))
              # We can add peakUI and download UI here, linked to scan_module_outputs
            )
          )
        } else {
          NULL # Don't show the card if no trait is selected for scanning
        }
      })

      # --- Trait Search Logic ---
      observeEvent(input[[ns_app_controller("trait_search_button")]], {
        shiny::req(main_selected_dataset_group(), main_selected_dataset_group()$group) # Ensure a dataset is selected
        
        searched_trait <- trimws(input[[ns_app_controller("trait_search_input")]])
        
        if (!nzchar(searched_trait)) {
          shiny::showNotification("Please enter a trait/gene ID to search.", type = "warning")
          return()
        }
        
        # Check if the searched trait is different from the current one to avoid re-triggering for no reason
        # Or if current is NULL, then definitely update.
        if (is.null(trait_for_lod_scan_rv()) || !identical(trait_for_lod_scan_rv(), searched_trait)) {
            message(paste("scanApp: Trait search triggered. Trait for LOD scan set to:", searched_trait, 
                          "for dataset:", main_selected_dataset_group()$group))
            trait_for_lod_scan_rv(searched_trait)
        } else {
            message(paste("scanApp: Trait search for already selected trait:", searched_trait, "- no change."))
        }
      })

      # --- Reactive for selected dataset group ---
      main_selected_dataset_group <- reactive({
        shiny::req(input[[ns_app_controller("specific_dataset_selector")]])
        selected_group <- input[[ns_app_controller("specific_dataset_selector")]]
        if (is.null(selected_group) || !nzchar(selected_group) || 
            selected_group %in% c("Select category first", "No datasets in category")) {
          return(NULL)
        }
        message(paste("Main selected dataset group:", selected_group))
        return(selected_group)
      })

      # Observer to clear the LOD scan trait when the main dataset changes
      observeEvent(main_selected_dataset_group(), {
        message("scanApp: Main dataset group changed. Clearing trait_for_lod_scan_rv.")
        trait_for_lod_scan_rv(NULL) # Clear any active LOD scan
        # Also clear the search input field
        updateTextInput(session, ns_app_controller("trait_search_input"), value = "")
      }, ignoreNULL = TRUE, ignoreInit = TRUE)

      # Observer for the "Back to Overview Plot" button
      observeEvent(input[[ns_app_controller("clear_lod_scan_btn")]], {
        message("scanApp: clear_lod_scan_btn clicked. Clearing trait_for_lod_scan_rv.")
        trait_for_lod_scan_rv(NULL)
        # Also clear the search input field
        updateTextInput(session, ns_app_controller("trait_search_input"), value = "")
      })
  }
  shiny::shinyApp(ui = ui, server = server)
}

# Server logic for scan related inputs and plot
scanServer <- function(id, trait_to_scan, selected_dataset_group, import_reactives, main_par_inputs) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_width_rv <- shiny::reactiveVal(1200)
    plot_height_rv <- shiny::reactiveVal(600)
    use_alternating_colors_rv <- shiny::reactiveVal(TRUE)
    clicked_plotly_point_details_lod_scan_rv <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$plot_width, { plot_width_rv(input$plot_width) }, ignoreNULL = TRUE)
    shiny::observeEvent(input$plot_height, { plot_height_rv(input$plot_height) }, ignoreNULL = TRUE)

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
      if(is.logical(input$color_toggle)){
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
      if (!is.data.frame(file_dir_val) || !("group" %in% names(file_dir_val))){
          stop("scanServer: import_reactives()$file_directory is not a valid data frame or missing 'group' column before calling trait_scan.")
      }

      result <- tryCatch({
        trait_scan(
          file_dir = file_dir_val,
          selected_dataset = dataset_group_val, 
          selected_trait = trait_val,
          cache_env = NULL
        )
      }, error = function(e) {
        message(paste0("scanServer (within scans reactive): ERROR DURING trait_scan CALL. Trait: '", trait_val, "', Dataset Group: '", dataset_group_val, "'. Error: ", e$message))
        return(NULL) 
      })
      
      message(paste0("scanServer (within scans reactive): RETURNED FROM trait_scan. Trait: '", trait_val, "', Dataset Group: '", dataset_group_val, "'. Result class: ", class(result), ", Result nrows: ", if(!is.null(result) && (is.data.frame(result) || is.data.table(result))) nrow(result) else "N/A"))
      
      # Additional check: if result is NULL due to error, or if it's an empty data frame, handle appropriately.
      if (is.null(result) || ( (is.data.frame(result) || is.data.table(result)) && nrow(result) == 0) ) {
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
        if(selected_chromosome == "X") sel_chr_num <- 20
        if(selected_chromosome == "Y") sel_chr_num <- 21
        if(selected_chromosome == "M") sel_chr_num <- 22
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
                        height = paste0(plot_height_rv(), "px")) |>
        shinycssloaders::withSpinner(type = 8, color = "#3498db")
    })
    
    output$render_plotly_plot <- plotly::renderPlotly({
      shiny::req(current_scan_plot_gg()) # Ensure ggplot object is ready
      plt <- plotly::ggplotly(current_scan_plot_gg(), 
                              source = ns("qtl_scan_plotly"), 
                              tooltip = c("x", "y", "chr"))
      
      plt <- plt %>% # Use the new pipe operator
        plotly::layout(
          dragmode = "zoom",
          hovermode = "closest",
          title = list(
            text = NULL 
          ),
          xaxis = list(title = current_scan_plot_gg()$labels$x), 
          yaxis = list(title = current_scan_plot_gg()$labels$y)  
        ) %>% # Use the new pipe operator
        plotly::config(
          displaylogo = FALSE,
          modeBarButtonsToRemove = c("select2d", "lasso2d", "hoverClosestCartesian", 
                                     "hoverCompareCartesian", "toggleSpikelines", "sendDataToCloud")
        )
      
      # Register the plotly_click event here for the LOD scan plot
      plt <- plotly::event_register(plt, 'plotly_click')
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
        selected_point_df_raw <- dplyr::filter(current_scan_data, markers == clicked_marker_id)
        
        if(nrow(selected_point_df_raw) > 0){
          selected_point_df_raw <- selected_point_df_raw[1, , drop = FALSE]
          cols_to_select <- c("markers", "chr", "position", "LOD")
          cols_to_select_present <- cols_to_select[cols_to_select %in% names(selected_point_df_raw)]
          
          if(length(cols_to_select_present) > 0){
            selected_point_df_processed <- dplyr::select(selected_point_df_raw, dplyr::all_of(cols_to_select_present))
            if("chr" %in% names(selected_point_df_processed)){
                selected_point_df_processed$chr <- chr_XYM(selected_point_df_processed$chr)
            }
            selected_point_df_processed <- dplyr::mutate(selected_point_df_processed, dplyr::across(dplyr::where(is.numeric), \(x) signif(x, 4)))
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
                             options = list(dom = 't'), rownames = FALSE))
      }
      DT::datatable(details, options = list(dom = 't', paging = FALSE), rownames = FALSE, selection = 'none')
    })

    output$download_qtl_plot_png <- shiny::downloadHandler(
      filename = function() {
        main_par_list <- main_par_inputs()
        trait_name <- current_trait_for_scan() %||% "plot" # Use the current scanned trait
        chr_suffix <- if(main_par_list$selected_chr() != "All") paste0("_chr", main_par_list$selected_chr()) else ""
        paste0("lod_plot_", trait_name, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".png")
      },
      content = function(file) {
        shiny::req(current_scan_plot_gg())
        ggplot2::ggsave(file, plot = current_scan_plot_gg(), 
                        width = plot_width_rv()/96, 
                        height = plot_height_rv()/96, 
                        dpi = 300, units = "in")
      }
    )

    output$download_qtl_plot_pdf <- shiny::downloadHandler(
      filename = function() {
        main_par_list <- main_par_inputs()
        trait_name <- current_trait_for_scan() %||% "plot"
        chr_suffix <- if(main_par_list$selected_chr() != "All") paste0("_chr", main_par_list$selected_chr()) else ""
        paste0("lod_plot_", trait_name, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".pdf")
      },
      content = function(file) {
        shiny::req(current_scan_plot_gg())
        ggplot2::ggsave(file, plot = current_scan_plot_gg(), 
                        width = plot_width_rv()/96, 
                        height = plot_height_rv()/96, 
                        device = cairo_pdf, units = "in")
      }
    )
    
    file_name_reactive <- shiny::reactive({
      main_par_list <- main_par_inputs()
      trait_name <- current_trait_for_scan() %||% "scan"
      instanceID <- trait_name
      selected_chromosome <- main_par_list$selected_chr()
      if(!is.null(selected_chromosome) && selected_chromosome != "All") {
        instanceID <- paste0(instanceID, "_chr", selected_chromosome)
      }
      paste("scan", instanceID, sep = "_")
    })
    
    `%||%` <- function(a, b) if (!is.null(a)) a else b

    # Return reactive values that might be useful for other modules (e.g., download peak info)
    # This is currently not used by the main app, but good practice for a module.
    return(shiny::reactiveValues(
      filename = file_name_reactive,
      tables = shiny::reactiveValues(scan = scan_table_chr),
      plots  = shiny::reactiveValues(scan = current_scan_plot_gg)
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