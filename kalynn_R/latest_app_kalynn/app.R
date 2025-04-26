# Set maximum file upload size for Shiny (20GB)
# options(shiny.maxRequestSize = 20000*1024^2)

message("Starting minimal app.R...")
# Load required data files using the import function
# imported_data <- import_data()
# message("Data loading complete. Assigning variables...")

# Assign imported data to variables (optional, but mirrors original structure)
# file_directory <- imported_data$file_directory
# gene_symbols <- imported_data$gene_symbols
# chr_breaks <- imported_data$chr_breaks
# annotation_list <- imported_data$annotation_list
# markers <- imported_data$markers

# Load necessary libraries after data import function might create caches
# require("rstudioapi")  
# require("dplyr")      
# require("stringr")     
# require("tidyverse")   
# require("BiocManager") 
# require("ggplot2")     
# require("qtl2")        
# require("grid")        
# require("ggrepel")     
# require("gridGraphics")
# require("ggpubr")     
require("shiny")      
# require("shinyFiles") 
# require("bslib")      
# require("spsComps")   
# require("DT")         
# require("shinyjs")    
# require("shinycssloaders")
# require("data.table") 
# require("reshape2")   
# require("plotly")     
# require("ggiraph")    
# require("writexl")   
# require("fontawesome")
# require("debounce")
# require("fst") # Ensure fst is loaded

# Source modules (assuming they are in the R/ directory relative to app)
# Function to source all R files in a directory
# source_dir <- function(path, trace = TRUE, ...) {
#   for (nm in list.files(path, pattern = "\\.[RrSsQq]")) {
#     if(trace) cat(nm,":")
#     source(file.path(path, nm), ...)
#     if(trace) cat("\n")
#   }
# }
# source_dir("R")

# Remove cache environment creations here, they are in import_data or individual functions now
# file_path_cache <- new.env(parent = emptyenv())
# trait_cache <- new.env(parent = emptyenv())
# peaks_cache <- new.env(parent = emptyenv())

# Sort gene symbols
# gene_symbols <- sort(gene_symbols)
# message("Loaded ", length(gene_symbols), " gene symbols")

# Create a group identifier by combining diet, trait compartment, trait type, and scan type
# Include 'sexes' to differentiate datasets like male-only vs both sexes, but only if not 'Both'
# file_directory$group <- paste0(file_directory$diet, " ", file_directory$trait_compartment, " ",
#                            file_directory$trait_type, 
#                            # Conditionally add sex info if it's not 'Both'
#                            ifelse(file_directory$sexes == "Both", "", paste0(" (", file_directory$sexes, ")")),
#                            ", ", file_directory$scan_type,
#                            ifelse(file_directory$scan_type == "interactive",
#                                   paste0(" (", file_directory$covars_interactive, ")"),
#                                   ""))

# set microfunctions==========================================================

# set UI=======================================================================
ui <- fluidPage(
  # useShinyjs(),
  titlePanel("Minimal Test App"),
  mainPanel(
    h2("App Started!")
  )
  # Add custom CSS
#   tags$head(
#       tags$style(HTML("
#   .well {
#     background-color: #ffffff;
#     border: 1px solid #e3e3e3;
#     border-radius: 8px;
#     box-shadow: 0 2px 4px rgba(0,0,0,0.1);
#     margin-bottom: 20px;
#   }
#   /* Add lever switch styling */
#   .lever-switch {
#     position: relative;
#     display: inline-block;
#     width: 60px;
#     height: 30px;
#     background-color: #2c3e50;
#     border-radius: 15px;
#     cursor: pointer;
#     transition: all 0.3s ease;
#     overflow: hidden;
#   }
#   .lever-switch:before {
#     content: '';
#     position: absolute;
#     width: 26px;
#     height: 26px;
#     left: 2px;
#     bottom: 2px;
#     background-color: white;
#     border-radius: 50%;
#     transition: all 0.3s ease;
#   }
#   .lever-switch.active {
#     background-color: #3498db;
#   }
#   .lever-switch.active:before {
#     transform: translateX(30px);
#   }
#   .lever-switch:hover {
#     box-shadow: 0 0 5px rgba(0,0,0,0.2);
#   }
#   .control-label {
#     color: #2c3e50;
#     font-weight: 500;
#     font-size: 14px;
#   }
#   .form-control {
#     border-radius: 6px;
#     border: 1px solid #dce4ec;
#   }
#   .btn {
#     border-radius: 6px;
#     text-transform: uppercase;
#     font-size: 12px;
#     font-weight: 600;
#     padding: 8px 16px;
#   }
#   .btn-default {
#     background-color: #3498db;
#     color: white;
#     border: none;
#   }
#   .btn-default:hover {
#     background-color: #2980b9;
#     color: white;
#   }
#   .selectize-input {
#     border-radius: 6px;
#     border: 1px solid #dce4ec;
#     position: relative;
#     z-index: 1;
#   }
#   .selectize-dropdown {
#     z-index: 10000 !important;
#   }
#   /* Ensure the strain effects selectize appears above the spinner */
#   #which_peak + .selectize-control .selectize-dropdown {
#     z-index: 99999 !important;
#   }
#   /* Ensure all selectize dropdowns remain on top */
#   .selectize-control {
#     position: relative;
#     z-index: 10;
#   }
#   .selectize-dropdown {
#     z-index: 99999 !important;
#   }
#   /* Fix positioning of the spinner relative to the plot */
#   .shiny-spinner-output-container {
#     position: relative;
#     z-index: 0;
#   }
#   .title-panel {
#     background-color: #2c3e50;
#     color: white;
#     padding: 20px;
#     margin-bottom: 30px;
#     border-radius: 8px;
#   }
#   .plot-container {
#     background-color: white;
#     padding: 20px;
#     border-radius: 8px;
#     box-shadow: 0 2px 4px rgba(0,0,0,0.1);
#     margin-bottom: 20px;
#   }
#   .datatable-container {
#     background-color: white;
#     padding: 20px;
#     border-radius: 8px;
#     box-shadow: 0 2px 4px rgba(0,0,0,0.1);
#   }
#   #error_message {
#     color: #e74c3c;
#     padding: 10px;
#     margin-top: 10px;
#     font-size: 14px;
#   }
#    /* Add animation for hover effects */
#   .well:hover {
#     transform: translateY(-2px);
#     transition: transform 0.2s ease;
#   }
#    .btn {
#     transition: all 0.3s ease;
#   }
#    .btn:hover {
#     transform: translateY(-1px);
#     box-shadow: 0 4px 8px rgba(0,0,0,0.1);
#   }
#    /* Add pulsing effect to the spinner */
#   .spinner-border {
#     animation: spinner-pulse 1s ease infinite;
#   }
#    @keyframes spinner-pulse {
#     0% { transform: scale(1); }
#     50% { transform: scale(1.1); }
#     100% { transform: scale(1); }
#   }
#    /* Add tooltip styling */
#   .tooltip-inner {
#     background-color: #2c3e50;
#     color: white;
#     border-radius: 4px;
#     padding: 8px 12px;
#   }
#    /* Add success message styling */
#   .success-message {
#     color: #18bc9c;
#     padding: 10px;
#     margin-top: 10px;
#     font-size: 14px;
#     display: none;
#   }
# "))
#   ),
  # Use bslib theme with custom settings
#   theme = bs_theme(
#       version = 4,
#       bootswatch = "flatly",
#       primary = "#3498db",
#       secondary = "#2c3e50",
#       success = "#18bc9c",
#       info = "#3498db",
#       warning = "#f39c12",
#       danger = "#e74c3c"
#   ),
  # Title Panel with modern styling
#   div(class = "title-panel",
#       h1("Pre-scanned QTL Visualizer for Diet DO Study",
#           style = "font-size: 28px; margin: 0;"),
#       p("Interactive visualization tool for QTL analysis",
#           style = "margin: 10px 0 0 0; opacity: 0.8;")
#   ),
#   fluidRow(
#       # Left column with inputs and allele effects
#       column(3,
#           wellPanel(
#               style = "padding: 20px;",
#               mainParInput("main_par"),
#               mainParUI("main_par"),
#               div(id = "error_message_container", style = "margin-bottom: 15px; color: #e74c3c; font-weight: bold; display: none;",
#                   textOutput("error_message")
#               )
#           ),
#           # Allele effects panel - REMOVE THIS ENTIRE wellPanel
#           # wellPanel(
#           #     style = "padding: 20px; position: relative; overflow: visible;",
#           #     h4("Strain Effects", ...),
#           #     div(style = "color: #7f8c8d; ...",
#           #         "Select a peak to see strain effects..."
#           #     ),
#           #     div(style = "margin-bottom: 20px; ...",
#           #         selectizeInput("which_peak", ...)
#           #     ),
#           #     div(style = "margin-bottom: 20px;",
#           #         downloadButton("download_effects_plot_png", ...),
#           #         downloadButton("download_effects_plot_pdf", ...)
#           #     ),
#           #     div(style = "margin-top: 20px; ...",
#           #         plotOutput("allele_effects", ...) %>%
#           #             withSpinner(...)
#           #     )
#           # )
#           # --- Add Peak Module UI Calls Here ---
#           wellPanel(
#               style = "padding: 20px; position: relative; overflow: visible;", # Keep wellPanel wrapper?
#               h4("Strain Effects", style = "color: #2c3e50; margin-bottom: 15px;"),
#               div(style = "color: #7f8c8d; margin-bottom: 15px;",
#                   "Select a peak to see strain effects (Only for additive scans)."
#               ),
#               peakInput("peak"), # Renders which_peak dropdown
#               # Add download buttons here? Or move to separate download module?
#               div(style = "margin-bottom: 20px;",
#                    downloadAlleleUI("download") # <-- Use download module UI
#               ),
#               peakUI("peak")    # Renders allele_effects plot
#           )
#       ),
#       # Right column with tabbed interface
#       column(9,
#           tabsetPanel(
#               # LOD Plot Tab
#               tabPanel("LOD Plot",
#                   div(class = "plot-container",
#                       div(style = "display: flex; flex-direction: column; gap: 15px; margin-bottom: 20px;",
#                           div(style = "display: flex; justify-content: space-between; align-items: center;",
#                               h3("LOD Score Plot", 
#                                  style = "margin: 0; color: #2c3e50; font-weight: 600; font-family: 'Montserrat', 'Helvetica Neue', sans-serif; font-size: 28px; letter-spacing: 0.5px;"),
#                               div(style = "display: flex; gap: 10px;",
#                                   # downloadButton("download_qtl_plot_png", "Download PNG", class = "btn-sm"),
#                                   # downloadButton("download_qtl_plot_pdf", "Download PDF", class = "btn-sm")
#                                   downloadQtlUI("download") # <-- Use download module UI
#                               )
#                           ),
#                           div(style = "display: flex; gap: 20px; align-items: center;",
#                               div(style = "display: flex; align-items: center; gap: 15px;",
#                                   div(style = "display: flex; align-items: center; gap: 10px;",
#                                       numericInput("plot_width", "Width:",
#                                           value = 1000, min = 400, max = 2000, step = 50,
#                                           width = "100px"),
#                                       numericInput("plot_height", "Height:",
#                                           value = 600, min = 300, max = 1200, step = 50,
#                                           width = "100px")
#                                   ),
#                                   div(style = "display: flex; gap: 5px;",
#                                       actionButton("preset_1to1", "1:1", class = "btn-sm"),
#                                       actionButton("preset_3to2", "3:2", class = "btn-sm"),
#                                       actionButton("preset_16to9", "16:9", class = "btn-sm")
#                                   ),
#                                   div(style = "display: flex; align-items: center; gap: 10px;",
#                                       div(class = "lever-switch", id = "color_toggle", 
#                                           style = "margin-left: 10px;",
#                                           onclick = "Shiny.setInputValue('toggle_colors', Date.now())")
#                                   )
#                               )
#                           )
#                       )
#                   ),
#                   scanlyOutput("scanly"),
#                   div(style = "margin-top: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 4px; border: 1px solid #e9ecef;",
#                       scanlyUI("scanly")
#                   )
#               ),
#               # Mediation Analysis Tab (Placeholder)
#               tabPanel("Mediation Analysis",
#                   div(style = "padding: 20px; text-align: center; color: #7f8c8d;",
#                       h3("Mediation Analysis", style = "color: #2c3e50;"),
#                       p("This feature will be implemented in a future update.",
#                         style = "font-size: 16px; margin-top: 20px;")
#                   )
#               ),
#               # Regression Analysis Tab (Placeholder)
#               tabPanel("Regression Analysis",
#                   div(style = "padding: 20px; text-align: center; color: #7f8c8d;",
#                       h3("Regression Analysis", style = "color: #2c3e50;"),
#                       p("This feature will be implemented in a future update.",
#                         style = "font-size: 16px; margin-top: 20px;")
#                   )
#               )
#           )
#       )
#   )
)


# Server=======================================================================
server <- function(input, output, session) {
    message("Minimal server function started.")
    # Call import data function
    # imported_data <- import_data()
    
    # Initialize the main parameters module
    # main_par <- mainParServer("main_par", reactive(imported_data))
    
    # Initialize the peak module
    # peak_mod_reactives <- peakServer("peak", main_par, reactive(imported_data))

    # scan_data reactiveVal
    # scan_data <- reactiveVal(NULL)
    
    # clicked_data reactiveVal
    # clicked_data <- reactiveVal(NULL)
    
    # current_trait reactiveVal
    # current_trait <- reactiveVal("")
    
    # trait_debounced reactive
    # trait_debounced <- reactive({
    #   main_par()$which_trait # Use the reactive value from the module
    # }) %>% debounce(800)
    
    # observeEvent for trait_debounced - Uses main_par()
    # observeEvent(trait_debounced(), {
    #   req(main_par()$selected_dataset) # Use main_par() reactive now
    #   trait <- trait_debounced() # Use main_par()$which_trait
    #   
    #   # Only search if the trait is not empty and different from the current one
    #   if (nchar(trait) > 0 && trait != current_trait()) {
    #     # Update current trait
    #     current_trait(trait)
    #     
    #     # Clear error message and hide the container
    #     output$error_message <- renderText("")
    #     shinyjs::hide("error_message_container")
    #     
    #     # Clear clicked point data when a new trait is searched
    #     clicked_data(NULL)
    #     
    #     # Try to load the scan data for the specific trait
    #     tryCatch({
    #       message("Searching for trait: ", trait, " in dataset: ", main_par()$selected_dataset)
    #       # Need to source the trait_scan function if not automatically done by package structure
    #       # source("R/trait_scan.R") 
    #       data <- trait_scan(file_directory, main_par()$selected_dataset, trait)
    #       scan_data(data)
    #     }, error = function(e) {
    #       # Show error message
    #       output$error_message <- renderText({
    #         paste("Error:", e$message)
    #       })
    #       shinyjs::show("error_message_container")
    #       scan_data(NULL)
    #     })
    #   }
    # })
    
    # observeEvent for main_par()$which_trait - Uses main_par()
    # observeEvent(main_par()$which_trait, {
    #   clicked_data(NULL)
    # })
 
    # plot_obj reactive - Uses main_par()
    # plot_obj <- reactive({
    #     req(scan_data())
    #     # Need to source QTL_plot_visualizer if not loaded
    #     # source("R/QTL_plot_visualizer.R")
    #     QTL_plot_visualizer(scan_data(), current_trait(), main_par()$LOD_thr, markers)
    # })
 
    # official_gene_symbol reactiveVal
    # official_gene_symbol <- reactiveVal("")
    
    # use_alternating_colors reactiveVal
    # use_alternating_colors <- reactiveVal(TRUE)
    
    # observeEvent for toggle_colors
    # observeEvent(input$toggle_colors, {
    #   use_alternating_colors(!use_alternating_colors())
    #   # Toggle the active class on the lever switch based on the new state
    #   runjs(sprintf("
    #     const toggle = document.getElementById('color_toggle');
    #     if (%s) {
    #       toggle.classList.add('active');
    #     } else {
    #       toggle.classList.remove('active');
    #     }
    #   ", tolower(use_alternating_colors())))
    # })
    
    # observeEvents for preset buttons
    # observeEvent(input$preset_1to1, {
    #   updateNumericInput(session, "plot_width", value = 800)
    #   updateNumericInput(session, "plot_height", value = 800)
    # })
    
    # observeEvent(input$preset_3to2, {
    #   updateNumericInput(session, "plot_width", value = 900)
    #   updateNumericInput(session, "plot_height", value = 600)
    # })
    
    # observeEvent(input$preset_16to9, {
    #   updateNumericInput(session, "plot_width", value = 1600)
    #   updateNumericInput(session, "plot_height", value = 900)
    # })
    
    # observeEvent for scan_data - Uses main_par() and peak_mod_reactives() (implicitly via plot_obj?)
    # observeEvent(scan_data(), {
    #   if (!is.null(scan_data()) && "Phenotype" %in% colnames(scan_data())) {
    #     # Set official gene symbol
    #     official_gene_symbol(unique(scan_data()$Phenotype)[1])
    #     
    #     # Automatically find and select the highest peak when a trait is loaded
    #     if (nrow(scan_data()) > 0) {
    #       message("Finding highest peak for trait: ", official_gene_symbol())
    #       
    #       # Get peaks data for the selected trait
    #       peaks_data <- peak_finder(file_directory, main_par()$selected_dataset, main_par()$which_trait)
    #       
    #       # Make sure we have data and the required columns
    #       if (nrow(peaks_data) > 0 && all(c("marker", "lod") %in% colnames(peaks_data))) {
    #         # Filter by LOD threshold and sort to find highest
    #         filtered_peaks <- peaks_data %>% 
    #           filter(lod >= main_par()$LOD_thr) %>%
    #           arrange(desc(lod))
    #         
    #         if (nrow(filtered_peaks) > 0) {
    #           # Get the highest peak
    #           highest_peak <- filtered_peaks$marker[1]
    #           message("Found highest peak: ", highest_peak, " with LOD: ", round(filtered_peaks$lod[1], 2))
    #           
    #           # Update dropdown to select this peak
    #           updateSelectizeInput(session, "which_peak", selected = highest_peak)
    #           
    #           # Find this peak in the plot data
    #           plot_data <- plot_obj()
    #           peak_point <- plot_data[[2]] %>% filter(markers == highest_peak)
    #           
    #           if (nrow(peak_point) > 0) {
    #             # Create a fake click event at this peak's position
    #             fake_click <- list(
    #               x = if(main_par()$selected_chr == "All") peak_point$BPcum[1] else peak_point$position[1],
    #               y = peak_point$LOD[1],
    #               curveNumber = 0,
    #               pointNumber = which(plot_data[[2]]$markers == highest_peak)[1]
    #             )
    #             
    #             # Update clicked data to automatically show point info for highest peak
    #             clicked_data(fake_click)
    #             message("Auto-selected highest peak with LOD: ", round(filtered_peaks$lod[1], 2))
    #           } else {
    #             message("Could not find highest peak marker in plot data")
    #           }
    #         } else {
    #           message("No peaks found above threshold (", main_par()$LOD_thr, ") for this trait")
    #         }
    #       } else {
    #         if (nrow(peaks_data) == 0) {
    #           message("No peaks found for this trait")
    #         } else {
    #           message("Required columns missing from peaks data. Have: ", 
    #                   paste(colnames(peaks_data), collapse=", "))
    #         }
    #       }
    #     }
    #   }
    # })
    
    # observeEvent for plotly_click
    # observeEvent(event_data("plotly_click", source = "scan_plot"), {
    #     clicked_data(event_data("plotly_click", source = "scan_plot"))
    # })
    
    # plot_base reactive - Uses main_par()
    # plot_base <- reactive({
    #     req(plot_obj())
    #     plot_data <- plot_obj()[[2]]
    #     
    #     # Filter data based on selected chromosome
    #     if (main_par()$selected_chr != "All") {
    #         chr_num <- switch(main_par()$selected_chr,
    #             "X" = 20,
    #             "Y" = 21,
    #             "M" = 22,
    #             as.numeric(main_par()$selected_chr)
    #         )
    #         plot_data <- plot_data %>% filter(chr == chr_num)
    #     }
    #     
    #     # Calculate chromosome axis positions
    #     if (main_par()$selected_chr == "All") {
    #         # For all chromosomes view
    #         axisdf <- plot_data %>%
    #             group_by(chr) %>%
    #             summarise(center = mean(BPcum))
    #         
    #         # Convert chromosome numbers to proper labels
    #         chr_labels <- as.character(axisdf$chr)
    #         chr_labels[chr_labels == "20"] <- "X"
    #         chr_labels[chr_labels == "21"] <- "Y"
    #         chr_labels[chr_labels == "22"] <- "M"
    #         
    #         # Create the base ggplot object for all chromosomes
    #         p <- ggplot(plot_data, aes(x = BPcum, y = LOD)) +
    #             geom_line(aes(color = if(use_alternating_colors()) as.factor(chr) else NULL), size = 0.75) +
    #             geom_hline(yintercept = main_par()$LOD_thr, color = "#e74c3c",
    #                       linetype = "dashed", size = 0.8) +
    #             scale_color_manual(values = if(use_alternating_colors()) rep(c("#3498db", "#2c3e50"), 22) else "#3498db") +
    #             scale_x_continuous(
    #                 breaks = axisdf$center,
    #                 labels = chr_labels,
    #                 expand = expansion(mult = c(0.01, 0.01))
    #             )
    #     } else {
    #         # For single chromosome view
    #         # Create the base ggplot object for a single chromosome
    #         p <- ggplot(plot_data, aes(x = position, y = LOD)) +
    #             geom_line(color = "#3498db", size = 0.75) +
    #             geom_hline(yintercept = main_par()$LOD_thr, color = "#e74c3c",
    #                       linetype = "dashed", size = 0.8) +
    #             scale_x_continuous(
    #                 expand = expansion(mult = c(0.02, 0.02))
    #             )
    #     }
    #     
    #     # Common theme and y-axis settings for both views
    #     p <- p +
    #         scale_y_continuous(
    #             expand = expansion(mult = c(0.02, 0.1))
    #         ) +
    #         theme_minimal() +
    #         theme(
    #             panel.grid.minor = element_blank(),
    #             panel.grid.major.x = element_blank(),
    #             panel.grid.major.y = element_line(color = "#ecf0f1", size = 0.2),
    #             axis.line = element_line(color = "#2c3e50", size = 0.5),
    #             axis.text = element_text(size = 11, color = "#2c3e50"),
    #             axis.title = element_text(size = 12, face = "bold", color = "#2c3e50"),
    #             axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    #             legend.position = "none",
    #             plot.title = element_text(size = 14, face = "bold", color = "#2c3e50"),
    #             plot.subtitle = element_text(size = 11, color = "#7f8c8d"),
    #             plot.background = element_rect(fill = "white", color = NA),
    #             panel.background = element_rect(fill = "white", color = NA),
    #             plot.margin = margin(t = 20, r = 20, b = 40, l = 40)
    #         ) +
    #         labs(
    #             x = if(main_par()$selected_chr == "All") "Chromosome" else paste("Position on Chromosome", main_par()$selected_chr, "(Mb)"),
    #             y = "LOD Score"
    #         )
    #     
    #     return(list(p = p, data = plot_data))
    # })
 
    # renderPlotly for scan_plot - REMOVE THIS
    # output$scan_plot <- renderPlotly({ ... })
 
    # observeEvent for plotly_doubleclick - REMOVE THIS (handled in scanly)
    # observeEvent(event_data("plotly_doubleclick", source = "scan_plot"), { ... })
 
    # renderDT for clicked_point_info - REMOVE THIS
    # output$clicked_point_info <- renderDT({ ... })
 
    # --- Scanly Module --- 
    # scan_object reactive (for scanly)
    # scan_object_reactive <- reactive({
    #     req(plot_base())
    #     # plot_base returns list(p = plot, data = data)
    #     list(plot = plot_base()$p, table = plot_base()$data)
    # })
    
    # scanlyServer("scanly", 
    #              main_par = main_par, 
    #              scan_object = scan_object_reactive, 
    #              peak_mod_reactives = peak_mod_reactives, 
    #              official_trait_symbol = official_gene_symbol)

    # --- Download Module --- 
    # Reactive for plot dimensions
    # plot_dims_reactive <- reactive({
    #     list(width = input$plot_width, height = input$plot_height)
    # })

    # Call download server
    # downloadServer("download",
    #                main_par = main_par,
    #                plot_base_obj = plot_base,
    #                peak_mod_reactives = peak_mod_reactives,
    #                official_gene_symbol = official_gene_symbol,
    #                plot_dims = plot_dims_reactive)

    # --- REMOVE OLD Download Handlers --- 
    # output$download_effects_plot_png <- downloadHandler({ ... })
    # output$download_effects_plot_pdf <- downloadHandler({ ... })
    # output$download_qtl_plot_png <- downloadHandler({ ... })
    # output$download_qtl_plot_pdf <- downloadHandler({ ... })
}

# Launch the Shiny application
shinyApp(ui = ui, server = server)

message("Reached end of minimal app.R")


