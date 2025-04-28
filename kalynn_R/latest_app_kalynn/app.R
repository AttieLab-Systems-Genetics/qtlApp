# Load required R packages for the application
require("rstudioapi")  
require("dplyr")      
require("stringr")     
require("tidyverse")   
require("BiocManager") 
require("ggplot2")     
require("qtl2")        
require("grid")        
require("ggrepel")     
require("gridGraphics")
require("ggpubr")     
require("shiny")      
require("shinyFiles") 
require("bslib")      
require("spsComps")   
require("DT")         
require("shinyjs")    
require("shinycssloaders")
require("data.table") 
require("reshape2")   
require("plotly")     
require("ggiraph")    
require("writexl")   
require("fontawesome")
# require("debounce") # Commented out as it's not available for R 4.4.2

# Source modules (assuming they are in the R/ directory relative to app)
# Function to source all R files in a directory
source_dir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
source_dir("R")

# Set maximum file upload size for Shiny (20GB)
options(shiny.maxRequestSize = 20000*1024^2)

# Load required data files
# Read the file index containing paths to all data files
file_directory <- read.csv("/data/dev/miniViewer_3.0/file_index.csv")

# Try to load gene symbols file
gene_symbols_path <- "/data/dev/miniViewer_3.0/gene_symbols.csv"

# Default gene symbols in case loading fails
gene_symbols <- c("Gnai3", "Cdc45", "Slc4a1", "Abca12", "Nadk", "Tfpi", "Scnn1b", "Cdc20", "Gpr89", "Cdc73")

if (file.exists(gene_symbols_path)) {
  message("Loading gene symbols from: ", gene_symbols_path)
  gene_symbols <- tryCatch({
    as.character(fread(gene_symbols_path)$gene_symbol)
  }, error = function(e) {
    warning("Error reading gene symbols file: ", e$message)
    c("Gnai3", "Cdc45", "Slc4a1", "Abca12", "Nadk", "Tfpi", "Scnn1b", "Cdc20", "Gpr89", "Cdc73")
  })
} else {
  warning("Gene symbols file not found at: ", gene_symbols_path, ". Using default symbols.")
}

# Sort gene symbols
gene_symbols <- sort(gene_symbols)
message("Loaded ", length(gene_symbols), " gene symbols")

# Load chromosome break points for mm11 genome
chr_breaks <- read.csv("/data/dev/miniViewer_3.0/chromosomal_sep_mm11.csv")

# Load gene annotations
annotation_list <- readRDS("/data/dev/miniViewer_3.0/annotation_list.rds")

# Load marker information for the diet DO study
markers <- readRDS(file.path("/data/dev/miniViewer_3.0/CHTC_dietDO_markers_RDSgrcm39.rds"))

# Create a group identifier by combining diet, trait compartment, trait type, and scan type
# Include 'sexes' to differentiate datasets like male-only vs both sexes, but only if not 'Both'
file_directory$group <- paste0(file_directory$diet, " ", file_directory$trait_compartment, " ",
                           file_directory$trait_type, 
                           # Conditionally add sex info if it's not 'Both'
                           ifelse(file_directory$sexes == "Both", "", paste0(" (", file_directory$sexes, ")")),
                           ", ", file_directory$scan_type,
                           ifelse(file_directory$scan_type == "interactive",
                                  paste0(" (", file_directory$covars_interactive, ")"),
                                  ""))

# Create a new environment for caching file paths to improve performance
file_path_cache <- new.env(parent = emptyenv())

trait_cache <- new.env(parent = emptyenv())

# Trait scan function is now expected to be sourced from R/trait_scan.R
# trait_scan <- function(...) { ... }

# Also add caching for peak_finder
peaks_cache <- new.env(parent = emptyenv())

# Peak finder function is now expected to be sourced from R/peak_finder.R
# peak_finder <- function(...) { ... }

# set microfunctions==========================================================
# QTL plot visualizer is now expected to be sourced from R/QTL_plot_visualizer.R
# QTL_plot_visualizer <- function(...) { ... }

# Function to convert CSV files to FST format is now expected to be sourced from R/csv2fst.R
# csv2fst <- function(...) { ... }

# Function to create row indices is now expected to be sourced from R/fst_rows.R
# fst_rows <- function(...) { ... }

# set UI=======================================================================
ui <- fluidPage(
  useShinyjs(),
  # Add custom CSS
  tags$head(
      tags$style(HTML("
  .well {
    background-color: #ffffff;
    border: 1px solid #e3e3e3;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    margin-bottom: 20px;
  }
  /* Add lever switch styling */
  .lever-switch {
    position: relative;
    display: inline-block;
    width: 60px;
    height: 30px;
    background-color: #2c3e50;
    border-radius: 15px;
    cursor: pointer;
    transition: all 0.3s ease;
    overflow: hidden;
  }
  .lever-switch:before {
    content: '';
    position: absolute;
    width: 26px;
    height: 26px;
    left: 2px;
    bottom: 2px;
    background-color: white;
    border-radius: 50%;
    transition: all 0.3s ease;
  }
  .lever-switch.active {
    background-color: #3498db;
  }
  .lever-switch.active:before {
    transform: translateX(30px);
  }
  .lever-switch:hover {
    box-shadow: 0 0 5px rgba(0,0,0,0.2);
  }
  .control-label {
    color: #2c3e50;
    font-weight: 500;
    font-size: 14px;
  }
  .form-control {
    border-radius: 6px;
    border: 1px solid #dce4ec;
  }
  .btn {
    border-radius: 6px;
    text-transform: uppercase;
    font-size: 12px;
    font-weight: 600;
    padding: 8px 16px;
  }
  .btn-default {
    background-color: #3498db;
    color: white;
    border: none;
  }
  .btn-default:hover {
    background-color: #2980b9;
    color: white;
  }
  .selectize-input {
    border-radius: 6px;
    border: 1px solid #dce4ec;
    position: relative;
    z-index: 1;
  }
  .selectize-dropdown {
    z-index: 10000 !important;
  }
  /* Ensure the strain effects selectize appears above the spinner */
  #peak-which_peak + .selectize-control .selectize-dropdown { /* Adjusted ID based on module */
    z-index: 99999 !important;
  }
  /* Ensure all selectize dropdowns remain on top */
  .selectize-control {
    position: relative;
    z-index: 10;
  }
  .selectize-dropdown {
    z-index: 99999 !important;
  }
  /* Fix positioning of the spinner relative to the plot */
  .shiny-spinner-output-container {
    position: relative;
    z-index: 0;
  }
  .title-panel {
    background-color: #2c3e50;
    color: white;
    padding: 20px;
    margin-bottom: 30px;
    border-radius: 8px;
  }
  .plot-container {
    background-color: white;
    padding: 20px;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    margin-bottom: 20px;
  }
  .datatable-container {
    background-color: white;
    padding: 20px;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
  }
  #error_message {
    color: #e74c3c;
    padding: 10px;
    margin-top: 10px;
    font-size: 14px;
  }
   /* Add animation for hover effects */
  .well:hover {
    transform: translateY(-2px);
    transition: transform 0.2s ease;
  }
   .btn {
    transition: all 0.3s ease;
  }
   .btn:hover {
    transform: translateY(-1px);
    box-shadow: 0 4px 8px rgba(0,0,0,0.1);
  }
   /* Add pulsing effect to the spinner */
  .spinner-border {
    animation: spinner-pulse 1s ease infinite;
  }
   @keyframes spinner-pulse {
    0% { transform: scale(1); }
    50% { transform: scale(1.1); }
    100% { transform: scale(1); }
  }
   /* Add tooltip styling */
  .tooltip-inner {
    background-color: #2c3e50;
    color: white;
    border-radius: 4px;
    padding: 8px 12px;
  }
   /* Add success message styling */
  .success-message {
    color: #18bc9c;
    padding: 10px;
    margin-top: 10px;
    font-size: 14px;
    display: none;
  }
"))
  ),
  # Use bslib theme with custom settings
  theme = bs_theme(
      version = 4,
      bootswatch = "flatly",
      primary = "#3498db",
      secondary = "#2c3e50",
      success = "#18bc9c",
      info = "#3498db",
      warning = "#f39c12",
      danger = "#e74c3c"
  ),
  # Title Panel with modern styling
  div(class = "title-panel",
      h1("Pre-scanned QTL Visualizer for Diet DO Study",
          style = "font-size: 28px; margin: 0;"),
      p("Interactive visualization tool for QTL analysis",
          style = "margin: 10px 0 0 0; opacity: 0.8;")
  ),
  fluidRow(
      # Left column: Main parameters and Allele Effects
      column(3,
          wellPanel(
              style = "padding: 20px;",
              # --- Main Parameters Module --- 
              mainParInput("main_par"), # UI for dataset/trait selection
              mainParUI("main_par"),    # Placeholder for potential outputs from mainPar
              div(id = "error_message_container", style = "margin-bottom: 15px; color: #e74c3c; font-weight: bold; display: none;",
                  textOutput("error_message") # Keep global error message display here
              )
          ),
          wellPanel(
              style = "padding: 20px; position: relative; overflow: visible;",
              h4("Strain Effects", style = "color: #2c3e50; margin-bottom: 15px;"),
              div(style = "color: #7f8c8d; margin-bottom: 15px;",
                  "Select a peak to see strain effects (Only for additive scans)."
              ),
              # --- Peak Module --- 
              peakInput("peak"),       # UI for peak selection dropdown
              downloadAlleleUI("download"), # UI for allele plot download buttons
              peakUI("peak")           # UI for displaying the allele effects plot
          )
      ),
      # Right column: LOD Plot and other tabs
      column(9,
          tabsetPanel(
              # LOD Plot Tab
              tabPanel("LOD Plot",
                  # --- Scanly Module (Handles Plot and Interactions) --- 
                  scanlyOutput("scanly"), # Output for the plotly object
                  div(style = "margin-top: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 4px; border: 1px solid #e9ecef;",
                     scanlyUI("scanly")    # UI elements from scanly (like clicked point info)
                  ),
                  # Add download buttons for QTL plot (using download module)
                  div(style="margin-top: 15px;", downloadQtlUI("download"))
              ),
              # Mediation Analysis Tab (Placeholder)
              tabPanel("Mediation Analysis",
                  div(style = "padding: 20px; text-align: center; color: #7f8c8d;",
                      h3("Mediation Analysis", style = "color: #2c3e50;"),
                      p("This feature will be implemented in a future update.",
                        style = "font-size: 16px; margin-top: 20px;")
                  )
              ),
              # Regression Analysis Tab (Placeholder)
              tabPanel("Regression Analysis",
                  div(style = "padding: 20px; text-align: center; color: #7f8c8d;",
                      h3("Regression Analysis", style = "color: #2c3e50;"),
                      p("This feature will be implemented in a future update.",
                        style = "font-size: 16px; margin-top: 20px;")
                  )
              )
          )
      )
  )
)


# Server=======================================================================
server <- function(input, output, session) {
    message("Server function started.")
    
    # --- Initial Data & Reactive Values --- 
    # Make sure file_directory, markers, annotation_list, gene_symbols are available
    # These are loaded at the top level of app.R now
    
    # Create a combined reactive list for data needed by modules
    imported_data_reactive <- reactive({
        list(
            file_directory = file_directory,
            annotation_list = annotation_list,
            gene_symbols = gene_symbols,
            markers = markers,
            chr_breaks = chr_breaks
        )
    })
    
    # Reactive value for the loaded scan data for the current trait
    scan_data <- reactiveVal(NULL)
    # Reactive value for the currently selected trait's official symbol/name
    official_trait_symbol <- reactiveVal("")
    # Reactive value for the full peaks table for the current trait
    trait_peaks_data <- reactiveVal(NULL)
    # Reactive for error message
    error_message_reactive <- reactiveVal("") 
    
    # --- Module Server Calls --- 
    
    # Main Parameters Module
    # Pass the combined reactive data list
    main_par_reactives <- mainParServer("main_par", 
                                        imported_data_reactive) 
    
    # Peak Module (Allele Effects)
    # Pass main parameters and the imported data reactive
    peak_mod_reactives <- peakServer("peak", 
                                     main_par_reactives, 
                                     imported_data_reactive) # Corrected arguments
    
    # Scanly Module (LOD Plot)
    # Needs main parameters, the scan data, peak module reactives, and official trait symbol
    scanlyServer("scanly", 
                 main_par_reactives, 
                 scan_data, 
                 peak_mod_reactives,
                 official_trait_symbol)
                 
    # Download Module
    # Needs main parameters, scan_data, peak module reactives, and official trait symbol
    downloadServer("download", 
                   main_par_reactives, 
                   scan_data, 
                   peak_mod_reactives, 
                   official_trait_symbol)
                   
    # --- Core Logic connecting Dataset/Trait selection to Data Loading --- 
    
    # Debounce the trait input from the main_par module
    trait_debounced <- reactive({ 
      main_par_reactives()$which_trait # Access the reactive value correctly
    }) %>% debounce(800)
    
    # Observe the debounced trait input to load scan data
    observeEvent(trait_debounced(), { 
      # Use req() on the reactive value itself
      req(main_par_reactives()$selected_dataset, trait_debounced()) 
      selected_trait_input <- trait_debounced() # This is the raw Symbol_ID input
      selected_dataset_name <- main_par_reactives()$selected_dataset
      
      # Only proceed if trait is not empty
      if (nchar(selected_trait_input) > 0) {
        # --- Extract the actual trait symbol using the helper function ---
        selected_trait_symbol <- helpers::get_selected_trait(
            file_directory = file_directory, # Pass the file directory
            which_trait = selected_trait_input, 
            selected_dataset = selected_dataset_name
        )
        
        # Check if symbol extraction worked
        if(is.null(selected_trait_symbol)){
             error_msg <- paste("Could not determine trait symbol from input:", selected_trait_input)
             message("ERROR: ", error_msg) 
             error_message_reactive(error_msg)
             shinyjs::show("error_message_container")
             scan_data(NULL) 
             trait_peaks_data(NULL)
             official_trait_symbol("")
             return() # Stop processing if symbol extraction failed
        }
        # --- End symbol extraction ---

        message("DEBUG: Trait observer fired for input: ", selected_trait_input, ", using symbol: ", selected_trait_symbol, " in dataset: ", selected_dataset_name)
        
        # Clear previous results and error
        scan_data(NULL)
        trait_peaks_data(NULL)
        official_trait_symbol("")
        error_message_reactive("")
        shinyjs::hide("error_message_container")

        # Try to load scan data using the extracted symbol
        tryCatch({
          message("DEBUG: Calling trait_scan for symbol: ", selected_trait_symbol)
          loaded_scan_data <- trait_scan(file_directory, selected_dataset_name, selected_trait_symbol) # Use symbol here
          message("DEBUG: trait_scan returned ", nrow(loaded_scan_data), " rows.") # DEBUG
          scan_data(loaded_scan_data) # Update reactive value
          
          # Extract and store the official symbol from the loaded data
          if("Phenotype" %in% colnames(loaded_scan_data) && nrow(loaded_scan_data) > 0) {
             official_symbol <- unique(loaded_scan_data$Phenotype)[1] # Use the actual Phenotype from data
             official_trait_symbol(official_symbol)
             message("DEBUG: Set official trait symbol: ", official_symbol)
          } else {
            # Fallback if Phenotype column missing - use the extracted symbol
            official_trait_symbol(selected_trait_symbol) 
            message("DEBUG: Could not set official trait symbol from scan_data. Using extracted symbol: ", selected_trait_symbol) 
          }
          
          # Try to load peaks data for this trait using the extracted symbol
          message("DEBUG: Calling peak_finder for symbol: ", selected_trait_symbol)
          loaded_peaks_data <- peak_finder(file_directory, selected_dataset_name, selected_trait_symbol) # Use symbol here
          message("DEBUG: peak_finder returned ", nrow(loaded_peaks_data), " rows.") # DEBUG
          trait_peaks_data(loaded_peaks_data) # Update reactive value
          
        }, error = function(e) {
          error_msg <- paste("Error loading data for trait ", selected_trait_symbol, ":", e$message) # Report error with the symbol used
          message("ERROR: ", error_msg) # DEBUG
          error_message_reactive(error_msg)
          shinyjs::show("error_message_container")
          scan_data(NULL) # Clear data on error
          trait_peaks_data(NULL)
          official_trait_symbol("")
        })
      } else {
          message("DEBUG: Trait observer fired, but trait was empty.") # DEBUG
      }
    }, ignoreInit = TRUE)
    
    # Observe dataset changes to clear caches and potentially reload current trait
    observeEvent(main_par_reactives()$selected_dataset, { 
        # Ensure it's not NULL or empty initially and is different from previous
        req(main_par_reactives()$selected_dataset) 
        selected_dataset_name <- main_par_reactives()$selected_dataset
        # Use isolate to get current trait without creating dependency loop
        current_trait_value <- isolate(main_par_reactives()$which_trait) 
        
        message("Dataset changed to: ", selected_dataset_name, ". Clearing caches.")
        
        # Clear caches
        rm(list = ls(envir = trait_cache), envir = trait_cache)
        rm(list = ls(envir = peaks_cache), envir = peaks_cache)
        
        # Clear reactive values
        scan_data(NULL)
        trait_peaks_data(NULL)
        official_trait_symbol("")
        error_message_reactive("")
        shinyjs::hide("error_message_container")
        
        # If a trait was selected before dataset change, trigger the debounced observer 
        if (isolate(nchar(current_trait_value) > 0)) { # Use isolate here too
           message("Re-triggering search for trait: ", current_trait_value, " in new dataset.")
           # This should happen automatically via the trait_debounced observer firing
        }
    }, ignoreNULL = TRUE, ignoreInit = TRUE) # Important: ignore initial NULL state and first run

    # Render the error message text
    output$error_message <- renderText({ 
      error_message_reactive() 
    })
    
    # --- REMOVE OLD SERVER LOGIC HANDLED BY MODULES --- 
    # (e.g., plot rendering, peak dropdown updates, click handling, download handlers are now in modules)

    message("Server function setup complete.")
}

# Launch the Shiny application
shinyApp(ui = ui, server = server)


