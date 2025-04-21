# UI definition for the Shiny application

ui <- fluidPage(
  shinyjs::useShinyjs(), 
  
  tags$head(
      tags$style(HTML("
  .well {
    background-color: #ffffff;
    border: 1px solid #e3e3e3;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    margin-bottom: 20px;
  }
  .lever-switch {
    position: relative;
    display: inline-block;
    width: 60px;
    height: 30px;
    background-color: #ccc;
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
    background-color: #3498db; /* Use primary color when active */
  }
  .lever-switch.active:before {
    transform: translateX(30px);
  }
  .lever-switch:hover {
    box-shadow: 0 0 5px rgba(0,0,0,0.2);
  }
  .control-label {
    color: #2c3e50; /* Use secondary color */
    font-weight: 500;
    font-size: 14px;
  }
  .form-control,
  .selectize-input {
    border-radius: 6px;
    border: 1px solid #dce4ec;
  }
  .btn {
    border-radius: 6px;
    text-transform: uppercase;
    font-size: 12px;
    font-weight: 600;
    padding: 8px 16px;
    transition: all 0.3s ease;
  }
  .btn:hover {
    transform: translateY(-1px);
    box-shadow: 0 4px 8px rgba(0,0,0,0.1);
  }
  .btn-primary {
    background-color: #3498db; /* Use primary color */
    color: white;
    border: none;
  }
  .btn-primary:hover {
    background-color: #2980b9; /* Darker primary */
    color: white;
  }
  .selectize-dropdown {
    z-index: 10000 !important; /* Ensure dropdowns are on top */
  }
  .title-panel {
    background-color: #2c3e50; /* Use secondary color */
    color: white;
    padding: 20px;
    margin-bottom: 30px;
    border-radius: 8px;
  }
  .plot-container, .datatable-container {
    background-color: white;
    padding: 20px;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    margin-bottom: 20px;
  }
  #error_message_container {
    color: #e74c3c; /* Use danger color */
    padding: 10px;
    margin-top: 10px;
    font-size: 14px;
    border: 1px dashed #e74c3c;
    border-radius: 4px;
    background-color: #fdedec;
    display: none; /* Hide by default */
  }
  /* Spinner styling */
  .spinner-border {
    color: #3498db; /* Use primary color */
    width: 3rem; 
    height: 3rem;
  }
  /* Tooltip styling */
  .tooltip-inner {
    background-color: #2c3e50;
    color: white;
    border-radius: 4px;
    padding: 8px 12px;
  }
  "))
  ),
  
  theme = bs_theme(
      version = 4,
      bootswatch = "flatly",
      primary = "#3498db",
      secondary = "#2c3e50",
      success = "#18bc9c",
      info = "#3498db", 
      warning = "#f39c12",
      danger = "#e74c3c",
      bg = "#ecf0f1", # Light background
      fg = "#2c3e50", # Default text color
      base_font = font_google("Lato"),
      heading_font = font_google("Montserrat")
  ),
  
  # Title Panel
  div(class = "title-panel",
      h1("Pre-scanned QTL Visualizer for Diet DO Study", 
         style = "font-size: 28px; margin: 0; font-family: 'Montserrat', sans-serif;"),
      p("Interactive visualization tool for QTL analysis", 
        style = "margin: 10px 0 0 0; opacity: 0.8; font-family: 'Lato', sans-serif;")
  ),
  
  fluidRow(
      
      column(3,
          # Input Panel
          wellPanel(
              style = "padding: 20px;",
              div(style = "margin-bottom: 25px;",
                  selectizeInput("selected_dataset", 
                      h4("Dataset Selection", style = "color: #2c3e50; margin-bottom: 15px;"),
                      choices = NULL, 
                      multiple = FALSE,
                      options = list(
                          placeholder = 'Select a dataset...', 
                          onInitialize = I('function() { this.setValue(""); }')
                      ))
              ),
              
              div(style = "margin-bottom: 25px;",
                  sliderInput("LOD_thr", 
                      h4("LOD Threshold", style = "color: #2c3e50; margin-bottom: 15px;"),
                      min = 4, max = 120, value = 7, step = 0.5,
                      ticks = TRUE)
              ),
              
              div(style = "margin-bottom: 25px;",
                  textInput("which_trait", 
                      h4("Search Trait", style = "color: #2c3e50; margin-bottom: 15px;"),
                      value = "", 
                      placeholder = "Input gene symbol (e.g., Gnai3)")
              ),
              
              
              div(id = "error_message_container", 
                  textOutput("error_message")
              ),
              
              
              div(style = "display: none;", 
                  actionButton("search_trait", "Search", class = "btn-primary", style = "width: 100%;")
              )
          ),
          
          
          wellPanel(
              style = "padding: 20px;",
              h4("Strain Effects", style = "color: #2c3e50; margin-bottom: 15px;"),
              div(style = "color: #7f8c8d; margin-bottom: 15px; font-size: 13px;",
                  "Select a peak to see strain effects (Available for additive scans)."
              ),
              div(style = "margin-bottom: 20px;",
                  selectizeInput("which_peak", "Choose peak", 
                      choices = NULL, 
                      multiple = FALSE,
                      options = list(
                          placeholder = 'Select a peak...', 
                          onInitialize = I('function() { this.setValue(""); }')
                      ))
              ),
              div(style = "margin-bottom: 20px; display: flex; gap: 10px;",
                  downloadButton("download_effects_plot_png", "PNG", class = "btn-sm"),
                  downloadButton("download_effects_plot_pdf", "PDF", class = "btn-sm")
              ),
              div(style = "margin-top: 10px;", 
                  # Use a spinner while the plot is rendering
                  shinycssloaders::withSpinner(
                      plotOutput("allele_effects", height = "400px"), 
                      type = 8 
                  )
              )
          )
      ),
      
      
      column(9,
          # LOD Plot Section
          div(class = "plot-container",
              # Plot Controls and Title
              div(style = "display: flex; flex-direction: column; gap: 15px; margin-bottom: 20px;",
                  # Title and Download 
                  div(style = "display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 10px;",
                      h3("LOD Score Plot", style = "margin: 0; color: #2c3e50; font-weight: 600; font-family: 'Montserrat', sans-serif; font-size: 24px; letter-spacing: 0.5px;"),
                      div(style = "display: flex; gap: 10px;",
                          downloadButton("download_qtl_plot_png", "PNG", class = "btn-sm"),
                          downloadButton("download_qtl_plot_pdf", "PDF", class = "btn-sm")
                      )
                  ),
                  
                  
                  div(style = "display: flex; gap: 15px; align-items: center; flex-wrap: wrap;",
                      # Chromosome Zoom Selector
                      selectInput("selected_chr", "Zoom Chr:", 
                          choices = c("All", as.character(1:19), "X", "Y", "M"), 
                          selected = "All", 
                          width = "120px"
                      ),
                      
                      # Plot Dimensions Inputs
                      numericInput("plot_width", "Width:", value = 1000, min = 400, max = 2000, step = 50, width = "90px"),
                      numericInput("plot_height", "Height:", value = 600, min = 300, max = 1200, step = 50, width = "90px"),
                      
                      # Dimension Presets
                      div(style = "display: flex; gap: 5px;",
                          actionButton("preset_1to1", "1:1", class = "btn-xs btn-outline-secondary"), # Use smaller outline buttons
                          actionButton("preset_3to2", "3:2", class = "btn-xs btn-outline-secondary"),
                          actionButton("preset_16to9", "16:9", class = "btn-xs btn-outline-secondary")
                      ),
                      
                      # Color Toggle Switch
                      div(style = "display: flex; align-items: center; gap: 5px;",
                          tags$label("Alt Colors:", `for` = "color_toggle", class = "control-label", style="margin-bottom: 0;"),
                          div(class = "lever-switch active", id = "color_toggle", # Start as active
                              onclick = "$(this).toggleClass('active'); Shiny.setInputValue('toggle_colors', $(this).hasClass('active'))") 
                      )
                  )
              ),
              
              # The Plotly Output with Spinner
              shinycssloaders::withSpinner(
                  plotlyOutput("scan_plot", width = "100%", height = "auto"), 
                  type = 8 
              )
          ),
          
          # Clicked Point Info Table Section
          div(class = "datatable-container",
              style = "margin-top: 10px;",
              h5("Selected Point / Peak Information", style = "margin: 0 0 15px 0; color: #2c3e50; font-size: 16px; font-weight: 600;"),
              DTOutput("clicked_point_info")
          )
      )
  )
)

