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
  }
  .selectize-dropdown {
    z-index: 10000 !important;
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
      # Left column with inputs and allele effects
      column(3,
          wellPanel(
              style = "padding: 20px;",
              div(style = "margin-bottom: 25px;",
                  selectizeInput("selected_dataset",
                      h4("Dataset Selection", style = "color: #2c3e50; margin-bottom: 15px;"),
                      choices = unique(file_directory$group),
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
             
              # Add error message display
              div(id = "error_message_container", style = "margin-bottom: 15px; color: #e74c3c; font-weight: bold; display: none;",
                  textOutput("error_message")
              ),
             
              # Hide the search button since we'll auto-search
              div(style = "margin-bottom: 25px; display: none;",
                  actionButton("search_trait", "Search",
                      class = "btn-primary",
                      style = "width: 100%;")
              )
          ),
    
          # Allele effects panel
          wellPanel(
              style = "padding: 20px;",
              h4("Strain Effects", style = "color: #2c3e50; margin-bottom: 15px;"),
              div(style = "color: #7f8c8d; margin-bottom: 15px;",
                  "Select a peak to see strain effects (Only for additive scans)."
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
              div(style = "margin-bottom: 20px;",
                  downloadButton("download_effects_plot_png", "Download PNG"),
                  downloadButton("download_effects_plot_pdf", "Download PDF")
              ),
              div(style = "margin-top: 20px;",
                  plotOutput("allele_effects", height = "400px") %>%
                      withSpinner(type = 8, color = "#3498db")
              )
          )
      ),
      # Right column with LOD plot and peaks table
      column(9,
          # LOD plot section with clicked point info
          div(class = "plot-container",
              div(style = "display: flex; flex-direction: column; gap: 15px; margin-bottom: 20px;",
                  # Title and download buttons row
                  div(style = "display: flex; justify-content: space-between; align-items: center;",
                      h3("LOD Score Plot", 
                         style = "margin: 0; color: #2c3e50; font-weight: 600; font-family: 'Montserrat', 'Helvetica Neue', sans-serif; font-size: 28px; letter-spacing: 0.5px;"),
                      div(style = "display: flex; gap: 10px;",
                          downloadButton("download_qtl_plot_png", "Download PNG", class = "btn-sm"),
                          downloadButton("download_qtl_plot_pdf", "Download PDF", class = "btn-sm")
                      )
                  ),
                  
                  # Controls row
                  div(style = "display: flex; gap: 20px; align-items: center;",
                      # Chromosome selector
                      div(style = "display: flex; align-items: center; gap: 10px;",
                          selectInput("selected_chr", "Zoom to Chromosome:",
                              choices = c("All", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                        "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y", "M"),
                              selected = "All",
                              width = "150px"
                          )
                      ),
                      
                      # Plot dimensions
                      div(style = "display: flex; align-items: center; gap: 15px;",
                          div(style = "display: flex; align-items: center; gap: 10px;",
                              numericInput("plot_width", "Width:",
                                  value = 1000, min = 400, max = 2000, step = 50,
                                  width = "100px"),
                              numericInput("plot_height", "Height:",
                                  value = 600, min = 300, max = 1200, step = 50,
                                  width = "100px")
                          ),
                          div(style = "display: flex; gap: 5px;",
                              actionButton("preset_1to1", "1:1", class = "btn-sm"),
                              actionButton("preset_3to2", "3:2", class = "btn-sm"),
                              actionButton("preset_16to9", "16:9", class = "btn-sm")
                          ),
                          # Add color toggle button
                          div(style = "display: flex; align-items: center; gap: 10px;",
                              div(class = "lever-switch", id = "color_toggle", 
                                  style = "margin-left: 10px;",
                                  onclick = "Shiny.setInputValue('toggle_colors', Date.now())")
                          )
                      )
                  )
              )
          ),
          plotlyOutput("scan_plot", width = "100%", height = "auto") %>%
              withSpinner(type = 8, color = "#3498db"),
          # Add clicked point info directly below plot
          div(style = "margin-top: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 4px; border: 1px solid #e9ecef;",
              h5("Selected Point Information", style = "margin: 0 0 10px 0; color: #2c3e50; font-size: 14px;"),
              DTOutput("clicked_point_info")
          )
      )
  )
)
