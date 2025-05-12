# UI styling module for the QTL viewer app
# This module contains all the custom CSS and styling elements

# Custom CSS for modern UI
custom_css <- HTML("
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
  #which_peak + .selectize-control .selectize-dropdown {
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
")

# Function to create a modern title panel
create_title_panel <- function(title, subtitle = NULL) {
  div(class = "title-panel",
      h1(title, style = "font-size: 28px; margin: 0;"),
      if (!is.null(subtitle)) {
        p(subtitle, style = "margin: 10px 0 0 0; opacity: 0.8;")
      }
  )
}

# Function to create a modern well panel
create_well_panel <- function(..., style = NULL) {
  wellPanel(
    style = paste("padding: 20px;", style),
    ...
  )
}

# Function to create a modern plot container
create_plot_container <- function(..., style = NULL) {
  div(class = "plot-container",
      style = style,
      ...
  )
}

# Function to create a modern datatable container
create_datatable_container <- function(..., style = NULL) {
  div(class = "datatable-container",
      style = style,
      ...
  )
}

# Function to create a modern button
create_button <- function(inputId, label, class = "btn-default", ...) {
  actionButton(inputId, label, class = class, ...)
}

# Function to create a modern select input
create_select_input <- function(inputId, label, choices, selected = NULL, ...) {
  selectizeInput(inputId, label, choices = choices, selected = selected, ...)
}

# Function to create a modern numeric input
create_numeric_input <- function(inputId, label, value, min = NULL, max = NULL, step = NULL, ...) {
  numericInput(inputId, label, value = value, min = min, max = max, step = step, ...)
}

# Function to create a modern slider input
create_slider_input <- function(inputId, label, min, max, value, step = NULL, ...) {
  sliderInput(inputId, label, min = min, max = max, value = value, step = step, ...)
}

# Function to create a modern download button
create_download_button <- function(outputId, label, class = "btn-default", ...) {
  downloadButton(outputId, label, class = class, ...)
}

# Function to create a modern error message container
create_error_message <- function(outputId) {
  div(id = "error_message_container", 
      style = "margin-bottom: 15px; color: #e74c3c; font-weight: bold; display: none;",
      textOutput(outputId)
  )
}

# Function to create a modern success message container
create_success_message <- function(outputId) {
  div(id = "success_message_container",
      class = "success-message",
      textOutput(outputId)
  )
}

# Function to create a modern lever switch
create_lever_switch <- function(inputId, label = NULL) {
  div(style = "display: flex; align-items: center; gap: 10px;",
      if (!is.null(label)) {
        span(label, class = "control-label")
      },
      div(class = "lever-switch", id = inputId,
          onclick = sprintf("Shiny.setInputValue('%s', Date.now())", inputId))
  )
}

# Function to create a modern plot output with spinner
create_plot_output <- function(outputId, width = "100%", height = "auto", spinner_type = 8, spinner_color = "#3498db") {
  plotlyOutput(outputId, width = width, height = height) %>%
    withSpinner(type = spinner_type, color = spinner_color)
}

# Function to create a modern datatable output
create_datatable_output <- function(outputId, caption = NULL) {
  DTOutput(outputId) %>%
    tagAppendAttributes(
      class = "datatable-container",
      style = if (!is.null(caption)) {
        "caption-side: top; text-align: left; color: #2c3e50; font-weight: bold; font-size: 14px;"
      }
    )
}

# Function to create a modern tab panel
create_tab_panel <- function(title, ..., id = NULL) {
  tabPanel(title, ..., id = id)
}

# Function to create a modern tabset panel
create_tabset_panel <- function(..., id = NULL) {
  tabsetPanel(..., id = id)
}

# Function to create a modern fluid row
create_fluid_row <- function(...) {
  fluidRow(...)
}

# Function to create a modern column
create_column <- function(width, ...) {
  column(width, ...)
}

# Function to create a modern fluid page
create_fluid_page <- function(..., theme = bs_theme(
  version = 4,
  bootswatch = "flatly",
  primary = "#3498db",
  secondary = "#2c3e50",
  success = "#18bc9c",
  info = "#3498db",
  warning = "#f39c12",
  danger = "#e74c3c"
)) {
  fluidPage(
    tags$head(tags$style(custom_css)),
    theme = theme,
    ...
  )
} 