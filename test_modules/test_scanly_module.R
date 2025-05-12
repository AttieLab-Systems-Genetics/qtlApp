# Load required packages
library(shiny)
library(bslib)
library(dplyr)
library(purrr)
library(stringr)
library(qtl2)
library(ggplot2)
library(plotly)
library(shinyjs)
library(shinycssloaders)
library(data.table)
library(ggiraph)
library(writexl)
library(fontawesome)
library(lubridate)
library(BiocManager)
library(grid)
library(ggrepel)
library(gridGraphics)
library(ggpubr)
library(shinyFiles)
library(spsComps)
library(DT)
library(reshape2)

# Source the necessary R files
source("R/importApp.R")
source("R/mainParApp.R")
source("R/scanApp.R")
source("R/peakApp.R")
source("R/scanlyApp.R")

# Create a minimal test app
ui <- fluidPage(
  titlePanel("Test Scanly Module"),
  sidebarLayout(
    sidebarPanel(
      # Add any necessary inputs here
      selectInput("test_chr", "Test Chromosome", choices = c("1", "2", "3")),
      numericInput("test_lod", "Test LOD", value = 3.0)
    ),
    mainPanel(
      # Add the scanly module UI
      scanlyOutput("test_scanly")
    )
  )
)

server <- function(input, output, session) {
  # Create mock data and parameters
  mock_import <- reactive({
    list(
      data = data.frame(
        chr = c("1", "1", "2", "2"),
        pos = c(10, 20, 15, 25),
        lod = c(2.5, 3.1, 2.8, 3.5)
      )
    )
  })
  
  mock_main_par <- reactive({
    list(
      selected_chr = reactive(input$test_chr),
      LOD_thr = reactive(input$test_lod)
    )
  })
  
  # Initialize the scanly module
  scanlyServer("test_scanly", mock_main_par, mock_import, NULL)
}

# Run the app
shinyApp(ui, server) 