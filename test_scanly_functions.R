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
cat("Sourcing R files...\n")
source("R/importApp.R")
source("R/mainParApp.R")
source("R/scanApp.R")
source("R/peakApp.R")
source("R/scanlyApp.R")

# Create a test function to verify plot generation
test_scanly_plot <- function() {
  cat("Testing scanly plot generation...\n")
  
  # Create a minimal test environment
  test_env <- new.env()
  
  # Create mock data
  test_data <- data.frame(
    chr = c("1", "1", "2", "2"),
    pos = c(10, 20, 15, 25),
    lod = c(2.5, 3.1, 2.8, 3.5)
  )
  
  # Test plot creation
  tryCatch({
    p <- ggplot2::ggplot(test_data, ggplot2::aes(x = pos, y = lod)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~chr)
    
    # Try converting to plotly
    p_plotly <- plotly::ggplotly(p)
    
    cat("✓ Basic plot creation successful\n")
    cat("✓ Plotly conversion successful\n")
    
    # Print plot structure
    cat("\nPlot structure:\n")
    str(p)
    
  }, error = function(e) {
    cat("Error in plot generation:", conditionMessage(e), "\n")
  })
}

# Run the test
test_scanly_plot() 