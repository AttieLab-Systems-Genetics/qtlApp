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
source("R/helpers.R")
source("R/importApp.R")
source("R/mainParApp.R")
source("R/scanApp.R")
source("R/peakApp.R")
source("R/scanlyApp.R")

# Test function to verify core scanly functionality
test_scanly_core <- function() {
  cat("\nTesting scanly core functionality...\n")
  
  # Create test data
  test_data <- data.frame(
    chr = c("1", "1", "2", "2"),
    pos = c(10, 20, 15, 25),
    lod = c(2.5, 3.1, 2.8, 3.5)
  )
  
  # Test 1: Basic plot creation
  cat("\nTest 1: Basic plot creation\n")
  tryCatch({
    p <- ggplot2::ggplot(test_data, ggplot2::aes(x = pos, y = lod)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~chr)
    cat("✓ Basic plot creation successful\n")
  }, error = function(e) {
    cat("✗ Basic plot creation failed:", conditionMessage(e), "\n")
  })
  
  # Test 2: Plotly conversion
  cat("\nTest 2: Plotly conversion\n")
  tryCatch({
    p_plotly <- plotly::ggplotly(p)
    cat("✓ Plotly conversion successful\n")
  }, error = function(e) {
    cat("✗ Plotly conversion failed:", conditionMessage(e), "\n")
  })
  
  # Test 3: Event registration
  cat("\nTest 3: Event registration\n")
  tryCatch({
    p_plotly <- plotly::event_register(p_plotly, 'plotly_click')
    p_plotly <- plotly::event_register(p_plotly, 'plotly_doubleclick')
    cat("✓ Event registration successful\n")
  }, error = function(e) {
    cat("✗ Event registration failed:", conditionMessage(e), "\n")
  })
  
  # Test 4: Peak finding
  cat("\nTest 4: Peak finding\n")
  tryCatch({
    peaks <- data.frame(
      marker = c("marker1", "marker2"),
      chr = c("1", "2"),
      pos = c(15, 20),
      lod = c(3.1, 3.5)
    )
    ordered_peaks <- highest_peaks(peaks, 3.0)
    cat("✓ Peak finding successful\n")
    cat("  Found peaks:", paste(ordered_peaks$marker, collapse=", "), "\n")
  }, error = function(e) {
    cat("✗ Peak finding failed:", conditionMessage(e), "\n")
  })
}

# Run the tests
test_scanly_core() 