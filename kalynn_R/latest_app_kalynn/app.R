# Load required R packages
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
library(lubridate) # Added based on potential dependencies in sourced files
library(BiocManager) # Added based on potential dependencies
library(grid)        # Added based on potential dependencies
library(ggrepel)     # Added based on potential dependencies
library(gridGraphics)# Added based on potential dependencies
library(ggpubr)      # Added based on potential dependencies
library(shinyFiles)  # Added based on potential dependencies
library(spsComps)    # Added based on potential dependencies
library(DT)          # Added based on potential dependencies
library(reshape2)    # Added based on potential dependencies
# Function to source all R files in a directory recursively
source_dir_recursive <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$", recursive = TRUE, full.names = TRUE)) {
    if(trace) cat(basename(nm),":")
    source(nm, ...)
    if(trace) cat("\n")
  }
}
# Source shared R files from the R directory relative to the Docker build context
# (assuming Dockerfile copies the top-level R into the container's R directory)
source_dir_recursive("R", trace = FALSE) # Source from copied R directory
# Set maximum file upload size (adjust as needed)
options(shiny.maxRequestSize = 20000*1024^2) # 20 GB
# Launch the application defined in R/qtlApp.R
qtlApp()