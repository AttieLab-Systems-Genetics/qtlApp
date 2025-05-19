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

# Source files in specific order
source("R/helpers.R")
source("R/importApp.R")
source("R/mainParApp.R")
source("R/scanApp.R")
source("R/peakApp.R")
source("R/scanlyApp.R") # Might not be used if qtlApp calls scanApp, but good to keep
source("R/qtlApp.R") # Source qtlApp itself
source("R/mergeApp.R") # Source mergeApp
source("R/downloadApp.R") # Source downloadApp
source("R/cisTransPlotApp.R") # Source the new module
source("R/ui_styles.R") # Source styles
source("R/plot_enhancements.R") # Source plot enhancements
source("R/plot_null.R") # Source plot null
source("R/data_handling.R") # Source data handling
source("R/peak_finder.R") # Source peak finder
source("R/trait_scan.R") # Source trait scan
source("R/ggplot_alleles.R") # Source ggplot alleles
source("R/ggplot_qtl_scan.R") # Source ggplot qtl scan
source("R/peak_info.R") # Source peak info
source("R/QTL_plot_visualizer.R") # Source plot visualizer
source("R/csv2fst.R") # Source csv2fst
source("R/fst_rows.R") # Source fst_rows
source("R/traitApp.R") # Source traitApp
source("R/traitTypeApp.R") # Source traitTypeApp

# Set maximum file upload size
options(shiny.maxRequestSize = 20000*1024^2) # 20 GB

# Launch the application (assuming scanApp is the desired main app for now)
scanApp()