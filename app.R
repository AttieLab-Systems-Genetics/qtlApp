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
library(DT)
library(reshape2)

# Source files in specific order
source("R/helpers.R")
source("R/data_handling.R")
source("R/import_data.R")
source("R/importApp.R")
source("R/mainParApp.R")
source("R/scanApp_monolithic_backup.R")
source("R/peakApp.R")
source("R/scanlyApp.R")
source("R/qtlApp.R")
source("R/mergeApp.R")
source("R/downloadApp.R")
source("R/cisTransPlotApp.R")
source("R/manhattanPlotApp.R")
source("R/ui_styles.R")
source("R/plot_enhancements.R")
source("R/plot_null.R")
source("R/peak_finder.R")
source("R/trait_scan.R")
source("R/ggplot_alleles.R")
source("R/ggplot_qtl_scan.R")
source("R/ggplotly_qtl_scan.R")
source("R/peak_info.R")
source("R/QTL_plot_visualizer.R")
source("R/fst_rows.R")
source("R/traitApp.R")
source("R/traitTypeApp.R")
source("R/traitProcessingModule.R")

# Set maximum file upload size
options(shiny.maxRequestSize = 20000 * 1024^2) # 20 GB

# Launch the application
scanApp() # Comment out the old app
# testGenesApp() # Use the new genes analysis test app
