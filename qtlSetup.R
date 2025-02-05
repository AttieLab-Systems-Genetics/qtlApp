#' This is the viewer script to view previously made qtl2 scans for various traits
#' it's designed to look at additive, interactive, and difference (interactive-additive) plots
#' all the scan traces are pre-computed, so this just reads the table of relevant values
#' and plots them
#'
#' @author Chris Emfinger, PhD. Couldn't really view anything well
#'
#' In addition to the relevant package citations, parts of the datatable code
#' followed formats shown on
#' https://clarewest.github.io/blog/post/making-tables-shiny/x
#'
#' some help with shiny
#' https://shiny.rstudio.com/gallery
#' https://shiny.posit.co/r/gallery
#' posit.cloud/spaces/298214/join

# load libraries=============================================================
# Now done through explicit reference package::function()
# require("rstudioapi")
# require("dplyr")
# require("stringr")
# require("tidyverse")
# require("BiocManager")
# require("ggplot2")
# require("qtl2")
# require("grid")
# require("ggrepel")
# require("gridGraphics")
# require("ggpubr")
# require("shiny")
# require("shinyFiles")
# require("bslib")
# require("spsComps")
# require("DT")
# require("shinyjs")
# require("shinycssloaders")
# require("data.table")
# require("reshape2")

# set shiny options==========================================================================
options(shiny.maxRequestSize = 20000*1024^2)  # Increase to 20GB needed for Genoprobs, etc

# load data==================================================================================
# load the file directory
#setwd(selectDirectory())

# load data==================================================================================
# load the file directory
#setwd(selectDirectory())
file_directory <- read.csv("/data/dev/miniViewer_3.0/file_index.csv")

# load the chromosomal breaks (not used?)
chr_breaks <- read.csv("/data/dev/miniViewer_3.0/chromosomal_sep_mm11.csv")

# load the annotation data for the different traits
# for making the object
# annotation_list <- list()
# annotation_list$isoforms <- mediation_isoforms$RNA_info$Liver[c(1,3,ncol(mediation_isoforms$RNA_info$Liver))]
# annotation_list$genes <- mediation_genes$RNA_info$Liver[c(1,6)]
# saveRDS(annotation_list, file="annotation_list.rds")
annotation_list <- readRDS("/data/dev/miniViewer_3.0/annotation_list.rds")

# load the markers
markers <- readRDS(file.path("/data/dev/miniViewer_3.0/CHTC_dietDO_markers_RDSgrcm39.rds"))

file_directory$group <- paste0(
    file_directory$diet," ",
    file_directory$trait_compartment, " ",
    file_directory$trait_type, ", ",
    file_directory$scan_type)

# source functions
source("trait_scan.R")
source("peak_finder.R")
source("QTL_plot_visualizer.R")
