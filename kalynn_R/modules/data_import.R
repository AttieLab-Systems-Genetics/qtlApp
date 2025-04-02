# data_import.R
# Load required data files

file_directory <- read.csv("/data/dev/miniViewer_3.0/file_index.csv")
chr_breaks <- read.csv("/data/dev/miniViewer_3.0/chromosomal_sep_mm11.csv")
annotation_list <- readRDS("/data/dev/miniViewer_3.0/annotation_list.rds")
markers <- readRDS(file.path("/data/dev/miniViewer_3.0/CHTC_dietDO_markers_RDSgrcm39.rds"))

# Create a group identifier for later use
file_directory$group <- paste0(file_directory$diet, " ", 
                               file_directory$trait_compartment, " ",
                               file_directory$trait_type, ", ", 
                               file_directory$scan_type)
