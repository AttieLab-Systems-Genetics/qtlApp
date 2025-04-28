# install_packages.R

# Set CRAN mirror
options(repos = c(CRAN = 'http://cran.rstudio.com/'))

# Install BiocManager (needed even if some Bioc packages are pre-installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
library(BiocManager)

# List of required packages NOT expected in rocker/shiny-verse base
# Assumes shiny, tidyverse, dplyr, ggplot2, stringr, data.table, reshape2 are present
pkgs <- c(
    # 'lme4', 'pbkrtest', 'car', 'rstatix', 'ragg', # Dependencies likely covered by shiny-verse
    # 'tidyverse', # Included in shiny-verse
    'ggpubr', 'plumber', # Keep ggpubr, plumber might not be included
    # 'rstudioapi', # Likely included
    # 'dplyr', 'stringr', 'ggplot2', # Included in shiny-verse
    # 'grid', # Base R
    'ggrepel', 'gridGraphics', 
    # 'shiny', # Included in shiny-verse
    'shinyFiles', 'bslib', 
    'spsComps', 'DT', 'shinyjs', 'shinycssloaders', 
    # 'data.table', # Included in shiny-verse 
    # 'reshape2', # Included in shiny-verse
    'plotly', 'ggiraph', 'writexl', 'fontawesome', 
    'fst', 'R.utils', 
    'qtl2' # Bioconductor package
)

# Install packages using BiocManager
print("--- Installing additional required packages ---")
BiocManager::install(pkgs, update = FALSE, ask = FALSE)

# Final check (only check packages we explicitly tried to install)
print("--- Final package check ---")
installed_pkgs <- installed.packages()[, "Package"]
required_pkgs_unique <- unique(pkgs)
missing_pkgs <- required_pkgs_unique[!(required_pkgs_unique %in% installed_pkgs)]

if (length(missing_pkgs) > 0) {
    stop(paste("Failed to install required packages:", paste(missing_pkgs, collapse=", ")))
} else {
    print("--- Additional required packages installed successfully ---")
} 