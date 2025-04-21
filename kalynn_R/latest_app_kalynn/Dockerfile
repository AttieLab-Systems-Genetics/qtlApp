# Use a base image that allows specifying R version
#FROM mini-viewer-image_v1
FROM rocker/shiny-verse
#FROM rocker/r-ver:4.1.1
#FROM rocker/shiny:4.3.1
#FROM dckremfinger/shiny_mediation:latest
# Install necessary system packages
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    && apt-get clean

# Install required R packages from CRAN
RUN R -e "install.packages(c(\
    'rstudioapi', 'dplyr', 'stringr', 'tidyverse', 'ggplot2', \
    'grid', 'ggrepel', 'gridGraphics', 'ggpubr', 'shiny', \
    'shinyFiles', 'bslib', 'spsComps', 'DT', 'shinyjs', \
    'shinycssloaders', 'data.table', 'reshape2', 'plotly', \
    'ggiraph', 'writexl', 'fontawesome', 'fst', 'debounce', \
    'R.utils', 'plumber', 'shinyjs'), repos='http://cran.rstudio.com/')"

# Install Bioconductor packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('qtl2', update = FALSE, ask = FALSE)"

# Create the data directory with proper permissions
RUN mkdir -p /data
RUN chmod -R 777 /data
#COPY ./example_file_directory_forviewer/file_index.csv /data
#COPY ./example_file_directory_forviewer/chromosomal_sep_mm11.csv /data
#COPY ./example_file_directory_forviewer/revised_markers_v3.rds /data
#COPY ./example_file_directory_forviewer/annotation_list.rds /data
#COPY ./example_file_directory_forviewer/file_index.xlsx /data
#COPY ./example_file_directory_forviewer/HC_results_sexplusgenlit_withinfo.csv /data
#COPY ./example_file_directory_forviewer/HF_scans_sexgenlitadditivelivgenes.csv /data
#COPY ./example_file_directory_forviewer/small_HC_csv_fortest.csv /data
#COPY ./example_file_directory_forviewer/small_HF_csv_fortest.csv /data

# Copy your Shiny app to the image
#COPY ./brian_app.R /srv/shiny-server/app.R
#COPY ./brian_setup.R /srv/shiny-server/
#COPY ./kalynn_qtlviz.R /srv/shiny-server/
#COPY ./chris_qtlviz.R /srv/shiny-server/
# Ensure the necessary permissions
COPY app.R /srv/shiny-server/app.R
RUN chown -R shiny:shiny /srv/shiny-server

# Expose the port
EXPOSE 3838

# Run the Shiny server
#CMD ["/usr/bin/shiny-server"]

