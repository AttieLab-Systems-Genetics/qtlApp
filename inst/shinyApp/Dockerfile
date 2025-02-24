# Use a base image that allows specifying R version
#FROM mini-viewer-image_v1
FROM rocker/shiny-verse
#FROM rocker/r-ver:4.1.1
#FROM rocker/shiny:4.3.1
#FROM dckremfinger/shiny_mediation:latest
# Install necessary system packages
RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libjpeg-dev \
    libxml2-dev \
    libxt-dev \
    libgit2-dev \
    curl \
    vim \
    wget \
    && apt-get clean

# Install Shiny Server from the correct URL
#RUN wget https://download3.rstudio.org/ubuntu-20.04/x86_64/shiny-server-1.5.20.1002-amd64.deb -O shiny-server.deb || \
#    wget https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-1.5.20.1002-amd64.deb -O shiny-server.deb || \
#    wget https://download3.rstudio.org/ubuntu-16.04/x86_64/shiny-server-1.5.20.1002-amd64.deb -O shiny-server.deb
#RUN gdebi -n shiny-server.deb
#RUN rm shiny-server.deb

# Install CRAN packages directly using install.packages()
RUN R -e "install.packages(c( \
    'rstudioapi', 'dplyr', 'stringr', 'tidyverse', 'ggplot2', 'cowplot', 'grid', 'ggrepel', \
    'gridGraphics', 'ggpubr', 'shiny', 'shinyFiles', 'bslib', 'spsComps', 'DT', 'shinyjs', \
    'shinycssloaders', 'data.table', 'reshape2', 'ggpubr', 'R.utils', \
    'fst', 'stringr' \
), repos = 'https://cloud.r-project.org/')"


# Install Bioconductor packages using BiocManager
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos = 'https://cloud.r-project.org/'); \
    BiocManager::install('qtl2', update = FALSE, ask = FALSE)"

RUN R -e "devtools::install_github('AttieLab-Systems-Genetics/qtlApp')"

# Create a shiny user and group if they don't exist
#RUN getent group shiny || groupadd shiny && \
#    getent passwd shiny || useradd -g shiny shiny

# Create the Shiny server directory
RUN mkdir -p /srv/shiny-server

# Create the data directory and copy data files
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
COPY qtlSetup.R /srv/shiny-server/qtlSetup.R
RUN chown -R shiny:shiny /srv/shiny-server

# Expose the port
EXPOSE 3838

# Run the Shiny server
#CMD ["/usr/bin/shiny-server"]

