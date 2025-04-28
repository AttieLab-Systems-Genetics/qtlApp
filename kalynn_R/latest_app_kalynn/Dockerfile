# Use a base image that allows specifying R version
#FROM mini-viewer-image_v1
FROM rocker/shiny-verse:latest
#FROM rocker/r-ver:4.1.1
#FROM rocker/shiny:4.3.1
#FROM dckremfinger/shiny_mediation:latest

USER root

# Install necessary *additional* system packages (shiny-verse base has many already)
# Keep Fortran/BLAS/LAPACK just in case qtl2 or others need them for compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
    gfortran \
    libblas-dev \
    liblapack-dev \
    cmake \
    libsodium-dev \
    libcurl4-gnutls-dev \
    # libcairo2-dev # Likely included in shiny-verse
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    # libpng-dev \ # Likely included
    # libjpeg-dev \ # Likely included
    # libtiff5-dev \ # Likely included
    # libfreetype6-dev \ # Likely included
    # libharfbuzz-dev \ # Likely included
    # libfribidi-dev \ # Likely included
    # libfontconfig1-dev \ # Likely included
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages using the script (will modify script next)
COPY install_packages.R /tmp/install_packages.R
RUN Rscript /tmp/install_packages.R && rm /tmp/install_packages.R

# Removed as qtl2 is now included above

# Force app logs to stderr (hopefully)
ENV SHINY_LOG_STDERR=1

# Create the data directory with proper permissions
RUN mkdir -p /data
RUN chmod -R 777 /data

# Copy your Shiny app to the image
COPY app.R /srv/shiny-server/app.R
COPY R/ /srv/shiny-server/R/
RUN chown -R shiny:shiny /srv/shiny-server

# Expose the port
EXPOSE 3838

USER shiny

# Run the Shiny server
CMD ["/usr/bin/shiny-server"]
