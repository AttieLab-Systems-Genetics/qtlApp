# Use a base image that allows specifying R version
FROM rocker/shiny-verse:latest

USER root

# Install necessary *additional* system packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    gfortran \
    libblas-dev \
    liblapack-dev \
    cmake \
    libsodium-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages using the script
COPY install_packages.R /tmp/install_packages.R
RUN Rscript /tmp/install_packages.R && rm /tmp/install_packages.R

# Force app logs to stderr
ENV SHINY_LOG_STDERR=1

# Create the data directory with proper permissions
RUN mkdir -p /data
RUN chmod -R 777 /data

# Copy the application files
COPY app.R /srv/shiny-server/
COPY R /srv/shiny-server/R/

# Set ownership of the app directory
RUN chown -R shiny:shiny /srv/shiny-server

# Expose the port Shiny Server runs on
EXPOSE 3838

USER shiny

# Run the Shiny server
CMD ["/usr/bin/shiny-server"]