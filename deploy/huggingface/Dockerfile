FROM rocker/shiny:4.4.1

# Install R packages using Posit Package Manager pre-built binaries (fast)
RUN R -e "install.packages( \
    c('geoR', 'ggplot2', 'dplyr'), \
    repos = 'https://packagemanager.posit.co/cran/__linux__/jammy/latest' \
  )"

# Configure shiny-server to listen on port 7860 (required by Hugging Face)
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

# Copy the Shiny app
COPY app/ /srv/shiny-server/

EXPOSE 7860

CMD ["/usr/bin/shiny-server"]
