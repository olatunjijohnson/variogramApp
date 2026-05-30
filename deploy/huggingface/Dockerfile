FROM eddelbuettel/r2u:22.04

# install2.r uses BSPM: installs as apt .deb where available, CRAN binary otherwise
# System dependencies are handled automatically
RUN install2.r --error shiny geoR ggplot2 dplyr

RUN R -e "library(shiny);   cat('shiny OK\n')"
RUN R -e "library(geoR);    cat('geoR OK\n')"
RUN R -e "library(ggplot2); cat('ggplot2 OK\n')"
RUN R -e "library(dplyr);   cat('dplyr OK\n')"

COPY app/ /app/

EXPOSE 7860

CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=7860, launch.browser=FALSE)"]
