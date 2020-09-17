FROM rocker/r-ver:3.6.3

RUN apt-get update && apt-get install -yq \
    libcurl4-openssl-dev \
    libjpeg-dev \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages("tidyverse")'
RUN Rscript -e 'install.packages("flashClust")'
RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'BiocManager::install("rhdf5")'
RUN Rscript -e 'install.packages("data.table")'
RUN Rscript -e 'BiocManager::install("WGCNA")'
RUN Rscript -e 'install.packages("arrow")'

ADD xmsannotator /xmsannotator
RUN R CMD INSTALL /xmsannotator \
 && rm -rf /xmsannotator
