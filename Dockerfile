FROM rocker/r-ver:3.6.3

RUN apt-get update && apt-get install -yq \
    libcurl4-openssl-dev \
    libjpeg-dev \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
 
# RUN Rscript \
#     -e 'install.packages("remotes")' \
#     -e 'remotes::install_version("doSNOW", version = "1.0.18")' \
#     -e 'remotes::install_version("flashClust", version = "1.01-2")' \
#     -e 'remotes::install_version("foreach", version = "1.5.0")' \
#     -e 'remotes::install_version("plyr", version = "1.8.6")' \
#     -e 'remotes::install_version("png", version = "0.1-7")' \
#     -e 'remotes::install_version("R2HTML", version = "2.3.2")' \
#     -e 'remotes::install_version("RCurl", version = "1.98-1.2")' \
#     -e 'remotes::install_version("rjson", version = "0.2.20")' \
#     -e 'remotes::install_version("snow", version = "0.4-3")' \
#     -e 'remotes::install_version("XML", version = "3.99-0.3")' \
#     -e 'remotes::install_version("BiocManager", version = "1.30.10")' \
#     -e 'BiocManager::install(version = "3.10")' \
#     -e 'BiocManager::install("KEGGREST")' \
#     -e 'BiocManager::install("pcaMethods")' \
#     -e 'BiocManager::install("Rdisop")' \
#     -e 'remotes::install_version("WGCNA", version = "1.69", repos = BiocManager::repositories())' \
#     -e 'remotes::install_version("SSOAP", version = "0.9-1", repos = c(getOption("repos"), "http://www.omegahat.net/R"))'

# RUN Rscript \
#     -e 'remotes::install_version("arrow", version = "0.17.0")' \
#     -e 'install.packages("purrr")' \
#     -e 'install.packages("furrr")' \
#     -e 'BiocManager::install("rhdf5 ")'

RUN Rscript \
    -e 'install.packages("dplyr")' \
    -e 'install.packages("flashClust")' \
    -e 'install.packages("WGCNA")' \
    -e 'install.packages("tidyr")' \
    -e 'install.packages("purr")' \
    -e 'install.packages("tibble")' \
    -e 'install.packages("arrow")' \
    -e 'install.packages("BiocManager")'
    -e 'BiocManager::install("rhdf5 ")'

ADD xmsannotator /xmsannotator
RUN R CMD INSTALL /xmsannotator \
 && rm -rf /xmsannotator
