FROM rocker/r-ver:4.0.0

RUN apt-get update && apt-get install -yq \
    libcurl4-openssl-dev \
    libhdf5-dev \
    libssl-dev \
    zlib1g-dev

RUN Rscript -e 'install.packages("devtools")'
# BH is the cran equivalent of boost-cpp from conda
RUN Rscript -e 'devtools::install_version("BH", version = "1.75.0-0", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'devtools::install_version("tibble", version = "3.1.5", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'devtools::install_version("tidyr", version = "1.1.4", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'devtools::install_version("testthat", version = "3.1.0", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'devtools::install_version("dplyr", version = "1.0.7", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'devtools::install_version("readr", version = "2.0.2", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'devtools::install_version("purrr", version = "0.3.4", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'devtools::install_version("flashClust", version = "1.01-2", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'devtools::install_version("data.table", version = "1.14.2", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'devtools::install_version("RANN", version = "2.6.1", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'devtools::install_version("pastecs", version = "1.3.21", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'devtools::install_version("gplots", version = "3.1.1", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'devtools::install_version("entropy", version = "1.3.1", repos = "http://cran.us.r-project.org")'
# rlist errors when installing through devtools
RUN Rscript -e 'install.packages("rlist")'
# WGCNA errors when installing through devtools
RUN Rscript -e 'devtools::install_version("BiocManager", version = "1.30.16", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install(version = "3.11", ask = FALSE)'
RUN Rscript -e 'BiocManager::install("WGCNA")'
# rcpp errors when installing through devtools
RUN Rscript -e 'install.packages("rcpp")'
RUN Rscript -e 'devtools::install_version("arrow", version = "5.0.0.2", repos = "http://cran.us.r-project.org")'
# rcdk requires java and also errors when installing through devtools
RUN apt-get install -yq default-jdk
RUN Rscript -e 'install.packages("rcdk")'
RUN Rscript -e 'devtools::install_version("patrick", version = "0.1.0", repos = "http://cran.us.r-project.org")'

ADD xmsannotator /xmsannotator
RUN R CMD INSTALL /xmsannotator

# ONLY image size reduction commands below
RUN find / -name '*.log' -delete
RUN find / -name '.cache' -delete
RUN rm -rf /var/lib/apt/lists/*
RUN rm -rf /var/cache/*

# https://askubuntu.com/questions/266738/how-to-truncate-all-logfiles
RUN truncate -s 0 /var/log/*log || true
RUN truncate -s 0 /var/log/**/*log || true

RUN apt -y autoremove
RUN apt -y clean
