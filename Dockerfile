FROM rocker/r-base:4.0.0

RUN apt-get update
RUN apt-get install libz-dev libxml2 libboost-all
ADD . /xmsannotator

RUN R -e "install.packages(c('BiocManager','data.table','digest', 'remotes'))"
RUN R -e "remotes::install_github('omegahat/XMLSchema')"
RUN R -e "remotes::install_github('cran/SSOAP')"
RUN R -e "BiocManager::install(c('KEGGREST','pcaMethods','Rdisop','GO.db','matrixStats','WGCNA'))"
RUN R CMD INSTALL /xmsannotator