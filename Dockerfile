FROM continuumio/miniconda3:4.10.3
ADD conda /conda
RUN conda init bash
RUN conda env update -n base -f /conda/environment-dev.yaml
ADD xmsannotator /xmsannotator
RUN R CMD INSTALL /xmsannotator
