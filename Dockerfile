FROM rocker/r-ver:3.6.2

RUN apt-get update \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*
