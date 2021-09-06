# Developer Documentation

## Setup

Use the provided conda environment files to setup the project with the followng command: `conda env create -f conda/environment.yaml`.
This will install all dependencies except for `dataCompareR` and `patrick`. Those can be installed in R using `install.packages(c('dataCompareR', 'patrick')).