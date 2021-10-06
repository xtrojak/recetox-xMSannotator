# Developer Documentation

## Setup

Use the provided conda environment files to setup the project with the followng command: `conda env create -f conda/environment.yaml`.
This will install all dependencies except for `dataCompareR` and `patrick`. Those can be installed in R using `install.packages(c('dataCompareR', 'patrick')).

## Code Style
This package is developed according to the `tidyverse` code formatting style. In order to have a uniform formatting of the source code, please use the [styler](https://styler.r-lib.org/reference/style_pkg.html) package for auto-formatting, as different IDEs might give different results. The `style_pkg(...)` function can be used to reformat a whole package.
