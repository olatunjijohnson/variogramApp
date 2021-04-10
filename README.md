
<!-- README.md is generated from README.Rmd. Please edit that file -->

# variogramApp

<!-- badges: start -->

<!-- badges: end -->

The goal of variogramApp is to allows user to understand variogram and
perform simulation in one and two
dimensions

## Installation

<!-- You can install the released version of variogramApp from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("variogramApp") -->

<!-- ``` -->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("olatunjijohnson/variogramApp", ref="master")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(variogramApp)
## basic example code
## run the App
run_variog_app()  # use the code
```

## Alternative way to run in R

You can also run the following line of code to run in
R

``` r
shiny::runGitHub(repo="variogramApp", username= "olatunjijohnson", ref="master", subdir = "inst/variogramApp")
```
