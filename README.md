brmsish: miscellaneous functions to customize brms regression models
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![lifecycle](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Dependencies](https://tinyverse.netlify.com/badge/brmsish)](https://cran.r-project.org/package=brmsish)
[![Project Status: Active The project has reached a stable, usable state
and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Licence](https://img.shields.io/badge/https://img.shields.io/badge/licence-GPL--2-blue.svg.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-%3E=%203.5.0-6666ff.svg)](https://cran.r-project.org/)[![packageversion](https://img.shields.io/badge/Package%20version-1.0.0-orange.svg?style=flat-square)](commits/develop)[![Last-changedate](https://img.shields.io/badge/last%20change-2023--02--24-yellowgreen.svg)](/commits/master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/brmsish)](https://cran.r-project.org/package=brmsish)
[![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/brmsish)](https://cranlogs.r-pkg.org/badges/grand-total/brmsish)

This package offers functions to visualize and/or modify
[brms](https://paul-buerkner.github.io/brms/index.html) models. It’s a
collection of functions that I have been accumulating as I have been
trying to customize the output of
[brms](https://paul-buerkner.github.io/brms/index.html) models. So far
the package has functions for:

- Running multiple comparisons (e.g. contrasts) between levels of a
  categorical variable
- Plotting posterior distributions, trace plots and diagnostics plots in
  a Rmarkdown friendly format
- Running models on a population of phylogenetic trees to account for
  phylogenetic uncertainty
- Combining a bunch of models saved as RDS files

<!-- Unlike other packages for setting up research compendiums, `brmsish` has very simple functionality. Hence, users can focus on the research project itself rather than on learning how to use a new R package. -->

To install the latest developmental version from
[github](https://github.com/) you will need the R package
[remotes](https://cran.r-project.org/package=remotes) (although it may
not install at all):

``` r

# From github
remotes::install_github("maRce10/brmsish")

# load package
library(brmsish)
```

The package documentation can be found
[here](https://marce10.github.io/brmsish/).

Please cite [brmsish](https://marce10.github.io/brmsish/) as follows:

Araya-Salas (2022), *brmsish: Miscellaneous functions to customize brms
bayesian regression models*. R package version 1.0.0.

# References

1.  Buerkner, P. C. (2017). brms: An R package for Bayesian multilevel
    models using Stan. Journal of statistical software, 80, 1-28.

2.  Buerkner, P. C. (2017). Advanced Bayesian multilevel modeling with
    the R package brms. arXiv preprint arXiv:1705.11123.

3.  Nalborczyk, L., Batailler, C., Lœvenbruck, H., Vilain, A., &
    Buerkner, P. C. (2019). An introduction to Bayesian multilevel
    models using brms: A case study of gender effects on vowel
    variability in standard Indonesian. Journal of Speech, Language, and
    Hearing Research, 62(5), 1225-1242.
