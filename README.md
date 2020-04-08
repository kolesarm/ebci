[![Travis build status](https://travis-ci.org/kolesarm/ebci.svg?branch=master)](https://travis-ci.org/kolesarm/ebci) [![Coverage status](https://codecov.io/gh/kolesarm/ebci/branch/master/graph/badge.svg)](https://codecov.io/github/kolesarm/ebci?branch=master)

# ebci

This R package implements robust empirical Bayes confidence intervals from
[Armstrong, Kolesár, and Plagborg-Møller
(2020)](https://arxiv.org/abs/2004.03448). See the package
[manual](doc/manual.pdf) for documentation of the package functions, an the
[package vignette](doc/ebci.pdf) for a description of the package and an example
of its usage (available through `vignette("ebci")` once the package is
installed).

## Installation

You can get the current development version from GitHub:

``` r
install.packages("remotes") # if not installed
remotes::install_github("kolesarm/ebci")
```

## Example

``` r
library("ebci")
## Only second moment constraints, B=sqrt(m_2)=2.
## Note that B is the square root of the second moment of the normalized bias.
cva(B=2, kappa=Inf, alpha=0.05)
## Add a constraint that kurtosis equals 3
cva(B=2, kappa=3, alpha=0.05)
```
