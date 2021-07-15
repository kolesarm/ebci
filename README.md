[![R-CMD-check](https://github.com/kolesarm/ebci/workflows/R-CMD-check/badge.svg)](https://github.com/kolesarm/ebci/actions) [![Coverage status](https://codecov.io/gh/kolesarm/ebci/branch/master/graph/badge.svg)](https://codecov.io/github/kolesarm/ebci?branch=master)

# ebci

This R package implements robust empirical Bayes confidence intervals from
[Armstrong, Kolesár, and Plagborg-Møller
(2020)](https://arxiv.org/abs/2004.03448). See the package
[manual](doc/manual.pdf) for documentation of the package functions, and the
[package vignette](doc/ebci.pdf) for a description of the package and an example
of its usage (available through `vignette("ebci")` once the package is
installed). See [ebci_matlab](https://github.com/mikkelpm/ebci_matlab) for a
Matlab version of this package.

## Installation

You can get the current development version from GitHub:

``` r
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("kolesarm/ebci")
```

## Example

Calculation of the critical values used to construct robust empirical Bayes
confidence intervals:

``` r
library("ebci")
## Only second moment constraint on normalized bias, m2=4.
cva(m2=4, kappa=Inf, alpha=0.05)
## Add a constraint that kurtosis equals 3
cva(m2=4, kappa=3, alpha=0.05)
```
