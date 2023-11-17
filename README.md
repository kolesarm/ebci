[![R-CMD-check](https://github.com/kolesarm/ebci/workflows/R-CMD-check/badge.svg)](https://github.com/kolesarm/ebci/actions) [![Coverage status](https://codecov.io/gh/kolesarm/ebci/branch/master/graph/badge.svg)](https://codecov.io/github/kolesarm/ebci?branch=master) [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ebci)](https://cran.r-project.org/package=ebci)

# ebci

This R package implements robust empirical Bayes confidence intervals from
[Armstrong, Kolesár, and Plagborg-Møller
(2022)](https://doi.org/10.3982/ECTA18597). See the package
[manual](doc/manual.pdf) for documentation of the package functions, and the
[package vignette](doc/ebci.pdf) for a description of the package and an example
of its usage (available through `vignette("ebci")` once the package is
installed). See [ebci_matlab](https://github.com/mikkelpm/ebci_matlab) for a
Matlab version of this package, and
[ebreg](https://github.com/kolesarm/ebciStata) for a Stata version.

This software package is based upon work supported by the National Science
Foundation under grant numbers SES-2049765 (Armstrong), SES-22049356 (Kolesár),
and SES-1851665 (Plagborg-Møller), and by work supported by the Alfred P. Sloan
Research Fellowship (Kolesár).

## Installation

You can install the released version of `ebci` from
[CRAN](https://CRAN.R-project.org/package=ebci) with:

``` r
install.packages("ebci")
```

Alternatively, you can get the current development version from GitHub:
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
## Usual critical value
cva(m2=0, kappa=Inf, alpha=0.05)
## Larger critical value that takes bias into account. Only uses second moment
## constraint on normalized bias, m2=4.
cva(m2=4, kappa=Inf, alpha=0.05)
## Add a constraint that kurtosis equals 3. This tightens the critical value
cva(m2=4, kappa=3, alpha=0.05)
```

Estimates and robust EBCIs for neighborhood effects, as in the empirical
application in [Armstrong, Kolesár, and Plagborg-Møller
(2020)](https://arxiv.org/abs/2004.03448). Shrink fixed-effect estimates of the
neighborhood effects, for children with parents at the 25th percentile of the
income distribution (`theta25`) toward average outcome for permanent residents
(stayers) at the 25th percentile of the income distribution. Use precision
weights proportional to the inverse of the squared standard error of the
fixed-effect estimates.

``` r
r <- ebci(theta25 ~ stayer25, data=cz, se=se25, weights=1/se25^2)
```
