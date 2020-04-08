#' Table of pre-computed critical values
#'
#' @format A data frame with 4 variables:
#' \describe{
#'
#' \item{\code{B}}{B Bound on the square root of the average squared normalized
#' bias}
#'
#' \item{\code{kappa}}{Bound on the kurtosis of the normalized bias}
#'
#' \item{\code{cv}}{Critical value}
#'
#' \item{\code{alpha}}{Determines confidence level, \eqn{1-\alpha}{1-alpha}.}
#'
#' }
#' @source Computed using the \code{cva} function in this package
"cva_tbl"

#' Neighborhood effects data from Chetty and Hendren (2018)
#'
#' This dataset contains a subset of the publicly available data from Chetty and
#' Hendren (2018). It contains raw estimates and standard errors of neighborhood
#' effects at the commuting zone level
#' @format A data frame with 741 rows corresponding to commuting zones (CZ) and
#'     10 columns corresponding to the variables:
#'
#' \describe{
#'
#' \item{cz}{Commuting zone ID}
#'
#' \item{czname}{Name of CZ}
#'
#' \item{state}{2-digit state code}
#'
#' \item{pop}{Population according to the year 2000 Census}
#'
#' \item{theta25}{Fixed-effect estimate of the causal effect of living in the CZ
#' for one year on children's percentile rank in the national distribution of
#' household earnings at age 26 relative to others in the same birth cohort for
#' children growing up with parents at the 25th percentile of national income
#' distribution}
#'
#' \item{theta75}{Fixed-effect estimate of the causal effect of living in the CZ
#' for one year on children's percentile rank in the national distribution of
#' household earnings at age 26 relative to others in the same birth cohort for
#' children growing up with parents at the 75th percentile of national income
#' distribution}
#'
#' \item{se25}{Standard error of \code{theta25}}
#'
#' \item{se75}{Standard error of \code{theta75}}
#'
#' \item{stayer25}{Average percentile rank in the national distribution of
#' household earnings at age 26 relative to others in the same birth cohort for
#' stayers (children who grew up in the CZ and did not move) with parents at the
#' 25th percentile of national income distribution.}
#'
#' \item{stayer75}{Average percentile rank in the national distribution of
#' household earnings at age 26 relative to others in the same birth cohort for
#' stayers (children who grew up in the CZ and did not move) with parents at the
#' 75th percentile of national income distribution.}
#'
#' }
#' @source \url{https://opportunityinsights.org/data/?paper_id=599}

#' @references{
#'
#' \cite{Chetty, R., & Hendren, N. (2018). The Impacts of Neighborhoods on
#' Intergenerational Mobility II: County-Level Estimates. The Quarterly Journal
#' of Economics, 133(3), 1163–1228. \url{https://doi.org/10.1093/qje/qjy006}}
#'
#' }
"cz"

#' Neighborhood effects data from Chetty and Hendren (2018)
#'
#' This dataset contains a subset of the publicly available data from Chetty and
#' Hendren (2018). It contains raw estimates and standard errors of neighborhood
#' effects at the county level
#' @format A data frame with 741 rows corresponding to counties and 12 columns
#'     corresponding to the variables:
#'
#' \describe{
#'
#' \item{ct}{County ID}
#'
#' \item{cz}{ID of commuting zone to which the county belongs.}
#'
#' \item{ctname}{County name}
#'
#' \item{state}{2-digit state code}
#'
#' \item{pop}{Population of the county according to the year 2000 Census}
#'
#' \item{pop_cz}{Population of the CZ to which the county belongs according to
#' the year 2000 Census}
#'
#' \item{theta25}{Fixed-effect estimate of the causal effect of living in the
#' county for one year on children's percentile rank in the national
#' distribution of household earnings at age 26 relative to others in the same
#' birth cohort for children growing up with parents at the 25th percentile of
#' national income distribution}
#'
#' \item{theta75}{Fixed-effect estimate of the causal effect of living in the
#' county for one year on children's percentile rank in the national
#' distribution of household earnings at age 26 relative to others in the same
#' birth cohort for children growing up with parents at the 75th percentile of
#' national income distribution}
#'
#' \item{se25}{Standard error of \code{theta25}}
#'
#' \item{se75}{Standard error of \code{theta75}}
#'
#' \item{stayer25}{Average percentile rank in the national distribution of
#' household earnings at age 26 relative to others in the same birth cohort for
#' stayers (children who grew up in the county and did not move) with parents at
#' the 25th percentile of national income distribution.}
#'
#' \item{stayer75}{Average percentile rank in the national distribution of
#' household earnings at age 26 relative to others in the same birth cohort for
#' stayers (children who grew up in the county and did not move) with parents at
#' the 75th percentile of national income distribution.}
#'
#' }
#' @source \url{https://opportunityinsights.org/data/?paper_id=599}

#' @references{
#'
#' \cite{Chetty, R., & Hendren, N. (2018). The Impacts of Neighborhoods on
#' Intergenerational Mobility II: County-Level Estimates. The Quarterly Journal
#' of Economics, 133(3), 1163–1228. \url{https://doi.org/10.1093/qje/qjy006}}
#'
#' }
"ct"
