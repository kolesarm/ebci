#' Table of pre-computed critical values
#'
#' @format A data frame with 18 rows and 19 variables:
#' \describe{
#'   \item{B}{B Bound on the square root of the average squared bias}
#'   \item{kurt}{Bound on the kurtosis of the bias}
#'   \item{cv}{Critical value}
#'   \item{alpha}{Determines CI level, \eqn{1-\alpha}{1-alpha}.}
#'  }
#' @source Computed using the \code{cva} function in this package
"cva_tbl"
