#' Optimal shrinkage for Empirical Bayes confidence intervals
#'
#' Compute linear shrinkage factor \eqn{w_{opt}}{w_opt} to minimize the length of the
#' resulting Empirical Bayes confidence interval (EBCI).
#' @param S Square root of the signal-to-noise ratio
#'     \eqn{\sqrt{\mu_{2}}/\sigma}{sqrt(mu_2)/sigma}, where \sqrt{\mu_{2}}{mu_2}
#'     is the variance of \eqn{\theta}{theta}
#' @param kappa Kurtosis of \eqn{\theta}{theta}
#' @param alpha determines CI level
#' @param cv_tbl Optionally, supply a data frame of critical values. (TODO)
#' @return (TODO)
#' @examples
#' w_opt(1, 3, cv_tbl=cva_tbl[cva_tbl$kurt==3 & cva_tbl$alpha==0.05, ])
#' @export
w_opt <- function(S, kappa, alpha=0.05, cv_tbl=NULL) {
    if (!is.null(cv_tbl)) {
        ## Linearly interpolate
        cv <- function(B) {
            idx <- which.max(cv_tbl$B>B)
            cv_tbl$cv[idx-1]+(B-cv_tbl$B[idx-1])/
                (cv_tbl$B[idx]-cv_tbl$B[idx-1]) *
                (cv_tbl$cv[idx]-cv_tbl$cv[idx-1])
        }
        maxbias <- max(cv_tbl$B)
    } else {
        cv <- function(B) cva(B, kurt=kappa, alpha=alpha, check=FALSE)$cv
        ## Maxium bias (i.e B) is 100 to prevent numerical issues
        maxbias <- 100
    }
    ci_length <- function(w) cv((1/w-1)*S)*w
    r <- stats::optimize(ci_length, c(1/(maxbias/S+1), 1), tol=tol)
    B <- (1/r$minimum-1)*S
    ## Recompute critical value, checking solution accuracy. If we're using
    ## cv_tbl, then optimum will be at one of the values of B in the table, so
    ## solution should always be accurate
    len <- cva((1/r$minimum-1)*S, kurt=kappa, alpha=alpha,
                  check=TRUE)$cv*r$minimum
    if (B > maxbias-1e-4)
        warning("Corner solution: optimum reached  maximal bias/sd ratio")
    list(w=r$minimum, length=len, B=B)
}

#' Empirical Bayes estimator and confidence intervals
#'
#' Compute linear shrinkage factor \eqn{w_{eb}}{w_eb} and normalize length of of
#' the associated robust Empirical Bayes confidence interval (EBCI).
#' @inheritParams w_opt
#' @return (TODO)
#' @export
w_eb <- function(S, kappa=Inf, alpha=0.05) {
    w <- S^2/(1+S^2)
    list(w=w, length=cva(1/S, kurt=kappa, alpha=alpha, check=TRUE)$cv*w, B=1/S)
}
