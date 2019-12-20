## TODO: pre-compute cva on a grid from 0 to 5, say, for the given kappa, then
## pass it as an argument.

#' Optimal shrinkage for Empirical Bayes confidence intervals
#'
#' @param S average signal, estimate of standard deviation of
#'     \eqn{\theta}{theta}
#' @param kappa kurtosis of \eqn{\theta}{theta}
#' @param alpha determines CI level
#' @param maxbias maximum bias of estimator, determines an upper bound on the
#'     amount of shrinkage
#' @param cv_tbl optionally, supply a data frame of critical values.
#' @return TODO
#' @examples
#' wopt(1, 3, cv_tbl=cva_tbl[cva_tbl$kurt==3 & cva_tbl$alpha==0.05, ])
#' @export
wopt <- function(S, kappa, alpha=0.05, cv_tbl=NULL, maxbias=10) {

    if (!is.null(cv_tbl)) {
        cv <- function(B) {
            idx <- which.max(cv_tbl$B>B)
            cv_tbl$cv[idx-1]+(B-cv_tbl$B[idx-1])/
                (cv_tbl$B[idx]-cv_tbl$B[idx-1]) *
                (cv_tbl$cv[idx]-cv_tbl$cv[idx-1])
        }
        maxbias <- max(cv_tbl$B)
    } else {
        cv <- function(B) cva(B, kurt=kappa, alpha=alpha, check=FALSE)$cv
    }
    le <- function(w) cv((1/w-1)*S)*w
    r <- stats::optimize(le, c(1/(maxbias/S+1), 1), tol=tol)
    B <- (1/r$minimum-1)*S
    if (B==maxbias)
        warning("Optimum reacher maximal bias/sd ratio")
    list(w=r$minimum, length=r$objective, B=B)
}

#' Optimal shrinkage for Empirical Bayes estimator
#' @inheritParams wopt
#' @export
wl2 <- function(S, kappa=Inf, alpha=0.05) {
    w <- S^2/(1+S^2)
    c(w=w, length=cva(B, kurt=kappa, alpha=alpha, check=TRUE)$cv*w, B=1/S)
}
