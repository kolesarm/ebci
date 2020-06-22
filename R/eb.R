#' Optimal shrinkage for Empirical Bayes confidence intervals
#'
#' Compute linear shrinkage factor \eqn{w_{opt}}{w_opt} to minimize the length
#' of the resulting Empirical Bayes confidence interval (EBCI).
#' @param S Signal-to-noise ratio \eqn{\mu_{2}/\sigma^2}{mu_2/sigma^2}, where
#'     \eqn{\mu_{2}}{mu_2} is the second moment of
#'     \eqn{\theta-X'\delta}{theta-X*delta}, and \eqn{\sigma^2}{sigma^2} is the
#'     variance of the preliminary estimator.
#' @param kappa Kurtosis of \eqn{\theta-X'\delta}{theta-X*delta}.
#' @param alpha Determines confidence level, \eqn{1-\alpha}{1-alpha}.
#' @param cv_tbl Optionally, supply a data frame of critical values.
#'     \code{cva(m2, kappa, alpha)}, for different values of \code{m2}, such as
#'     a subset of the data frame \code{\link{cva_tbl}}, that matches the
#'     supplied values of \code{alpha} and \code{kappa}. The data frame needs to
#'     contain two variables, \code{m2}, corresponding to the value of second
#'     moment of the normalized bias, and \code{cv}, with the corresponding
#'     value of \code{cva(m2, kappa, alpha)}. If non \code{NULL}, for the
#'     purposes of optimizing the shrinkage factor, compute the critical value
#'     \code{cva} by interpolating between the critical values in this data
#'     frame, instead of computing them from scratch. This can speed up the
#'     calculations.
#' @return Returns a list with 3 components: \describe{
#'
#' \item{\code{w}}{Optimal shrinkage factor \eqn{w_{opt}}{w_opt}}
#'
#' \item{\code{length}}{Normalized half-length of the corresponding confidence
#' interval, so that the interval obtains by taking the estimator based on
#' shrinkage given by \code{w}, and adding and subtracting \code{length} times
#' the standard error \eqn{\sigma}{sigma} of the preliminary estimator.}
#'
#' \item{\code{m2}}{The second moment of the normalized bias,
#' \eqn{(1/w-1)^{2}S}{(1/w-1)^2*S}}}
#' @references{
#'
#' \cite{Armstrong, Timothy B., Kolesár, Michal, and Plagborg-Møller, Mikkel
#' (2020): Robust Empirical Bayes Confidence Intervals,
#' \url{https://arxiv.org/abs/2004.03448}}
#'
#' }
#' @examples
#' w_opt(4, 3)
#' ## Use precomputed critical value table
#' w_opt(4, 3, cv_tbl=cva_tbl[cva_tbl$kappa==3 & cva_tbl$alpha==0.05, ])
#' @export
w_opt <- function(S, kappa, alpha=0.05, cv_tbl=NULL) {
    if (!is.null(cv_tbl)) {
        ## Linearly interpolate
        cv <- function(m2) {
            idx <- which.max(cv_tbl$m2>m2)
            idx0 <- max(idx-1, 1)
            cv_tbl$cv[idx0]+(m2-cv_tbl$m2[idx0])/
                (cv_tbl$m2[idx]-cv_tbl$m2[idx0]) *
                (cv_tbl$cv[idx]-cv_tbl$cv[idx0])
        }
        maxbias <- max(cv_tbl$m2)
        minbias <- min(cv_tbl$m2)
    } else {
        cv <- function(m2) cva(m2, kappa=kappa, alpha=alpha, check=FALSE)$cv
        ## Maxium bias (i.e m2) is 100^2 to prevent numerical issues
        maxbias <- 100^2
        minbias <- 0
    }
    ci_length <- function(w) cv((1/w-1)^2*S)*w
    if (S > tol) {
        r <- stats::optimize(ci_length, c(1/(sqrt(maxbias/S)+1),
                                          1/(sqrt(minbias/S)+1)), tol=tol)
        m2 <- (1/r$minimum-1)^2*S
    } else {
        r <- list(minimum=0, objective=NA)
        m2 <- Inf
    }
    ## Recompute critical value, checking solution accuracy. If we're using
    ## cv_tbl, then optimum will be at one of the values of m2 in the table, so
    ## solution should always be accurate
    len <- cva(m2, kappa=kappa, alpha=alpha, check=TRUE)$cv*r$minimum
    if (m2 > maxbias-1e-4)
        warning("Corner solution: optimum reached at maximal normalized bias")
    list(w=r$minimum, length=len, m2=m2)
}

#' Empirical Bayes estimator and confidence intervals
#'
#' Compute empirical Bayes shrinkage factor \eqn{w_{eb}}{w_eb} and the
#' normalized half-length of of the associated robust Empirical Bayes confidence
#' interval (EBCI).
#' @inheritParams w_opt
#' @return Returns a list with 3 components:
#' \describe{
#'
#' \item{\code{w}}{Empirical Bayes shrinkage factor \eqn{w_{eb}}{w_eb}}
#'
#' \item{\code{length}}{Normalized half-length of the corresponding confidence
#' interval, so that the interval obtains by taking the estimator based on
#' shrinkage given by \code{w}, and adding and subtracting \code{length} times
#' the standard error \eqn{\sigma}{sigma} of the preliminary estimator.}
#'
#' \item{\code{m2}}{Second moment of the normalized bias, \code{1/m2}.}
#'
#' }
#' @references{
#'
#' \cite{Armstrong, Timothy B., Kolesár, Michal, and Plagborg-Møller, Mikkel
#' (2020): Robust Empirical Bayes Confidence Intervals,
#' \url{https://arxiv.org/abs/2004.03448}}
#'
#' }
#' @examples
#' w_eb(2, 3)
#' ## No constraint on kurtosis yields doesn't affect shrinkage, but yields
#' ## larger half-length
#' w_eb(2, Inf)
#' @export
w_eb <- function(S, kappa=Inf, alpha=0.05) {
    w <- S/(1+S)
    list(w=w, length=cva(1/S, kappa=kappa, alpha=alpha, check=TRUE)$cv*w,
         m2=1/S)
}


#' Compute empirical Bayes confidence intervals by shrinking toward regression
#'
#' Computes empirical Bayes estimators based on shrinking towards a regression,
#' and associated robust empirical Bayes confidence intervals (EBCIs), as well
#' as length-optimal robust EBCIs.
#' @param formula object of class \code{"formula"} (or one that can be coerced
#'     to that class) of the form \code{Y ~ predictors}, where \code{Y} is a
#'     preliminary unbiased estimator, and \code{predictors} are predictors
#'     \eqn{X} that guide the direction of shrinkage. For shrinking toward the
#'     grand mean, use \code{Y ~ 1}, and for shrinking toward \code{0} use
#'     \code{Y ~ 0}
#' @param data optional data frame, list or environment (or object coercible by
#'     \code{as.data.frame} to a data frame) containing the preliminary
#'     estimator \code{Y} and the predictors. If not found in \code{data}, these
#'     variables are taken from \code{environment(formula)}, typically the
#'     environment from which the function is called.
#' @param se Standard errors \eqn{\sigma}{sigma} associated with the preliminary
#'     estimates \code{Y}
#' @param weights An optional vector of weights to be used in the fitting
#'     process in computing \eqn{\delta}{delta}, \eqn{\mu_2}{mu_2} and
#'     \eqn{\kappa}{kappa}. Should be \code{NULL} or a numeric vector.
#' @param alpha Determines confidence level, \eqn{1-\alpha}{1-alpha}.
#' @param kappa If non-\code{NULL}, use pre-specified value for the kurtosis
#'     \eqn{\kappa}{kappe} of \eqn{\theta-X'\delta}{theta-X*delta} (such as
#'     \code{Inf}), instead of computing it.
#' @param fs_correction Finite-sample correction method used to compute
#'     \eqn{\mu_2}{mu_2} and \eqn{\kappa}{kappa}. These corrections ensure that
#'     we do not shrink the preliminary estimates \code{Y} all the way to zero.
#'     If \code{"PMT"}, use posterior mean truncation, if \code{"FPLIB"} use
#'     limited information Bayesian approach with a flat prior, and if
#'     \code{"none"}, truncate the estimates at \code{0} for \eqn{\mu_2}{mu_2}
#'     and \code{1} for \eqn{\kappa}{kappa}.
#' @param tstat If \code{TRUE}, shrink the t-statistics \code{Y/se} rather than
#'     the preliminary estimates \code{Y}.
#' @param wopt If \code{TRUE}, also compute length-optimal robust EBCIs. These
#'     are robust EBCIs centered at length-optimal shrinkage factors
#'     \code{\link{w_opt}}.
#' @param cores Number of cores to use when computing length-optimal robust
#'     EBCIs. By default, the computation of the length-optimal shrinkage
#'     factors \code{\link{w_opt}} is parallelized if there are more than 30
#'     observations to speed up the calculations.
#' @return Returns a list with the following components: \describe{
#'
#' \item{\code{sqrt_mu2}}{Square root of the estimated second moment of
#' \eqn{\theta-X'\delta}{theta-X*delta}, \eqn{\sqrt{\mu_2}}{mu_2^(1/2)}. Vector
#' of length 2, the first element corresponds to the estimate after the
#' finite-sample correction as specified by \code{fs_correction}, the second
#' element is the uncorrected estimate.}
#'
#' \item{\code{kappa}}{Estimated kurtosis \eqn{\kappa}{kappa} of
#' \eqn{\theta-X'\delta}{theta-X*delta}. Vector of length 2 with the same
#' structure as \code{sqrt_mu2}.}
#'
#' \item{\code{delta}}{Estimated regression coefficients \eqn{\delta}{delta}}
#'
#' \item{\code{df}}{Data frame with components described below.}
#' }
#'
#' \code{df} has the following components:
#'
#' \describe{
#'
#' \item{\code{w_eb}}{EB shrinkage factors,
#'    \eqn{\mu_{2}/(\mu_{2}+\sigma^2_i)}{mu_2/(mu_2+sigma^2_i)}}
#'
#' \item{\code{w_opt}}{Optimal shrinkage factors \code{\link{w_opt}}}
#'
#' \item{\code{ncov_pa}}{Maximal non-coverage of parametric EBCIs}
#'
#' \item{\code{len_eb}}{Half-length of robust EBCIs based on EB shrinkage, so
#' that the intervals take the form \code{cbind(th_eb-len_eb, th_eb+len_eb)}}
#'
#' \item{\code{len_op}}{Half-length of robust EBCIs based on length-optimal
#' shrinkage, so that the intervals take the form \code{cbind(th_op-len_op,
#' th_op+len_op)}}
#'
#' \item{\code{len_pa}}{Half-length of parametric EBCIs, which take the form
#' \code{cbind(th_eb-len_pa, th_eb+len_a)}}
#'
#' \item{\code{len_us}}{Half-length of unshrunk CIs, which take the form
#' \code{cbind(th_us-len_us, th_us+len_us)}}
#'
#' \item{\code{th_us}}{Unshrunk estimate \eqn{Y}}
#'
#' \item{\code{th_eb}}{EB estimate.}
#'
#' \item{\code{th_eb}}{Estimate based on length-optimal shrinkage.}
#'
#' \item{\code{se}}{Standard error \eqn{\sigma}{sigma}, as supplied by the
#' argument \code{se}.}
#' }
#' @references{
#'
#' \cite{Armstrong, Timothy B., Kolesár, Michal, and Plagborg-Møller, Mikkel
#' (2020): Robust Empirical Bayes Confidence Intervals,
#' \url{https://arxiv.org/abs/2004.03448}}
#'
#' }
#' @examples
#' ## Using only commuting zones in California
#' ebci(theta25 ~ stayer25, cz[cz$state=="NY", ],
#'      se25, 1/se25^2, tstat=FALSE, cores=1)
#' @export
ebci <- function(formula, data, se, weights=NULL, alpha=0.1, kappa=NULL,
                 wopt=FALSE, fs_correction="PMT", tstat=FALSE,
                 cores=max(parallel::detectCores()-1L, 1L)) {
    ## Construct model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "se", "weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$formula <- stats::as.formula(formula)
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    Y <- stats::model.response(mf, "numeric")
    wgt <- as.vector(stats::model.weights(mf))
    ## Do not weight
    if (is.null(wgt))
        wgt <- rep(1L, length(Y))
    se <- mf$"(se)"
    X <- stats::model.matrix(attr(mf, "terms"), mf)
    za <- stats::qnorm(1-alpha/2)

    ## Compute delta, and moments
    Yt <- if(tstat) Y/se else Y
    set <- if (tstat) 1L else se
    r1 <- stats::lm.wfit(X, Yt, wgt)
    mu1 <- unname(r1$fitted.values)
    delta <- r1$coefficients
    w2 <- (Yt-mu1)^2 - set^2
    w4 <- (Yt-mu1)^4 - 6*set^2*(Yt-mu1)^2 + 3*set^4

    ## Weighted sample variance and mean of truncated normal
    wgtV <- function(Z, wgt)
        sum(wgt^2*(Z^2-stats::weighted.mean(Z, wgt)^2))/(sum(wgt)^2-sum(wgt^2))
    tmean <- function(m, V)  m + sqrt(V)*stats::dnorm(m/sqrt(V))/
                                 stats::pnorm(m/sqrt(V))
    ## \tilde{mu}_2, \tilde{mu}_4)
    tmu <- c(stats::weighted.mean(w2, wgt), stats::weighted.mean(w4, wgt))
    pmt_trim <- c(2 *mean(wgt^2*set^4) / (sum(wgt)*mean(wgt*set^2)),
                  32*mean(wgt^2*set^8) / (sum(wgt)*mean(wgt*set^4)))
    mu2 <- switch(fs_correction,
                  "none"=max(tmu[1], 0),
                  "PMT"=max(tmu[1], pmt_trim[1]),
                  "FPLIB"=tmean(tmu[1], wgtV(w2, wgt)),
                  stop("Unknown fs_correction method"))
    if (mu2 == 0) {
        warning("mu2 estimate is 0")
        df <- data.frame(w_eb=se*0, w_opt=se*0, ncov_pa=se*NA, len_eb=se*NA,
                         len_op=se*NA, len_pa=se*NA, len_us=za*se, th_us=Y,
                         th_eb=Yt-mu1, th_op=Yt-mu1, se=se)
        return(list(sqrt_mu2=sqrt(mu2), kappa=kappa, delta=delta, df=df))
    }
    if (is.null(kappa)) {
        kappa <- switch(fs_correction,
                  "none"=max(tmu[2]/mu2^2, 1),
                  "PMT"=max(tmu[2]/mu2^2, 1 + pmt_trim[2]/mu2^2),
                  "FPLIB"= 1+tmean(tmu[2]-tmu[1]^2,
                                 wgtV(w4-2*mu2*w2, wgt))/mu2^2)
    }

    if (tstat) {
        eb <- w_eb(mu2, kappa, alpha)
        th_eb <- se*(mu1+eb$w*(Yt-mu1))
        op <- w_opt(mu2, kappa, alpha)
        th_op <- se*(mu1+op$w*(Yt-mu1))
        ncov_pa <- rho(1/eb$w-1, kappa, za/sqrt(eb$w), check=TRUE)$alpha
    } else {
        ## EB shrinkage
        eb <- vapply(se, function(se) unlist(w_eb(mu2/se^2, kappa, alpha)),
                     numeric(3))
        eb <- as.data.frame(t(eb))
        th_eb <- mu1+eb$w*(Yt-mu1)

        ## Optimal shrinkage. First pre-compute cva for this kurtosis
        if (length(se) > 30 & wopt) {
            mmin <- w_opt(mu2/min(se)^2, kappa, alpha=alpha)$m2
            mmax <- w_opt(mu2/max(se)^2, kappa, alpha=alpha)$m2
            df <- data.frame(m2=seq(mmin, mmax, length.out=501)^2)
            cvj <- function(j)
                ebci::cva(df$m2[j], kappa=kappa, alpha, check=TRUE)$cv
            df$cv <- if (cores>1)
                         unlist(parallel::mclapply(seq_len(nrow(df)), cvj,
                                                   mc.cores=cores))
                     else
                         vapply(seq_len(nrow(df)), cvj, numeric(1))
        } else {
            df <- NULL
        }
        if (wopt) {
            op <- vapply(se, function(se)
                unlist(w_opt(mu2/se^2, kappa, cv_tbl=df,
                             alpha=alpha)), numeric(3))
            op <- as.data.frame(t(op))
            th_op <- mu1+op$w*(Yt-mu1)
        } else {
            op <- eb*NA
            th_op <- th_eb*NA
        }

        par_cov <- function(se) rho(se^2/mu2, kappa,
                                    za*sqrt(1+se^2/mu2), check=TRUE)$alpha
        ncov_pa <- vapply(se, par_cov, numeric(1))
    }
    df <- data.frame(w_eb=eb$w, w_opt=op$w, ncov_pa=ncov_pa,
                     len_eb=eb$length*se, len_op=op$length*se,
                     len_pa=sqrt(eb$w)*za*se, len_us=za*se, th_us=Y,
                     th_eb=th_eb, th_op=th_op, se=se)
    list(mu2=c("estimate"=mu2, "uncorrected_estimate"=tmu[1]),
         kappa=c("estimate"=kappa, "uncorrected_estimate"=tmu[2]/tmu[1]^2),
         delta=delta, df=df)
}
