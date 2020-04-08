#' Optimal shrinkage for Empirical Bayes confidence intervals
#'
#' Compute linear shrinkage factor \eqn{w_{opt}}{w_opt} to minimize the length
#' of the resulting Empirical Bayes confidence interval (EBCI).
#' @param S Square root of the signal-to-noise ratio
#'     \eqn{\sqrt{\mu_{2}}/\sigma}{sqrt(mu_2)/sigma}, where
#'     \eqn{\sqrt{\mu_{2}}}{mu_2} is the variance of
#'     \eqn{\theta-X'\delta}{theta-X*delta}, and \eqn{\sigma^2}{sigma^2} is the
#'     variance of the preliminary estimator.
#' @param kappa Kurtosis of \eqn{\theta-X'\delta}{theta-X*delta}.
#' @param alpha Determines confidence level, \eqn{1-\alpha}{1-alpha}.
#' @param cv_tbl Optionally, supply a data frame of critical values.
#'     \code{cva(B, kappa, alpha)}, for different values of \code{B}, such as a
#'     subset of the data frame \code{\link{cva_tbl}}, that matches the supplied
#'     values of \code{alpha} and \code{kappa}. The data frame needs to contain
#'     two variables, \code{B}, corresponding to the value of average squared
#'     normalized bias, and \code{cv}, with the corresponding value of
#'     \code{cva(B, kappa, alpha)}. If non \code{NULL}, for the purposes of
#'     optimizing the shrinkage factor, compute the critical value \code{cva} by
#'     interpolating between the critical values in this data frame, instead of
#'     computing them from scratch. This can speed up the calculations.
#' @return Returns a list with 3 components:
#' \describe{
#'
#' \item{\code{w}}{Optimal shrinkage factor \eqn{w_{opt}}{w_opt}}
#'
#' \item{\code{length}}{Normalized half-length of the corresponding confidence
#' interval, so that the interval obtains by taking the estimator based on
#' shrinkage given by \code{w}, and adding and subtracting \code{length} times
#' the standard error \eqn{\sigma}{sigma} of the preliminary estimator.}
#'
#' \item{\code{B}}{Square root of the normalized bias,
#' \eqn{(1/w-1)S}{(1/w-1)*S}}}
#' @references{
#'
#' \cite{Armstrong, Timothy B., Kolesár, Michal, and Plagborg-Møller, Mikkel
#' (2020): Robust Empirical Bayes Confidence Intervals,
#' \url{https://arxiv.org/abs/2004.03448}}
#'
#' }
#' @examples
#' w_opt(1, 3)
#' ## Use precomputed critical value table
#' w_opt(1, 3, cv_tbl=cva_tbl[cva_tbl$kappa==3 & cva_tbl$alpha==0.05, ])
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
        cv <- function(B) cva(B, kappa=kappa, alpha=alpha, check=FALSE)$cv
        ## Maxium bias (i.e B) is 100 to prevent numerical issues
        maxbias <- 100
    }
    ci_length <- function(w) cv((1/w-1)*S)*w
    r <- stats::optimize(ci_length, c(1/(maxbias/S+1), 1), tol=tol)
    B <- (1/r$minimum-1)*S
    ## Recompute critical value, checking solution accuracy. If we're using
    ## cv_tbl, then optimum will be at one of the values of B in the table, so
    ## solution should always be accurate
    len <- cva((1/r$minimum-1)*S, kappa=kappa, alpha=alpha,
                  check=TRUE)$cv*r$minimum
    if (B > maxbias-1e-4)
        warning("Corner solution: optimum reached  maximal bias/sd ratio")
    list(w=r$minimum, length=len, B=B)
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
#' \item{\code{B}}{Square root of the normalized bias, \eqn{1/S}.}
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
#' w_opt(1, 3)
#' ## No constraint on kurtosis yields doesn't affect shrinkage, but yields
#' ## larger half-length
#' w_opt(1, Inf)
#' @export
w_eb <- function(S, kappa=Inf, alpha=0.05) {
    w <- S^2/(1+S^2)
    list(w=w, length=cva(1/S, kappa=kappa, alpha=alpha, check=TRUE)$cv*w, B=1/S)
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
#' @param tstat If \code{TRUE}, shrink the t-statistics \code{Y/se} rather than
#'     the preliminary estimates \code{Y}.
#' @param cores Number of cores to use. By default, the computation of the
#'     length-optimal shrinkage factors \code{\link{w_opt}} is parallelized to
#'     speed up the calculations.
#' @return Returns a list with the following components: \describe{
#'
#' \item{\code{sqrt_mu2}}{Square root of the estimated second moment of
#' \eqn{\theta-X'\delta}{theta-X*delta}, \eqn{\sqrt{\mu_2}}{mu_2^(1/2)}}
#'
#' \item{\code{kappa}}{Estimated kurtosis \eqn{\kappa}{kappa} of
#' \eqn{\theta-X'\delta}{theta-X*delta}}
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
#' \item{\code{len_eb}}{Half-length of EBCIs based on EB shrinkage, so that the
#' intervals take the form
#' \code{cbind(th_eb-len_eb, th_eb+len_eb)}}
#'
#' \item{\code{len_op}}{Half-length of EBCIs based on length-optimal shrinkage,
#' so that the intervals take the form \code{cbind(th_op-len_op, th_op+len_op)}}
#'
#' \item{\code{len_pa}}{Half-length of parametric EBCIs, which take the form
#' \code{cbind(th_eb-len_pa, th_eb+len_a)}}
#'
#' \item{\code{len_us}}{Half-length of unshrunk CIs that take the form
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
#' ebci(theta25 ~ stayer25, cz, se25, pop/pop, tstat=TRUE)
#' ebci(theta25 ~ 0, cz, se25, pop/pop, tstat=TRUE)
#' @export
ebci <- function(formula, data, se, weights, alpha=0.1, kappa=NULL,
                     tstat=FALSE, cores=max(parallel::detectCores()-1L, 1L)) {
    ## Construct model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "se", "weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    formula <- stats::as.formula(formula)
    mf$formula <- formula
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    Y <- stats::model.response(mf, "numeric")
    wgt <- as.vector(stats::model.weights(mf))
    ## Do not weight
    if (is.null(wgt))
        wgt <- rep(1L, length(Y))
    se <- mf$"(se)"
    X <- stats::model.matrix(stats::terms(formula, data = data), mf)

    ## Compute delta, and moments
    Yt <- if(tstat) Y/se else Y
    set <- if (tstat) 1L else se
    r1 <- stats::lm.wfit(X, Yt, wgt)
    mu1 <- unname(r1$fitted.values)
    ## delta <- r1$coefficients[!is.na(r1$coefficients)]
    delta <- r1$coefficients
    mu2 <- stats::weighted.mean((Yt-mu1)^2-set^2, wgt)
    if (is.null(kappa)) {
        ## Formula without assuming epsilon_i and sigma_i are uncorrelated,
        ## alternative is weighted.mean((Yt-mu1)^4 - 6*set^2*mu2 - 3*set^4, wgt)
        kappa <- stats::weighted.mean((Yt-mu1)^4 - 6*set^2*(Yt-mu1)^2 +
                                      3*set^4, wgt) / mu2^2
        if (kappa < 1) {
            message("Kurtosis estimate is smaller than 1, setting it to Inf")
            kappa <- Inf
        }
    }
    za <- stats::qnorm(1-alpha/2)
    if (tstat) {
        op <- w_opt(sqrt(mu2), kappa, alpha)
        eb <- w_eb(sqrt(mu2), kappa, alpha)

        th_eb <- se*(mu1+eb$w*(Yt-mu1))
        th_op <- se*(mu1+op$w*(Yt-mu1))
        ncov_pa <- rho(1/eb$w-1, kappa, za/sqrt(eb$w),
                       check=TRUE)$size
    } else {
        ## Pre-compute cva for this kurtosis
        df <- data.frame(B=seq(0, 20, length.out=2001), kappa=kappa)
        cvj <- function(j)
            ebci::cva(df$B[j], kappa=kappa, alpha, check=TRUE)$cv
        df$cv <- if (cores>1)
                     unlist(parallel::mclapply(seq_len(nrow(df)), cvj,
                                               mc.cores=cores))
                 else
                     vapply(seq_len(nrow(df)), cvj, numeric(1))
        df$alpha <- alpha

        ## Optimal shrinkage for each individual
        op <- vapply(se, function(se)
            unlist(w_opt(sqrt(mu2)/se, kappa, cv_tbl=df,
                              alpha=alpha)), numeric(3))
        op <- as.data.frame(t(op))
        ## EB shrinkage
        eb <- vapply(se, function(se)
            unlist(w_eb(sqrt(mu2)/se, kappa, alpha)), numeric(3))
        eb <- as.data.frame(t(eb))

        th_eb <- mu1+eb$w*(Yt-mu1)
        th_op <- mu1+op$w*(Yt-mu1)
        par_cov <- function(se)
            rho(se^2/mu2, kappa,
                 za*sqrt(1+se^2/mu2), check=TRUE)$size
        ncov_pa <- vapply(se, par_cov, numeric(1))
    }
    df <- data.frame(w_eb=eb$w,
                     w_opt=op$w,
                     ncov_pa=ncov_pa,
                     len_eb=eb$length*se,
                     len_op=op$length*se,
                     len_pa=sqrt(eb$w)*za*se,
                     len_us=za*se,
                     th_us=Y,
                     th_eb=th_eb,
                     th_op=th_op,
                     se=se)
    list(sqrt_mu2=sqrt(mu2), kappa=kappa, delta=delta, df=df)
}
