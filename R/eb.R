#' Optimal shrinkage for Empirical Bayes confidence intervals
#'
#' Compute linear shrinkage factor \eqn{w_{opt}}{w_opt} to minimize the length
#' of the resulting Empirical Bayes confidence interval (EBCI).
#' @param S Square root of the signal-to-noise ratio
#'     \eqn{\sqrt{\mu_{2}}/\sigma}{sqrt(mu_2)/sigma}, where
#'     \eqn{\sqrt{\mu_{2}}}{mu_2} is the variance of \eqn{\theta}{theta}
#' @param kappa Kurtosis of \eqn{\theta}{theta}
#' @param alpha determines CI level
#' @param cv_tbl Optionally, supply a data frame of critical values. (TODO)
#' @return (TODO)
#' @examples
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
#' Compute linear shrinkage factor \eqn{w_{eb}}{w_eb} and normalize length of of
#' the associated robust Empirical Bayes confidence interval (EBCI).
#' @inheritParams w_opt
#' @return (TODO)
#' @export
w_eb <- function(S, kappa=Inf, alpha=0.05) {
    w <- S^2/(1+S^2)
    list(w=w, length=cva(1/S, kappa=kappa, alpha=alpha, check=TRUE)$cv*w, B=1/S)
}


weighted.var <- function(y, w, na.rm=FALSE) {
    if (na.rm) {
        idx <- !is.na(y) & !is.na(w)
        y <- y[idx]
        w <- w[idx]
    }
    length(y)/(length(y)-1)*sum(w*(y-weighted.mean(y, w))^2)/sum(w)
}

#' EBCIs in applications
#' @param kappa Use pre-specified value for kurtosis (such as Inf). If NULL,
#'     then compute it
ebci <- function(formula, data, se, weights, alpha=0.1, kappa=NULL,
                     tstat=FALSE) {
    ## Construct model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "se", "weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    formula <- as.formula(formula)
    mf$formula <- formula
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    Y <- stats::model.response(mf, "numeric")
    wgt <- as.vector(stats::model.weights(mf))
    se <- mf$"(se)"
    X <- model.matrix(stats::terms(formula, data = data), mf)

    ## Compute delta, and moments
    Yt <- if(tstat) Y/se else Y
    set <- if (tstat) 1L else se
    r1 <- lm.wfit(X, Yt, wgt)
    mu1 <- unname(r1$fitted.values)
    ## delta <- r1$coefficients[!is.na(r1$coefficients)]
    delta <- r1$coefficients
    mu2 <- weighted.mean((Yt-mu1)^2-set^2, wgt)
    if (is.null(kappa)) {
        ## Formula without assuming epsilon_i and sigma_i are uncorrelated,
        ## alternative is weighted.mean((Yt-mu1)^4 - 6*set^2*mu2 - 3*set^4, wgt)
        kappa <- weighted.mean((Yt-mu1)^4 - 6*set^2*(Yt-mu1)^2 +
                               3*set^4, wgt) / mu2^2
        if (kappa < 1) {
            message("Kurtosis estimate is smaller than 1, setting it to Inf")
            kappa <- Inf
        }
    }
    if (tstat) {
        op <- w_opt(sqrt(mu2), kappa, alpha)
        eb <- w_eb(sqrt(mu2), kappa, alpha)

        th_eb <- se*(mu1+eb$w*(Yt-mu1))
        th_op <- se*(mu1+op$w*(Yt-mu1))
        ncov_pa <- rho2(1/eb$w-1, kappa, qnorm(1-alpha/2)/sqrt(eb$w),
                       check=TRUE)$size
    } else {
        ## Pre-compute cva for this kurtosis
        df <- data.frame(B=seq(0, 20, length.out=2001), kappa=kappa)
        cvj <- function(j)
            ebci::cva(df$B[j], kappa=kappa, alpha, check=TRUE)$cv
        cores <- max(parallel::detectCores()-1L, 1L)
        df$cv <- unlist(parallel::mclapply(seq_len(nrow(df)), cvj, mc.cores=cores))
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
            rho2(se^2/mu2, kappa,
                 qnorm(1-alpha/2)*sqrt(1+se^2/mu2), check=TRUE)$size
        ncov_pa <- vapply(se, par_cov, numeric(1))
    }
    df <- data.frame(w_eb=eb$w,
                     w_opt=op$w,
                     ncov_pa=ncov_pa,
                     len_eb=eb$length*se,
                     len_op=op$length*se,
                     len_pa=sqrt(eb$w)*qnorm(1-alpha/2)*se,
                     len_us=qnorm(1-alpha/2)*se,
                     th_us=Y,
                     th_eb=th_eb,
                     th_op=th_op,
                     se=se)

    list(sqrt_mu2=sqrt(mu2), kappa=kappa, delta=delta, df=df)
}
