tol <- 1e-12

## Function called r_0 in the paper
r <- function(t, chi) {
    idx <- sqrt(t)-chi <= 5
    r <- rep(1, length(t))
    ## Use normal rather than stats::pchisq(..., lower.tail=FALSE), since
    ## non-central chi^2 is not accurate for extreme values. It's also much
    ## faster.
    r[idx] <- stats::pnorm(-sqrt(t[idx])-chi)+stats::pnorm(sqrt(t[idx])-chi)
    r
}

## Derivative of r
r1 <- function(t, chi) {
    ## Apply L'Hospital's rule. This gives maximum absolute error 1e-9
    ifelse(t >= 1e-8,
    (stats::dnorm(sqrt(t)-chi)-stats::dnorm(sqrt(t)+chi)) / (2*sqrt(t)),
    chi*stats::dnorm(chi))
}

## Second Derivative of r
r2 <- function(t, chi) {
    ## Apply L'Hospital's rule 3x. This gives maximum abs error about 8e-8
    ifelse(t >= 2e-6,
           (stats::dnorm(sqrt(t)+chi)*(chi*sqrt(t)+t+1)+
            stats::dnorm(sqrt(t)-chi)*(chi*sqrt(t)-t-1)) / (4*t^(3/2)),
           stats::dnorm(chi)*chi*(chi^2-3)/6)
}

## Third Derivative of r
r3 <- function(t, chi) {
    ## Apply L'Hospital's rule 5x. This gives maximum abs error about 3e-6
    ifelse(t >= 2e-4,
    (stats::dnorm(chi-sqrt(t))*(t^2-2*chi*t^(3/2)+(2+chi^2)*t-3*chi*sqrt(t)+3)
    -stats::dnorm(chi+sqrt(t))*(t^2+2*chi*t^(3/2)+(2+chi^2)*t+3*chi*sqrt(t)+3))
    /(8*t^(5/2)),
    stats::dnorm(chi)*(chi^5-10*chi^3+15*chi)/60)
}

## find t0 and inflection point, called t1 in the paper.
rt0 <- function(chi) {
    ## Find point where we touch origin
    f0 <- function(t) r(t, chi)-t*r1(t, chi)-r(0, chi)
    if (chi^2 < 3) {
        t0 <- ip <- 0
    } else {
        ## Avoid issues when chi is numerically very large
        if (abs(r2(chi^2-3/2, chi=chi)) < tol | (chi^2-3)-chi^2==0L)
            ip <- chi^2-3/2
        else
            ip <- stats::uniroot(r2, c(chi^2-3, chi^2), tol=tol, chi=chi)$root

        ## Make sure upper endpoint of interval is positive; it always is for
        ## chi< 100,000, so we should never enter the while loop
        up <- 2*chi^2
        lo <- ip
        while (f0(up) < 0) {
            lo <- up
            up <- 2*up
        }
        t0 <- lo
        if (f0(lo) < 0) {
            t0 <- stats::uniroot(f0, c(lo, up), tol=tol)$root
        } else if (f0(lo) > tol) {
            warning("Failed to solve for t0 using rt0 at chi=", chi)
        }
    }
    list(t0=t0, ip=ip)
}

## \rho(m_{2}, \chi)
rho0 <- function(t, chi) {
    t0 <- rt0(chi)$t0
    idx <- (t >= t0)
    res <- r(t0, chi) + (t-t0)*r1(t0, chi)
    if (any(idx))
        res[idx] <- r(t[idx], chi)
    res
}


## \delta(x; x_{0})
delta <- function(x, x0, chi) {
    ## Apply L'Hospital's rule 2x. This gives maximum abs error about 1e-7
    ifelse(abs(x-x0)<1e-4,
           r2(x0, chi)/2,
           (r(x, chi) - r(x0, chi) - r1(x0, chi)*(x-x0)) / (x-x0)^2)
}

delta1 <- function(x, x0, chi) {
    ## Apply L'Hospital's rule 2x. This gives maximum abs error about 7e-6
    ifelse(abs(x-x0)<1e-3,
           r3(x0, chi)/6,
           ((r1(x, chi)+r1(x0, chi))-2*(r(x, chi)-r(x0, chi))/(x-x0))/(x-x0)^2)
}

## maximize delta(x, x0, chi) over x.
lam <- function(x0, chi) {
    ## Check derivatives at 0, inflection point, t0, and x0 If we're above
    ## inflection point, then maximum is below it, and it's at zero if
    ## derivative at zero is negative.
    xs <- sort(unlist(rt0(chi)))
    if (x0 >= xs[1]) {
        xs <- c(0, xs[1])
    } else {
        xs <- c(0, x0, xs)
    }
    ## We want >= since derivative at zero may be numerically close to zero if
    ## chi is large
    der <- delta1(xs, x0, chi) >= 0

    if (all(der<=0) | max(abs(delta1(xs, x0, chi))) < tol^2)
        return(list(lam=delta(0, x0, chi), x0=0))
    else if (all(diff(der)<=0)) {
        ## Function first increasing, then decreasing
        start <- xs[which.min(der)-1]
        end <- xs[which.min(der)]
    } else {
        stop(paste0("There are multiple optima in the function delta(x, x0=",
                    x0, ", chi=", chi, ")"))
    }

    rr1 <- stats::optimize(function(x) -delta(x, x0, chi),
                           c(start, end), tol=tol)

    list(lam=-rr1$objective, x0=unname(rr1$minimum))
}

## \rho(m_{2}, mu_{4}, \chi)
rho <- function(m2, kappa, chi, check=TRUE, len=5000) {
    r0 <- rho0(m2, chi)
    t0 <- rt0(chi)$t0

    if (kappa == 1) {
        list(alpha=r(m2, chi), x=c(0, m2), p=c(0, 1))
    } else if (m2 >= t0) {
        ## Concave part of parameter space
        list(alpha=r0, x=c(0, m2), p=c(0, 1))
    } else if (kappa==Inf || m2*kappa >= t0) {
        ## LF under rho: (0, t0) wp (1-m2/t0, m2/t0), E[t^2]=m2*t0
        ## So here kappa doesn't bind
        list(alpha=r0, x=c(0, t0), p=c(1-m2/t0, m2/t0))
    } else {
        ## First determine where delta(x, 0) is maximized
        tbar <- lam(0, chi)$x0
        lammax <- function(x0) {
            ## delta(x, x0) maximized at 0
            ifelse(x0 >= tbar, delta(0, x0, chi), max(lam(x0, chi)$lam, 0))
        }
        obj <-  function(x0)
            r(x0, chi)+r1(x0, chi)*(m2-x0)+
                lammax(x0)*(kappa*m2^2-2*x0*m2+x0^2)
        ## Optimize separately below and above \bar{t}, since there are
        ## typically be multiple local minima
        rb <- if (tbar >0)
            stats::optimize(obj, interval=c(0, tbar), tol=10^3*tol)
              else
            list(minimum=0, objective=obj(0))
        ra <- stats::optimize(obj, interval=c(tbar, t0), tol=10^3*tol)
        rr <- if (rb$objective<ra$objective) rb else ra
        ## LF points
        xs0 <- sort(c(rr$minimum, lam(rr$minimum, chi)$x0))
        p <- (m2-xs0[2])/(xs0[1]-xs0[2])
        ps0 <- c(p, 1-p)
        ## Rejection rates, m2, and kappa at LF solution
        primal <- c(sum(c(r(xs0[1], chi), r(xs0[2], chi))*ps0),
                   sum(xs0*ps0),
                   sum(xs0^2*ps0)/sum(xs0*ps0)^2)
        ## If LF solution close to dual, no need to check linear program
        if (max(abs(primal-c(rr$objective, m2, kappa))) > 1e-4 && check) {
            ## Add m2 here for cases where it's very small, so we can satisfy
            ## the constraint
            xs <- sort(unique(c(m2, xs0, seq(0, t0, length.out=len))))
            ## Find the optimal solution; use <= for last constraint
            opt <-  lpSolve::lp(direction="max",
                                objective.in = r(xs, chi),
                                const.mat = rbind(xs^0, xs, xs^2),
                                const.dir = c("==", "==", "<="),
                                const.rhs = c(1, m2, m2^2*kappa))
            ## 0: success, 2: infeasible
            if (opt$status!=0 | abs(rr$objective-opt$objval)>=1e-4) {
                msg <- paste0("Linear program finds rejection ", opt$objval,
                              ". Direct approach finds rejection ",
                              rr$objective,
                              ". Difference>0.001. This happened for ",
                              "chi = ", chi, ", m2 = ", m2, ", kappa = ", kappa)
                warning(msg)
                return(list(alpha=NA, x=NA, p=NA))
            }
        }
        list(alpha=unname(rr$objective), x=xs0, p=ps0)
    }
}

## Critical value from Armstrong and Kolesar
CVb <- function(B, alpha=0.05) {
    idx <- B<10
    r <- B+stats::qnorm(1-alpha)
    r[idx] <- sqrt(stats::qchisq(1-alpha, df = 1, ncp = B[idx]^2))
    r
}

#' Compute average coverage critical value under moment constraints.
#'
#' Computes the critical value \eqn{cva_{\alpha}(m_{2}, \kappa)}{cva_alpha(m_2,
#' kappa)} from Armstrong, Kolesár, and Plagborg-Møller (2020).
#' @param m2 Bound on second moment of the normalized bias, \eqn{m_{2}}{m_2}
#' @param kappa Bound on the kurtosis of the normalized bias,
#'     \eqn{\kappa}{kappa}
#' @param alpha Determines confidence level, \eqn{1-\alpha}{1-alpha}.
#' @param check If \code{TRUE}, verify accuracy of the solution by checking that
#'     the implied least favorable distribution satisfies the \code{m2} and
#'     \code{kappa} constraints and yields the same non-coverage rate. If this
#'     fails (perhaps due to numerical accuracy issues), solve a finite-grid
#'     approximation (by discretizing the support of the normalized bias) to the
#'     primal linear programing problem, and check that it agrees with the dual
#'     solution.
#' @return Returns a list with 4 components: \describe{
#'
#' \item{\code{cv}}{Critical value for constructing two-sided confidence
#' intervals.}
#'
#' \item{\code{alpha}}{The argument \code{alpha}.}
#'
#' \item{\code{x}}{Support points for the least favorable distribution for the
#' squared normalized bias, \eqn{b^2}.}
#'
#' \item{\code{p}}{Probabilities associated with the support points.}
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
#' ## Critical value without imposing a constraint on kurtosis
#' cva(1, kappa=Inf)
#' ## With a constraint
#' cva(1, kappa=3)
#' @export
cva <- function(m2, kappa=Inf, alpha=0.05, check=TRUE) {
    if (kappa==1 | m2==0) {
        list(cv=CVb(sqrt(m2), alpha), alpha=alpha, x=c(0, m2), p=c(0, 1))
    } else {
        ## For very large values of m2, assume kappa=Inf
        if (1/m2 < tol & kappa!=Inf) {
            res <- cva(m2, kappa=Inf, alpha=alpha, check=FALSE)
            if (kappa < sum(res$x^2*res$p)/m2^2)
                warning("Value of m2 (", m2, ") is too large to reliably",
                        " compute the critical value. ",
                        "Assuming kappa constraint not binding")
            return(res)
        }
        ## limits: critical values under kappa=1 and kappa=Inf to get bounds on
        ## cv
        lo <- CVb(sqrt(m2), alpha)-0.01
        limits <- c(lo, NA)
        ## Use Chebyshev inequality
        up <- sqrt((1+m2)/alpha)
        ## If upper rejection rate close to alpha, keep it
        if (abs(rho0(m2, up)-alpha) < 9e-6)
            limits[2] <- up
        else
            limits[2] <- stats::uniroot(function(chi) rho0(m2, chi)-alpha,
                                        c(lo, up), tol=tol)$root

        ## If rejection rate is already close to alpha, keep cv under kappa=Inf
        if (rho(m2, kappa, limits[2])$alpha-alpha < -1e-5)
            cv <- stats::uniroot(function(chi) rho(m2, kappa, chi,
                                            check=FALSE)$alpha-alpha,
                          limits, tol=tol)$root
        else
            cv <- limits[2]
        c(cv=cv, rho(m2, kappa, cv, check=check, len=5000))
    }
}
