tol <- 1e-12

r <- function(t, c0) {
    idx <- sqrt(t)-c0 <= 5
    r <- rep(1, length(t))
    ## Use normal rather than stats::pchisq(..., lower.tail=FALSE), since
    ## non-central chi^2 is not accurate for extreme values. It's also much
    ## faster.
    r[idx] <- stats::pnorm(-sqrt(t[idx])-c0)+stats::pnorm(sqrt(t[idx])-c0)
    r
}

## Derivative of r
r1 <- function(t, c0) {
    ## Apply L'Hospital's rule. This gives maximum absolute error 1e-9
    ifelse(t >= 1e-8,
    (stats::dnorm(sqrt(t)-c0)-stats::dnorm(sqrt(t)+c0)) / (2*sqrt(t)),
    c0*stats::dnorm(c0))
}


## Second Derivative of r
r2 <- function(t, c0) {
    ## Apply L'Hospital's rule 3x. This gives maximum abs error about 8e-8
    ifelse(t >= 2e-6,
           (stats::dnorm(sqrt(t)+c0)*(c0*sqrt(t)+t+1)+
            stats::dnorm(sqrt(t)-c0)*(c0*sqrt(t)-t-1)) / (4*t^(3/2)),
           stats::dnorm(c0)*c0*(c0^2-3)/6)
}

## Third Derivative of r
r3 <- function(t, c0) {
    ## Apply L'Hospital's rule 5x. This gives maximum abs error about 3e-6
    ifelse(t >= 2e-4,
    (stats::dnorm(c0-sqrt(t))*(t^2-2*c0*t^(3/2)+(2+c0^2)*t-3*c0*sqrt(t)+3)
    -stats::dnorm(c0+sqrt(t))*(t^2+2*c0*t^(3/2)+(2+c0^2)*t+3*c0*sqrt(t)+3))
    /(8*t^(5/2)),
    stats::dnorm(c0)*(c0^5-10*c0^3+15*c0)/60)
}


## find t0 and inflection point
rt0 <- function(c0) {
    ## Find point where we touch origin
    f0 <- function(t) r(t, c0)-t*r1(t, c0)-r(0, c0)
    if (c0^2<3) {
        t0 <- ip <- 0
    } else {
        t0 <- stats::uniroot(f0, c(c0^2-3, 5*c0^2), tol=tol)$root
        ip <- stats::uniroot(r2, c(c0^2-3, c0^2), tol=tol, c0=c0)$root
    }
    list(t0=t0, ip=ip)
}

rho <- function(t, c0) {
    t0 <- rt0(c0)$t0
    idx <- (t >= t0)
    res <- r(t0, c0) + (t-t0)*r1(t0, c0)
    if (any(idx))
        res[idx] <- r(t[idx], c0)
    res
}


delta <- function(x, x0, c0) {
    ## Apply L'Hospital's rule 2x. This gives maximum abs error about 1e-7
    ifelse(abs(x-x0)<1e-4,
           r2(x0, c0)/2,
           (r(x, c0) - r(x0, c0) - r1(x0, c0)*(x-x0)) / (x-x0)^2)
}

delta1 <- function(x, x0, c0) {
    ## Apply L'Hospital's rule 2x. This gives maximum abs error about 7e-6
    ifelse(abs(x-x0)<1e-3,
           r3(x0, c0)/6,
           ((r1(x, c0)+r1(x0, c0))-2*(r(x, c0)-r(x0, c0))/(x-x0))/(x-x0)^2)
}

lam <- function(x0, c0) {
    ## Check derivatives at 0, inflection point, t0, and x0 If we're above
    ## inflection point, then maximum is below it, and it's at zero if
    ## derivative at zero is negative. Otherwise between 0 and inflection point
    xs <- sort(unlist(rt0(c0)))
    if (x0 >= xs[1]) {
        xs <- c(0, xs[1])
    } else {
        xs <- c(0, x0, xs)
    }
    ## We want >= since derivative at zero may be zero if c0 is large
    der <- delta1(xs, x0, c0) >= 0

    if (all(der<=0))
        return(list(lam=delta(0, x0, c0), x0=0))
    else if (all(diff(der)<=0)) {
        ## Function first increasing, then decreasing
        start <- xs[which.min(der)-1]
        end <- xs[which.min(der)]
    } else {
        stop("There seem to be multiple optima in lam")
    }

    rr1 <- stats::optimize(function(x) -delta(x, x0, c0),
                           c(start, end), tol=tol)

    list(lam=-rr1$objective, x0=unname(rr1$minimum))
}

rho2 <- function(mu, c0, check=TRUE, len=5000) {
    r0 <- rho(mu[1], c0)
    t0 <- rt0(c0)$t0
    if (mu[2] == mu[1]^2) {
        list(size=r(mu[1], c0), x=c(0, mu[1]), p=c(0, 1))
    } else if (mu[1] >= t0) {
        ## Concave part of parameter space
        list(size=r0, x=c(0, mu[1]), p=c(0, 1))
    } else if (mu[2] >= mu[1]*t0) {
        ## LF under rho: (0, t0) wp (1-mu[1]/t0, mu[1]/t0), E[t^2]=mu[1]*t0
        ## So here mu4 doesn't bind
        list(size=r0, x=c(0, t0), p=c(1-mu[1]/t0, mu[1]/t0))
    } else {
        ## First determine t1
        t1 <- lam(0, c0)$x0
        lammax <- function(x0) {
            ## Delta maximized at 0
            ifelse(x0 >= t1, delta(0, x0, c0), lam(x0, c0)$lam)
        }
        obj <-  function(x0)
            r(x0, c0)+r1(x0, c0)*(mu[1]-x0)+
                lammax(x0)*(mu[2]-2*x0*mu[1]+x0^2)
        rr <- stats::optimize(obj, interval=c(0, t0), tol=1e-8)
        ## LF points
        xs0 <- sort(c(rr$minimum, lam(rr$minimum, c0)$x0))

        ## Now double-check we found the optimum by solving the primal
        if (check) {
            ## Add mu[1] here for cases where it's very small, so we can satisfy
            ## the constraint
            xs <- sort(unique(c(mu, xs0, seq(0, t0, length.out=len))))
            ## Find the optimal solution; use <= for last constraint
            opt <-  lpSolve::lp(direction="max",
                                objective.in = r(xs, c0),
                                const.mat = rbind(xs^0, xs, xs^2),
                                const.dir = c("==", "==", "<="),
                                const.rhs = c(1, mu))
            ## 0: success, 2: infeasible
            if (opt$status!=0 | abs(rr$objective-opt$objval)>=1e-4) {
                msg <- paste0("Linear program finds rejection", opt$objval,
                              "direct approach finds rejection", rr$objective,
                              "Difference>0.001. This happened for",
                              "c0=", c0, "mu[1]=", mu[1], "mu[2]=", mu[2])
                warning(msg)
                return(list(size=NA, x=NA, p=NA))
            }
        }
        p <- (mu[1]-xs0[2])/(xs0[1]-xs0[2])
        list(size=unname(rr$objective), x=xs0, p=c(p, 1-p))
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
#' @param B Bound on the square root of the average squared bias
#' @param kurt Bound on the kurtosis of the bias
#' @param alpha Determines CI level, \eqn{1-\alpha}{1-alpha}.
#' @param check Check critical value by solving linear program
#' @return Critical value for constructing two-sided confidence intervals. TODO:
#'     also returns least favorable info
#' @examples
#' ## Without imposing a constraint on kurtosis
#' cva(1, kurt=Inf)
#' ## With a constraint
#' cva(1, kurt=3)
#' @export
cva <- function(B, kurt=Inf, alpha=0.05, check=TRUE) {
    if (kurt==1 | B==0) {
        list(cv=CVb(B, alpha), size=alpha, x=c(0, B), p=c(0, 1))
    } else {
        ## Critical value under kappa=Inf and kappa=1 to get bounds on cv
        limits <- c(CVb(B, alpha)-0.01, stats::qnorm(1-alpha/2)+10*B+0.01)
        limits[2] <- stats::uniroot(function(c0) rho(B^2, c0)-alpha, limits,
                             tol=tol)$root
        mu <- c(B^2, kurt*ifelse(B>0, B, 1)^4)
        if (rho2(mu, limits[2])$size-alpha < -1e-5)
            cv <- stats::uniroot(function(c0) rho2(mu, c0,
                                            check=FALSE)$size-alpha,
                          limits, tol=tol)$root
        else
            cv <- limits[2]
        c(cv=cv,
             rho2(c(B^2, kurt*ifelse(B>0, B, 1)^4), cv, check=check, len=5000))
    }
}
