context("Check critical value calculations")

test_that("Numerical issues with rt0", {

    ## Convert return value from list to vector
    rt00 <- function(chi) unlist(rt0(chi))

    chis <- sqrt(3)+ exp(seq(log(1e-12), log(1e-3), length.out=200))
    f0 <- vapply(chis, rt00, numeric(2))
    expect_equal(sum(f0[1, ]<f0[2, ]), 0L)

    chis <- 10^seq(0, 5, length.out=200)
    f1 <- vapply(chis, rt00, numeric(2))
    expect_equal(sum(f1[1, ]<f1[2, ]), 0L)

    ## For large values of m2, Chebyshev inequality should be sharp
    m2 <- 1/1.11022e-16
    expect_equal(cva(m2)$cv/sqrt((1+m2)/0.05), 1L)
    expect_equal(cva(m2/100)$cv/sqrt((1+m2/100)/0.05), 1L)

    ## Assume kappa constraint not binding for large m2
    expect_warning(cv0 <- cva(m2, kappa=3, check=FALSE)$cv)
    expect_warning(rho(m2, kappa=3, chi=cv0))

    m2s <- seq(1/1000, 100, length.out=10)*m2
    vals <- vapply(m2s,
           function(r) cva(r, kappa=400, check=FALSE, alpha=0.1)$cv,
           numeric(1))
    expect_equal(vals/sqrt((1+m2s)/0.1), rep(1, length(m2s)))

    ## This is increasing, maximum at 4.8518e-26, close to 1
    expect_lt(abs(lam(x0=180144474255572736,
                     chi=424434295.36931574345)$lam/4.851805e-26 - 1L), 1e-6)
    ## This one is decreasing, max should be at 0
    expect_equal(lam(x0=12577.803781, chi=109.857326)$x0, 0L)
    ## First derivative evaluates to negative, but should probably be possitive;
    ## this shouldn't throw us off
    expect_equal(lam(x0 = 0.000759676644983789, chi = 2.85710241617874)$lam,
                 0.00827903)


})

test_that("Simple critical value sanity checks", {
    expect_equal(cva(m2=0, kappa=1)$cv, stats::qnorm(0.975))
    expect_equal(cva(m2=0, kappa=4)$cv, stats::qnorm(0.975))
    expect_equal(cva(m2=0, kappa=4, alpha=0.2)$cv, stats::qnorm(0.9))
    expect_equal(cva(m2=0, kappa=Inf)$cv, stats::qnorm(0.975))
    expect_equal(cva(m2=1, kappa=1)$cv, CVb(B=1))
    expect_equal(cva(m2=100^2, kappa=1)$cv, CVb(B=100))
    expect_equal(cva(m2=0, kappa=1)$cv, CVb(B=0))
    expect_equal(rho(m2=4, kappa=1,
                     chi=cva(m2=4, kappa=1)$cv)$alpha, 0.05)

    expect_equal(cva(m2=1, kappa=10000, check=TRUE, alpha=0.2)$cv,
                 cva(m2=1, kappa=Inf, check=TRUE, alpha=0.2)$cv)
    ## Here under kappa=Inf, solution has kurtosis (sol$x[2]^2*sol$p[2])/16=17.4
    expect_equal(cva(m2=4, kappa=20, check=TRUE, alpha=0.05)$cv,
                 cva(m2=4, kappa=Inf, check=TRUE, alpha=0.05)$cv)
    ## Here under kappa=Inf, solution has kurtosis (sol$x[2]^2*sol$p[2])/16=16.6
    expect_lt(cva(m2=0.25, kappa=15, check=TRUE, alpha=0.05)$cv,
                 cva(m2=0.25, kappa=Inf, check=TRUE, alpha=0.05)$cv)

    expect_equal(cva(m2=25, kappa=3, check=TRUE)$cv, 11.88358367)
    expect_equal(cva(m2=1, kappa=3, check=TRUE, alpha=0.2)$cv, 1.8494683)
    expect_lt(abs(cva(m2=25, kappa=1.000001)$cv-CVb(B=5)), 1e-5)
    expect_lt(abs(cva(m2=0.01^2, kappa=15)$cv-CVb(B=0.01)), 1e-5)
    expect_lt(abs(cva(m2=0.01^2,
                      kappa=15, alpha=0.2)$cv-CVb(B=0.01, alpha=0.2)),
              1e-5)

    ## Test LF distribution
    r1 <- cva(m2=4, kappa=1.000001, alpha=0.1)
    r2 <- cva(m2=4, kappa=1L, alpha=0.1)
    expect_lt(r1$x[which.max(r1$p)]-r2$x[which.max(r2$p)], 1e-3)
    expect_lt(abs(r1$p[which.max(r1$p)]-1), 1e-5)
})

test_that("Check large values", {
    expect_equal(lam(0, 42.2201)$x0, 1932.78377442)
    expect_equal(cva(10^2, 3)$cv, 24.86172217)
    expect_equal(rho(m2=1, kappa=3, chi=11.9699639845401)$alpha, 6.0216e-05)
    ## Previously this was failing
    expect_equal(cva(m2=2.8^2, alpha=0.1, kappa=8)$cv, 7.0779747158)
    expect_equal(rho(2.8^2, 8, 7.1324290147839279896)$alpha, 0.09805099)

    store_cva <- function(m2, kappa, alpha, check=TRUE) {
        df <- expand.grid(m2=m2, kappa=kappa, alpha=alpha)
        cvj <- function(j)
            cva(df$m2[j], kappa=df$kappa[j], df$alpha[j], check)$cv
        df$cv <- vapply(seq_len(nrow(df)), cvj, numeric(1))
        df
    }
    m2 <-  seq(0, 5, by=1)^2
    ## Takes 30s on X1 Carbon
    kappa <- c(1:10, Inf)
    expect_silent(store_cva(m2, kappa, alpha=0.1))
})

test_that("Check delta1 has crosses zero at most once", {
    checkdelta <- function(chi) {
        ts <- sort(unlist(rt0(chi)))
        xs <-  c(seq(0, ts[1], length.out=100),
                 seq(ts[1], ts[2], length.out=100))
        x0s <- xs
        ## If chi <=2.856, then derivative should always be negative
        if (chi <=  2.856) {
            expect_equal(sum(vapply(x0s, function(x0)
                sum(delta1(xs, x0, chi) >= 0), numeric(1))),
                0L)
        } else {
            ## Derivative should be first positive, then negative. Only need
            ## this to hold on [0, t1] if x>t1
            der <- function(x0) {
                    delta1(xs[1:(100*(x0 >= ts[1])+(x0<ts[1])*200)], x0, chi)
            }
            expect_equal(sum(vapply(x0s, function(x0)
                is.unsorted(!(der(x0) >= 0)), logical(1))),
                0L)
        }
    }
    chis <- seq(sqrt(3), 50, length.out=300)
    vapply(chis, checkdelta, numeric(1))
})
