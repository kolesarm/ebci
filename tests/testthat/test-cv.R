context("Check critical value calculations")

test_that("Simple critical value sanity checks", {
    expect_equal(cva(B=0, kurt=1)$cv, stats::qnorm(0.975))
    expect_equal(cva(B=0, kurt=4)$cv, stats::qnorm(0.975))
    expect_equal(cva(B=0, kurt=4, alpha=0.2)$cv, stats::qnorm(0.9))
    expect_equal(cva(B=0, kurt=Inf)$cv, stats::qnorm(0.975))
    expect_equal(cva(B=1, kurt=1)$cv, CVb(B=1))
    expect_equal(cva(B=100, kurt=1)$cv, CVb(B=100))
    expect_equal(cva(B=0, kurt=1)$cv, CVb(B=0))

    expect_equal(cva(B=1, kurt=10000, check=TRUE, alpha=0.2)$cv,
                 cva(B=1, kurt=Inf, check=TRUE, alpha=0.2)$cv)
    ## Here under kurt=Inf, solution has kurtosis (sol$x[2]^2*sol$p[2])/16=17.4
    expect_equal(cva(B=2, kurt=20, check=TRUE, alpha=0.05)$cv,
                 cva(B=2, kurt=Inf, check=TRUE, alpha=0.05)$cv)
    ## Here under kurt=Inf, solution has kurtosis (sol$x[2]^2*sol$p[2])/16=16.6
    expect_lt(cva(B=0.5, kurt=15, check=TRUE, alpha=0.05)$cv,
                 cva(B=0.5, kurt=Inf, check=TRUE, alpha=0.05)$cv)

    expect_equal(cva(B=5, kurt=3, check=TRUE)$cv, 11.88358367)
    expect_equal(cva(B=1, kurt=3, check=TRUE, alpha=0.2)$cv, 1.8494683)
    expect_lt(abs(cva(B=5, kurt=1.000001)$cv-CVb(B=5)), 1e-5)
    expect_lt(abs(cva(B=0.01, kurt=15)$cv-CVb(B=0.01)), 1e-5)
    expect_lt(abs(cva(B=0.01, kurt=15, alpha=0.2)$cv-CVb(B=0.01, alpha=0.2)), 1e-5)
})

test_that("Check large values", {
    expect_equal(lam(0, 42.2201)$x0, 1932.78377442)
    expect_equal(cva(10, 3)$cv, 24.86172217)
    expect_equal(rho2(mu=c(1, 3), chi=11.9699639845401)$size, 6.0216e-05)
})

test_that("Check delta1 has crosses zero at most once", {
    checkdelta <- function(chi) {
        ts <- sort(unlist(rt0(chi)))
        xs <-  c(seq(0, ts[1], length.out=100), seq(ts[1], ts[2], length.out=100))
        x0s <- xs
        ## If chi <=2.856, then derivative should always be negative
        if (chi <=  2.856) {
            expect_equal(sum(sapply(x0s, function(x0)
                sum(delta1(xs, x0, chi) >= 0))),
                0L)
        } else {
            ## Derivative should be first positive, then negative. Only need
            ## this to hold on [0, t1] if x>t1
            der <- function(x0) {
                    delta1(xs[1:(100*(x0 >= ts[1])+(x0<ts[1])*200)], x0, chi)
            }
            expect_equal(sum(sapply(x0s, function(x0)
                is.unsorted(!(der(x0) >= 0)))),
                0L)
        }
    }
    chis <- seq(sqrt(3), 50, length.out=300)
    vapply(chis, checkdelta, numeric(1))
})
