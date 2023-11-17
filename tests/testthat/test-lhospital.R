context("Check L'Hospital rule approximations")

test_that("Test L'Hospital rule approximations are accurate", {

    xx <- seq(0, 10, length.out=10000)

    ff <- function(c0) r1(1e-8, c0)-c0*dnorm(c0)
    maxx <- max(abs(vapply(xx, ff, numeric(1))))
    expect_lt(maxx, 1e-8)
    expect_gt(maxx, 1e-10)

    ff <- function(c0) r2(2e-6, c0)-dnorm(c0)*c0 * (c0^2-3)/6
    maxx <- max(abs(vapply(xx, ff, numeric(1))))
    expect_lt(maxx, 1e-7)
    expect_gt(maxx, 1e-8)

    ff <- function(c0) r3(2e-4, c0)-dnorm(c0) * (c0^5-10*c0^3+15*c0)/60
    maxx <- max(abs(vapply(xx, ff, numeric(1))))
    expect_lt(maxx, 4e-6)
    expect_gt(maxx, 1e-6)

    xx <- seq(0, 20, length.out=1000)
    ff <- function(c0) max(abs(delta(xx, xx+1e-4, c0) -r2(xx, c0)/2))
    maxx1 <- max(abs(vapply(seq(1, 20, length.out=100), ff, numeric(1))))
    ff <- function(c0) max(abs(delta(xx+1e-4, xx, c0) -r2(xx+1e-4, c0)/2))
    maxx2 <- max(abs(vapply(seq(1, 20, length.out=100), ff, numeric(1))))
    expect_lt(max(maxx1, maxx2), 4e-6)
    expect_gt(max(maxx1, maxx2), 1e-6)

    ff <- function(c0) max(abs(delta1(xx, xx+1e-3, c0) -r3(xx, c0)/6))
    maxx1 <- max(abs(vapply(seq(1, 20, length.out=100), ff, numeric(1))))
    ff <- function(c0) max(abs(delta1(xx+1e-3, xx, c0) -r3(xx+1e-3, c0)/6))
    maxx2 <- max(abs(vapply(seq(1, 20, length.out=100), ff, numeric(1))))
    expect_lt(max(maxx1, maxx2), 4e-6)
    expect_gt(max(maxx1, maxx2), 1e-8)
})
