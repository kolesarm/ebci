context("Check critical value calculations")

test_that("Simple critical value sanity checks", {

    expect_equal(cva(B=0, kurt=1)$cv, stats::qnorm(0.975))
    expect_equal(cva(B=0, kurt=4)$cv, stats::qnorm(0.975))
    expect_equal(cva(B=0, kurt=Inf)$cv, stats::qnorm(0.975))
    expect_equal(cva(B=1, kurt=1)$cv, CVb(B=1))
    expect_equal(cva(B=100, kurt=1)$cv, CVb(B=100))
    expect_equal(cva(B=0, kurt=1)$cv, CVb(B=0))

    expect_equal(cva(B=5, kurt=3, check=TRUE)$cv, 11.88358367)
    expect_lt(abs(cva(B=5, kurt=1.000001)$cv-CVb(B=5)), 1e-5)
})

test_that("Check large values", {
    expect_equal(lam(0, 42.2201)$x0, 1932.78377442)
    expect_equal(cva(10, 3)$cv, 24.86172217)
})
