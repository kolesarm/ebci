context("Check optimal and EB shrinkage")

test_that("Check shrinkage with cv_tbl ", {
    ## cva_tbl too coarse
    cv_tbl <- cva_tbl[cva_tbl$kappa==3 & cva_tbl$alpha==0.05, ]
    cv_tbl2 <- cv_tbl[cv_tbl$m2<=1, ]
    expect_warning(w_opt(0.1^2, 3, cv_tbl=cv_tbl2))
    expect_lt(max(abs(unlist(w_opt(1, 3, cv_tbl=cv_tbl))-unlist(w_opt(1, 3)))),
              3e-4)
    expect_lt(max(abs(unlist(w_opt(0.01, 3, cv_tbl=cv_tbl))-
                      unlist(w_opt(0.01, 3)))), 4e-4)
    cv_tbl <- cva_tbl[cva_tbl$kappa==Inf & cva_tbl$alpha==0.05, ]
    expect_lt(max(abs(unlist(w_opt(0.01, Inf, cv_tbl=cv_tbl))-
                      unlist(w_opt(0.01, Inf)))), 3e-4)
    ## eb shrinkage should be close to optimal
    expect_lt(max(abs(unlist(w_opt(20^2, 4))-unlist(w_eb(20^2, 4)))),
              5e-5)

})

test_that("Table 1 in paper", {
    ts <- cz[!is.na(cz$theta25), ]
    eb_table <- function(res, wgt) {
        res$df$sn <- res$mu2[1]/res$df$se^2
        l1 <- vapply(res$df, weighted.mean, FUN.VALUE=numeric(1), w=wgt)
        ret <- c("sqrt{mu_2}"=sqrt(res$mu2[1]),
                 "kappa"=res$kappa[1],
                 l1[length(l1)],
                 res$delta, l1[1:7], l1[4]/l1[5:7])
        ret
    }
    res2 <- ebci(theta25 ~ stayer25, ts, se25, pop, tstat=TRUE)
    testthat::expect_equal(unname(round(eb_table(res2, ts$pop), 6)),
                           c(1.195109,   3.200384, 193.144423,  -7.458751,
                             0.157330,   0.588187,   0.613905,   0.102815,
                             0.221920,   0.221464,   0.220051,   0.286923,
                             1.002058,   1.008492,   0.773446))
    ## Only large CZ to  make it faster
    df <- cz[sort(cz$pop, index.return=TRUE, decreasing=TRUE)$ix[1:10], ]
    i2 <- ebci(formula=theta25~stayer25, data=df, se=se25, alpha=0.1, wopt=TRUE,
               fs_correction="none")
    i3 <- ebci(formula=theta25~0, data=df, se=se25, alpha=0.05, kappa=Inf,
               wopt=TRUE)

    testthat::expect_equal(unname(round(eb_table(i2), 6)),
                           c(0.060852, 5.959142, 0.956639, -0.913932, 0.020049,
                             0.456490, 0.526469, 0.115947, 0.078879, 0.075404,
                             0.073284, 0.113034, 1.046083, 1.076356, 0.697836))
    testthat::expect_equal(unname(round(eb_table(i3), 6)),
                           c(0.085168, Inf, 1.873921, 0.612815, 0.695995,
                             0.063068, 0.114249, 0.107078, 0.102651, 0.134688,
                             1.066973, 1.112989, 0.848247))
})


test_that("ebci function problems", {
    ## negative mu2 estimate
    df <- cz[sort(cz$pop, index.return=TRUE, decreasing=TRUE)$ix[1:100], ]
    expect_warning(r <- ebci(formula=theta25~stayer25, data=df, kappa=Inf,
                             se=se25, alpha=0.1, cores=1, fs_correction="none"))
})

test_that("Zero bias estimate", {
    expect_identical(cva(m2=Inf, kappa=1)$cv, NA)
    expect_warning(r <- w_opt(S=0, kappa=3))
    expect_true(is.na(r$length))
    expect_true(is.na(w_eb(S=0, kappa=1)$length))

})

test_that("fs corrections", {
    df <- cz[sort(cz$pop, index.return=TRUE, decreasing=TRUE)$ix[1:40], ]
    r <- ebci(formula=theta25~stayer25, data=df, se=se25,
              alpha=0.1, fs_correction="FPLIB", cores=1, wopt=TRUE)
    r2 <- ebci(formula=theta25~stayer25, data=df, se=se25,
              alpha=0.1, wopt=FALSE, fs_correction="PMT")
    expect_equal(unname(c(sqrt(r$mu2[1]), r$kappa[1])),
                 c(0.04386307, 34.60536787))
    expect_equal(as.vector(summary(r$df$len_op)),
                 c(0.0510135245, 0.06513867, 0.071050167, 0.069886713,
                   0.076510241, 0.080472608))
    expect_equal(unname(sqrt(r2$mu2)), c(0.0247533358, 0.016287771))
    expect_equal(unname(r2$kappa), c(475.964641182, 1719.44225761))

    ## When noise is small, FPLIB should be the similar
    set.seed(42)
    th <- 42+rnorm(500)/4
    se <- sqrt(seq(1, 10, length.out=500))/10
    Y <- th + rnorm(500)*se
    rt <- ebci(formula=Y~1, se=se, weights=1/se^2, wopt=FALSE)
    rn <- ebci(formula=Y~1, se=se, weights=1/se^2, wopt=FALSE,
               fs_correction = "none")
    rb <- ebci(formula=Y~1, se=se, weights=1/se^2, wopt=FALSE,
               fs_correction = "FPLIB")
    expect_identical(c(rt$kappa, rt$mu2), c(rn$kappa, rn$mu2))
    expect_lt(max(abs(c(rb$kappa, rb$mu2)-c(rn$kappa, rn$mu2))), 1e-2)
    est <- colSums(abs(rt$df-rb$df))
    expect_lt(max(est[!is.na(est)]), 1e-2)
    ## Test format output while we're at it
    expect_identical(dim(rt$df), c(500L, 11L))
    expect_identical(dim(rb$df), c(500L, 11L))
    rtt <- ebci(formula=Y~1, se=se, weights=1/se^2, wopt=TRUE,
               fs_correction = "FPLIB", tstat=TRUE)
    expect_identical(dim(rt$df), dim(rtt$df))
})
