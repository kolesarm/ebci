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
        res$df$sn <- res$sqrt_mu2^2/res$df$se^2
        l1 <- vapply(res$df, weighted.mean, FUN.VALUE=numeric(1), w=wgt)
        ret <- c("sqrt{mu_2}"=res$sqrt_mu2,
                 "kappa"=res$kappa,
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
    ## i2 <- ebci(theta25 ~ stayer25, ts[ts$pop> 300000, ],
    ##            se25, 1/se25^2, tstat=FALSE, kappa=Inf)
    i2 <- ebci(theta25 ~ stayer25, ts[ts$pop> 300000, ],
               se25, 1/se25^2, tstat=FALSE, kappa=Inf, cores=1)
    testthat::expect_equal(unname(round(eb_table(i2,
                                                 ts$pop[ts$pop>300000]), 6)),
                           c(0.063182, Inf,  0.641786, -1.081536,  0.023759,
                             0.315798, 0.411483, 0.147695, 0.107307, 0.091945,
                             0.084915, 0.194525, 1.167075, 1.263687, 0.551635))
})
