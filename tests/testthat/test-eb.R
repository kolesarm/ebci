context("Check optimal and EB shrinkage")

test_that("Check shrinkage with cv_tbl ", {
    ## cva_tbl too coarse
    cv_tbl <- cva_tbl[cva_tbl$kurt==3 & cva_tbl$alpha==0.05, ]
    cv_tbl2 <- cv_tbl[cv_tbl$B<=1, ]
    expect_warning(w_opt(0.1, 3, cv_tbl=cv_tbl2))
    expect_lt(max(abs(unlist(w_opt(1, 3, cv_tbl=cv_tbl))-unlist(w_opt(1, 3)))),
              1e-4)
    expect_lt(max(abs(unlist(w_opt(0.1, 3, cv_tbl=cv_tbl))-unlist(w_opt(0.1, 3)))),
              4e-4)
    cv_tbl <- cva_tbl[cva_tbl$kurt==Inf & cva_tbl$alpha==0.05, ]
    expect_lt(max(abs(unlist(w_opt(0.1, Inf, cv_tbl=cv_tbl))-unlist(w_opt(0.1, Inf)))),
              3e-4)
    ## eb shrinkage should be close to optimal
    expect_lt(max(abs(unlist(w_opt(20, 4))-unlist(w_eb(20, 4)))),
              5e-5)

})
