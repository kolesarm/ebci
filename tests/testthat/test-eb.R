context("Check optimal and EB shrinkage")

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
                             se=se25, alpha=0.1, fs_correction="none"))
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
              alpha=0.1, fs_correction="FPLIB", wopt=TRUE)
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
    expect_identical(dim(rt$df), c(500L, 13L))
    expect_identical(dim(rb$df), c(500L, 13L))
    rtt <- ebci(formula=Y~1, se=se, weights=1/se^2, wopt=TRUE,
                fs_correction = "FPLIB")
    expect_identical(dim(rt$df), dim(rtt$df))
})

test_that("Pre-computing critical values", {
    y <- c(-0.0912394115425668, 0.270659035981612, 0.0505183828184417,
           0.401285011472904, -0.116494855854479, -0.334314741812537,
           -0.18928733418011, 0.165181263555627, 0.132184880945018,
           -0.131267800385369, 0.584675521401188, 0.105172038899963,
           0.470647126771099, 0.289695479626726, 0.202471324044271,
           0.0146653486070987, 0.568678856504733, 0.216916114930276,
           -0.259804655342694, -0.582460577829165, -0.449075030332067,
           -0.251549682686536, -0.331692925639294, -0.116896950204574,
           -0.175038163494127, 0.178108950152456, 0.211683997890061,
           -0.911066040123521, -0.164931392459329, -0.161697036212087,
           -0.0662609931760586, 0.0316248074589969)
    se <- c(0.224386445769098, 0.257751786695356, 0.276711829264707,
            0.403355134534814, 0.343841207609935, 0.260088725292071,
            0.239258750261168, 0.177225564240537, 0.165070998037012,
            0.215345806604199, 0.329392513934246, 0.166990070498265,
            0.356719957739564, 0.338361696727126, 0.26596125123086,
            0.278943346050078, 0.286066216054302, 0.211277847211765,
            0.197890061668623, 0.254834949524809, 0.161320618248102,
            0.159717394888581, 0.190422698118735, 0.11198110341991,
            0.128750097633114, 0.137084993056956, 0.144134717528643,
            0.198405131349222, 0.235877391582392, 0.129194532604923,
            0.0592319350677268, 0.45161460229795)

    ## Should be no warning, previously we did get them
    expect_silent(c2 <- ebci(y ~ 0, data.frame(y=y, se=se),
                             se, wopt=TRUE, alpha=0.05, kappa=Inf))
    expect_equal(dim(c2$X), c(32L, 0L))
})
