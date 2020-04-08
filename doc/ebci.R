## ---- include=FALSE, cache=FALSE----------------------------------------------
library("knitr")
knitr::opts_knit$set(self.contained = FALSE)
knitr::opts_chunk$set(tidy = TRUE, collapse=TRUE, comment = "#>",
                      tidy.opts=list(blank=FALSE, width.cutoff=55))

## -----------------------------------------------------------------------------
library("ebci")
## If B=sqrt(m_2)=0, then we get the usual critical value
cva(B=0, kappa=Inf, alpha=0.05)$cv
## Otherwise the critical value is larger. Suppose m_2=4, so that B=2:
cva(B=2, kappa=Inf, alpha=0.05)$cv
## Imposing a constraint on kurtosis tightens it
cva(B=2, kappa=3, alpha=0.05)$cv

## -----------------------------------------------------------------------------
## For illustration, only use 20 largest CZ
df <- cz[sort(cz$pop, index.return=TRUE, decreasing=TRUE)$ix[1:20], ]

## As Y_i, use fixed effect estimate theta25 of causal effect of neighborhood for children with parents at the 25th percentile of income distribution. The standard error for this estimate is se25. As predictors use average outcome for permanent residents (stayers), stayer25. Let us use 90% CIs, and let us set kappa=Inf (so we only don't use a constraint on kurtosos) in order to decrease the computational time. Given this, we do not need to parallelize the computation, and can set cores=1:
r <- ebci(formula=theta25~stayer25, data=df, se=se25, alpha=0.1, cores=1, kappa=Inf)

## -----------------------------------------------------------------------------
r$delta

## -----------------------------------------------------------------------------
c(r$sqrt_mu2, r$kappa)

## -----------------------------------------------------------------------------
names(r$df)

## -----------------------------------------------------------------------------
cva(B=((1-1/r$df$w_eb[1])*r$sqrt_mu2/r$df$se[1]), Inf, alpha=0.1)$cv*
r$df$w_eb[1]*r$df$se[1]
r$df$len_eb[1]

## -----------------------------------------------------------------------------
knitr::kable(data.frame(name=paste0(df$czname, ", ", df$state), estimate=r$df$th_eb,
lower_ci=r$df$th_eb-r$df$len_eb, upper_ci=r$df$th_eb+r$df$len_eb), digits=3)

## -----------------------------------------------------------------------------
knitr::kable(data.frame(name=paste0(df$czname, ", ", df$state), estimate=r$df$th_op,
lower_ci=r$df$th_op-r$df$len_op, upper_ci=r$df$th_op+r$df$len_op), digits=3)

## -----------------------------------------------------------------------------
mean(r$df$len_op)/mean(r$df$len_eb)

## -----------------------------------------------------------------------------
mean(r$df$ncov_pa)

