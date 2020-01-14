## The original data, online_table_3.dta and online_table_4.dta, has been
## downloaded from
## https://opportunityinsights.org/wp-content/uploads/2018/04/replicate-4.zip,
## subdirectoty data (not online tables, which contains the source files for
## these two .dta files)

## 1. CZ dataset: 741 obs by 459 variables
dt <- readstata13::read.dta13("online_table_3.dta")
## cc2: specification
## kr26: HH income at age 26
cz <- data.frame(cz=dt$cz,
                 czname=dt$czname,
                 state=dt$stateabbrv,
                 pop=dt$pop2000,
                 theta25=dt$Bj_p25_czkr26_cc2,
                 theta75=dt$Bj_p75_czkr26_cc2,
                 se25=dt$Bj_p25_czkr26_cc2_bs_se,
                 se75=dt$Bj_p75_czkr26_cc2_bs_se,
                 stayer25=dt$e_rank_b_kr26_p25,
                 stayer75=dt$e_rank_b_kr26_p75)

## There should be no missing data for large counties
testthat::expect_equal(sum(is.na(cz[cz$population>=25000, ])), 0L)
testthat::expect_equal(sum(!is.na(cz[cz$population<25000, 5:6])), 0L)

## Estimates should have population-weighted mean equal to zero
testthat::expect_lt(abs(weighted.mean(cz$theta25,
                                      w=cz$pop, na.rm=TRUE)), 1e-8)
testthat::expect_lt(abs(weighted.mean(cz$theta75,
                                      w=cz$pop, na.rm=TRUE)), 1e-8)

## 2. County dataset, 3138 obs by 611 vars
dt <- readstata13::read.dta13("online_table_4.dta")

ct <- data.frame(ct=dt$cty,
                 cz=dt$cz,
                 ctname=dt$county_name,
                 state=dt$stateabbrv,
                 pop=dt$cty_pop2000,
                 pop_cz=dt$cz_pop2000,
                 theta25=dt$Bj_p25_cz_cty_kr26_cc2,
                 theta75=dt$Bj_p75_cz_cty_kr26_cc2,
                 se25=dt$Bj_p25_cz_cty_kr26_cc2_bs_se,
                 se75=dt$Bj_p75_cz_cty_kr26_cc2_bs_se,
                 stayer25=dt$e_rank_b_kr26_p25,
                 stayer75=dt$e_rank_b_kr26_p75,
                 se_cc25=dt$Bj_p25_cty_kr26_cc2_se,
                 se_cc75=dt$Bj_p75_cty_kr26_cc2_se,
                 theta_cc25=dt$Bj_p25_cty_kr26_cc2,
                 theta_cc75=dt$Bj_p75_cty_kr26_cc2)
## Last four columns are counties within CZ, we only use these for sanity checks

## In the 747 counties with dc$pop < 10000 | dc$pop_cz < 25000, then we don't
## have theta_ct. Focus on remaining 2391 counties
testthat::expect_equal(0L,
                       sum(!is.na(ct[ct$pop<10000 | ct$pop_cz<25000, 7:8])))
## Missing data: 33 counties. 9 don't have stayers outcomes. 12 are in a CZ with
## only 1 county, so no county-within-CZ esitmate (which is fine). 12 don't have
## such estimate for unknown reason, probably data limitations.

## Check county-within-CZ means are zero
for (cz0 in unique(ct$cz)) {
    idx <- ct$cz %in% cz0
    avg25 <- weighted.mean(ct$theta_cc25[idx], w=ct$pop[idx])
    avg75 <- weighted.mean(ct$theta_cc75[idx], w=ct$pop[idx])
    testthat::expect_true(is.na(avg25) | avg25 < 1e-7)
    testthat::expect_true(is.na(avg75) | avg75 < 1e-6)
}

## Check standard errors for counties line up
df2 <- merge(ct, cz[, c("cz", "theta25", "theta75", "se25", "se75")], by="cz")
names(df2)[c(7:10, 17:20)] <- c("theta_ct25", "theta_ct75", "se_ct25",
                                      "se_ct75", "theta_cz25", "theta_cz75",
                                      "se_cz25", "se_cz75")
df2$se_t25 <- sqrt(df2$se_cz25^2+ df2$se_cc25^2)
df2$se_t75 <- sqrt(df2$se_cz75^2+ df2$se_cc75^2)
## for 12 cases with no cc, use cz standard error
idx <- is.na(df2$se_cc25) & !is.na(df2$theta_ct25)
df2$se_t25[idx] <- df2$se_cz25[idx]
df2$se_t75[idx] <- df2$se_cz75[idx]

testthat::expect_lt(max(abs(df2$se_t25-df2$se_ct25), na.rm=TRUE), 1e-6)
testthat::expect_lt(max(abs(df2$se_t75-df2$se_ct75), na.rm=TRUE), 2e-6)
testthat::expect_equal(is.na(df2$se_t25), is.na(df2$se_ct25))
testthat::expect_equal(is.na(df2$se_t75), is.na(df2$se_ct75))

## Check theta_cc+theta_cz=theta_ct
th25 <- df2$theta_cz25+df2$theta_cc25 - df2$theta_ct25
th75 <- df2$theta_cz75+df2$theta_cc75 - df2$theta_ct75
testthat::expect_true(max(th25, na.rm=TRUE)<min(th25, na.rm=TRUE)+1e-6)
testthat::expect_true(max(th75, na.rm=TRUE)<min(th75, na.rm=TRUE)+5e-6)


usethis::use_data(cz, overwrite=TRUE, internal=FALSE)
ct <- ct[, 1:12]
usethis::use_data(ct, overwrite=TRUE, internal=FALSE)
