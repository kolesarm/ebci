## The original data, online_table_3.dta, has been
## downloaded from
## https://opportunityinsights.org/wp-content/uploads/2018/04/replicate-4.zip,
## subdirectory data (not online tables, which contains the source files for
## this .dta file)

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

usethis::use_data(cz, overwrite=TRUE, internal=FALSE)
