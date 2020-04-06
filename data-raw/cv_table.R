store_cva <- function(B, kappa, alpha, check=TRUE) {
    df <- expand.grid(B=B, kappa=kappa, alpha=alpha)
    cvj <- function(j) cva(df$B[j], kappa=df$kappa[j], df$alpha[j], check)$cv
    cores <- max(parallel::detectCores()-1L, 1L)
    df$cv <- unlist(parallel::mclapply(seq_len(nrow(df)), cvj, mc.cores=cores))
    df
}

Bs <-  seq(0, 5, by=0.001)
## Takes 3.5hrs on X1 Carbon
kappa <- c(1:10, Inf)

cv_tbl5 <- data.frame(store_cva(Bs, kappa, alpha=0.05))
cv_tbl1 <- data.frame(store_cva(Bs, kappa, alpha=0.1))
cva_tbl <- rbind(cv_tbl5, cv_tbl1)

usethis::use_data(cva_tbl, overwrite=TRUE, internal=FALSE)


## B=2.8, alpha=0.1, kappa=8
## cva(B=2.8, alpha=0.1, kappa=8)
## workds for kappa=7 or 9
## rho2(2.8^2, 8, 7.1324290147839279896)

if (FALSE) {
    r1 <- rho2(m2, kappa=10, chi, check=FALSE)
    ## This should give the same rejection rate
    sum(r(r1$x, chi)*r1$p)
}
