store_cva <- function(m2, kappa, alpha, check=TRUE) {
    df <- expand.grid(m2=m2, kappa=kappa, alpha=alpha)
    cvj <- function(j) cva(df$m2[j], kappa=df$kappa[j], df$alpha[j], check)$cv
    cores <- max(parallel::detectCores()-1L, 1L)
    df$cv <- unlist(parallel::mclapply(seq_len(nrow(df)), cvj, mc.cores=cores))
    df
}

m2 <-  seq(0, 5, by=0.001)^2
## Takes 4 hrs on X1 Carbon
kappa <- c(1:10, Inf)
## Takes 2 hrs on X1 Carbon
## kappa <- c(1:4, Inf)

cv_tbl5 <- data.frame(store_cva(m2, kappa, alpha=0.05))
cv_tbl1 <- data.frame(store_cva(m2, kappa, alpha=0.1))
cva_tbl <- rbind(cv_tbl5, cv_tbl1)

usethis::use_data(cva_tbl, overwrite=TRUE, internal=FALSE)
