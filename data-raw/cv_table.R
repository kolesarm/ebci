store_cva <- function(B, kurt, alpha=0.05, check=TRUE) {
    df <- expand.grid(B=B, kurt=c(kurt, Inf))
    cvj <- function(j) cva(df$B[j], kurt=df$kurt[j], alpha, check)$cv
    cores <- max(parallel::detectCores()-1L, 1L)
    df$cv <- unlist(parallel::mclapply(seq_len(nrow(df)), cvj, mc.cores=cores))
    df
}

Bs <-  seq(0, 5, by=0.01)
kurt <- c(seq(1, 10, length.out=181), seq(10, 100, by=5))
kurt <- 1:4

cv_tbl5 <- data.frame(store_cva(Bs, kurt, alpha=0.05), alpha=0.05)
cv_tbl1 <- data.frame(store_cva(Bs, kurt, alpha=0.1), alpha=0.1)
cva_tbl <- rbind(cv_tbl5, cv_tbl1)

usethis::use_data(cva_tbl, overwrite=TRUE, internal=FALSE)
