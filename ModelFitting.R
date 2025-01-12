source("~/GitHub/BayesRegDTR/datagen.R")
library(gtools)
compute_Zt <- function(A, At_max, X, t) {
    st <- function(A, key, t) {
        return(apply(A[,1:t-1] == rep(key, each = n), 1, all))
    }
    
    Z_tilde <- X[,1:t-1]
    
    perms <- permutations(At_max, t-1, repeats.allowed = TRUE)
    Zt <- matrix(0, nrow = n, ncol = (t-1) * nrow(perms))
    for (i in 1:nrow(perms)) {
        colnum <- (i-1) * (t-1) + 1
        Zt[, colnum:(colnum + t - 2)] <- st(A, perms[i,], t) * Z_tilde
    }
    if ((ans <- ncol(Zt)) != (proper <- (t-1) * At_max^(t-1)))
        stop("Zt does not have the right amount of columns: ", ans, " != ", proper)
    
    return(Zt)
}

# Summary Stats
compute_thetat_hat <- function(Zt, Xt) {
    solve(t(Zt) %*% Zt) %*% t(Zt) %*% Xt
}

compute_RSSt <- function(Zt, thetat_hat, Xt) {
    sum((Xt - Zt %*% thetat_hat) ^ 2)
}

compute_sigmat <- function(Zt, tau) {
    t(Zt) %*% Zt + diag(tau, nrow = ncol(Zt))
}

compute_ct <- function(Zt, thetat_hat, sigmat) {
    drop(t(thetat_hat) %*% t(Zt) %*% Zt %*% solve(sigmat) %*% t(Zt) %*% Zt %*% thetat_hat -
        t(thetat_hat) %*% t(Zt) %*% Zt %*% thetat_hat)
}

compute_mt <- function(Zt, thetat_hat, sigmat) {
    solve(sigmat) %*% t(Zt) %*% Zt %*% thetat_hat
}

t <- 3
tau <- 10
Zt          <- compute_Zt(A, At_max, X, t)
thetat_hat  <- compute_thetat_hat(Zt, X[,t])
RSSt        <- compute_RSSt(Zt, thetat_hat, X[,t])
sigmat      <- compute_sigmat(Zt, tau)
ct          <- compute_ct(Zt, thetat_hat, sigmat)
mt          <- compute_mt(Zt, thetat_hat, sigmat)
