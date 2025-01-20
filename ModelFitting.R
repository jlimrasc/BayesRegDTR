source("~/GitHub/BayesRegDTR/datagen.R")
library(gtools)
compute_Zt <- function(A, At_len, X, t, n) {
    st <- function(A, key, t) {
        return(
            if(is.vector(tempRes <- A[,1:(t-1)] == rep(key, each = n)))
                tempRes
            else
                apply(tempRes, 1, all))
    }
    
    Z_tilde <- matrix(X[,,1:(t-1)], nrow = n, ncol = (t-1) * p_t)

    # Calculate permutations of a1, ..., an
    perms <- permutations(At_len, t-1, repeats.allowed = TRUE)
    
    # Preallocate Zt
    Zt <- matrix(0, nrow = n, ncol = (t-1) * p_t * nrow(perms))

    # Compute Zt
    for (i in 1:nrow(perms)) {
        colnum <- (i-1) * (t-1) * p_t + 1
        Zt[, colnum:(colnum + (t-1) * p_t - 1)] <- st(A, perms[i,], t) * Z_tilde
    }

    # Check size
    if ((ans <- ncol(Zt)) != (proper <- (t-1) * p_t * At_len^(t-1)))
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

compute_omegat <- function(Zt, tau) {
    t(Zt) %*% Zt + diag(tau, nrow = ncol(Zt))
}

compute_ct <- function(Zt, thetat_hat, omegat) {
    drop(t(thetat_hat) %*% t(Zt) %*% Zt %*% solve(omegat) %*% t(Zt) %*% Zt %*% thetat_hat -
        t(thetat_hat) %*% t(Zt) %*% Zt %*% thetat_hat)
}

compute_mt <- function(Zt, thetat_hat, omegat) {
    solve(omegat) %*% t(Zt) %*% Zt %*% thetat_hat
}

t <- 3
tau <- 0.01
Zt          <- compute_Zt(A, At_len, X, t, n)
thetat_hat  <- compute_thetat_hat(Zt, X[,,t])
RSSt        <- compute_RSSt(Zt, thetat_hat, X[,,t])
omegat      <- compute_omegat(Zt, tau)
ct          <- compute_ct(Zt, thetat_hat, omegat)
mt          <- compute_mt(Zt, thetat_hat, omegat)
