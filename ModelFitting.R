source("~/GitHub/DTR/datagen.R")
library(gtools)
compute_Zt <- function(A, At_max, X, t) {
    st <- function(A, key, t) {
        return(apply(A[,1:t] == rep(key[1:t], each = n), 1, all))
    }
    
    Z_tilde <- X[,1:t]
    
    
    perms <- permutations(At_max, t, repeats.allowed = TRUE)
    Z_list <- vector(mode = "list", length = nrow(perms))
    for (i in 1:nrow(perms)) {
        Z_list[[i]] <- st(A, perms[i,], t) * Z_tilde
    }
    names(Z_list) <- apply(perms, 1, paste, collapse = "") # Name each element as the list of `a` values
    # Each name is actually modified ternary. To calculate index, -1 from all 
    # bits other than LSB, then convert from ternary to decimal
    return(Z_list)
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

t <- 2
tau <- 5
Zt_list     <- compute_Zt(A, At_max, X, t)
thetat_hat  <- lapply(X = Zt_list, FUN = compute_thetat_hat, Xt = X[,t])
RSSt        <- mapply(compute_RSSt, Zt = Zt_list, thetat_hat = thetat_hat, MoreArgs = list(Xt = X[,t]), SIMPLIFY = FALSE)
sigmat      <- lapply(Zt_list, compute_sigmat, tau = tau)
ct          <- mapply(compute_ct, Zt = Zt_list, thetat_hat = thetat_hat, sigmat = sigmat, SIMPLIFY = FALSE)
mt          <- mapply(compute_mt, Zt = Zt_list, thetat_hat = thetat_hat, sigmat = sigmat, SIMPLIFY = FALSE)