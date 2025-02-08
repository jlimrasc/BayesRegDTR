library(gtools)
compute_Zt <- function(A, At_len, X, t, n, p_list) {
    st <- function(A, key, t) {
        return(
            if(is.vector(tempRes <- A[,1:(t-1)] == rep(key, each = n)))
                tempRes
            else
                apply(tempRes, 1, all)) # Makes sure all of row is true
    }
    
    p_sum   <- sum(p_list[1:(t-1)])
    Z_tilde <- matrix(unlist(X[1:(t-1)]), nrow = n, ncol = p_sum)

    # Calculate permutations of a1, ..., an
    perms <- permutations(At_len, t-1, repeats.allowed = TRUE)
    
    # Preallocate Zt
    Zt <- matrix(0, nrow = n, ncol = sum(p_list[1:(t-1)]) * nrow(perms))

    # Compute Zt
    for (i in 1:nrow(perms)) {
        colnum <- (i-1) * p_sum + 1
        Zt[, colnum:(colnum + p_sum - 1)] <- st(A, perms[i,], t) * Z_tilde
    }

    # Check size
    if ((ans <- ncol(Zt)) != (proper <- p_sum * At_len^(t-1)))
        stop("Zt does not have the right amount of columns: ", ans, " != ", proper)
    
    return(Zt)
}

# Summary Stats
compute_thetat_hat <- function(Zt, Xt) {
    solve(t(Zt) %*% Zt) %*% t(Zt) %*% Xt
}

compute_omegat <- function(Zt, tau) {
    t(Zt) %*% Zt + diag(tau, nrow = ncol(Zt))
}

compute_ct <- function(Zt, Xt, omegat, n) {
    drop(t(Xt) %*% (diag(n) - Zt %*% solve(omegat) %*% t(Zt)) %*% Xt)
}

compute_mt <- function(Zt, thetat_hat, omegat) {
    solve(omegat) %*% t(Zt) %*% Zt %*% thetat_hat
}

t <- 3
tau <- 0.01
Zt          <- compute_Zt(A, At_len, X, t, n, p_list)
thetat_hat  <- compute_thetat_hat(Zt, X[[t]])
omegat      <- compute_omegat(Zt, tau)
ct          <- compute_ct(Zt, X[[t]], omegat, n)
mt          <- compute_mt(Zt, thetat_hat, omegat)
