compute_Zt <- function(A, At_lens, X, t, n, p_list) {
    st <- function(A, key, t) {
        return(
            if(is.vector(tempRes <- A[, 1:(t - 1), drop = FALSE] == rep(key, each = n)))
                tempRes
            else
                apply(tempRes, 1, all)) # Makes sure all of row is true
    }

    p_sum   <- sum(p_list[1:(t-1)])

    # Fix X format
    if (is.list(X)) {
        Z_tilde <- matrix(unlist(X[1:(t-1)]), nrow = n, ncol = p_sum)
    } else if (is.matrix(X) &&
               NROW(X) >= n &&
               NCOL(X) >= sum(p_list[1:(t - 1)])) {
        Z_tilde <- X[1:n, 1:sum(p_list[1:(t - 1)])]
    } else {
        stop("X must be a list or matrix")
    }


    # Calculate permutations of a1, ..., an
    perms <- vec_permutations(At_lens[1:(t-1)])

    # Preallocate Zt
    Zt <- matrix(0, nrow = n, ncol = p_sum * nrow(perms))

    # Compute Zt
    for (i in 1:nrow(perms)) {
        colnum <- (i-1) * p_sum + 1
        Zt[, colnum:(colnum + p_sum - 1)] <- st(A, perms[i,], t) * Z_tilde
    }

    # Check size
    if ((ans <- ncol(Zt)) != (proper <- p_sum * prod(At_lens[1:(t-1)])))
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

compute_ct <- function(Zt, Xt, omegat_inv, n) {
    drop(t(Xt) %*% (diag(n) - Zt %*% omegat_inv %*% t(Zt)) %*% Xt)
}

compute_mt <- function(Zt, thetat_hat, omegat_inv) {
    omegat_inv %*% t(Zt) %*% Zt %*% thetat_hat
}

vec_permutations <- function(max_vals) {
    # Create a list of sequences for each position
    ranges <- lapply(max_vals, function(x) 1:x)

    # Use expand.grid to generate all combinations
    perms <- expand.grid(rev(ranges))

    # Reverse columns back to original order
    perms <- perms[, rev(seq_along(perms))]

    # Convert to matrix or leave as data frame
    perms <- as.matrix(perms)
    dimnames(perms) <- NULL

    return(perms)
}
# t <- 3
# tau <- 0.01
# Zt          <- compute_Zt(A, At_lens, X, t, n, p_list)
# thetat_hat  <- compute_thetat_hat(Zt, X[[t]])
# omegat      <- compute_omegat(Zt, tau)
# ct          <- compute_ct(Zt, X[[t]], omegat, n)
# mt          <- compute_mt(Zt, thetat_hat, solve(omegat))
