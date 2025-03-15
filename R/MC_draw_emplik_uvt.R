# h_t_wrap <- function(X_i1t, a_i1t, thetat, At_len, p_list, t) {
#     Z_ti <- compute_Zt(a_i1t, At_len, X_i1t, t, 1, p_list)
#     return(t(Z_ti) %*% (tail(x_i1t, n = 1) - Z_ti %*% thetat))
# }
#
# h_y_wrap <- function(yi, X_i1T, a_i1t, thetaT, At_len, p_list, t) {
#     Z_T1i <- compute_Zt(a_i1t, At_len, X_i1T, t, 1, p_list)
#     return(t(Z_T1i) %*% (tail(x_i1t, n = 1) - Z_T1i %*% thetaT))
# }

h_t <- function(Zti, thetat) {
    print("Zti:")
    print(Zti)
    print("thetat:")
    print(thetat)
    x_it <- Zti[length(Zti) / 2]
    return(t(Zti) %*% (x_it - Zti %*% thetat))
}

h_ty <- function(ZT1i, thetat) {
    yi <- ZT1i[1]
    return(t(ZT1i) %*% (yi - ZT1i %*% thetat)) # Assume Beta is meant 2 be theta?
}

delth_h <- function(Zti, th) {
    return(-t(Zti) %*% Zti)
}

delbet_hy <- function(ZT1i, th) {
    return(-t(ZT1i) %*% ZT1i)
}

delth_logpi <- function(thetat) {
    return(-thetat / 100)
}

delbet_logpi <- function(beta) {
    ret-bet / 100
}

library(VBel)
X_mat <- matrix(unlist(X), nrow = n, ncol = T)

t <- T
tau <- 0.01

Zt          <- compute_Zt(A, At_len, X, t, n, p_list)
omegat      <- compute_omegat(Zt, tau)
omegat_inv  <- solve(omegat)

Mnt <- omegat_inv %*% t(Zt) %*% X[[t]]
qt <- ncol(Zt)

thetat_b_list <- vector(mode = "list", length = T)

# Step 2 for (t in 2:T)
res <- compute_GVA(mu0 = Mnt, C0 = diag(qt), h = h_t,
            delthh = delth_h, delth_logpi = delth_logpi, lam0 = rep(0, qt),
            rho = 0.9, epsil = 1e-6, a = 1e-3,
            z = cbind(X_mat[,1:t], A[,1:(t - 1)]))

thetat_b <- rnorm(B, res$mu_FC, res$C_FC %*% t(res$C_FC))

thetat_b_list[t] <- thetat_b

# Step 3
resy <- compute_GVA(mu0 = MnT1, C0 = diag(qT1), h = h_ty,
            delthh = delbet_hy, delth_logpi = delbet_logpi, lam0 = rep(0, qT1),
            rho = 0.9, epsil = 1e-6, a = 1e-3,
            z = cbind(y, X_mat, A))

# Step 4
bet_b <- rnorm(B, resy$mu_FC, resy$C_FC %*% t(resy$C_FC))

