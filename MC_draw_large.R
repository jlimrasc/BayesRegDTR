library(mvtnorm)
library(MixMatrix)
B <- 10000
nu0 <- 3
V0 <- 4
# compute_MC_draws <- function(D, tau, At_mags, B, nu0, V0) {
draw_sigmat_b <- function(Zt, Xt, Mnt, nu0, V0, tau, n) {
    Vn <- V0 + t(Xt - Zt %*% Mnt) %*% (Xt - Zt %*%Mnt) + tau * t(Mnt) %*% Mnt
    return(1 / rWishart(B, df = n + nu0, Sigma = Vn))
}

draw_Wt_b <- function(Zt, omegat, Mnt, tau, sigmat_b) {
    # apply(rmatrixnorm(1, mean = Mnt, U = solve(omegat), sigmat_b))
}

# thetat_b_list <- matrix(0, nrow = B, ncol = (t-1) * p_t *At_len^(t-1))
sigmat_b_list   <- vector(mode = "list", length = T)
Wt_b_list       <- vector(mode = "list", length = T)
for (t in 2:T) {
    Zt          <- compute_Zt(A, At_len, X, t, n)
    omegat      <- compute_omegat(Zt, tau)

    Mnt <- solve(omegat) %*% t(Zt) %*% X[,,t]
    
    sigmat_b <- draw_sigmat_b(Zt, X[,,t], Mnt, nu0, V0, tau, n)
    browser()
    Wt_b <- draw_Wt_b(Zt, omegat, Mnt, tau, sigmat_b)
    

    sigmat_b_list[[t]]  <- sigmat_b
    Wt_b_list[[t]] <- Wt_b
}
    
#     return(list(sigmat_b_list, Wt_b_list))
# }
# 
# res <- compute_MC_draws(D = 1, tau = 0.01, At_mags = 1, B = 10000, nu0 = 3, V0 = 4)
