library(mvtnorm)
library(MixMatrix)
library(tictoc)
compute_MC_draws <- function(D, tau, At_mags, B, nu0, V0) {
    draw_sigmat_b <- function(Zt, Xt, Mnt, nu0, V0, tau, n, B, t) {
        Vn <- V0[[t]] + t(Xt - Zt %*% Mnt) %*% (Xt - Zt %*%Mnt) + tau * t(Mnt) %*% Mnt
        test <- rWishart(B, df = n + nu0, Sigma = solve(Vn))
        return(apply(test, 3, solve, simplify = FALSE))
    }
    
    draw_Wt_b <- function(omegat, Mnt, tau, sigmat_b) {
        lapply(sigmat_b, FUN = function(z) rmatrixnorm(1, mean = Mnt, U = solve(omegat), z))
    }
    
    X <- Data[-1]
    
    sigmat_b_list   <- vector(mode = "list", length = T)
    Wt_b_list       <- vector(mode = "list", length = T)
    for (t in 2:T) {
        # Compute summary stats
        Zt          <- compute_Zt(A, At_len, X, t, n, p_list)
        omegat      <- compute_omegat(Zt, tau)
    
        Mnt <- solve(omegat) %*% t(Zt) %*% X[[t]]
        
        # Draw
        sigmat_b <- draw_sigmat_b(Zt, X[[t]], Mnt, nu0, V0, tau, n, B, t)
        
        tic(t)
        Wt_b <- draw_Wt_b(omegat, Mnt, tau, sigmat_b)
        toc()
        
    
        # Store
        sigmat_b_list[[t]]  <- sigmat_b
        Wt_b_list[[t]] <- Wt_b
    }
        
        return(list(sigmat_b_list, Wt_b_list))
}

# source("~/GitHub/BayesRegDTR/datagen_multivariate.R")
# source("~/GitHub/BayesRegDTR/ModelFitting.R")
res <- compute_MC_draws(D = Data, tau = 0.01, At_mags = rep(3, T), B = 10000, nu0 = 3, 
                        V0 = mapply(diag, p_list, SIMPLIFY = FALSE))
