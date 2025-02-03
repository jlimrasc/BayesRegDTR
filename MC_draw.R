library(mvtnorm)
compute_MC_draws <- function(D, tau, At_mags, B, alph, bet) {
    draw_thetat_b <- function(RSSt, ct, mt, omegat, B, alph, bet, n) {
        t(rmvt(B, sigma = (ct + 2*bet) / (n + 2 * alph) * solve(omegat), 
             df = n+2*alph, 
             delta = mt,
             type = "shifted"))
    }
    
    draw_sigmat_2b <- function(thetat_b, Zt, Xt, tau, alph, bet, n) {
        1 / rgamma(B, shape = alph + (n + ncol(Zt)) / 2, 
               rate = bet + 1/2 * (sum((Xt - Zt %*% thetat_b)^2) + tau * sum(thetat_b^2)))
    }

    thetat_b <- draw_thetat_b(RSSt, ct, mt, omegat, B, alph, bet, n)
    sigmat_2b <- draw_sigmat_2b(thetat_b, Zt, X[,,t], tau, alph, bet, n)

    # thetat_b_list <- matrix(0, nrow = B, ncol = (t-1) * p_t *At_len^(t-1))
    thetat_b_list  <- vector(mode = "list", length = T)
    sigmat_2b_list <- vector(mode = "list", length = T)
    for (t in 2:T) {
        Zt          <- compute_Zt(A, At_len, X, t, n)
        thetat_hat  <- compute_thetat_hat(Zt, X[,,t])
        RSSt        <- compute_RSSt(Zt, thetat_hat, X[,,t])
        omegat      <- compute_omegat(Zt, tau)
        ct          <- compute_ct(Zt, X[,,t], omegat, n)
        mt          <- compute_mt(Zt, thetat_hat, omegat)

        thetat_b <- draw_thetat_b(RSSt, ct, mt, omegat, B, alph, bet, n)
        sigmat_2b <- draw_sigmat_2b(thetat_b, Zt, X[,,t], tau, alph, bet, n)

        thetat_b_list[[t]]  <- thetat_b
        sigmat_2b_list[[t]] <- sigmat_2b
    }

    return(list(thetat_b_list, sigmat_2b_list))
}

res <- compute_MC_draws(D = 1, tau = 0.01, At_mags = 1, B = 10000, alph = 3, bet = 4)
