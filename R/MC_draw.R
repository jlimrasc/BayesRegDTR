
# D <- Data
# tau <- 0.01
# At_mags <- 1
# B <- 10000
# alph <- 3
# bet <- 4
compute_MC_draws <- function(D, tau, At_mags, B, alph, bet, p_list) {
    library(mvtnorm)
    draw_thetat_B <- function(ct, mt, omegat, B, alph, bet, n) {
        t(rmvt(B, sigma = (ct + 2*bet) / (n + 2 * alph) * solve(omegat),
             df = n+2*alph,
             delta = mt,
             type = "shifted")) # Draw all B at once
    }

    draw_sigmat_2B <- function(thetat_b, Zt, Xt, tau, alph, bet, n) {
        draw_sigmat_2b_inner <- function(thetat_b, Zt, Xt, tau, alph, bet, n) {
            1 / rgamma(1, shape = alph + (n + ncol(Zt)) / 2,
                       rate = bet + 1/2 * (sum((Xt - Zt %*% thetat_b)^2) + tau * sum(thetat_b^2)))
        }

        # Run function for each b
        browser()
        return(apply(thetat_b, 2, draw_sigmat_2b_inner,
              Zt = Zt, Xt = Xt, tau = tau, alph = alph, bet = bet, n = n))
    }
    browser()

    X <- D[-1]
    T <- length(X)
    n <- nrow(X[[1]])

    thetat_b_list  <- vector(mode = "list", length = T)
    sigmat_2b_list <- vector(mode = "list", length = T)
    for (t in 2:T) {
        # Compute summary stats
        browser()
        Zt          <- compute_Zt(A, At_len, X, t, n, p_list)
        thetat_hat  <- compute_thetat_hat(Zt, X[[t]])
        omegat      <- compute_omegat(Zt, tau)
        ct          <- compute_ct(Zt, X[[t]], omegat, n)
        mt          <- compute_mt(Zt, thetat_hat, omegat)

        # Draw
        thetat_b <- draw_thetat_B(ct, mt, omegat, B, alph, bet, n)
        sigmat_2b <- draw_sigmat_2B(thetat_b, Zt, X[[t]], tau, alph, bet, n)

        # Store
        thetat_b_list[[t]]  <- thetat_b
        sigmat_2b_list[[t]] <- sigmat_2b
    }

    return(list(thetat_b_list, sigmat_2b_list))
}

# source("~/GitHub/BayesRegDTR/R/datagen_univariate.R")
# source("~/GitHub/BayesRegDTR/R/ModelFitting.R")
res <- compute_MC_draws(D = Data, tau = 0.01, At_mags = 1, B = 10000, alph = 3, bet = 4, p_list = rep(1, T))
