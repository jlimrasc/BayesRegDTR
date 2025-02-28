compute_MC_draws <- function(D, tau, At_lens, B, alph, gam, p_list) {
    library(mvtnorm)
    draw_thetat_B <- function(ct, mt, omegat_inv, B, alph, gam, n) {
        t(rmvt(B, sigma = (ct + 2*gam) / (n + 2 * alph) * omegat_inv,
             df = n+2*alph,
             delta = mt,
             type = "shifted")) # Draw all B at once
    }

    draw_sigmat_2B <- function(thetat_B, Zt, Xt, tau, alph, gam, n) {
        draw_sigmat_2B_inner <- function(thetat_B, Zt, Xt, tau, alph, gam, n) {
            1 / rgamma(1, shape = alph + (n + ncol(Zt)) / 2,
                       rate = gam + 1/2 * (sum((Xt - Zt %*% thetat_B)^2) + tau * sum(thetat_B^2)))
        }

        # Run function for each b
        return(apply(thetat_B, 2, draw_sigmat_2B_inner,
              Zt = Zt, Xt = Xt, tau = tau, alph = alph, gam = gam, n = n))
    }

    draw_beta_B <- function(cT1, mT1, omegaT1_inv, B, alph, gam, n) {
        return(draw_thetat_B(cT1, mT1, omegaT1_inv, B, alph, gam, n))
    }

    draw_sigmay_2B <- function(beta_b, ZT1, y, tau, alph, gam, n) {
        return(draw_sigmat_2B(beta_b, ZT1, y, tau, alph, gam, n))
    }

    T <- length(D) - 2
    X <- D[2:(T+1)]
    y <- D[[1]]
    A <- D[[T+2]]
    n <- nrow(X[[1]])

    thetat_B_list  <- vector(mode = "list", length = T)
    sigmat_2B_list <- vector(mode = "list", length = T)
    for (t in 2:T) {
        tic(paste(t, "summary"))
        # Compute summary stats
        Zt          <- compute_Zt(A, At_len, X, t, n, p_list)
        thetat_hat  <- compute_thetat_hat(Zt, X[[t]])
        omegat      <- compute_omegat(Zt, tau)
        omegat_inv  <- solve(omegat)
        ct          <- compute_ct(Zt, X[[t]], omegat_inv, n)
        mt          <- compute_mt(Zt, thetat_hat, omegat_inv)
        toc()

        # Draw
        tic('thetat_B')
        thetat_B <- draw_thetat_B(ct, mt, omegat, B, alph, gam, n)
        toc()
        tic('sigmat_2B')
        sigmat_2B <- draw_sigmat_2B(thetat_B, Zt, X[[t]], tau, alph, gam, n)
        toc()

        # Store
        thetat_B_list[[t]]  <- thetat_B
        sigmat_2B_list[[t]] <- sigmat_2B
    }

    tic(paste(T+1, "summary"))
    ZT1         <- compute_Zt(A, At_len, X, T+1, n, p_list)
    thetaT1_hat <- compute_thetat_hat(ZT1, y)
    omegaT1     <- compute_omegat(ZT1, tau)
    omegaT1_inv <- solve(omegaT1)
    cT1         <- compute_ct(ZT1, y, omegaT1_inv, n)
    mT1         <- compute_mt(ZT1, thetaT1_hat, omegaT1_inv)
    toc()

    tic('thetat_B')
    beta_B <- draw_beta_B(cT1, mT1, omegaT1_inv, B, alph, gam, n)
    toc()
    tic('sigmat_2B')
    sigmay_2B <- draw_sigmay_2B(beta_B, ZT1, y, tau, alph, gam, n)
    toc()

    return(list(thetat_B_list = thetat_B_list, sigmat_2B_list = sigmat_2B_list,
                beta_B = beta_B, sigmay_2B = sigmay_2B))
}

# source("~/GitHub/BayesRegDTR/R/datagen_univariate.R")
# source("~/GitHub/BayesRegDTR/R/ModelFitting.R")
# D <- Data
# tau <- 0.01
# At_lens <- 1
# B <- 10000
# alph <- 3
# gam <- 4
# res <- compute_MC_draws(D = Data, tau = 0.01, At_lens = 1, B = 10000, alph = 3, gam = 4, p_list = rep(1, T))
