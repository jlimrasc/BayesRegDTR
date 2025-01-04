library(mvtnorm)
compute_MC_draws <- function(D, tau, At_mags, B, alph, bet) {
    draw_thetat_b <- function(RSSt, ct, mt, sigmat, B, alph, bet, n) {
        t(rmvt(B, sigma = (RSSt + ct + 2*bet) /
                 (n + 2 * alph) * solve(sigmat), 
             df = n+2*alph, 
             delta = mt))
    }
    
    draw_sigmat_2b <- function(thetat_b, Zt, Xt, tau, alph, bet, n) {
        1 / rgamma(B, shape = alph + (n + ncol(Zt)) / 2, 
               rate = bet + 1/2 * (sum((Xt - Zt %*% thetat_b)^2) + tau * thetat_b^2))
    }
    thetat_b <- mapply(draw_thetat_b, RSSt, ct, mt, sigmat,
                       MoreArgs = list(B = B, alph = alph, bet = bet, n = n),
                       SIMPLIFY = FALSE)
    sigmat_2b <- mapply(draw_sigmat_2b, thetat_b, Zt_list, 
                        MoreArgs = list(Xt = X[,t], tau = tau, alph = alph,
                                        bet = bet, n = n),
                        SIMPLIFY = FALSE)
    browser()
    thetat_b_list <- vector(mode = "list", length = T)
    sigmat_2b_list <- vector(mode = "list", length = T)
    for (t in 1:T) {
        
    }
}