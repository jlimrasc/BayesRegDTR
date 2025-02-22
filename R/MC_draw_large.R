library(mvtnorm)
library(MixMatrix)
library(tictoc)
library(expm)
compute_MC_draws_mvt <- function(D, tau, At_lens, B, nu0,
                                 V0 = mapply(diag, p_list, SIMPLIFY = FALSE),
                                 p_list) {
    draw_sigmat_b <- function(Zt, Xt, Mnt, nu0, V0, tau, n, B, t) {
        Vn <- V0[[t]] + t(Xt - Zt %*% Mnt) %*% (Xt - Zt %*%Mnt) + tau * t(Mnt) %*% Mnt
        temp_rwish <- rWishart(B, df = n + nu0, Sigma = solve(Vn))
        return(apply(temp_rwish, 3, solve, simplify = FALSE))
    }

    draw_Wt_b <- function(omegat_inv, Mnt, tau, sigmat_b, t, pt, qt) {
        # Matrix norm
        # lapply(sigmat_b, FUN = function(z) rmatrixnorm(1, mean = Mnt, U = solve(omegat), V = z))

        # Use sqrtm
        # Eigen version:
        # ev <- eigen(solve(omegat))
        # sqrtm_omegat <- ev$vectors %*% diag(sqrt(ev$values)) %*% solve(ev$vectors)

        # Sqrtm version:
        R_list <- array(rnorm(B * pt * qt), dim = c(qt, pt, B))

        omegat_inv_sqrtm <- sqrtm(omegat_inv)
        omegR_list <- apply(R_list, 3, FUN = function(R) omegat_inv_sqrtm %*% R, simplify = FALSE)
        return(mapply(FUN = function(omegR, sigmatb) Mnt + omegR %*% sqrtm(sigmatb), omegR_list, sigmat_b, SIMPLIFY = FALSE))
    }

    # Unpack data
    X <- D[2:(T+1)]
    A <- D[[T+2]]
    T <- length(X)
    n <- nrow(X[[1]])

    # Input validation
    if (length(At_lens) == 1) At_lens <- rep(At_lens, T)
    if (length(p_list) == 1) {
        if ((is.list(V0) && length(V0) == 1) || is.atomic(V0))
            V0 <- replicate(T, matrix(unlist(V0), ncol = p_list), simplify = FALSE)
        p_list <- rep(p_list, T)
    }

    if (!all(c(length(At_lens), length(V0), length(p_list)) == T))
        stop("Length of At_lens, V0 and p_list must be equal and of length 1 or T")

    if (tau < 0) stop("Value of tau must be positive")

    sigmat_b_list   <- vector(mode = "list", length = T)
    Wt_b_list       <- vector(mode = "list", length = T)
    for (t in 2:T) {
        # Compute summary stats
        Zt          <- compute_Zt(A, At_lens[t], X, t, n, p_list)
        omegat      <- compute_omegat(Zt, tau)
        omegat_inv  <- solve(omegat)

        Mnt <- omegat_inv %*% t(Zt) %*% X[[t]]

        # Draw
        sigmat_b <- draw_sigmat_b(Zt, X[[t]], Mnt, nu0, V0, tau, n, B, t)

        tic(t)
        Wt_b <- draw_Wt_b(omegat_inv, Mnt, tau, sigmat_b, t, p_list[t], ncol(Zt))
        toc()


        # Store
        sigmat_b_list[[t]]  <- sigmat_b
        Wt_b_list[[t]] <- Wt_b
    }

        return(list(sigmat_b_list = sigmat_b_list, Wt_b_list = Wt_b_list))
}

# source("~/GitHub/BayesRegDTR/R/datagen_multivariate.R")
# source("~/GitHub/BayesRegDTR/R/ModelFitting.R")
res6 <- compute_MC_draws_mvt(D = Data, tau = 0.01, At_lens = 3, B = 10000, nu0 = 3,
                        V0 = diag(2), p_list = 2)
