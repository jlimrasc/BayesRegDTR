#' Title
#'
#' @param D
#' @param tau
#' @param At_lens
#' @param B
#' @param nu0
#' @param V0
#' @param p_list
#'
#' @returns
#' @useDynLib BayesRegDTR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
#'
#' @examples
compute_MC_draws_mvt <- function(D, tau, At_lens, B, nu0,
                                 V0 = mapply(diag, p_list, SIMPLIFY = FALSE),
                                 p_list) {
    library(mvtnorm)
    library(MixMatrix)
    library(tictoc)
    library(expm)
    draw_sigmat_b <- function(Zt, Xt, Mnt, nu0, V0, tau, n, B, t) {
        Vn <- V0[[t]] + t(Xt - Zt %*% Mnt) %*% (Xt - Zt %*%Mnt) + tau * t(Mnt) %*% Mnt
        temp_rwish <- rWishart(B, df = n + nu0, Sigma = solve(Vn))

        return(apply(temp_rwish, 3, solve, simplify = FALSE))
    }

    # Unpack data
    T <- length(D) - 2
    X <- D[2:(T+1)]
    A <- D[[T+2]]
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
        tic(cat('\n', 't = ', t, ' sigmat_b: ', sep = ''))
        sigmat_b <- draw_sigmat_b(Zt, X[[t]], Mnt, nu0, V0, tau, n, B, t)
        toc(log = TRUE)

        tic(cat('\n', 't = ', t, ' Wtb: ', sep = ''))
        Wt_b <- draw_Wt_b_cpp(omegat_inv, Mnt, tau, sigmat_b, p_list[t], ncol(Zt), B)
        toc(log = TRUE)


        # Store
        sigmat_b_list[[t]]  <- sigmat_b
        # Wt_b_list[[t]] <- Wt_b
    }

    return(list(sigmat_b_list = sigmat_b_list, Wt_b_list = Wt_b_list))
}

# res3 <- compute_MC_draws_mvt(D = Data, tau = 0.01, At_lens = 3, B = 10000, nu0 = 3,
#                              V0 = diag(2), p_list = 2)
