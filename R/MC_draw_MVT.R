#' Compute Monte Carlo Draws from Multivariate Dataset
#'
#' @description
#' Obtain Monte Carlo draws from posterior distribution of stagewise regression parameters
#'
#' @param Data      Observed data organised as a list of {y, X, A} where y is a vector of the final outcomes,
#' X is a list of matrices of the intermediate covariates
#' and A is a matrix of the assigned treatments
#' @param tau       Prior precision scale. Should be specified with a small value
#' @param At_lens   Vector of number of treatment options at each stage
#' @param B         Number of MC draws
#' @param nu0       Inverse-Wishart degres of freedom
#' @param V0        Inverse-Wishart scale matrix
#' @param p_list    Vector of dimension for each stage
#'
#' @returns Monte Carlo draws??? A list containing:  \enumerate{
#'              \item sigmat_B_list: Desc. A list of length T with each element a vector of size B x p_t
#'              \item Wt_B_list: Desc. A list of length T with each element a matrix of size B x p_t
#'              \item beta_B: Desc. A list of length B
#'              \item sigmay_2B: Desc. A list of length B
#'              }
#' @useDynLib BayesRegDTR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
#'
#' @examples
#' # -----------------------------
#' # Initialise Inputs
#' # -----------------------------
#' n           <- 5000
#' num_treats  <- 3
#' p_list      <- rep(2, num_treats)
#' At_lens     <- rep(3, num_treats)
#' Data        <- generate_dataset_mvt(n, num_treats, p_list, At_lens)
#'
#' # -----------------------------
#' # Main
#' # -----------------------------
#' res_mvt <- compute_MC_draws_mvt(Data = Data, tau = 0.01, At_lens = 3, B = 100, nu0 = 3, V0 = diag(2), alph = 3, gam = 4, p_list = 2)
compute_MC_draws_mvt <- function(Data, tau, At_lens, B, nu0,
                                 V0 = mapply(diag, p_list, SIMPLIFY = FALSE),
                                 alph, gam, p_list ) {
    library(mvtnorm)
    library(MixMatrix)
    library(tictoc)
    library(expm)
    draw_sigmat_B <- function(Zt, Xt, Mnt, nu0, V0, tau, n, B, t) {
        Vn <- V0[[t]] + t(Xt - Zt %*% Mnt) %*% (Xt - Zt %*%Mnt) + tau * t(Mnt) %*% Mnt
        temp_rwish <- rWishart(B, df = n + nu0, Sigma = solve(Vn))

        return(apply(temp_rwish, 3, solve, simplify = FALSE))
    }

    draw_thetat_B <- function(ct, mt, omegat_inv, B, alph, gam, n) {
        t(rmvt(B, sigma = (ct + 2*gam) / (n + 2 * alph) * omegat_inv,
               df = n+2*alph,
               delta = mt,
               type = "shifted")) # Draw all B at once
    }

    draw_sigmat_2B <- function(thetat_b, Zt, Xt, tau, alph, gam, n) {
        draw_sigmat_2b_inner <- function(thetat_b, Zt, Xt, tau, alph, gam, n) {
            1 / rgamma(1, shape = alph + (n + ncol(Zt)) / 2,
                       rate = gam + 1/2 * (sum((Xt - Zt %*% thetat_b)^2) + tau * sum(thetat_b^2)))
        }

        # Run function for each b
        return(apply(thetat_b, 2, draw_sigmat_2b_inner,
                     Zt = Zt, Xt = Xt, tau = tau, alph = alph, gam = gam, n = n))
    }

    draw_beta_B <- function(cT1, mT1, omegaT1_inv, B, alph, gam, n) {
        return(draw_thetat_B(cT1, mT1, omegaT1_inv, B, alph, gam, n))
    }

    draw_sigmay_2B <- function(beta_b, ZT1, y, tau, alph, gam, n) {
        return(draw_sigmat_2B(beta_b, ZT1, y, tau, alph, gam, n))
    }

    # Unpack data
    T <- length(Data) - 2
    X <- Data[2:(T+1)]
    y <- Data[[1]]
    A <- Data[[T+2]]
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

    sigmat_B_list   <- vector(mode = "list", length = T)
    Wt_B_list       <- vector(mode = "list", length = T)
    for (t in 2:T) {
        # Compute summary stats
        Zt          <- compute_Zt(A, At_lens, X, t, n, p_list)
        omegat      <- compute_omegat(Zt, tau)
        omegat_inv  <- solve(omegat)

        Mnt <- omegat_inv %*% t(Zt) %*% X[[t]]

        # Draw
        tic(paste('t =', t, 'sigmat_B: '))
        sigmat_B <- draw_sigmat_B(Zt, X[[t]], Mnt, nu0, V0, tau, n, B, t)
        toc(log = TRUE)

        tic(paste('t =', t, 'Wtb: '))
        Wt_B <- draw_Wt_B_cpp(omegat_inv, Mnt, tau, sigmat_B, p_list[t], ncol(Zt), B)
        toc(log = TRUE)


        # Store
        sigmat_B_list[[t]]  <- sigmat_B
        Wt_B_list[[t]] <- Wt_B
    }

    tic(paste(T+1, "summary"))
    ZT1         <- compute_Zt(A, At_lens, X, T+1, n, p_list)
    thetaT1_hat <- compute_thetat_hat(ZT1, y)
    omegaT1     <- compute_omegat(ZT1, tau)
    omegaT1_inv <- solve(omegaT1)
    cT1         <- compute_ct(ZT1, y, omegaT1_inv, n)
    mT1         <- compute_mt(ZT1, thetaT1_hat, omegaT1_inv)
    toc()

    tic('thetat_b')
    beta_B <- draw_beta_B(cT1, mT1, omegaT1_inv, B, alph, gam, n)
    toc()
    tic('sigmat_2b')
    sigmay_2B <- draw_sigmay_2B(beta_B, ZT1, y, tau, alph, gam, n)
    toc()

    return(list(sigmat_B_list = sigmat_B_list[-1], Wt_B_list = Wt_B_list[-1],
                beta_B = beta_B, sigmay_2B = sigmay_2B))
}
