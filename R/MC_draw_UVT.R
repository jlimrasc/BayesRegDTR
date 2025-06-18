#' Compute Monte Carlo Draws from Univariate Dataset
#'
#' @description
#' Obtain Monte Carlo draws from posterior distribution of stagewise regression parameters
#'
#' @param Data          Observed data organised as a list of \eqn{\{y, X, A\}} where y is a vector of the final outcomes,
#' X is a list of matrices of the intermediate covariates and A is a matrix of the assigned treatments
#' @param tau           Prior precision scale. Should be specified with a small value
#' @param num_treats    Vector of number of treatment options at each stage
#' @param B             Number of MC draws
#' @param alph          Inverse-Gamma prior shape parameter for regression error variance of y. default:  1
#' @param gam           Inverse-Gamma prior rate parameter for regression error variance of y. default:  1
#' @param p_list        Vector of dimension for each stage
#' @param showBar       Whether to show a progress bar. Uses bar from \link[progress]{progress_bar} deafult: TRUE
#'
#' @returns Monte Carlo draws??? A list containing:  \enumerate{
#'              \item thetat_B_list: Desc. A list of length num_stages with each element a vector of length B
#'              \item sigmat_2B_list: Desc. A list of length num_stages with each element a vector of length B
#'              \item beta_B: Desc. A list of length B
#'              \item sigmay_2B: Desc. A list of length B
#'              }
#' @keywords internal
compute_MC_draws_uvt <- function(Data, tau, num_treats, B, alph, gam, p_list, showBar = TRUE) {
    draw_thetat_B <- function(ct, mt, omegat_inv, B, alph, gam, n) {
        t(mvtnorm::rmvt(B, sigma = (ct + 2*gam) / (n + 2 * alph) * omegat_inv,
             df = n+2*alph,
             delta = mt,
             type = "shifted")) # Draw all B at once
    }

    draw_sigmat_2B <- function(thetat_B, Zt, Xt, tau, alph, gam, n) {
        draw_sigmat_2B_inner <- function(thetat_B, Zt, Xt, tau, alph, gam, n) {
            1 / stats::rgamma(1, shape = alph + (n + ncol(Zt)) / 2,
                              rate = gam + 1/2 * (sum((Xt - Zt %*% thetat_B)^2) +
                                                      tau * sum(thetat_B^2)))
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

    num_stages <- length(Data) - 2
    X <- Data[2:(num_stages+1)]
    y <- Data[[1]]
    A <- Data[[num_stages+2]]
    n <- nrow(X[[1]])

    # Create bar
    progressr::with_progress({
        p <- progressr::progressor(steps = num_stages, enable = showBar,
                                   message = "Computing MC Draws (t=2)")

        thetat_B_list  <- vector(mode = "list", length = num_stages)
        sigmat_2B_list <- vector(mode = "list", length = num_stages)
        for (t in 2:num_stages) {
            # Compute summary stats
            Zt          <- compute_Zt(A, num_treats, X, t, n, p_list)
            thetat_hat  <- compute_thetat_hat(Zt, X[[t]])
            omegat      <- compute_omegat(Zt, tau)
            omegat_inv  <- solve(omegat)
            ct          <- compute_ct(Zt, X[[t]], omegat_inv, n)
            mt          <- compute_mt(Zt, thetat_hat, omegat_inv)

            # Draw
            thetat_B <- draw_thetat_B(ct, mt, omegat_inv, B, alph, gam, n)
            sigmat_2B <- draw_sigmat_2B(thetat_B, Zt, X[[t]], tau, alph, gam, n)

            # Store
            thetat_B_list[[t]]  <- thetat_B
            sigmat_2B_list[[t]] <- sigmat_2B

            p(message = sprintf("Computing MC Draws (t=%d)", t + 1))
        }

        ZT1         <- compute_Zt(A, num_treats, X, num_stages+1, n, p_list)
        thetaT1_hat <- compute_thetat_hat(ZT1, y)
        omegaT1     <- compute_omegat(ZT1, tau)
        omegaT1_inv <- solve(omegaT1)
        cT1         <- compute_ct(ZT1, y, omegaT1_inv, n)
        mT1         <- compute_mt(ZT1, thetaT1_hat, omegaT1_inv)

        beta_B <- draw_beta_B(cT1, mT1, omegaT1_inv, B, alph, gam, n)
        sigmay_2B <- draw_sigmay_2B(beta_B, ZT1, y, tau, alph, gam, n)

        # Complete bar
        p(message = sprintf("Computing MC Draws (t=%d)", num_stages + 1))
    })

    return(list(thetat_B_list = thetat_B_list[-1], sigmat_2B_list = sigmat_2B_list[-1],
                beta_B = beta_B, sigmay_2B = sigmay_2B))
}
