#' Title
#'
#' @param Dat.train     Training data in format returned by `generate_dataset`:
#' organised as a list of \eqn{\{y, X_1, X_2..., X_{num\_states}, A\}} where y is a
#' vector of the final outcomes, \eqn{X_1, X_2..., X_{num\_states}} is a list of matrices
#' of the intermediate covariates and A is an \eqn{n \times num\_stages}{n x num_stages} matrix of the
#' assigned treatments
#' @param Dat.pred      Prediction data in format returned by `generate_dataset`:
#' organised as a list of \eqn{\{X_1, X_2..., X_t, A\}} where
#' \eqn{X_1, X_2..., X_t} is a list of matrices of the intermediate
#' covariates and A is an \eqn{n \times (t-1)}{n x (t-1)} matrix of the assigned treatments
#' @param n.train       Number of samples/individuals in the training data
#' @param n.pred        Number of samples/individuals in the prediction data
#' @param num_stages    Total number of stages
#' @param num_treats    Vector of number of treatment options at each stage
#' @param p_list        Vector of intermediate covariate dimensions for each stage
#' @param t             Prediction stage t, where \eqn{\leq}{<=} num_stages
#' @param R             Draw size from posterior distribution of intermediate covariates. default:  30
#' @param numCores      Number of cores in the system to use. default uses parallel::detectCores() - 1
#' @param tau           Normal prior scale parameter for regression coefficients. Should be specified with a small value. default:  0.01
#' @param B             Number of MC draws from posterior of regression parameters. default 10000
#' @param nu0           Inverse-Wishart prior degrees of freedom for regression error Vcov matrix. default: 3
#' @param V0            List of Inverse-Wishart prior scale matrix for regression error Vcov matrix. default: list of identity matrices
#' @param alph          Inverse-Gamma prior shape parameter for regression error variance of y. default:  1
#' @param gam           Inverse-Gamma prior rate parameter for regression error variance of y. default:  1
#'
#' @returns list of {GCV_results, freqs} where GCV_results is an array of size
        #' n.pred x p_list\[t\] x B, and freqs is a n.pred x p_list\[t\] matrix
        #' of the frequency of each treatment type at stage t
#' @export
#'
#' @examples
#' # -----------------------------
#' # Initialise Inputs
#' # -----------------------------
#' num_stages  <- 3
#' p_list      <- rep(2, num_stages)
#' num_treats  <- rep(2, num_stages)
#' n.train     <- 5000
#' n.pred      <- 100
#'
#' tau     <- 0.01
#' B       <- 100
#' alph    <- 3
#' gam     <- 4
#' nu0     <- 3
#' V0      <- mapply(diag, p_list, SIMPLIFY = FALSE)
#' R       <- 30
#' numCores<- parallel::detectCores()
#'
#' # -----------------------------
#' # Generate Dataset
#' # -----------------------------
#' Dat.train <- generate_dataset_mvt(n.train, num_stages, p_list, num_treats)
#' Dat.pred  <- generate_dataset_mvt(n.pred,  num_stages, p_list, num_treats)
#' Dat.pred  <- Dat.pred[-1]
#' Dat.pred[[num_stages+1]]  <- Dat.pred[[num_stages+1]][1:n.pred, 1:(t-1), drop = FALSE]
#'
#' # -----------------------------
#' # Main
#' # -----------------------------
#' gcv_res <- testParallelGCV(Dat.train, Dat.pred, n.train, n.pred, num_stages, num_treats, p_list, t)
testParallelGCV <- function(Dat.train, Dat.pred, n.train, n.pred, num_stages, num_treats,
                            p_list, t, R = 30, numCores,
                            tau = 0.01, B = 10000, nu0 = 3,
                            V0 = mapply(diag, p_list, SIMPLIFY = FALSE),
                            alph = 1, gam = 1
                            ) {

    # Verify inputs
    stopifnot("t must be less than or equal to num_stages" = t <= num_stages)

    # Retrieve training data and train model
    # Dat.train <- c(list(Dat[[1]][1:n.train]), lapply(Dat[-1], function(x) x[1:n.train,]))
    res_mc <- compute_MC_draws_mvt(Data = Dat.train, tau = tau, num_treats = num_treats, B = B,
                                    nu0 = nu0, V0 = V0, alph = alph, gam = gam, p_list = p_list)

    # Set up parallel processing
    if (missing(numCores)) numCores <- parallel::detectCores() - 1
    doParallel::registerDoParallel(numCores)  # use multicore, set to the number of our cores

    # Inner loop
    inner_b_GCV <- function(i, p_t, p_list) {
        res_GCV_1B <- matrix(0, nrow = B, ncol = p_t)
        for (b in 1:B) {
            # histDat <- c(lapply(Dat[2:t], function(x) x[i,,drop = FALSE]), list(Dat[[num_stages + 2]][i,1,drop = FALSE]))
            # currDat <- Dat[[t]][i,,drop = FALSE]
            histDat <- c(lapply(Dat.pred[1:(t-1)], function(x) x[i,,drop = FALSE]), list(Dat.pred[[num_stages + 1]][i,1:(t-1),drop = FALSE]))
            currDat <- Dat.pred[[t]][i,,drop = FALSE]
            Wt      <- lapply(res_mc$Wt_B_list, function(x) x[[b]])
            Sigmat  <- lapply(res_mc$sigmat_B_list, function(x) x[[b]])

            res_GCV_1B[b,] <-
                GiveChoiceValue(Wt = Wt, Sigmat = Sigmat, bet = res_mc$beta_B[,b],
                                sigmay = res_mc$sigmay_2B[b], t = t,
                                num_stages = num_stages, p_list = p_list,
                                histDat = histDat, currDat = currDat, R = R,
                                num_treats = num_treats)
        }
        return(res_GCV_1B)
    }

    # Calculate all GCVs
    p_t <- p_list[t]
    res_GCV <- foreach(i=1:n.pred, .inorder = TRUE, .packages = "mvtnorm") %dopar% inner_b_GCV(i, p_t, p_list)
    res_GCV <- array(unlist(res_GCV), dim = c(B, p_t, n.train)) # Reformat into array

    # Calculate frequencies
    freqs <- matrix(0, ncol = 2, nrow = n.train)
    for (i in 1:n.train){
             temp <- apply(res_GCV[,,i], 1, which.max)
             freqs[i,1] <- sum(temp==1)
             freqs[i,2] <- sum(temp==2)
    }
    # all(apply(freqs, 1, sum) == n) # No need check?
    return(list("GCV_results" = res_GCV, "freqs" = freqs))
}
