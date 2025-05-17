#' Title
#'
#' @param Dat           Data in format returned by `generate_dataset_mvt`:
#' organised as a list of {y, X, A} where y is a vector of the final outcomes,
#' X is a list of matrices of the intermediate covariates and A is a matrix of
#' the assigned treatments
#' @param n             Number of samples/individuals
#' @param n.train       Train model using MC Draws from 1:n.train
#' @param num_stages    Total number of stages per individual
#' @param num_treats    Vector of number of treatment options at each stage
#' @param p_list        Vector of dimension for each stage
#' @param t             Prediction stage t, where \eqn{\leq}{<=} num_stages
#' @param R             Draw size. default = 30
#' @param numCores      Number of cores in the system to use. default uses parallel::detectCores()
#' @param tau           Prior precision scale. Should be specified with a small value. default = 0.01
#' @param B             Number of MC draws
#' @param nu0           Inverse-Wishart degres of freedom. default: 3
#' @param V0            Inverse-Wishart scale matrix. default: diagonalisation of p_list
#' @param alph          ???
#' @param gam           ???
#'
#' @returns list of {GCV_results, freqs} where GCV_results is an array of size
        #' n-n.train x p_list[t] x B, and freqs is a n-n.train x p_list[t] matrix
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
#' n           <- 1000
#' n.train     <- 500
#'
#' tau     <- 0.01
#' B       <- 5
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
#' Dat <- generate_dataset_mvt(n, num_stages, p_list, num_treats)
#'
#' # -----------------------------
#' # Main
#' # -----------------------------
#' gcv_res <- TestParallelGCV(Dat, n, n.train, num_stages, p_list, t, R, numCores,
#'                            tau, B, nu0, V0, alph, gam)
testParallelGCV <- function(Dat, n, n.train, num_stages, num_treats, p_list, t, R = 30, numCores,
                            tau = 0.01, B, nu0 = 3, V0 = mapply(diag, p_list, SIMPLIFY = FALSE), alph, gam
                            ) {

    # Verify inputs
    stopifnot("n.train must be less than n" = n.train < n)
    stopifnot("t must be less than or equal to num_stages" = t <= num_stages)

    # Retrieve training data and train model
    Dat.train <- c(list(Dat[[1]][1:n.train]), lapply(Dat[-1], function(x) x[1:n.train,]))
    res_mc <- compute_MC_draws_mvt(Data = Dat.train, tau = tau, num_treats = num_treats, B = B,
                                    nu0 = nu0, V0 = V0, alph = alph, gam = gam, p_list = p_list)

    # Set up parallel processing
    if (missing(numCores)) numCores <- parallel::detectCores()
    doParallel::registerDoParallel(numCores)  # use multicore, set to the number of our cores

    # Inner loop
    inner_b_GCV <- function(i, p_t, p_list) {
        res_GCV_1B <- matrix(0, nrow = B, ncol = p_t)
        for (b in 1:B) {
            histDat <- c(lapply(Dat[2:t], function(x) x[i,,drop = FALSE]), list(Dat[[num_stages + 2]][i,1,drop = FALSE]))
            currDat <- Dat[[t]][i,,drop = FALSE]
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
    res_GCV <- foreach(i=(n.train + 1):n, .inorder = TRUE, .packages = "mvtnorm") %dopar% inner_b_GCV(i, p_t, p_list)
    res_GCV <- array(unlist(res_GCV), dim = c(B, p_t, (n - n.train))) # Reformat into array

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
