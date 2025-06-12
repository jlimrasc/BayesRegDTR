#' Main function for fitting a Bayesian likelihood-based linear regression model
#'
#' Fits the Bayesian likelihood-based linear model to obtain an estimated posterior
#' distribution of the optimal treatment option at a user-specified prediction stage.
#' Uses backward induction and dynamic programming theory for computing
#' expected values.
#'
#' Utilises a \link[future]{future} framework, so to enable
#' parallel processing and register a parallel backend, \link[future]{plan} and
#' \link[doFuture]{registerDoParallel} must be called first. Additionally,
#' progress bars use \link[progressr]{progressr} API, and a non-default progress
#' bar (e.g. cli) is recommended. See below or \link[doFuture]{registerDoParallel} and
#' \link[progressr]{handlers} for examples. Note that to have a progress bar for
#' the parallel sections, future must be used.
#'
#'
#' @param Dat.train     Training data in format returned by `generate_dataset`:
#' organised as a list of \eqn{\{y, X_1, X_2..., X_{num\_stages}, A\}} where y is a
#' vector of the final outcomes, \eqn{X_1, X_2..., X_{num\_stages}} is a list of matrices
#' of the intermediate covariates and A is an \eqn{n.train \times num\_stages}{n.train x num_stages} matrix of the
#' assigned treatments, where num_stages is the total number of stages
#' @param Dat.pred      Prediction data in format returned by `generate_dataset`:
#' organised as a list of \eqn{\{X_1, X_2..., X_t, A\}} where
#' \eqn{X_1, X_2..., X_t} is a list of matrices of the intermediate
#' covariates and A is an \eqn{n.pred \times (t-1)}{n.pred x (t-1)} matrix of the assigned treatments,
#' where t is the prediction stage
#' @param n.train       Number of samples/individuals in the training data
#' @param n.pred        Number of samples/individuals in the prediction data
#' @param num_stages    Total number of stages
#' @param num_treats    Vector of number of treatment options at each stage
#' @param p_list        Vector of intermediate covariate dimensions for each stage
#' @param t             Prediction stage t, where t \eqn{\leq}{<=} num_stages
#' @param R             Draw size from distribution of intermediate covariates. default:  30
#' @param numCores      Number of cores in the system to use. default uses parallel::detectCores() - 1
#' @param tau           Normal prior scale parameter for regression coefficients. Should be specified with a small value. default:  0.01
#' @param B             Number of MC draws from posterior of regression parameters. default 10000
#' @param nu0           Inverse-Wishart prior degrees of freedom for regression error Vcov matrix. Ignored if using a univariate dataset. default: 3
#' @param V0            List of Inverse-Wishart prior scale matrix for regression error Vcov matrix. Ignored if using a univariate dataset. default: list of identity matrices
#' @param alph          Inverse-Gamma prior shape parameter for regression error variance of y. default:  1
#' @param gam           Inverse-Gamma prior rate parameter for regression error variance of y. default:  1
#' @param showBar       Whether to show a progress bar. Uses bar from \link[progress]{progress_bar}
#' and \link[doParabar]{doParabar} for parallel integration deafult: TRUE
#'
#' @returns
#' \item{GCV_results}{An array of dimension
#' \eqn{n.pred \times num\_treats[t] \times B}{n.pred x num_treats[t] x B},
#' indicating the expected value under each treatment option at stage t.}
#' \item{post.prob}{An \eqn{n.pred \times num\_treats[t]}{n.pred x num_treats[t]}
#' matrix of the posterior probability that each treatment type at stage t is optimal}
#' \item{MC_draws.train}{A list of Monte Carlo draws containing:\itemize{
#' \item{\emph{sigmat_B_list} - A list of length num_stages with each element a
#' vector of size \eqn{B \times p\_list[t]}{B x p_list[t]}}
#' \item{\emph{Wt_B_list} - A list of length num_stages with each element a
#' matrix of size \eqn{B \times p\_list[t]}{B x p_list[t]}}
#' \item{\emph{beta_B} - A list of length B}
#' \item{\emph{sigmay_2B} - A list of length B}
#' }}
#' @export
#'
#' @examples
#' # Code does not run within 10 seconds, so don't run
#' \dontrun{
#' progressr::handlers("cli")          # Set handler to something with title/text
#' future::plan(future::multisession)  # or plan(multicore) on Unix
#' doFuture::registerDoFuture()        # or doParallel::registerDoParallel()
#'                                     # if no progress bar is needed and future
#'                                     # is unwanted
#'
#' ## UVT
#' # -----------------------------
#' # Initialise Inputs
#' # -----------------------------
#' num_stages  <- 5
#' t           <- 3
#' p_list      <- rep(1, num_stages)
#' num_treats  <- rep(2, num_stages)
#' n.train     <- 5000
#' n.pred      <- 10
#' numCores    <- parallel::detectCores()
#'
#' # -----------------------------
#' # Generate Dataset
#' # -----------------------------
#' Dat.train  <- generate_dataset(n.train,  num_stages, p_list, num_treats)
#' Dat.pred  <- generate_dataset(n.pred,  num_stages, p_list, num_treats)
#' Dat.pred  <- Dat.pred[-1]
#' Dat.pred[[num_stages+1]]  <- Dat.pred[[num_stages+1]][1:n.pred, 1:(t-1), drop = FALSE]
#'
#' # -----------------------------
#' # Main
#' # -----------------------------
#' gcv_uvt <- BayesLinRegDTR.model.fit(Dat.train, Dat.pred, n.train, n.pred,
#'                                     num_stages, num_treats,
#'                                     p_list, t, R = 30, numCores,
#'                                     tau = 0.01, B = 500, nu0 = NULL,
#'                                     V0 = NULL, alph = 3, gam = 4)
#'
#' ## MVT
#' # -----------------------------
#' # Initialise Inputs
#' # -----------------------------
#' num_stages  <- 3
#' t           <- 2
#' p_list      <- rep(2, num_stages)
#' num_treats  <- rep(2, num_stages)
#' n.train     <- 5000
#' n.pred      <- 10
#' numCores    <- parallel::detectCores()
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
#' gcv_res <- BayesLinRegDTR.model.fit(Dat.train, Dat.pred, n.train, n.pred,
#'                                     num_stages, num_treats,
#'                                     p_list, t, R = 30, numCores,
#'                                     tau = 0.01, B = 500, nu0 = 3,
#'                                     V0 = mapply(diag, p_list, SIMPLIFY = FALSE),
#'                                     alph = 3, gam = 4)
#' }
BayesLinRegDTR.model.fit <- function(Dat.train, Dat.pred, n.train, n.pred, num_stages, num_treats,
                                     p_list, t, R = 30, numCores,
                            tau = 0.01, B = 10000, nu0 = 3,
                            V0 = mapply(diag, p_list, SIMPLIFY = FALSE),
                            alph = 1, gam = 1, showBar = TRUE
                            ) {

    # Verify inputs
    stopifnot("t must be less than or equal to num_stages" = t <= num_stages)

    # Check parallelism
    current_plan <- future::plan()
    doParName <- foreach::getDoParName()
    # Check if the current plan is sequential (no parallelism)
    if (doParName == "doFuture" && inherits(current_plan, "sequential")) {
        warning("No parallel backend detected: future plan is 'sequential'. ",
                "For better performance, consider setting a parallel plan, e.g., plan(multisession).")
    }

    # Check progress bar
    if (showBar && any(vapply(progressr::handlers(), identical, logical(1), progressr::handler_txtprogressbar())))
        reporting <- TRUE
    else
        reporting <- FALSE


    # Retrieve training data and train model
    if (reporting) print("=== Computing MC Draws ===")
    if (any(p_list > 1)) {
        res_mc <- compute_MC_draws_mvt(Data = Dat.train, tau = tau, num_treats = num_treats, B = B,
                                        nu0 = nu0, V0 = V0, alph = alph, gam = gam, p_list = p_list,
                                       showBar = showBar)
        res_mc_f <- list("Wt_B_list" = res_mc$Wt_B_list, "sigmat_B_list" = res_mc$sigmat_B_list) # Formatted copy without extra data
    }
    else {
        res_mc <- compute_MC_draws_uvt(Data = Dat.train, tau = tau, num_treats = num_treats, B = B,
                                       alph = alph, gam = gam, p_list = p_list, showBar = showBar)
        res_mc_f <- list("thetat_B_list" = res_mc$thetat_B_list, "sigmat_2B_list" = res_mc$sigmat_2B_list) # Formatted copy without extra data
    }
    if (reporting) print("=== MC Draws Completed ===")


    # Inner loop
    # MVT Version
    inner_b_GCV_MVT <- function(i, ntreats_t, p_list, histDat, currDat) {
        res_GCV_1B <- matrix(0, nrow = B, ncol = ntreats_t)
        histDat <- c(lapply(Dat.pred[1:(t-1)], function(x) x[i,,drop = FALSE]), list(Dat.pred[[num_stages + 1]][i,1:(t-1),drop = FALSE]))
        currDat <- Dat.pred[[t]][i,,drop = FALSE]
        for (b in 1:B) {
            Wt      <- lapply(res_mc_f$Wt_B_list, function(x) x[[b]])
            Sigmat  <- lapply(res_mc_f$sigmat_B_list, function(x) x[[b]])

            res_GCV_1B[b,] <-
                GiveChoiceValue(Wt = Wt, Sigmat = Sigmat, bet = res_mc$beta_B[,b],
                                sigmay = res_mc$sigmay_2B[b], t = t,
                                num_stages = num_stages, p_list = p_list,
                                histDat = histDat, currDat = currDat, R = R,
                                num_treats = num_treats)
        }
        return(res_GCV_1B)
    }

    # UVT Version
    inner_b_GCV_UVT <- function(i, ntreats_t, p_list) {
        res_GCV_1B <- matrix(0, nrow = B, ncol = ntreats_t)
        histDat <- c(lapply(Dat.pred[1:(t-1)], function(x) x[i,,drop = FALSE]), list(Dat.pred[[num_stages + 1]][i,1:(t-1),drop = FALSE]))
        currDat <- Dat.pred[[t]][i,,drop = FALSE]
        for (b in 1:B) {
            thetat  <- lapply(res_mc_f$thetat_B_list, function(x) matrix(x[,b]))
            Sigmat  <- lapply(res_mc_f$sigmat_2B_list, function(x) matrix(x[b]))

            res_GCV_1B[b,] <-
                GiveChoiceValue(Wt = thetat, Sigmat = Sigmat, bet = res_mc$beta_B[,b],
                                sigmay = res_mc$sigmay_2B[b], t = t,
                                num_stages = num_stages, p_list = p_list,
                                histDat = histDat, currDat = currDat, R = R,
                                num_treats = num_treats)
        }
        return(res_GCV_1B)
    }
    # Set up parallel processing
    # if (missing(numCores)) numCores <- parallel::detectCores() - 1
    # doParallel::registerDoParallel(numCores)  # use multicore, set to the number of our cores
    # with(doFuture::registerDoFuture(), local = TRUE) # register future backend with foreach
    # backend <- parabar::start_backend(cores = numCores, cluster_type = "psock", backend_type = "async") # Backend for progress bar
    # doParabar::registerDoParabar(backend) # Register it with the `foreach` package.

    # Calculate all GCVs
    if (reporting) print("=== Predicting Data ===")
    ntreats_t <- num_treats[t]
    progressr::with_progress({
        p <- progressr::progressor(steps = n.pred + 1, message = "Predicting data", enable = showBar)
        if (any(p_list > 1)) {
            # res_GCV <- foreach(i=1:n.pred, .inorder = TRUE, .packages = "mvtnorm",
            #                    .export = c("inner_b_GCV_MVT", "ntreats_t", "p_list",
            #                                "Dat.pred", "res_mc_f", "GiveChoiceValue")) %dopar% inner_b_GCV_MVT(i, ntreats_t, p_list)
            res_GCV <- foreach(i = 1:n.pred, .inorder = TRUE, .packages = "mvtnorm", .export = c("GiveChoiceValue")) %dorng% {
                p()
                inner_b_GCV_MVT(i, ntreats_t, p_list)
            }
        }

        else
            res_GCV <- foreach(i = 1:n.pred, .inorder = TRUE, .packages = "mvtnorm") %dorng% {
                p()
                inner_b_GCV_UVT(i, ntreats_t, p_list)
            }
        p()
        res_GCV <- array(unlist(res_GCV), dim = c(B, ntreats_t, n.train)) # Reformat into array

    })


    if (reporting) print("=== Prediction Completed ===")

    # Stop the backend
    # parabar::stop_backend(backend)

    # Calculate frequencies
    if (reporting) print("=== Computing Frequencies ===")
    # Create bar
    # showBarOld <- getOption("progress_enabled")
    # options(progress_enabled = showBar)
    # pb <- progress::progress_bar$new( format = " Calculating Freqs [:bar] :percent eta: :eta", total = n.train)
    progressr::with_progress({
        p <- progressr::progressor(steps = floor(n.train/10), enable = showBar,
                                   message = "Computing Frequencies")
        freqs <- matrix(0, ncol = ntreats_t, nrow = n.train)
        foreach (i = 1:n.train) %dopar% {
                 temp <- apply(res_GCV[,,i], 1, which.max)
                 freqs[i, 1:ntreats_t] <- tabulate(temp, nbins = ntreats_t)
                 # cli::cli_progress_update()
                 # pb$tick()
                 if (i%%10 == 0) p()
        }
        post.prob <- freqs/B
    })

    # Complete bar
    # cli::cli_progress_update()
    # cli::cli_progress_done()
    # pb$tick(tokens = list(stg = num_stages + 1))
    # pb$terminate()
    # Reset bar settings
    # options(progress_enabled = showBarOld)
    if (reporting) print("=== Frequencies Completed ===")


    return(list("GCV_results" = res_GCV, "post.prob" = post.prob, "MC_draws.train" = res_mc))
}
