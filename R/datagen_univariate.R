#' Generate Univariate Dataset
#'
#' @param n             Number of samples/individuals to generate
#' @param num_stages    Total number of stages per individual
#' @param num_treats    Vector of number of treatment options at each stage
#'
#' @returns Observed data organised as a list of \eqn{\{y, X, A\}} where y is a
#' vector of the final outcomes, X is a list of matrices of the intermediate
#' covariates and A is a matrix of the assigned treatments
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # -----------------------------
#' # Initialise Inputs
#' # -----------------------------
#' n           <- 500
#' num_stages  <- 5
#' num_treats  <- rep(3, num_stages)
#'
#' # -----------------------------
#' # Main
#' # -----------------------------
#' Data <- generate_dataset_uvt(n, num_stages, num_treats)
#' }
generate_dataset_uvt <- function(n, num_stages, num_treats) {

    stopifnot("num_treats length must equal num_stages" = length(force(num_treats)) == num_stages)

    p_list      <- rep(1, num_stages)

    # Step 1
    x_i1 <- matrix(stats::rt(n * p_list[1], df = 10), nrow = n, ncol = p_list[1])

    # Step 2
    A <- sapply(num_treats, function(a_max) sample(a_max, n, replace = TRUE))

    # Step 3
    X <- vector(mode = 'list', length = num_stages) # Preallocate

    # Run first 3 lines fist bc can't be bothered making checks for numeric(0)s
    if (num_stages > 0) {
        X[[1]] <- x_i1
    }

    if (num_stages > 1) {
    t <- 2
        X[[t]] <- mvtnorm::rmvnorm(n, sigma = diag(0.5^2, p_list[t]), mean = rep(0, p_list[t])) +
            (A[,t-1] == 2) * (t * X[[t-1]]) + # If a==2
            (A[,t-1] == 3) * (-t * X[[t-1]]) # Elif a==3
    }

    if (num_stages > 2) {
    t <- 3
        X[[t]] <- mvtnorm::rmvnorm(n, sigma = diag(0.5^2, p_list[t]), mean = rep(0, p_list[t])) +
            (A[,t-1] == 2) * (t * X[[t-1]] - (t-1) * X[[t-2]]) + # If a==2
            (A[,t-1] == 3) * (-t * X[[t-1]] + sqrt(t-1) * X[[t-2]]) # Elif a==3
    }

    if (num_stages > 3) {
        for (t in 4:num_stages) {
            # Could use dplyr instead
            X[[t]] <- mvtnorm::rmvnorm(n, sigma = diag(0.5^2, p_list[t]), mean = rep(0, p_list[t])) +
                (A[,t-1] == 2) * (t * X[[t-1]] - (t-1) * X[[t-2]] + (t-2) * X[[t-3]]) + # If a==2
                (A[,t-1] == 3) * (-t * X[[t-1]] + sqrt(t-1) * X[[t-2]] + sqrt(t-2) * X[[t-3]]) # Elif a==3
        }
    }

    # Normalise X
    # if (num_stages > 0) {
    #     X <- lapply(X, function(Xt) apply(Xt, 2, function(z) (z-mean(z))/sd(z)))
    # }

    # Step 4
    mi <- 3

    if (num_stages > 0) {
        t <- 1
        mi <- mi + (A[,t] == 2) * (sin(10*t) * X[[t]]) + (A[,t] == 3) * (cos(10*t) * X[[t]])
    }

    if (num_stages > 1) {
        t <- 2
        mi <- mi + (A[,t] == 2) * (sin(10*t) * X[[t]] - sin(10*t - 10) * X[[t-1]]) +
            (A[,t] == 3) * (cos(10*t) * X[[t]] - cos(10*t - 10) * X[[t-1]])
    }

    if (num_stages > 2) {
        for (t in 3:num_stages) {
            mi <- mi +
                (A[,t] == 2) * (sin(10*t) * X[[t]] - sin(10*t - 10) * X[[t-1]] + sin(10*t - 20) * X[[t-2]]) +
                (A[,t] == 3) * (cos(10*t) * X[[t]] - cos(10*t - 10) * X[[t-1]] + sqrt(abs(cos(10*t - 20))) * X[[t-2]])
        }
    }

    yi <- stats::rnorm(n, mean = mi, sd = 1)

    return(c(list(yi), X, list(A)))
}

