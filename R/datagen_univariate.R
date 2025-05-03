#' Generate Univariate Dataset
#'
#' @param n             Number of samples to generate
#' @param num_treats    Total number of stages per sample
#' @param At_lens        Vector of number of treatment options at each stage
#'
#' @returns Observed data organised as a list of {y, X, A} where y is a vector of the final outcomes,
#' X is a list of matrices of the intermediate covariates
#' and A is a matrix of the assigned treatments
#' @export
#' @examples
#' # -----------------------------
#' # Initialise Inputs
#' # -----------------------------
#' n           <- 500
#' num_treats  <- 5
#' At_lens     <- rep(3, num_treats)
#'
#' # -----------------------------
#' # Main
#' # -----------------------------
#' Data <- generate_dataset_uvt(n, num_treats, At_lens)
#'
generate_dataset_uvt <- function(n, num_treats, At_lens) {
    library(mvtnorm)

    stopifnot("At_lens length must equal num_treats" = length(force(At_lens)) == num_treats)

    p_list      <- rep(1, num_treats)

    # Step 1
    x_i1 <- matrix(rt(n * p_list[1], df = 10), nrow = n, ncol = p_list[1])

    # Step 2
    A <- sapply(At_lens, function(a_max) sample(a_max, n, replace = TRUE))

    # Step 3
    X <- vector(mode = 'list', length = num_treats) # Preallocate

    # Run first 3 lines fist bc can't be bothered making checks for numeric(0)s
    if (num_treats > 0) {
        X[[1]] <- x_i1
    }

    if (num_treats > 1) {
    t <- 2
        X[[t]] <- rmvnorm(n, sigma = diag(0.5^2, p_list[t]), mean = rep(0, p_list[t])) +
            (A[,t-1] == 2) * (t * X[[t-1]]) + # If a==2
            (A[,t-1] == 3) * (-t * X[[t-1]]) # Elif a==3
    }

    if (num_treats > 2) {
    t <- 3
        X[[t]] <- rmvnorm(n, sigma = diag(0.5^2, p_list[t]), mean = rep(0, p_list[t])) +
            (A[,t-1] == 2) * (t * X[[t-1]] - (t-1) * X[[t-2]]) + # If a==2
            (A[,t-1] == 3) * (-t * X[[t-1]] + sqrt(t-1) * X[[t-2]]) # Elif a==3
    }

    if (num_treats > 3) {
        for (t in 4:num_treats) {
            # Could use dplyr instead
            X[[t]] <- rmvnorm(n, sigma = diag(0.5^2, p_list[t]), mean = rep(0, p_list[t])) +
                (A[,t-1] == 2) * (t * X[[t-1]] - (t-1) * X[[t-2]] + (t-2) * X[[t-3]]) + # If a==2
                (A[,t-1] == 3) * (-t * X[[t-1]] + sqrt(t-1) * X[[t-2]] + sqrt(t-2) * X[[t-3]]) # Elif a==3
        }
    }

    # Normalise X
    # if (num_treats > 0) {
    #     X <- lapply(X, function(Xt) apply(Xt, 2, function(z) (z-mean(z))/sd(z)))
    # }

    # Step 4
    mi <- 3

    if (num_treats > 0) {
        t <- 1
        mi <- mi + (A[,t] == 2) * (sin(10*t) * X[[t]]) + (A[,t] == 3) * (cos(10*t) * X[[t]])
    }

    if (num_treats > 1) {
        t <- 2
        mi <- mi + (A[,t] == 2) * (sin(10*t) * X[[t]] - sin(10*t - 10) * X[[t-1]]) +
            (A[,t] == 3) * (cos(10*t) * X[[t]] - cos(10*t - 10) * X[[t-1]])
    }

    if (num_treats > 2) {
        for (t in 3:num_treats) {
            mi <- mi +
                (A[,t] == 2) * (sin(10*t) * X[[t]] - sin(10*t - 10) * X[[t-1]] + sin(10*t - 20) * X[[t-2]]) +
                (A[,t] == 3) * (cos(10*t) * X[[t]] - cos(10*t - 10) * X[[t-1]] + sqrt(abs(cos(10*t - 20))) * X[[t-2]])
        }
    }

    yi <- rnorm(n, mean = mi, sd = 1)

    return(c(list(yi), X, list(A)))
}

