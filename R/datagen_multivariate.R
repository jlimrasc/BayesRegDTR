#' Generate Multivariate dataset
#'
#' @param n             Number of samples/individuals to generate
#' @param num_stages    Total number of stages per individual
#' @param p_list        Vector of dimension for each stage
#' @param num_treats    Vector of number of treatment options at each stage
#'
#' @returns Observed data organised as a list of {y, X, A} where y is a vector of the final outcomes,
         #' X is a list of matrices of the intermediate covariates
         #' and A is a matrix of the assigned treatments
#' @export
#'
#' @examples
#' # -----------------------------
#' # Initialise Inputs
#' # -----------------------------
#' n           <- 5000
#' num_stages  <- 3
#' p_list      <- rep(2, num_stages)
#' num_treats  <- rep(3, num_stages)
#'
#' # -----------------------------
#' # Main
#' # -----------------------------
#' Data        <- generate_dataset_mvt(n, num_stages, p_list, num_treats)
generate_dataset_mvt <- function(n, num_stages, p_list, num_treats) {

    stopifnot("num_treats length must equal num_stages" = length(force(num_treats)) == num_stages)

    # Step 1
    x_i1 <- mvtnorm::rmvt(n, sigma = diag(p_list[1]), df = 10, delta = rep(0, p_list[1]), type = c("shifted"))

    # Step 2
    A <- sapply(num_treats, function(a_max) sample(a_max, n, replace = TRUE))

    # Step 3
    X <- vector(mode = 'list', length = num_stages) # Preallocate

    # Run first 3 lines fist bc inefficient to check for numeric(0)s
    if (num_stages > 0) {
        X[[1]] <- x_i1
    }

    if (num_stages > 1) {
        t <- 2
        C_max <- matrix(0, nrow = max(p_list), ncol = max(p_list))
        C_max <- (-1) ^ (row(C_max) + col(C_max))

        C <- C_max[1:p_list[t], 1:p_list[t-1]]

        gen_xit <- function(X_it1, A2, A3, t) {
            A2 * ( t * C %*% X_it1) + # If a==2
            A3 * (-t * C %*% X_it1) # Elif a==3
        }

        X[[t]] <- mvtnorm::rmvnorm(n, sigma = diag(0.5^2, p_list[t]), mean = rep(0, p_list[t])) +
            t(unname(mapply(FUN = gen_xit,
                     split(X[[t-1]], row(X[[t-1]])),
                     split(A[,t-1] == 2, 1:n),
                     split(A[,t-1] == 3, 1:n),
                     t = t)))
    }

    if (num_stages > 2) {
        t <- 3
        C <- C_max[1:p_list[t], 1:p_list[t-1]]
        C2 <- C_max[1:p_list[t], 1:p_list[t-2]]
        gen_xit <- function(X_it1, X_it2, A2, A3, t, p_t) {
            A2 * ( t * C %*% X_it1 -     (t-1) * C2 %*% X_it2) + # If a==2
            A3 * (-t * C %*% X_it1 + sqrt(t-1) * C2 %*% X_it2) # Elif a==3
        }
        X[[t]] <- mvtnorm::rmvnorm(n, sigma = diag(0.5^2, p_list[t]), mean = rep(0, p_list[t])) +
            t(unname(mapply(FUN = gen_xit,
                     split(X[[t-1]], row(X[[t-1]])),
                     split(X[[t-2]], row(X[[t-2]])),
                     split(A[,t-1] == 2, 1:n),
                     split(A[,t-1] == 3, 1:n),
                     t = t)))
    }

    if (num_stages > 3) {
        for (t in 4:num_stages) {
            # Could use dplyr instead
            C <- C_max[1:p_list[t], 1:p_list[t-1]]
            C2 <- C_max[1:p_list[t], 1:p_list[t-2]]
            C3 <- C_max[1:p_list[t], 1:p_list[t-3]]

            gen_xit <- function(X_it1, X_it2, X_it3, A2, A3, t, p_t) {
                A2 * ( t * C %*% X_it1 -     (t-1) * C2 %*% X_it2 +     (t-2) * C3 %*% X_it3) + # If a==2
                A3 * (-t * C %*% X_it1 + sqrt(t-1) * C2 %*% X_it2 + sqrt(t-2) * C3 %*% X_it3) # Elif a==3
            }

            X[[t]] <- mvtnorm::rmvnorm(n, sigma = diag(0.5^2, p_list[t]), mean = rep(0, p_list[t])) +
                t(unname(mapply(FUN = gen_xit,
                         split(X[[t-1]], row(X[[t-1]])),
                         split(X[[t-2]], row(X[[t-2]])),
                         split(X[[t-3]], row(X[[t-3]])),
                         split(A[,t-1] == 2, 1:n),
                         split(A[,t-1] == 3, 1:n),
                         t = t)))
        }
    }

    # # Normalise X
    # if (num_stages > 0) {
    #     X <- lapply(X, function(Xt) apply(Xt, 2, function(z) (z-mean(z))/sd(z)))
    # }

    # Step 4
    mi <- 3

    if (num_stages > 0) {
        t <- 1
        # browser()
        mi <- mi + (A[,t] == 2) * (sin(10*t) * X[[t]]) %*% rep(1, p_list[t]) + (A[,t] == 3) * (cos(10*t) * X[[t]] %*% rep(1, p_list[t]))

    }

    if (num_stages > 1) {
        t <- 2
        mi <- mi + (A[,t] == 2) * (sin(10*t) * X[[t]] %*% rep(1, p_list[t]) - sin(10*t - 10) * X[[t-1]] %*% rep(1, p_list[t-1]))  +
                   (A[,t] == 3) * (cos(10*t) * X[[t]] %*% rep(1, p_list[t]) - cos(10*t - 10) * X[[t-1]] %*% rep(1, p_list[t-1]))
    }

    if (num_stages > 2) {
        for (t in 3:num_stages) {
            mi <- mi +
                (A[,t] == 2) * (sin(10*t) * X[[t]] %*% rep(1, p_list[t]) - sin(10*t - 10) * X[[t-1]] %*% rep(1, p_list[t-1]) + sin(10*t - 20) * X[[t-2]] %*% rep(1, p_list[t-2]) )+
                (A[,t] == 3) * (cos(10*t) * X[[t]] %*% rep(1, p_list[t]) - cos(10*t - 10) * X[[t-1]] %*% rep(1, p_list[t-1]) + sqrt(abs(cos(10*t - 20))) * X[[t-2]] %*% rep(1, p_list[t-2]))
        }
    }

    yi <- rnorm(n, mean = mi, sd = 1)
    return(c(list(yi), X, list(A)))
}
