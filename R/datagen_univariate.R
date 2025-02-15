generate_dataset_uvt <- function(n = 5000, T = 5, At_len = 3) {
    library(mvtnorm)

    # Step 1
    x_i1 <- matrix(rt(n * p_list[1], df = 10), nrow = n, ncol = p_list[1])

    # Step 2
    A <- matrix(sample(1:At_len, n*T, replace = TRUE), nrow = n, ncol = T)

    # Step 3
    X <- vector(mode = 'list', length = T) # Preallocate

    # Run first 3 lines fist bc can't be bothered making checks for numeric(0)s
    if (T > 0) {
        X[[1]] <- x_i1
    }

    if (T > 1) {
    t <- 2
        X[[t]] <- rmvnorm(n, sigma = diag(0.5^2, p_list[t]), mean = rep(0, p_list[t])) +
            (A[,t-1] == 2) * (t * X[[t-1]]) + # If a==2
            (A[,t-1] == 3) * (-t * X[[t-1]]) # Elif a==3
    }

    if (T > 2) {
    t <- 3
        X[[t]] <- rmvnorm(n, sigma = diag(0.5^2, p_list[t]), mean = rep(0, p_list[t])) +
            (A[,t-1] == 2) * (t * X[[t-1]] - (t-1) * X[[t-2]]) + # If a==2
            (A[,t-1] == 3) * (-t * X[[t-1]] + sqrt(t-1) * X[[t-2]]) # Elif a==3
    }

    if (T > 3) {
        for (t in 4:T) {
            # Could use dplyr instead
            X[[t]] <- rmvnorm(n, sigma = diag(0.5^2, p_list[t]), mean = rep(0, p_list[t])) +
                (A[,t-1] == 2) * (t * X[[t-1]] - (t-1) * X[[t-2]] + (t-2) * X[[t-3]]) + # If a==2
                (A[,t-1] == 3) * (-t * X[[t-1]] + sqrt(t-1) * X[[t-2]] + sqrt(t-2) * X[[t-3]]) # Elif a==3
        }
    }

    # Normalise X
    if (T > 0) {
        X <- lapply(X, function(Xt) apply(Xt, 2, function(z) (z-mean(z))/sd(z)))
    }

    # Step 4
    mi <- 3

    if (T > 0) {
        t <- 1
        mi <- mi + (A[,t] == 2) * (sin(10*t) * X[[t]]) + (A[,t] == 3) * (cos(10*t) * X[[t]])
    }

    if (T > 1) {
        t <- 2
        mi <- mi + (A[,t] == 2) * (sin(10*t) * X[[t]] - sin(10*t - 10) * X[[t-1]]) +
            (A[,t] == 3) * (cos(10*t) * X[[t]] - cos(10*t - 10) * X[[t-1]])
    }

    if (T > 2) {
        for (t in 3:T) {
            mi <- mi +
                (A[,t] == 2) * (sin(10*t) * X[[t]] - sin(10*t - 10) * X[[t-1]] + sin(10*t - 20) * X[[t-2]]) +
                (A[,t] == 3) * (cos(10*t) * X[[t]] - cos(10*t - 10) * X[[t-1]] + sqrt(abs(cos(10*t - 20))) * X[[t-2]])
        }
    }

    yi <- rnorm(n, mean = mi, sd = 1)

    return(c(list(yi), X, list(A)))
}

set.seed(1) #remove later

# Starting Scalar values
n       <- 5000
T       <- 5
p_list  <- rep(1, T)
At_len  <- 3

Data <- generate_dataset_uvt(n, T, At_len)
# Data <- res[1:(T+1)]
A <- Data[[T+2]]
