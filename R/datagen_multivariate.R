generate_dataset_mvt <- function(n, T, p_list, At_len) {
    library(mvtnorm)

    # Step 1
    x_i1 <- rmvt(n, sigma = diag(p_list[1]), df = 10, delta = rep(0, p_list[1]), type = c("shifted"))

    # Step 2
    A <- matrix(sample(1:At_len, n*T, replace = TRUE), nrow = n, ncol = T)

    # Step 3
    # X <- array(0, dim = c(n, p_t, T)) # Preallocate
    X <- vector(mode = 'list', length = T) # Preallocate

    # Run first 3 lines fist bc inefficient to check for numeric(0)s
    if (T > 0) {
        X[[1]] <- x_i1
    }

    if (T > 1) {
        t <- 2
        C_max <- matrix(0, nrow = max(p_list), ncol = max(p_list))
        C_max <- (-1) ^ (row(C_max) + col(C_max))

        C <- C_max[1:p_list[t], 1:p_list[t-1]]

        gen_xit <- function(X_it1, A2, A3, t) {
            A2 * ( t * C %*% X_it1) + # If a==2
            A3 * (-t * C %*% X_it1) # Elif a==3
        }

        X[[t]] <- rmvnorm(n, sigma = diag(0.5^2, p_list[t]), mean = rep(0, p_list[t])) +
            t(unname(mapply(FUN = gen_xit,
                     split(X[[t-1]], row(X[[t-1]])),
                     split(A[,t-1] == 2, 1:n),
                     split(A[,t-1] == 3, 1:n),
                     t = t)))
    }

    if (T > 2) {
        t <- 3
        C <- C_max[1:p_list[t], 1:p_list[t-1]]
        C2 <- C_max[1:p_list[t], 1:p_list[t-2]]
        gen_xit <- function(X_it1, X_it2, A2, A3, t, p_t) {
            A2 * ( t * C %*% X_it1 -     (t-1) * C2 %*% X_it2) + # If a==2
            A3 * (-t * C %*% X_it1 + sqrt(t-1) * C2 %*% X_it2) # Elif a==3
        }
        X[[t]] <- rmvnorm(n, sigma = diag(0.5^2, p_list[t]), mean = rep(0, p_list[t])) +
            t(unname(mapply(FUN = gen_xit,
                     split(X[[t-1]], row(X[[t-1]])),
                     split(X[[t-2]], row(X[[t-2]])),
                     split(A[,t-1] == 2, 1:n),
                     split(A[,t-1] == 3, 1:n),
                     t = t)))
    }

    if (T > 3) {
        for (t in 4:T) {
            # Could use dplyr instead
            C <- C_max[1:p_list[t], 1:p_list[t-1]]
            C2 <- C_max[1:p_list[t], 1:p_list[t-2]]
            C3 <- C_max[1:p_list[t], 1:p_list[t-3]]

            gen_xit <- function(X_it1, X_it2, X_it3, A2, A3, t, p_t) {
                A2 * ( t * C %*% X_it1 -     (t-1) * C2 %*% X_it2 +     (t-2) * C3 %*% X_it3) + # If a==2
                A3 * (-t * C %*% X_it1 + sqrt(t-1) * C2 %*% X_it2 + sqrt(t-2) * C3 %*% X_it3) # Elif a==3
            }

            X[[t]] <- rmvnorm(n, sigma = diag(0.5^2, p_list[t]), mean = rep(0, p_list[t])) +
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
    # if (T > 0) {
    #     X <- lapply(X, function(Xt) apply(Xt, 2, function(z) (z-mean(z))/sd(z)))
    # }

    # Step 4
    mi <- 3

    if (T > 0) {
        t <- 1
        # browser()
        mi <- mi + (A[,t] == 2) * (sin(10*t) * X[[t]]) %*% rep(1, p_list[t]) + (A[,t] == 3) * (cos(10*t) * X[[t]] %*% rep(1, p_list[t]))

    }

    if (T > 1) {
        t <- 2
        mi <- mi + (A[,t] == 2) * (sin(10*t) * X[[t]] %*% rep(1, p_list[t]) - sin(10*t - 10) * X[[t-1]] %*% rep(1, p_list[t-1]))  +
                   (A[,t] == 3) * (cos(10*t) * X[[t]] %*% rep(1, p_list[t]) - cos(10*t - 10) * X[[t-1]] %*% rep(1, p_list[t-1]))
    }

    if (T > 2) {
        for (t in 3:T) {
            mi <- mi +
                (A[,t] == 2) * (sin(10*t) * X[[t]] %*% rep(1, p_list[t]) - sin(10*t - 10) * X[[t-1]] %*% rep(1, p_list[t-1]) + sin(10*t - 20) * X[[t-2]] %*% rep(1, p_list[t-2]) )+
                (A[,t] == 3) * (cos(10*t) * X[[t]] %*% rep(1, p_list[t]) - cos(10*t - 10) * X[[t-1]] %*% rep(1, p_list[t-1]) + sqrt(abs(cos(10*t - 20))) * X[[t-2]] %*% rep(1, p_list[t-2]))
        }
    }

    yi <- rnorm(n, mean = mi, sd = 1)
    return(c(list(yi), X, list(A)))
}

set.seed(1) #remove later

# Starting Scalar values
n       <- 5000
T       <- 3
p_list  <- rep(2, T)
At_len  <- 3

Data <- generate_dataset_mvt(n, T, p_list, At_len)
# Data <- res[1:(T+1)]
# X <- Data[2:(T+1)]
# A <- Data[[T+2]]
# y <- Data[[1]]
