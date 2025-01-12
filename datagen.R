set.seed(1) #remove later
library(mvtnorm)
# Starting Scalar values
n       <- 5000
T       <- 5
p_t     <- 1
At_max  <- 3

# Step 1
x_i1 <- rt(n, 10) # rmvt(n, sigma = diag(pt), df = 10) # not shifted? so no delta?

# Step 2
A <- matrix(sample(1:At_max, n*T, replace = TRUE), nrow = n, ncol = T)

# Step 3
X <- matrix(0, nrow = n, ncol = T) # Preallocate

# Run first 3 lines fist bc can't be bothered making checks for numeric(0)s
X[,1] <- x_i1
t <- 2
X[,t] <- rmvnorm(n, sigma = matrix(0.5^2, p_t), mean = matrix(0, p_t)) +
    (A[,t-1] == 2) * (t * X[,t-1]) + # If a==2
    (A[,t-1] == 3) * (-t * X[,t-1]) # Elif a==3
t <- 3
X[,t] <- rmvnorm(n, sigma = matrix(0.5^2, p_t), mean = matrix(0, p_t)) +
    (A[,t-1] == 2) * (t * X[,t-1] - (t-1) * X[,t-2]) + # If a==2
    (A[,t-1] == 3) * (-t * X[,t-1] + sqrt(t-1) * X[,t-2]) # Elif a==3

for (t in 4:T) {
    # Could use dplyr instead
    X[,t] <- rmvnorm(n, sigma = matrix(0.5^2, p_t), mean = matrix(0, p_t)) +
        (A[,t-1] == 2) * (t * X[,t-1] - (t-1) * X[,t-2] + (t-2) * X[,t-3]) + # If a==2
        (A[,t-1] == 3) * (-t * X[,t-1] + sqrt(t-1) * X[,t-2] + sqrt(t-2) * X[,t-3]) # Elif a==3
}

# Normalise X
X <- apply(X, 2, function(z) (z-mean(z))/sd(z))

# gen_xit <- function(a, x_i, t) {
#     return(
#         ifelse (a == 2, t * x_i[t-1] - max(t-1, 0) * x_i[t-2] + max(t-2, 0) * x_i[t-3], # If a==2
#         ifelse (a == 3, -t * x_i[t-1] + sqrt(max(t-1, 0)) * x_i[t-2] + sqrt(max(t-2, 0)) * x_i[t-3]), # Elif a==3
#                 0) # Else
#     )
# }

# Step 4
mi <- 3
t <- 1
mi <- mi + (A[,t] == 2) * (sin(10*t) * X[,t]) + (A[,t] == 3) * (cos(10*t) * X[,t])
t <- 2
mi <- mi + (A[,t] == 2) * (sin(10*t) * X[,t] - sin(10*t - 10) * X[,t-1]) + 
    (A[,t] == 3) * (cos(10*t) * X[,t] - cos(10*t - 10) * X[,t-1])
for (t in 3:T) {
    mi <- mi + 
        (a<-(A[,t] == 2) * (sin(10*t) * X[,t] - sin(10*t - 10) * X[,t-1] + sin(10*t - 20) * X[,t-2])) +
        (b<-(A[,t] == 3) * (cos(10*t) * X[,t] - cos(10*t - 10) * X[,t-1] + sqrt(abs(cos(10*t - 20))) * X[,t-2]))
}

yi <- rnorm(n, mean = mi, sd = 1)