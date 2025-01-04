set.seed(1)
# Starting Scalar values
n       <- 5000
T       <- 5
p_t     <- 1
At_max  <- 3

# Step 1
x_i1 <- rt(n, 10)

# Step 2
A <- matrix(sample(1:At_max, n*T, replace = TRUE), nrow = n, ncol = T)

# Step 3
X <- matrix(0, nrow = n, ncol = T) # Preallocate

# Run first 3 lines fist bc can't be bothered making checks for numeric(0)s
X[,1] <- x_i1
t <- 2
X[,t] <- rnorm(n, mean = 0.5^2) +
    (A[,t-1] == 2) * (t * X[,t-1]) + # If a==2
    (A[,t-1] == 3) * (-t * X[,t-1]) # Elif a==3
t <- 3
X[,t] <- rmvnorm(n, mean = 0.5^2) +
    (A[,t-1] == 2) * (t * X[,t-1] - (t-1) * X[,t-2]) + # If a==2
    (A[,t-1] == 3) * (-t * X[,t-1] + sqrt(t-1) * X[,t-2]) # Elif a==3

for (t in 4:T) {
    xi_t <- rmvnorm(n, mean = 0.5^2)
    # Could use dplyr or matrix operations instead of ifelse
    X[,t] <- xi_t +
        (A[,t-1] == 2) * (t * X[,t-1] - (t-1) * X[,t-2] + max(t-2, 0) * X[,t-3]) + # If a==2
        (A[,t-1] == 3) * (-t * X[,t-1] + sqrt(t-1) * X[,t-2] + sqrt(max(t-2, 0)) * X[,t-3]) # Elif a==3
}
# gen_xit <- function(a, x_i, t) {
#     return(
#         ifelse (a == 2, t * x_i[t-1] - max(t-1, 0) * x_i[t-2] + max(t-2, 0) * x_i[t-3], # If a==2
#         ifelse (a == 3, -t * x_i[t-1] + sqrt(max(t-1, 0)) * x_i[t-2] + sqrt(max(t-2, 0)) * x_i[t-3]), # Elif a==3
#                 0) # Else
#     )
# }

# Step 4
sumcum <- 0
t <- 1
sumcum <- sumcum + (A[,t] == 2) * (sin(10*t) * X[,t])
sumcum <- sumcum + (A[,t] == 3) * (cos(10*t) * X[,t])
t <- 2
sumcum <- sumcum + (A[,t] == 2) * (sin(10*t) * X[,t] - sin(10*t - 10) * X[,t-1])
sumcum <- sumcum + (A[,t] == 3) * (cos(10*t) * X[,t] - cos(10*t - 10) * X[,t-1])
for (t in 3:T) {
    sumcum <- sumcum + (A[,t] == 2) * (sin(10*t) * X[,t] - sin(10*t - 10) * X[,t-1] + sin(10*t - 20) * X[,t-2])
    sumcum <- sumcum + (A[,t] == 3) * (cos(10*t) * X[,t] - cos(10*t - 10) * X[,t-1] + cos(10*t - 20) * X[,t-2])
}
mi <- 3 + sumcum

yi <- rnorm(n, mean = mi, sd = 1)