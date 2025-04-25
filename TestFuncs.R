# MVT
set.seed(1) #remove later

# Starting Scalar values
n           <- 500#0
num_treats  <- 3
p_list      <- rep(2, num_treats)
At_len      <- 3

Data <- generate_dataset_mvt(n, num_treats, p_list, At_len)
# X <- Data[2:(num_treats+1)]
# A <- Data[[num_treats+2]]
# y <- Data[[1]]
res_mvt <- compute_MC_draws_mvt(Data = Data, tau = 0.01, At_lens = 3, B = 100, nu0 = 3,
                             V0 = diag(2), alph = 3, gam = 4, p_list = 2)



# UVT
set.seed(1) #remove later

# Starting Scalar values
n           <- 100#500#0
num_treats  <- 2#5
p_list      <- rep(1, num_treats)
At_len      <- 3

Data <- generate_dataset_uvt(n, num_treats, At_len)
X <- Data[2:(num_treats+1)]
A <- Data[[num_treats+2]]
y <- Data[[1]]

tau <- 0.01
At_lens <- 1
B <- 10000
alph <- 3
gam <- 4
res_uvt <- compute_MC_draws_uvt(Data = Data, tau = 0.01, At_lens = 1, B = 10000,
                                alph = 3, gam = 4, p_list = rep(1, num_treats))

