# # MVT
# set.seed(1) #remove later
#
# # Starting Scalar values
# n           <- 5000
# num_treats  <- 3
# p_list      <- rep(2, num_treats)
# num_treats  <- rep(3, num_treats)
#
# Data <- generate_dataset_mvt(n, num_treats, p_list, num_treats)
# # X <- Data[2:(num_treats+1)]
# # A <- Data[[num_treats+2]]
# # y <- Data[[1]]
# res_mvt <- compute_MC_draws_mvt(Data = Data, tau = 0.01, num_treats = 3, B = 100, nu0 = 3,
#                              V0 = diag(2), alph = 3, gam = 4, p_list = 2)



# # UVT
# set.seed(1) #remove later
#
# # Starting Scalar values
# n           <- 100#500#0
# num_treats  <- 2#5
# p_list      <- rep(1, num_treats)
# num_treats  <- rep(3, num_treats)
#
# Data <- generate_dataset_uvt(n, num_treats, num_treats)
# X <- Data[2:(num_treats+1)]
# A <- Data[[num_treats+2]]
# y <- Data[[1]]
#
# tau <- 0.01
# B <- 10000
# alph <- 3
# gam <- 4
# res_uvt <- compute_MC_draws_uvt(Data = Data, tau = 0.01, num_treats = num_treats, B = 10000,
#                                 alph = 3, gam = 4, p_list = rep(1, num_treats))

vec_permutations <- function(max_vals) {
    # Create a list of sequences for each position
    ranges <- lapply(max_vals, function(x) 1:x)

    # Use expand.grid to generate all combinations
    perms <- expand.grid(rev(ranges))

    # Reverse columns back to original order
    perms <- perms[, rev(seq_along(perms))]

    # Convert to matrix or leave as data frame
    perms <- as.matrix(perms)
    dimnames(perms) <- NULL

    return(perms)
}

# GCV
set.seed(1)
numTreats   <- 3
p_list      <- rep(2, numTreats)
At_lens     <- rep(2, numTreats)
n           <- 1000

Dat <- generate_dataset_mvt(n, numTreats, p_list, At_lens)

tau     <- 0.01
B       <- 1000
alph    <- 3
gam     <- 4
nu0     <- 3

Dat_to_500 <- c(list(Dat[[1]][1:500]), lapply(Dat[-1], function(x) x[1:500,]))
res_GCV <- compute_MC_draws_mvt(Data = Dat_to_500, tau = tau, At_lens = At_lens, B = B,
                                nu0 = nu0, alph = alph, gam = gam, p_list = p_list)

# source("~/GitHub/BayesRegDTR/R/GiveChoiceValue.R")
# 
# tic("i = 501:n")
# res_GCV2 <- matrix(0, nrow = 500, ncol = 3)
# for (i in 501:n) {
#     if (i%%10 == 0) print(i)
#     histDat <- list(Dat[[2]][i,,drop = FALSE], Dat[[5]][i,1,drop = FALSE])
#     currDat <- Dat[[3]][i,,drop = FALSE]
#     Wt      <- lapply(res_GCV$Wt_B_list, function(x) x[[B]])
#     Sigmat  <- lapply(res_GCV$sigmat_B_list, function(x) x[[B]])
#     R       <- 30
# 
#     res_GCV2[i-500,] <-
#         GiveChoiceValue(Wt = Wt, Sigmat = Sigmat, bet = res_GCV$beta_B[,B],
#                         sigmay = res_GCV$sigmay_2B[B], t = 2, numTreats = numTreats,
#                         histDat = histDat, currDat = currDat, R = R,
#                         At_lens = At_lens)
# }
# res_GCV2 <- cbind(res_GCV2,0)
# res_GCV2[,3] <- apply(res_GCV2[,1:2], 1, which.max)
# toc(log = TRUE)
# 
# tic("b = 1:B")
# res_GCV3 <- matrix(0, nrow = B, ncol = 3)
# for (b in 1:B) {
#     if (b%%100 == 0) print(b)
#     histDat <- list(Dat[[2]][501,,drop = FALSE], Dat[[5]][501,1,drop = FALSE])
#     currDat <- Dat[[3]][501,,drop = FALSE]
#     Wt      <- lapply(res_GCV$Wt_B_list, function(x) x[[b]])
#     Sigmat  <- lapply(res_GCV$sigmat_B_list, function(x) x[[b]])
#     R       <- 30
# 
#     res_GCV3[b,] <-
#         GiveChoiceValue(Wt = Wt, Sigmat = Sigmat, bet = res_GCV$beta_B[,b],
#                         sigmay = res_GCV$sigmay_2B[b], t = 2, numTreats = numTreats,
#                         histDat = histDat, currDat = currDat, R = R,
#                         At_lens = At_lens)
# }
# res_GCV3 <- cbind(res_GCV3,0)
# res_GCV3[,3] <- apply(res_GCV3[,1:2], 1, which.max)
# toc(log = TRUE)
# 
# 
