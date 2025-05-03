# # MVT
# set.seed(1) #remove later
#
# # Starting Scalar values
# n           <- 5000
# num_treats  <- 3
# p_list      <- rep(2, num_treats)
# At_lens     <- rep(3, num_treats)
#
# Data <- generate_dataset_mvt(n, num_treats, p_list, At_lens)
# # X <- Data[2:(num_treats+1)]
# # A <- Data[[num_treats+2]]
# # y <- Data[[1]]
# res_mvt <- compute_MC_draws_mvt(Data = Data, tau = 0.01, At_lens = 3, B = 100, nu0 = 3,
#                              V0 = diag(2), alph = 3, gam = 4, p_list = 2)



# # UVT
# set.seed(1) #remove later
#
# # Starting Scalar values
# n           <- 100#500#0
# num_treats  <- 2#5
# p_list      <- rep(1, num_treats)
# At_lens     <- rep(3, num_treats)
#
# Data <- generate_dataset_uvt(n, num_treats, At_lens)
# X <- Data[2:(num_treats+1)]
# A <- Data[[num_treats+2]]
# y <- Data[[1]]
#
# tau <- 0.01
# B <- 10000
# alph <- 3
# gam <- 4
# res_uvt <- compute_MC_draws_uvt(Data = Data, tau = 0.01, At_lens = At_lens, B = 10000,
#                                 alph = 3, gam = 4, p_list = rep(1, num_treats))


# GCV
set.seed(1)
numTreats   <- 3
p_list      <- rep(2, numTreats)
At_lens     <- rep(2, numTreats)
n           <- 501

Dat <- generate_dataset_mvt(n, numTreats, p_list, At_lens)

tau     <- 0.01
B       <- 5
alph    <- 3
gam     <- 4
nu0     <- 3

res_GCV <- compute_MC_draws_mvt(Data = Dat, tau = tau, At_lens = At_lens, B = B,
                                nu0 = nu0, alph = alph, gam = gam, p_list = p_list)

# histDat <- c(lapply(Dat[2:(numTreats)], function(x) x[,-p_list[1]]), list(Dat[[numTreats+1]][,-p_list[1]]))
# currDat <- lapply(Dat[2:(numTreats)], function(x) x[,p_list[1]])
histDat <- list(Dat[[2]][501,,drop = FALSE], Dat[[5]][501,1,drop = FALSE])
currDat <- Dat[[3]][501,,drop = FALSE]
Wt      <- lapply(res_GCV$Wt_B_list, function(x) x[[B]])
Sigmat  <- lapply(res_GCV$sigmat_B_list, function(x) x[[B]])
R       <- 30
res_GCV2<- GiveChoiceValue(Wt = Wt, Sigmat = Sigmat, bet = res_GCV$beta_B[,B],
                           sigmay = res_GCV$sigmay_2B[B], t = 2, numTreats = numTreats,
                           histDat = histDat, currDat = currDat, R = R,
                           At_lens = At_lens)
