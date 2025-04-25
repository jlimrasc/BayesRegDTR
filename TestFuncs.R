# # MVT
# set.seed(1) #remove later
#
# # Starting Scalar values
# n           <- 500#0
# num_treats  <- 3
# p_list      <- rep(2, num_treats)
# At_len      <- 3
#
# Data <- generate_dataset_mvt(n, num_treats, p_list, At_len)
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
# At_len      <- 3
#
# Data <- generate_dataset_uvt(n, num_treats, At_len)
# X <- Data[2:(num_treats+1)]
# A <- Data[[num_treats+2]]
# y <- Data[[1]]
#
# tau <- 0.01
# At_lens <- 1
# B <- 10000
# alph <- 3
# gam <- 4
# res_uvt <- compute_MC_draws_uvt(Data = Data, tau = 0.01, At_lens = 1, B = 10000,
#                                 alph = 3, gam = 4, p_list = rep(1, num_treats))

# GCV
set.seed(1)
numTreats   <- 3
p_list      <- rep(2, numTreats)
At_len      <- 2
n           <- 501

Dat <- generate_dataset_mvt(n, numTreats, p_list, At_len)

tau     <- 0.01
B       <- 5
alph    <- 3
gam     <- 4
nu0     <- 3

res_GCV <- compute_MC_draws_mvt(Data=Dat, tau=tau, At_lens=At_len, B=B, nu0=nu0, alph=alph, gam=gam, p_list=p_list)

histDat <- c(lapply(Dat[2:(numTreats)], function(x) x[,-p_list[1]]), list(Dat[[numTreats+1]][,-p_list[1]]))
currDat <- lapply(Dat[2:(numTreats)], function(x) x[,p_list[1]])
Wt      <- lapply(res_GCV$Wt_B_list, function(x) x[[B]])
Sigmat  <- lapply(res_GCV$sigmat_B_list, function(x) x[[B]])
R       <- 30
GiveChoiceValue(Wt=Wt, Sigmat=Sigmat, bet=res_GCV$beta_B[,B], sigmay=res_GCV$sigmay_2B[B],
                t=2, numTreats=numTreats, histDat=histDat, currDat = currDat, R=R)
