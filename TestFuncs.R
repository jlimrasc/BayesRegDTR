# Parallel testing
test <- function(){
    options(progress_enabled = NULL)
    # cli::cli_progress_bar(name = "Monte-Carlo draws", status = "stat?", type = "iterator", total = 100)
    # p <- progressr::progressor(steps = 100,)
    pb <- progress_bar$new( format = " Computing MC Draws [:bar] :percent eta: :eta", total = 20)
    x<-function(){
    for (i in 1:20) {
        Sys.sleep(5/100)
        # cli::cli_progress_update()
        # p(sprintf("step %d", i))
        pb$tick()
    }}
    x()
    # cli_progress_done()
    pb$terminate()
}


test()

test <- function(){
    results <- foreach(i = 1:20) %dopar% {
        Sys.sleep(0.2)   # simulate work
        i^2               # your computation result
    }
}
backend <- parabar::start_backend(cores = 4, cluster_type = "psock", backend_type = "async")

# Register it with the `foreach` package.
doParabar::registerDoParabar(backend)

# Stop the backend.
parabar::stop_backend(backend)

future::plan(future::multisession)  # or plan(multicore) on Unix

doFuture::registerDoFuture()  # register future backend with foreach
test <- function(){
    test2<-function(i) i^2
    progressr::with_progress({
        p <- progressr::progressor(steps = 10)
        for(i in 1:10) {
            p(message = paste("Step", i, "of 10"))
            Sys.sleep(0.5)
        }
        p <- progressr::progressor(steps = 100, message = "starting", enable = TRUE)  # initialize the progressor

        results <- foreach(i = 1:100) %dopar% {
            p()  # update progress
            Sys.sleep(0.05)  # simulate some work

        }
    })
}
progressr::handlers("cli")
test()
test()


library(doFuture)    ## %dofuture%
plan(multisession, workers = 2)

library(progressr)
handlers(global = TRUE)
handlers("cli")

my_fcn <- function(xs) {
    p <- progressor(along = xs)
    foreach(x = xs) %dofuture% {
        Sys.sleep((10.0-x)/2)
        p(sprintf("x=%g", x))
        sqrt(x)
    }
}

y <- my_fcn(1:10)

# MVT
set.seed(1) #remove later

# Starting Scalar values
n           <- 5000
num_stages  <- 4
p_list      <- c(2, 2, rep(1, num_stages-2))
num_treats  <- rep(3, num_stages)

Data <- generate_dataset_mvt(n, num_stages, p_list, num_treats)
res_mvt <- compute_MC_draws_mvt(Data = Data, tau = 0.01, num_treats = num_treats, B = 500, nu0 = 3,
                                V0 = mapply(diag, p_list, SIMPLIFY = FALSE), alph = 3, gam = 4, p_list = p_list)

set.seed(1) #remove later
Data2 <- generate_dataset(n, num_stages, p_list, num_treats)
# X <- Data[2:(num_stages+1)]
# A <- Data[[num_stages+2]]
# y <- Data[[1]]
res_mvt2 <- compute_MC_draws_mvt(Data = Data2, tau = 0.01, num_treats = num_treats, B = 100, nu0 = 3,
                             V0 = mapply(diag, p_list, SIMPLIFY = FALSE), alph = 3, gam = 4, p_list = p_list)
equality <- c()
equality <- c(equality, all(unlist(Data) == unlist(Data2)), all(unlist(res_mvt) == unlist(res_mvt2)))

t<-3
b<-1
i<-1
Dat.pred  <- generate_dataset_mvt(n,  num_stages, p_list, num_treats)
Dat.pred  <- Dat.pred[-1]
Dat.pred[[num_stages+1]]  <- Dat.pred[[num_stages+1]][1:n, 1:(t-1), drop = FALSE]
histDat <- c(lapply(Dat.pred[1:(t-1)], function(x) x[i,,drop = FALSE]), list(Dat.pred[[num_stages + 1]][i,1:(t-1),drop = FALSE]))
currDat <- Dat.pred[[t]][i,,drop = FALSE]
Wt      <- lapply(res_mvt$Wt_B_list, function(x) x[[b]])
Sigmat  <- lapply(res_mvt$sigmat_B_list, function(x) x[[b]])

temp <-
    GiveChoiceValue(Wt = Wt, Sigmat = Sigmat, bet = res_mvt$beta_B[,b],
                    sigmay = res_mvt$sigmay_2B[b], t = t,
                    num_stages = num_stages, p_list = p_list,
                    histDat = histDat, currDat = currDat, R = 30,
                    num_treats = num_treats)

rm(list = setdiff(ls(), "equality"))

# UVT
set.seed(1) #remove later

# Starting Scalar values
n           <- 5000
num_stages  <- 5
p_list      <- rep(1, num_stages)
num_treats  <- rep(3, num_stages)

Data <- generate_dataset_uvt(n, num_stages, num_treats)
# X <- Data[2:(num_stages+1)]
# A <- Data[[num_stages+2]]
# y <- Data[[1]]
#
tau <- 0.01
B <- 100
alph <- 3
gam <- 4
res_uvt <- compute_MC_draws_uvt(Data = Data, tau = tau, num_treats = num_treats, B = B,
                                alph = alph, gam = gam, p_list = p_list)

set.seed(1) #remove later
Data2 <- generate_dataset(n, num_stages, p_list, num_treats)
res_uvt2 <- compute_MC_draws_uvt(Data = Data2, tau = tau, num_treats = num_treats, B = B,
                                alph = alph, gam = gam, p_list = p_list)

equality <- c(equality, all(unlist(Data) == unlist(Data2)), all(unlist(res_uvt) == unlist(res_uvt2)))
rm(list = setdiff(ls(), "equality"))

# UVT GCV Test
n           <- 5000
num_stages  <- 5
p_list      <- rep(1, num_stages)
num_treats  <- rep(2, num_stages)
tau <- 0.01
B <- 100
alph <- 3
gam <- 4
t<-3
b<-1
i<-1
tic()
Dat.train  <- generate_dataset(n,  num_stages, p_list, num_treats)
res_uvt <- compute_MC_draws_uvt(Data = Dat.train, tau = tau, num_treats = num_treats, B = B,
                                alph = alph, gam = gam, p_list = p_list, showBar = FALSE)
toc()
Dat.pred  <- generate_dataset(n,  num_stages, p_list, num_treats)
Dat.pred  <- Dat.pred[-1]
Dat.pred[[num_stages+1]]  <- Dat.pred[[num_stages+1]][1:n, 1:(t-1), drop = FALSE]
histDat <- c(lapply(Dat.pred[1:(t-1)], function(x) x[i,,drop = FALSE]), list(Dat.pred[[num_stages + 1]][i,1:(t-1),drop = FALSE]))
currDat <- Dat.pred[[t]][i,,drop = FALSE]
thetat  <- lapply(res_uvt$thetat_B_list, function(x) matrix(x[,b]))
Sigmat  <- lapply(res_uvt$sigmat_2B_list, function(x) matrix(x[b]))

tic()
gcv_uvt <-
    GiveChoiceValue(Wt = thetat, Sigmat = Sigmat, bet = res_uvt$beta_B[,b],
                    sigmay = res_uvt$sigmay_2B[b], t = t,
                    num_stages = num_stages, p_list = p_list,
                    histDat = histDat, currDat = currDat, R = 10,
                    num_treats = num_treats)
toc()



## BayesLinRegDTR.model.fit Tests
future::plan(future::multisession)  # or plan(multicore) on Unix
doFuture::registerDoFuture()
progressr::handlers("cli")

## UVT
num_stages  <- 5
t           <- 3
p_list      <- rep(1, num_stages)
num_treats  <- rep(3, num_stages)
n.train     <- 5000
n.pred      <- 15
numCores    <- parallel::detectCores()

Dat.train  <- generate_dataset(n.train,  num_stages, p_list, num_treats)
Dat.pred  <- generate_dataset(n.pred,  num_stages, p_list, num_treats)
Dat.pred  <- Dat.pred[-1]
Dat.pred[[num_stages+1]]  <- Dat.pred[[num_stages+1]][1:n.pred, 1:(t-1), drop = FALSE]

gcv_uvt <- BayesLinRegDTR.model.fit(Dat.train, Dat.pred, n.train, n.pred,
                                    num_stages, num_treats,
                                    p_list, t, R = 30,
                                    tau = 0.01, B = 500, nu0 = NULL,
                                    V0 = NULL, alph = 3, gam = 4)

## MVT
num_stages  <- 3
t           <- 2
p_list      <- rep(2, num_stages)
num_treats  <- rep(2, num_stages)
n.train     <- 5000
n.pred      <- 10
numCores    <- parallel::detectCores()

Dat.train <- generate_dataset_mvt(n.train, num_stages, p_list, num_treats)
Dat.pred  <- generate_dataset_mvt(n.pred,  num_stages, p_list, num_treats)
Dat.pred  <- Dat.pred[-1]
Dat.pred[[num_stages+1]]  <- Dat.pred[[num_stages+1]][1:n.pred, 1:(t-1), drop = FALSE]

future::plan(future::multisession)  # or plan(multicore) on Unix
doFuture::registerDoFuture()
# progressr::handlers("cli")
doParallel::registerDoParallel(numCores)  # use multicore, set to the number of our cores
gcv_res <- BayesLinRegDTR.model.fit(Dat.train, Dat.pred, n.train, n.pred,
                                    num_stages, num_treats,
                                    p_list, t, R = 30,
                                    tau = 0.01, B = 500, nu0 = 3,
                                    V0 = mapply(diag, p_list, SIMPLIFY = FALSE),
                                    alph = 3, gam = 4)

# vec_permutations <- function(max_vals) {
#     # Create a list of sequences for each position
#     ranges <- lapply(max_vals, function(x) 1:x)
#
#     # Use expand.grid to generate all combinations
#     perms <- expand.grid(rev(ranges))
#
#     # Reverse columns back to original order
#     perms <- perms[, rev(seq_along(perms))]
#
#     # Convert to matrix or leave as data frame
#     perms <- as.matrix(perms)
#     dimnames(perms) <- NULL
#
#     return(perms)
# }

# GCV
# set.seed(1)
# num_stages  <- 3
# t           <- 2
# p_list      <- rep(2, num_stages)
# num_treats  <- rep(2, num_stages)
# n.train     <- 500
# n.pred      <- 100
#
# Dat.train <- generate_dataset_mvt(n.train, num_stages, p_list, num_treats)
# Dat.pred  <- generate_dataset_mvt(n.pred,  num_stages, p_list, num_treats)
# Dat.pred  <- Dat.pred[-1]
# Dat.pred[[num_stages+1]]  <- Dat.pred[[num_stages+1]][1:n.pred, 1:(t-1), drop = FALSE]
#
# # tau     <- 0.01
# # B       <- 5
# # alph    <- 3
# # gam     <- 4
# # nu0     <- 3
# # V0      <- mapply(diag, p_list, SIMPLIFY = FALSE)
# # R       <- 30
# # numCores<- parallel::detectCores()
#
# gcv_res <- BayesLinRegDTR.model.fit(Dat.train, Dat.pred, n.train, n.pred, num_stages, num_treats, p_list, t)
# , R, numCores,
#                            tau, B, nu0, V0, alph, gam)

# Dat_to_500 <- c(list(Dat[[1]][1:500]), lapply(Dat[-1], function(x) x[1:500,]))
# res_GCV <- compute_MC_draws_mvt(Data = Dat_to_500, tau = tau, num_treats = num_treats, B = B,
#                                 nu0 = nu0, alph = alph, gam = gam, p_list = p_list)

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
#                         sigmay = res_GCV$sigmay_2B[B], t = 2,
#                         num_stages = num_stages, p_list = p_list
#                         histDat = histDat, currDat = currDat, R = R,
#                         num_treats = num_treats)
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
#                         sigmay = res_GCV$sigmay_2B[b], t = 2,
#                         num_stages = num_stages, p_list = p_list
#                         histDat = histDat, currDat = currDat, R = R,
#                         num_treats = num_treats)
# }
# res_GCV3 <- cbind(res_GCV3,0)
# res_GCV3[,3] <- apply(res_GCV3[,1:2], 1, which.max)
# toc(log = TRUE)
#
#
