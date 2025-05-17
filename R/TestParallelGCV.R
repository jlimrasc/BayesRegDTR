library(foreach)
library(doParallel)
library(mvtnorm)
library(tictoc)
source("~/Github/BayesRegDTR/R/GiveChoiceValue.R")

numCores <- detectCores()
registerDoParallel(numCores)  # use multicore, set to the number of our cores

inner_b_GCV <- function(i) {
    res_GCV3 <- matrix(0, nrow = B, ncol = 2)
    for (b in 1:B) {
        histDat <- list(Dat[[2]][i,,drop = FALSE], Dat[[5]][i,1,drop = FALSE])
        currDat <- Dat[[3]][i,,drop = FALSE]
        Wt      <- lapply(res_GCV$Wt_B_list, function(x) x[[b]])
        Sigmat  <- lapply(res_GCV$sigmat_B_list, function(x) x[[b]])
        R       <- 30

        res_GCV3[b,1:p_list[1]] <-
            GiveChoiceValue(Wt = Wt, Sigmat = Sigmat, bet = res_GCV$beta_B[,b],
                            sigmay = res_GCV$sigmay_2B[b], t = 2, numTreats = numTreats,
                            histDat = histDat, currDat = currDat, R = R,
                            At_lens = At_lens)
    }
    return(res_GCV3)
}
# res_GCV3[,3] <- apply(res_GCV3[,1:2], 1, which.max)

fori <- n
tic("b = 1:B")
result <- foreach(i=501:fori, .inorder = TRUE, .packages = "mvtnorm") %dopar% inner_b_GCV(i)
toc(log = TRUE)

result <- array(unlist(result), dim = c(B, p_list[1], (fori - 500)))

# result <- array(result, dim = c(B, p_list[1], (fori-500)))

# res_GCV3 <- cbind(res_GCV3,0)
# res_GCV3[,3] <- apply(res_GCV3[,1:2], 1, which.max)
freqs <- matrix(0, ncol = 2, nrow = 500)
for (i in 1:500){
         temp <- apply(result[,1:2,i], 1, which.max)
         freqs[i,1] <- sum(temp==1)
         freqs[i,2] <- sum(temp==2)
}
all(apply(freqs, 1, sum)==1000)