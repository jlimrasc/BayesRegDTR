GiveChoiceValue <- function(Wt, Sigmat, bet, sigmay2, t, numTreats, histDat, currDat, R, At_lens) {
    stopifnot(length(histDat) == t)

    histDatX <- histDat[1:(t-1)]
    histDatA <- histDat[[t]]

    stopifnot("length of historical data X must be t-1" = length(histDatX) == t-1,
              "length of historical data A must be t-1" = NCOL(histDatA) == t-1,
              "length of Wt must be numTreats-1" = length(Wt) == numTreats-1,
              "length of Sigmat must be numTreats-1" = length(Sigmat) == numTreats-1
              )

    if (t == numTreats) {
        n       <- NROW(histDatA)
        X_T     <- matrix(c(unlist(histDatX), currDat), nrow = n)

        # Calculate permutations of a1, ..., an
        perms   <- vec_permutations(At_lens[1:t])

        p_sum   <- sum(p_list[1:t])

        # Preallocations
        Q_T     <- matrix(0, nrow = n, ncol = p_sum * nrow(perms))
        Q_T_vec <- rep(0, At_lens[t])

        for (a in 1:At_lens[t]) {
            A_T     <- cbind(histDatA, a)

            for (i in 1:nrow(perms)) {
                colnum <- (i-1) * p_sum + 1

                # Check for treatments that match specific perm
                selector <- apply(A_T == rep(perms[i,], each=n), 1, all)

                Q_T[,colnum:(colnum + p_sum - 1)] <- rep(selector, each=n) * X_T
            }
            Q_T_vec[a] <- sum(t(Q_T) * bet) # Dot product
        }

        return(Q_T_vec)
    }

    else {
        n       <- NROW(histDatX[[1]])
        X_t     <- matrix(c(unlist(histDatX), currDat), nrow = n)

        # Calculate permutations of a1, ..., an
        perms   <- vec_permutations(At_lens[1:t])

        p_sum   <- sum(p_list[1:t])
        x_t1ra_mid  <- matrix(0, nrow = n, ncol = p_sum * nrow(perms))
        # x_t1ra_l<- vector(mode="list", length = At_len[t])

        X_ta    <- vector(mode = 'list', length = At_lens[t])
        for (a in 1:At_lens[t]) {
            X_ta[[a]] <- rep(0, R) # Preallocate
            A_t     <- cbind(histDatA, a)

            for (r in 1:R){

                for (i in 1:nrow(perms)) {
                    colnum <- (i-1) * p_sum + 1

                    # Check for treatments that match specific perm
                    selector <- apply(A_t == rep(perms[i,], each=n), 1, all)

                    x_t1ra_mid[,colnum:(colnum + p_sum - 1)] <- rep(selector, each=n) * X_t
                }

                x_t1ra  <- t(Wt[[t]]) %*% t(x_t1ra_mid) + t(rmvnorm(1, sigma = Sigmat[[t]]))

                X_ta[[a]][r] <- max(
                    GiveChoiceValue(
                        Wt=Wt, Sigmat=Sigmat, bet=bet, sigmay=sigmay,
                        t=t+1, numTreats=numTreats,
                        histDat=c(histDat[-t],list(currDat), list(matrix(c(histDat[t], a), ncol = t))),
                        currDat = x_t1ra, R=r, At_lens = At_lens)
                    )

            }
        }
    }
    return(sapply(X_ta, mean))
}
