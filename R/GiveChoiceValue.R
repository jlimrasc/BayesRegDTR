GiveChoiceValue <- function(Wt, Sigmat, bet, sigmay, t, num_stages, p_list, histDat, currDat, R, num_treats) {
    stopifnot(length(histDat) == t)

    histDatX <- histDat[1:(t-1)]
    histDatA <- histDat[[t]]

    stopifnot("length of historical data X must be t-1" = length(histDatX) == t-1,
              "length of historical data A must be t-1" = NCOL(histDatA) == t-1,
              "length of Wt must be num_stages-1" = length(Wt) == num_stages-1,
              "length of Sigmat must be num_stages-1" = length(Sigmat) == num_stages-1
              )

    if (t == num_stages) {
        n       <- NROW(histDatA)
        X_T     <- matrix(c(unlist(histDatX), currDat), nrow = n)

        # Calculate permutations of a1, ..., an
        perms   <- vec_permutations(num_treats[1:t])

        p_sum   <- sum(p_list[1:t])

        # Preallocations
        Q_T     <- matrix(0, nrow = n, ncol = p_sum * nrow(perms))
        Q_T_vec <- rep(0, num_treats[t]) # Single vector of Q_T

        for (a in 1:num_treats[t]) {
            A_T     <- cbind(histDatA, a) # Combine historical A with current a

            for (i in 1:nrow(perms)) {
                # Starting index
                colnum <- (i-1) * p_sum + 1

                # Check for treatments that match specific perm
                selector <- apply(A_T == rep(perms[i,], each=n), 1, all)

                # Dot product selection
                Q_T[,colnum:(colnum + p_sum - 1)] <- rep(selector, each=n) * X_T
            }
            Q_T_vec[a] <- sum(t(Q_T) * bet) # Dot product
        }

        return(Q_T_vec)
    }

    else {
        n       <- NROW(histDatX[[1]])
        X_t     <- matrix(c(unlist(histDatX), currDat), nrow = n) # Combine X from 1 to t into a single matrix

        # Calculate permutations of a1, ..., an
        perms   <- vec_permutations(num_treats[1:t])

        p_sum   <- sum(p_list[1:t])
        x_t1ra_mid  <- matrix(0, nrow = n, ncol = p_sum * nrow(perms)) # Preallocate

        X_ta    <- vector(mode = 'list', length = num_treats[t])
        for (a in 1:num_treats[t]) {
            X_ta[[a]] <- rep(0, R) # Preallocate
            A_t     <- cbind(histDatA, a)

            for (r in 1:R){

                for (i in 1:nrow(perms)) {
                    #Starting index
                    colnum <- (i-1) * p_sum + 1

                    # Check for treatments that match specific perm
                    selector <- apply(A_t == rep(perms[i,], each=n), 1, all)

                    # Dot product selection
                    x_t1ra_mid[,colnum:(colnum + p_sum - 1)] <- rep(selector, each=n) * X_t
                }

                # Compute x^r_{t+1}(a)
                x_t1ra  <- t(Wt[[t]]) %*% t(x_t1ra_mid) + t(mvtnorm::rmvnorm(1, sigma = Sigmat[[t]]))

                X_ta[[a]][r] <- max(
                    GiveChoiceValue(
                        Wt=Wt, Sigmat=Sigmat, bet=bet, sigmay=sigmay,
                        t=t+1, num_stages=num_stages, p_list = p_list,
                        histDat=c(histDat[-t],list(currDat), list(matrix(c(histDat[[t]], a), ncol = t))),
                        currDat = x_t1ra, R=r, num_treats = num_treats)
                    )

            }
        }
    }
    return(sapply(X_ta, mean)) # For each treatment, find the mean of R predictions
}
