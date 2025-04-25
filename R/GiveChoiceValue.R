GiveChoiceValue <- function(Wt, Sigmat, bet, sigmay2, t, numTreats, histDat, currDat, R) {
    stopifnot(length(histDat) == (t+1))

    histDatX <- histDat[1:(t-1)]
    histDatA <- histDat[[t]]

    stopifnot("length of historical data X must be t-1" = length(histDatX) == t-1,
              "length of historical data A must be t-1" = NCOL(histDatA) == t-1,
              "length of Wt must be numTreats-1" = length(Wt) == numTreats-1,
              "length of Sigmat must be numTreats-1" = length(Sigmat) == numTreats-1
              )


}
