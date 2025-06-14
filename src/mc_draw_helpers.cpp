// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <vector>
#include <iostream>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List draw_Wt_B_cpp(const arma::mat& omegat_inv, const arma::mat& Mnt,
                                 const double tau, std::vector<arma::mat>& sigmatb,
                                 const int pt, const int qt, const int B) {
    Rcpp::List omegR_list(B);
    arma::mat omegat_inv_sqrtm = arma::sqrtmat_sympd(omegat_inv); // Pre-calculate sqrtm
    for (int i = 0; i < B; i++) {
        omegR_list[i] = Mnt + omegat_inv_sqrtm * arma::randn(qt, pt) * arma::sqrtmat_sympd(sigmatb[i]);
    }

    return omegR_list;
}
