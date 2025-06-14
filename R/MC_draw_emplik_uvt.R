# # h_t_wrap <- function(X_i1t, a_i1t, thetat, At_len, p_list, t) {
# #     Z_ti <- compute_Zt(a_i1t, At_len, X_i1t, t, 1, p_list)
# #     return(t(Z_ti) %*% (tail(x_i1t, n = 1) - Z_ti %*% thetat))
# # }
# #
# # h_y_wrap <- function(yi, X_i1T, a_i1t, thetaT, At_len, p_list, t) {
# #     Z_T1i <- compute_Zt(a_i1t, At_len, X_i1T, t, 1, p_list)
# #     return(t(Z_T1i) %*% (tail(x_i1t, n = 1) - Z_T1i %*% thetaT))
# # }
# # compute_Zt <- function(A, At_len, X, t, n, p_list)
# h_t_wrap <- function(At_len, p_list, t) {
#     h_t <- function(data, thetat) {
#         mid     <- length(data) %/% 2 + 1
#         x_it    <- data[mid]
#         X       <- t(head(data, mid))
#         A       <- t(tail(data, mid - 1))
#         Zti <- compute_Zt(A, At_len, X, t, 1, p_list)
#
#         return(t(Zti) %*% (x_it - Zti %*% thetat))
#     }
#
#     return(h_t)
# }
#
# h_ty_wrap <- function(At_len, p_list, t) {
#     h_ty <- function(data, beta) {
#         mid     <- length(data) %/% 2 + 1
#         X       <- t(data[2:mid])
#         A       <- t(tail(data, mid - 1))
#         ZT1i <- compute_Zt(A, At_len, X, t, 1, p_list)
#
#         yi <- data[1]
#         return(t(ZT1i) %*% (yi - ZT1i %*% beta))
#     }
# }
#
# delth_h_wrap <- function(At_len, p_list, t) {
#     delth_h <- function(data, th) {
#         mid     <- length(data) %/% 2 + 1
#         X       <- t(head(data, mid))
#         A       <- t(tail(data, mid - 1))
#         Zti <- compute_Zt(A, At_len, X, t, 1, p_list)
#
#         return(-t(Zti) %*% Zti)
#     }
# }
#
# delbet_hy_wrap <- function(At_len, p_list, t) {
#     delbet_hy <- function(data, th) {
#         mid     <- length(data) %/% 2 + 1
#         X       <- t(data[2:mid])
#         A       <- t(tail(data, mid - 1))
#         ZT1i <- compute_Zt(A, At_len, X, t, 1, p_list)
#
#         return(-t(ZT1i) %*% ZT1i)
#     }
# }
#
# delth_logpi <- function(thetat) {
#     return(-thetat / 100)
# }
#
# delbet_logpi <- function(beta) {
#     return(-beta / 100)
# }
#
# library(VBel)
# X_mat <- matrix(unlist(X), nrow = n, ncol = T)
#
# t <- 2
# tau <- 0.01
#
# Zt          <- compute_Zt(A, At_len, X, t, n, p_list)
# omegat      <- compute_omegat(Zt, tau)
# omegat_inv  <- solve(omegat)
#
# Mnt <- omegat_inv %*% t(Zt) %*% X[[t]]
# qt <- ncol(Zt)
#
# thetat_b_list <- vector(mode = "list", length = T)
#
# # Step 2 for (t in 2:T)
# res1 <- compute_GVA(mu0 = Mnt, C0 = diag(qt), h = h_t_wrap(At_len, p_list, t),
#             delthh = delth_h_wrap(At_len, p_list, t),
#             delth_logpi = delth_logpi, lam0 = rep(0, qt),
#             rho = 0.9, epsil = 1e-6, a = 1e-3,
#             z = cbind(X_mat[,1:t], A[,1:(t - 1)]), verbosity = 100,
#             SGD_iters = 10000, AEL_iters = 500)
#
# thetat_b <- rnorm(B, res$mu_FC, res$C_FC %*% t(res$C_FC))
#
# thetat_b_list[t] <- thetat_b
#
# # Step 3
# t <- 3
# B <- 100
# ZT          <- compute_Zt(A, At_len, X, t, n, p_list)
# omegaT      <- compute_omegat(ZT, tau)
# omegaT_inv  <- solve(omegaT)
#
# MnT <- omegaT_inv %*% t(ZT) %*% y
# qT <- ncol(ZT)
# library(tictoc)
# tic("resy")
# resy <- compute_GVA(mu0 = MnT, C0 = diag(qT), h = h_ty_wrap(At_len, p_list, t),
#             delthh = delbet_hy_wrap(At_len, p_list, t),
#             delth_logpi = delbet_logpi, lam0 = rep(0, qT),
#             rho = 0.9, epsil = 1e-6, a = 1e-3,
#             z = cbind(y, X_mat, A), verbosity = 100,
#             SGD_iters = 10000, AEL_iters = 500)
# toc()
#
# # Step 4
# bet_b <- rnorm(B, resy$mu_FC, resy$C_FC %*% t(resy$C_FC))
# #
# # testfun <- function(th, h, lam0, a, z, iters = 500) {
# #     p <- ncol(z)
# #     n <- nrow(z) + 1
# #     h_sum <- 0
# #     H_Zth <- c()
# #
# #     for (i in 1:(n - 1)) {
# #         zi <- matrix(z[i, ], nrow = p) # Row of z as vertical vector
# #         h_zith <- h(zi, th)
# #
# #         h_sum <- h_sum + h_zith # For h(zn,th)
# #         H_Zth <- rbind(H_Zth, t(h_zith)) # Build up H(Z,th)
# #     }
# #
# #     h_znth <- -a / (n - 1) * h_sum
# #     H_Zth <- rbind(H_Zth, t(h_znth)) # Last row of H is h(zn,th)
# # }
