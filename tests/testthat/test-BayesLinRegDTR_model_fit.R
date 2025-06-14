test_that("Error enforced", {
    n           <- 10
    num_stages  <- 3
    p_list      <- rep(2, num_stages)
    num_treats  <- rep(2, num_stages)

    Data        <- generate_dataset(n, num_stages, p_list, num_treats)
    expect_error(BayesLinRegDTR.model.fit(Data, Data, n, n,
                                          num_stages, num_treats,
                                          p_list, t = 10, R = 30,
                                          tau = 0.01, B = 500, nu0 = 3,
                                          V0 = mapply(diag, p_list, SIMPLIFY = FALSE),
                                          alph = 3, gam = 4),
                 "t must be less than or equal to num_stages")
    expect_error(
        expect_warning(
            BayesLinRegDTR.model.fit(Data, Data, n, n,
                                              num_stages, num_treats,
                                              p_list, t = 2, R = 30,
                                              tau = 0.01, B = 500, nu0 = 3,
                                              V0 = 2,
                                              alph = 3, gam = 4),
            "No parallel backend detected: future plan is 'sequential'. For better performance, consider setting a parallel plan, e.g., plan(multisession)."),
        "Length of num_treats, V0 and p_list must be equal and of length 1 or num_stages")

})
