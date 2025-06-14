test_that("Univariate output sizes expected", {
    n           <- 10
    num_stages  <- 3
    p_list      <- rep(1, num_stages)
    num_treats  <- rep(2, num_stages)

    Data        <- generate_dataset(n, num_stages, p_list, num_treats)
    expect_length(Data, num_stages + 2)
    expect_length(Data[[1]], n)
    expect_equal(dim(Data[[2]]), c(n, p_list[1]))
    expect_equal(dim(Data[[num_stages+2]]), c(n, num_stages))
})

test_that("Multivariate output sizes expected", {
    n           <- 10
    num_stages  <- 3
    p_list      <- c(rep(2, num_stages-1), 1)
    num_treats  <- rep(2, num_stages)

    Data        <- generate_dataset(n, num_stages, p_list, num_treats)
    expect_length(Data, num_stages + 2)
    expect_length(Data[[1]], n)
    expect_equal(dim(Data[[2]]), c(n, p_list[1]))
    expect_equal(dim(Data[[num_stages+1]]), c(n, p_list[num_stages]))
    expect_equal(dim(Data[[num_stages+2]]), c(n, num_stages))
})
