test_that("KS rexp no CI", {
  # Arrange
  nreps_outer <- 50 # 50 tests that we get precisely the same p-values as KScorrect
  nreps_inner <- 100 # with 100 data points in each test
  set.seed(1)
  x_mtx <- matrix(rexp(nreps_outer * nreps_inner, rate = 2), nrow = nreps_outer)

  # Act
  actual_expected <- vapply(seq_len(nreps_outer), function(i) {
    set.seed(2)
    actual <- gof_test_sim(x_mtx[i,], "KS", dist = "exp", fn_estimate_params = estimate_exp_ol, nreps = 999)$p_value

    set.seed(2)
    expected <- KScorrect::LcKS(x_mtx[i,], "pexp", nreps = 999)$p.value

    return(c(actual, expected))
  }, numeric(2))


  # Assert
  expect_equal(actual_expected[1,], actual_expected[2,])
})

test_that("KS rnorm", {
  # Arrange
  nreps_outer <- 50 # 50 tests that we get precisely the same p-values as KScorrect
  nreps_inner <- 100 # with 100 data points in each test
  set.seed(1)
  x_mtx <- matrix(rnorm(nreps_outer * nreps_inner), nrow = nreps_outer)

  # Act
  actual_expected <- vapply(seq_len(nreps_outer), function(i) {
    set.seed(2)
    actual <- gof_test_sim(x_mtx[i,], "KS", dist = "norm", nreps = 999)$p_value

    set.seed(2)
    expected <- KScorrect::LcKS(x_mtx[i,], "pnorm", nreps = 999)$p.value

    return(c(actual, expected))
  }, numeric(2))


  # Assert
  expect_equal(actual_expected[1,], actual_expected[2,])
})

test_that("KS rnorm - pass function not 'KS' string", {
  # Arrange
  nreps_outer <- 50 # 50 tests that we get precisely the same p-values as KScorrect
  nreps_inner <- 100 # with 100 data points in each test
  set.seed(1)
  x_mtx <- matrix(rnorm(nreps_outer * nreps_inner), nrow = nreps_outer)

  # Act
  actual_expected <- vapply(seq_len(nreps_outer), function(i) {
    set.seed(2)
    actual <- gof_test_sim(x_mtx[i,], calc_ks_test_stat, dist = "norm", nreps = 999)$p_value

    set.seed(2)
    expected <- KScorrect::LcKS(x_mtx[i,], "pnorm", nreps = 999)$p.value

    return(c(actual, expected))
  }, numeric(2))


  # Assert
  expect_equal(actual_expected[1,], actual_expected[2,])
})

test_that("KS rnorm - pass function not 'KS' string", {
  # Arrange - nothing to do

  # Act and Assert
  expect_error(gof_test_sim(1:10, test_type = list(), dist = "norm", nreps = 999), "either be a character variable .* or a function")
})

test_that("AD rnorm", {
  # Arrange
  nreps_outer <- 50 # 50 tests that we get p-values from nortest that are within the 99% confidence intervals
  nreps_inner <- 100 # with 100 data points in each test
  set.seed(1)
  x_mtx <- matrix(rnorm(nreps_outer * nreps_inner), nrow = nreps_outer)

  # Act
  actual_expected <- vapply(seq_len(nreps_outer), function(i) {
    actual <- gof_test_sim(x_mtx[i,], "AD", dist = "norm", nreps = 1999, bs_ci = 0.99)

    expected <- nortest::ad.test(x_mtx[i,])$p.value

    return(c(actual$p_value_lower, actual$p_value, actual$p_value_upper, expected))
  }, numeric(4)) # 1 is lower bound on p-value, 2 is p-value, 3 is upper bound, 4 is expected from nortest

  # Assert - need a small amount of headroom on the downside, via the 0.97
  # - not unexpected given we are working with simulation and bootstrapped confidence intervals
  expect_true(all((actual_expected[1,] * 0.97 < actual_expected[4,]) & (actual_expected[3,] > actual_expected[4,])))
})

test_that("AD rnorm - parallelised (NB 2 cores only otherwise rejected by R CMD check)", {
  # Arrange
  nreps_outer <- 50 # 50 tests that we get p-values from nortest that are within the 99% confidence intervals
  nreps_inner <- 100 # with 100 data points in each test
  set.seed(1)
  x_mtx <- matrix(rnorm(nreps_outer * nreps_inner), nrow = nreps_outer)

  # Act
  actual_expected <- vapply(seq_len(nreps_outer), function(i) {
    actual <- gof_test_sim(x_mtx[i,], "AD", dist = "norm", nreps = 1999, bs_ci = 0.99, parallelise = TRUE, ncores = 2)

    expected <- nortest::ad.test(x_mtx[i,])$p.value

    return(c(actual$p_value_lower, actual$p_value, actual$p_value_upper, expected))
  }, numeric(4)) # 1 is lower bound on p-value, 2 is p-value, 3 is upper bound, 4 is expected from nortest

  # Assert - need a small amount of headroom on the downside, via the 0.97
  # - not unexpected given we are working with simulation and bootstrapped confidence intervals
  expect_true(all((actual_expected[1,] * 0.97 < actual_expected[4,]) & (actual_expected[3,] > actual_expected[4,])))
})

test_that("AD rnorm - parallelised values are reproducible", {
  # Arrange
  set.seed(1)
  x <- rnorm(100)

  # Act
  set.seed(2)
  actual_1 <- gof_test_sim(x, "AD", dist = "norm", nreps = 1999, bs_ci = 0.99, parallelise = TRUE, ncores = 2)

  set.seed(2)
  actual_2 <- gof_test_sim(x, "AD", dist = "norm", nreps = 1999, bs_ci = 0.99, parallelise = TRUE, ncores = 2)

  # Assert
  expect_equal(actual_1, actual_2)
})

