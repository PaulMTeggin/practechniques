test_that("rnorm no CI", {
  # Arrange
  nreps_outer <- 50
  nreps_inner <- 100
  x_mtx <- matrix(rnorm(nreps_outer * nreps_inner), nrow = nreps_outer)

  # Act
  actual_expected <- vapply(seq_len(nreps_outer), function(i) {
    set.seed(1)
    actual <- gof_test_sim(x_mtx[i,], "KS", dist = "norm", nreps = 999)$p_value

    set.seed(1)
    expected <- KScorrect::LcKS(x_mtx[i,], "pnorm", nreps = 999)$p.value

    return(c(actual, expected))
  }, numeric(2))


  # Assert
  expect_equal(actual_expected[1,], actual_expected[2,])
})

test_that("rnorm with CI - CI doesn't change p values", {
  # Arrange
  nreps_outer <- 50
  nreps_inner <- 100
  x_mtx <- matrix(rnorm(nreps_outer * nreps_inner), nrow = nreps_outer)

  # Act
  actual_expected <- vapply(seq_len(nreps_outer), function(i) {
    set.seed(1)
    actual <- gof_test_sim(x_mtx[i,], "KS", dist = "norm", nreps = 999, bs_ci = 0.95)$p_value

    set.seed(1)
    expected <- KScorrect::LcKS(x_mtx[i,], "pnorm", nreps = 999)$p.value

    return(c(actual, expected))
  }, numeric(2))


  # Assert
  expect_equal(actual_expected[1,], actual_expected[2,])
})

