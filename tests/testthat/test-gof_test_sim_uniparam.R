test_that("Estimation fails, get NA back", {
  # Arrange
  fn_estimate_params <- function(x, noverlap = 1) {
    stop("Error for test purposes")
  }
  fn_p <- function(x, params) pnorm(x, params$mean, params$sd)
  fn_test_statistic <- function(x, est_params) calc_ks_test_stat(x, est_params, fn_p)
  fn_simulate <- function(N, est_params) stats::rnorm(N, est_params$mean, est_params$sd)

  # Act
  actual <- gof_test_sim_uniparam(rnorm(100), fn_estimate_params, fn_test_statistic, fn_simulate)

  # Assert
  expect_true(is.na(actual$ts))
  expect_true(is.na(actual$p_value))
  expect_true(is.na(actual$count_NA))
  expect_true(is.na(actual$p_value_lower))
  expect_true(is.na(actual$p_value_upper))
})

test_that("error with parallelise / ncores validation", {
  # Arrange
  fn_p <- function(x, params) pnorm(x, params$mean, params$sd)
  fn_test_statistic <- function(x, est_params) calc_ks_test_stat(x, est_params, fn_p)
  fn_simulate <- function(N, est_params) stats::rnorm(N, est_params$mean, est_params$sd)

  # Act and Assert
  expect_error(gof_test_sim_uniparam(rnorm(100), estimate_mean_sd_ol, fn_test_statistic, fn_simulate,
                                     parallelise = FALSE, ncores = 2),
               "multiple cores but not to run in parallel")

  expect_error(gof_test_sim_uniparam(rnorm(100), estimate_mean_sd_ol, fn_test_statistic, fn_simulate,
                                     parallelise = TRUE, ncores = 1),
               "ncores to an integer greater than 1")

  expect_error(gof_test_sim_uniparam(rnorm(100), estimate_mean_sd_ol, fn_test_statistic, fn_simulate,
                                     parallelise = TRUE, ncores = list()),
               "ncores to an integer greater than 1")

  skip("R CMD check rejects more than 2 cores")
  expect_warning(gof_test_sim_uniparam(rnorm(100), estimate_mean_sd_ol, fn_test_statistic, fn_simulate,
                                       parallelise = TRUE, ncores = parallel::detectCores() + 1),
                 "more cores than your computer contains")
})

