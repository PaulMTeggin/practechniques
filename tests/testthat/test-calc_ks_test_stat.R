test_that("rnorm", {
  # Arrange
  set.seed(1)
  x <- rnorm(100)
  fn_p <- function(x, params) pnorm(x, params$mean, params$sd)
  expected <- unname(stats::ks.test(x, "pnorm", mean(x), sd(x))$statistic)

  # Act
  actual <- calc_ks_test_stat(x, list(mean = mean(x), sd = sd(x)), fn_p)

  # Assert
  expect_equal(actual, expected)
})

test_that("rt", {
  # Arrange
  set.seed(1)
  x <- rnorm(100)
  df <- 5
  sd_fac <- df / (df - 2)
  fn_p <- function(x, params) pt((x - params$mean) / (params$sd * sd_fac), df = params$df)
  expected <- unname(stats::ks.test((x - mean(x)) / (sd(x) * sd_fac), "pt", df = 5)$statistic)

  # Act
  actual <- calc_ks_test_stat(x, list(mean = mean(x), sd = sd(x), df = 5), fn_p)

  # Assert
  expect_equal(actual, expected)
})

test_that("rlogis", {
  # Arrange
  set.seed(1)
  x <- rlogis(100)
  sd_fac <- sqrt(3) / pi
  fn_p <- function(x, params) plogis(x, location = params$location, scale = params$scale)
  expected <- unname(stats::ks.test(x, "plogis", location = mean(x), scale = sd(x) * sd_fac)$statistic)

  # Act
  actual <- calc_ks_test_stat(x, list(location = mean(x), scale = sd(x) * sd_fac), fn_p)

  # Assert
  expect_equal(actual, expected)
})
