test_that("rnorm", {
  # Arrange
  set.seed(1)
  fn_p <- function(x, params) pnorm(x, params$mean, params$sd)

  # Act - compare 1000 KS statistics from the code in this package to the value produced by ADGofTest::ad.test
  actual_expected <- replicate(1000, {
    x <- rnorm(100)
    expected <- unname(ADGofTest::ad.test(x, pnorm, mean(x), sd(x))$statistic)
    actual <- calc_ad_test_stat(x, list(mean = mean(x), sd = sd(x)), fn_p)
    return(c(actual, expected))
  })

  # Assert
  expect_equal(actual_expected[1,], actual_expected[2,])
})

test_that("rt(5)", {
  # Arrange
  set.seed(1)
  df <- 5
  sd_fac <- df / (df - 2)
  fn_p <- function(x, params) pt((x - params$mean) / (params$sd * sd_fac), df = params$df)

  # Act - compare 1000 KS statistics from the code in this package to the value produced by ADGofTest::ad.test
  actual_expected <- replicate(1000, {
    x <- rnorm(100)
    expected <- unname(ADGofTest::ad.test((x - mean(x)) / (sd(x) * sd_fac), pt, df = 5)$statistic)
    actual <- calc_ad_test_stat(x, list(mean = mean(x), sd = sd(x), df = 5), fn_p)
    return(c(actual, expected))
  })

  # Assert
  expect_equal(actual_expected[1,], actual_expected[2,])
})

test_that("rlogis", {
  # Arrange
  set.seed(1)
  sd_fac <- sqrt(3) / pi
  fn_p <- function(x, params) plogis(x, location = params$location, scale = params$scale)

  # Act - compare 1000 KS statistics from the code in this package to the value produced by ADGofTest::ad.test
  actual_expected <- replicate(1000, {
    x <- rnorm(100)
    expected <- unname(ADGofTest::ad.test(x, plogis, location = mean(x), scale = sd(x) * sd_fac)$statistic)
    actual <- calc_ad_test_stat(x, list(location = mean(x), scale = sd(x) * sd_fac), fn_p)
    return(c(actual, expected))
  })

  # Assert
  expect_equal(actual_expected[1,], actual_expected[2,])
})
