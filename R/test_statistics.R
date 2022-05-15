#' Calculate the one-sample KS test statistic of data against a fitted distribution
#'
#' The purpose of this function is to allow the KS test statistic to be calculated
#' for a wide range of distributions, i.e. to abstract the calculation of the test
#' statistic from the CDF (p) function.
#' This requires the estimated parameters to be specified as a single object.
#'
#' @param x The data.
#' @param params The parameters of the distribution, generally estimated from the data,
#'   as a single object.
#' @param fn_p The cumulative distribution function of the distribution,
#'   in the form of a function that takes the data and estimated parameters as a single object,
#'   and returns cumulative probability values.
#'
#' @return The KS test statistic.
#' @family Test statistic functions
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' params <- list(mean = mean(x), sd = sd(x))
#' calc_ks_test_stat(x, params, function(x, params) pnorm(x, params$mean, params$sd))
#' stats::ks.test(x, "pnorm", mean(x), sd(x))$statistic
calc_ks_test_stat <- function(x, params, fn_p)
{
  N <- length(x)
  p <- fn_p(sort(x), params)

  diffs <- p - (0:(N - 1)) / N
  D <- max(c(diffs, 1 / N - diffs))
  return(D)
}

#' Calculate the AD test statistic of data against a fitted distribution
#'
#' The purpose of this function is to allow the AD test statistic to be calculated
#' for a wide range of distributions, i.e. to abstract the calculation of the test
#' statistic from the CDF (p) function.
#' This requires the estimated parameters to be specified as a single object.
#'
#' @param x The data.
#' @param params The parameters of the distribution, generally estimated from the data,
#'   as a single object.
#' @param fn_p The cumulative distribution function of the distribution,
#'   in the form of a function that takes the data and estimated parameters as a single object,
#'   and returns cumulative probability values.
#'
#' @return The AD test statistic.
#' @family Test statistic functions
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' params <- list(mean = mean(x), sd = sd(x))
#' calc_ad_test_stat(x, params, function(x, params) pnorm(x, params$mean, params$sd))
#' if (require(ADGofTest)) ADGofTest::ad.test(x, pnorm, mean(x), sd(x))$statistic
calc_ad_test_stat <- function(x, params, fn_p) {
  N <- length(x)
  p <- fn_p(sort(x), params)

  A <- p * (1 - rev(p))
  A <- (2 * seq(p) - 1) * log(A)
  return(-mean(A) - N)
}
