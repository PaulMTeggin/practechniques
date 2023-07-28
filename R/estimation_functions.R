#' Estimate the exponential distribution rate parameter by method of moments
#'
#' This function has the right signature to be used in \code{\link{gof_test_sim}},
#' and estimates the rate parameter as the reciprocal of the sample mean.
#'
#' @param x The data.
#' @param noverlap The extent to which the data is overlapped. This is not used,
#'   but required for compatibility with \code{\link{gof_test_sim}}.
#'
#' @return A list with one item:
#' \describe{
#'   \item{rate}{the rate parameter, calculated as the reciprocal of the mean of `x`}
#' }
#' @family Estimation functions
#' @export
#'
#' @examples
#' estimate_exp_ol(rexp(100))
estimate_exp_ol <- function(x, noverlap = 1) {
  list(rate = 1 / mean(x))
}

#' Estimate mean and standard deviation by method of moments allowing for potentially overlapping data
#'
#' This function serves two purposes: (1) collate the mean and standard deviation into a single list object,
#' and (2) allow for bias in the standard deviation arising from the overlapping data
#' by using \code{\link{calc_sd_ol_bias_fac}}.
#'
#' @param x The data.
#' @param noverlap The extent to which the data is overlapped. \code{noverlap == 1} (the default) means no overlap,
#'   and in this case \code{\link{calc_sd_ol_bias_fac}} returns 1, so no adjustment is applied.
#'
#' @return A list with two items:
#' \describe{
#'   \item{mean}{the mean of `x`}
#'   \item{sd}{the standard deviation of \code{x} adjusted for overlapping data bias, if any}
#' }
#' @family Estimation functions
#' @export
#'
#' @examples
#' estimate_mean_sd_ol(rnorm(100))
estimate_mean_sd_ol <- function(x, noverlap = 1) {
  list(mean = mean(x), sd = stats::sd(x) * calc_sd_ol_bias_fac(length(x), noverlap))
}

#' Estimate location and scale parameters for the logistic distribution by method of moments allowing for potentially overlapping data
#'
#' This function serves two purposes: (1) collate the location and scale parameters into a single list object,
#' and (2) allow for bias in the standard deviation arising from the overlapping data
#' by using \code{\link{calc_sd_ol_bias_fac}}.
#'
#' The scale parameter is derived from the sample standard deviation by multiplying by `sqrt(3) / pi`.
#' See [the formula for logistic variance on Wikipedia](https://en.wikipedia.org/wiki/Logistic_distribution).
#' Note that the adjustment for overlapping data bias is multiplicative to the sample standard deviation
#' and so applies to the logistic scale parameter in the same way that it does to the Normal sample standard deviation.
#'
#' @param x The data whose mean and standard deviation are required.
#' @param noverlap The extent to which the data is overlapped. \code{noverlap == 1} (the default) means no overlap,
#'   and in this case \code{\link{calc_sd_ol_bias_fac}} returns 1, so no adjustment is applied.
#'
#' @return
#' \describe{
#'   \item{location}{the mean of `x`}
#'   \item{scale}{the scale parameter derived from standard deviation of \code{x} adjusted for overlapping data bias, if any}
#' }

#' @family Estimation functions
#' @export
#'
#' @examples
#' estimate_logis_ol(rlogis(100))
estimate_logis_ol <- function(x, noverlap = 1) {
  mean_sd_ol <- estimate_mean_sd_ol(x, noverlap)
  sd_fac_logis <- sqrt(3) / pi
  list(location = mean_sd_ol$mean, scale = mean_sd_ol$sd * sd_fac_logis)
}
