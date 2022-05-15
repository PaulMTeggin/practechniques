#' Calculate theoretical serial correlation at a given lag from overlapped data
#'
#' Serial correlation arises from overlapped data because some of the information is 'shared'
#' between successive overlapping observations. The amount of serial correlation
#' is a function of the lag between observations and the extent of the overlap.
#' For example the serial correlation between Jan-Dec and Feb-Jan observations is higher
#' than between Jan-Dec and Sep-Aug-next-year observations.
#'
#' The formula follows from writing each overlapping observation as the sum of
#' `noverlap` independent observations and counting the extent to which these two sets overlap.
#' Where the lag is higher than the overlap no information is shared and the theoretical
#' serial correlation becomes zero.
#'
#' @param noverlap The extent of the overlap, e.g. 12 means annual overlaps from monthly data.
#'   1 means no overlap.
#' @param lag The lag in time between two overlapped observations.
#'   The serial correlation at lags equal to or higher than `noverlap` is zero.
#'
#' @return The theoretical serial correlation.
#' @family Overlapping data functions
#' @export
#'
#' @examples
#' calc_theo_sercor(12, 0:12)
calc_theo_sercor <- function(noverlap, lag) {
  return(pmax(noverlap - abs(lag), 0) / noverlap)
}

#' Build a matrix of theoretical serial correlations
#'
#' The rows and columns of the matrix represent indexes to the overlapped data.
#' All values on the leading diagonal will be 1, since a variable is always perfectly correlated
#' with itself. Values off the diagonal will depend on the size of the overlap (`noverlap`).
#'
#' @param N The size of the matrix (number of rows and columns)
#' @inheritParams calc_theo_sercor
#'
#' @return The matrix of theoretical serial correlations.
#' @family Overlapping data functions
#' @export
#'
#' @examples
#' build_theo_sercor_mtx(24, 12)
build_theo_sercor_mtx <- function(N, noverlap) {
  return(outer(seq_len(N), seq_len(N), function(i, j) {
    calc_theo_sercor(noverlap, abs(i - j))
  }))
}

#' Calculate the overlapped bias factor for sample standard deviation
#'
#' When overlapped data is used, the sample variance and standard deviation are biased.
#' This function calculates a factor that corrects for the bias in the sample standard deviation.
#' It assumes [Bessel's correction](https://en.wikipedia.org/wiki/Bessel%27s_correction)
#' has already been applied in the calculation of the sample
#' standard deviation and removes this to avoid double-counting.
#' In practice this means this factor can be applied to the results of \code{\link[stats]{sd}}.
#'
#' The factor is derived in a
#' [2009 Risk.net article](https://www.risk.net/risk-management/1509219/error-var-overlapping-intervals)
#' by Sun, Nelken et al. Where `noverlap == 1` the factor is 1 and has no numerical effect.
#'
#' @param N The number of overlapped data points.
#'   There will be `N + noverlap - 1` in the data prior to taking overlapped samples,
#'   which reduces to `N` after overlapping by `noverlap`.
#'   E.g. we need 11 more observations than `N` in order to have `N` annual monthly-overlapped observations.
#' @inheritParams calc_theo_sercor
#'
#' @return A factor to multiply the sample standard deviation by to give an unbiased estimate.
#' @references \url{https://www.risk.net/risk-management/1509219/error-var-overlapping-intervals}

#' @export
#'
#' @examples
#' calc_sd_ol_bias_fac(c(24, 60, 120, 180, 240), 12)
calc_sd_ol_bias_fac <- function(N, noverlap) {
  var_bias <- (N - (noverlap - (noverlap * noverlap - 1) / (3 * N))) / (N - 1)
  return(1/sqrt(var_bias))
}

#' Generate a matrix of overlapped Normal variates
#'
#' The function operates by calculating a matrix of 'pthly' variates
#' (where `p` is the extent of the overlap, i.e. the parameter `noverlap`)
#' having 1/pth the mean and variance, and then iterating over rows and
#' summing the row values in overlapping sets of `noverlap`.
#'
#' @param nsims The number of simulations (rows in the matrix).
#' @param nsteps The number of steps (columns in the matrix).
#' @inheritParams calc_theo_sercor
#' @param mu The mean of the Normal distribution.
#' @param sigma The standard deviation of the Normal distribution.
#'
#' @return A \code{nsims} by \code{steps} matrix with overlapped values.
#' @export
#'
#' @examples
#' generate_ol_rnorm_mtx(100, 10, noverlap = 12)
generate_ol_rnorm_mtx <- function(nsims, nsteps, noverlap, mu = 0, sigma = 1) {
  steps_noverlap <- nsteps + noverlap - 1
  pthly_mtx <- matrix(stats::rnorm(nsims * steps_noverlap, mean = mu / noverlap, sd = sigma / sqrt(noverlap)),
                      nrow = nsims, ncol = steps_noverlap)

  # Sum columns in sets of noverlap to create the overlapping matrix
  mtx <- vapply(seq_len(nsteps), function(i) {
    rowSums(pthly_mtx[,i:(i + noverlap - 1), drop = FALSE])
  }, numeric(nsims))

  return(mtx)
}
