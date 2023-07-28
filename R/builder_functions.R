#' Build a uniparameter-friendly list of arguments
#'
#' In general params will be a list returned by an estimation function, so just tack the x argument
#' on the front for do.call purposes later. However sometimes the params may be a bespoke object
#' of some kind, e.g. the result of calling ghyp() from the ghyp package, and needs to be wrapped in a list.
#'
#' @param x The non-parameter element of the arguments (e.g. `p`, for `q`-functions).
#' @param params The parameters.
#'
#' @return A list consisting of `x` and `params` in a format suitable for passing to \code{\link{do.call}}.
#' @export
#'
#' @examples
#' build_uniparam_args(1:10, list(mean = 0, sd = 1))
build_uniparam_args <- function(x, params) {
  if (is.list(params)) {
    arg_list <- c(list(x), params)
  } else {
    arg_list <- c(list(x), list(params))
  }

  return(arg_list)
}

#' Get a prefixed distribution function
#'
#' Helper function to get, say, the function \code{\link[stats]{pnorm}}
#' by passing `"p"` and `"norm"`.
#'
#' @param prefix The prefix (generally one of `p`, `q`, `r` or `d`).
#' @param dist_name The name of the distribution.
#'
#' @return The function corresponding to the prefix and distribution_name.
#' @export
#'
#' @examples
#' get_prefixed_dist_fn("p", "norm")
get_prefixed_dist_fn <- function(prefix, dist_name) {
  get(paste0(prefix, dist_name), mode = "function")
}

#' Build a uniparameter-friendly function
#'
#' If `fn` is `NULL`, [get_prefixed_dist_fn] will be called to get an function.
#' This is then converted to a uniparameter-friendly version that calls [build_uniparam_args]
#' internally so it can then call [do.call].
#'
#' @inheritParams get_prefixed_dist_fn
#' @param fn If not `NULL`, will override the result of calling [get_prefixed_dist_fn]
#'   with `prefix` and `dist`
#'
#' @return A uniparameter version of the input function.
#' @export
#'
#' @examples
#' uniparam_pnorm <- build_uniparam_fn("p", "norm")
#' uniparam_pnorm(0, list(mean = 0, sd = 1))
build_uniparam_fn <- function(prefix, dist_name, fn = NULL) {
  if (is.null(fn)) fn <- get_prefixed_dist_fn(prefix, dist_name)

  fn_uniparam <- function(x, params) {
    arg_list <- build_uniparam_args(x, params)

    return(do.call(fn, arg_list))
  }

  return(fn_uniparam)
}

#' Build an estimated function that captures the estimated parameters
#'
#' @inheritParams get_prefixed_dist_fn
#' @param fn_estimate_params A function that takes the data and the extent of the overlap in the data,
#'   and returns a single object holding estimated parameters of the distribution being fitted.
#'   The method of estimation should be unbiased. Note that for many distributions, MLE
#'   only gives *asymptotically* unbiased parameters. Users should validate that their estimation
#'   functions are unbiased and if necessary adjust the threshold p-value to compensate for this.
#' @param x_data The data from which the parameters will be estimated.
#' @param noverlap The extent of any overlap in the data. `1` means no overlap.
#'
#' @return A function with captured estimated parameters that is called with a single argument
#'   (e.g. `q` for a `p` function).
#' @export
#'
#' @examples
#' x_data <- rnorm(100)
#' pnorm_est <- build_estimated_fn("p", "norm", estimate_mean_sd_ol, x_data)
#' pnorm_est(0)
#' pnorm(0, mean(x_data), sd(x_data))
build_estimated_fn <- function(prefix, dist_name, fn_estimate_params, x_data, noverlap = 1) {
  fn_uniparam <- build_uniparam_fn(prefix, dist_name, fn = NULL)
  estimated_params <- fn_estimate_params(x_data, noverlap)

  fn_uniparam_est <- function(x) {
    fn_uniparam(x, estimated_params)
  }

  return(fn_uniparam_est)
}

# Internal function that builds a function to calculate test statistics from data and estimated parameters
# by capturing the CDF (p) function. This gives a unified signature for use in gof_test_sim_uniparam.
.build_calc_test_stat <- function(fn_calc, fn_p) {
  captured_fn_p <- fn_p
  function(x, params) fn_calc(x, params, captured_fn_p)
}
