#' Perform a goodness of fit test using simulation with uniparameter plug-in functions
#'
#' Many statistical tests have null hypotheses that assume a distribution is fully specified
#' (with its parameters known in advance). It is common to estimate parameters from data,
#' and in this case a general method for adapting the statistical test is to use
#' Monte Carlo to produce a simulated distribution of the test statistic, and derive the p-value from this distribution.
#' This approach is used in the \code{\link[KScorrect]{LcKS}} function of the \code{KScorrect} package.
#' However, the implementation in \code{LcKS} only supports the KS test and a closed list of distributions,
#' because it has bespoke code for each supported distribution
#' for estimating parameters and simulating values using the estimated parameters.
#' This function generalises the approach in \code{\link[KScorrect]{LcKS}} by
#' adopting the underlying `LcKS` algorithm and allowing general estimation, test statistic and simulation
#' functions to be plugged into that algorithm.
#'
#' This function uses the same general approach as \code{LcKS}, which is to:
#'
#' * Estimate parameters from the input data `x`
#' * Calculate a test statistic for `x` against the specified distribution function with these parameters
#' * Use Monte Carlo simulation to produce a simulated distribution of potential alternative values
#' for the test statistic.
#' * Derive a p-value by comparing the test statistic of `x` against the simulated distribution.
#' The p-value is calculated as the proportion of Monte Carlo samples with test statistics at least as extreme
#' as the test statistic of `x`. A value of 1 is added to both the numerator and denominator for the same reasons as
#' `KScorrect`, which among other reasons has the benefit of avoiding estimated p-values that are precisely zero.
#'
#' However this function is more generic:
#' * General distributions are supported, rather than the closed list used by `KScorrect`.
#' * Multiple statistical tests are supported, not just [KS](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test).
#' * Testing can be performed against distributions fitted to overlapping data, not just IID data,
#' using the idea of a Gaussian copula to induce autocorrelation consistent with overlapping data suggested in section 4.2 of the
#' [2019 paper](https://www.cambridge.org/core/journals/british-actuarial-journal/article/calibration-of-var-models-with-overlapping-data/B20D66D81DB918AFD3BBDF9EDAC20863)
#' by the [Extreme Events Working Party](https://www.actuaries.org.uk/practice-areas/life/research-working-parties/extreme-events)
#' of the UK [Institute and Faculty of Actuaries](https://www.actuaries.org.uk/).
#'
#' The genericity is achieved by requiring all statistical functions involved
#' to be 'uniparameter', i.e. to have all their parameters put into a single object.
#' This entails wrapping (say) \code{\link[stats]{pnorm}} so the wrapper function takes a list containing the
#' `mean` and `sd` parameters, and passes them on.
#'
#' By making all functions take their parameters as single objects, the algorithm
#' used in the `KScorrect` package can be abstracted from the functions for estimating parameters
#' (`fn_estimate_params`), calculating test statistics (`fn_calc_test_stat`),
#' and simulating values using those estimated parameters (`fn_simulate`),
#' These functions are 'plugged in' to the algorithm and called at the appropriate points.
#' They must be mutually compatible with each other.
#'
#' For simplicity and to ensure compatibility,
#' the function \code{\link{gof_test_sim}} sets up the plug-in functions automatically,
#' based on the un-prefixed name of the distribution (e.g. `"norm"`).
#' This has a slight performance hit as it uses \code{\link[base]{do.call}},
#' but this can be avoided if performance is key, by hand-writing the wrapper function.
#'
#' Similarly, adapting to overlapping data requires the simulation to be done in a way that induces the autocorrelation
#' consistent with overlapping data. This function can perform testing on overlapping data
#' by suitable choice of the plug-in function \code{fn_simulate}.
#' In this case the estimation function \code{fn_estimate_params} should also allow for bias
#' in parameter estimation induced by the overlap. There is no need to adapt the test statistic function
#' `fn_calc_test_stat` to overlapping data.
#'
#' For some distributions the estimation of parameters may occasionally fail within the simulation.
#' In this case the test statistic is set to `NA` and disregarded when calculating p-values.
#' Warnings produced in parameter estimation are suppressed as (e.g. when using `MASS::fitdistr`)
#' these often arise from estimating the uncertainty around the estimated parameters, which is not used here.
#'
#' The framework here can in principle also be used where parameters are known in advance
#' rather than estimated from the data (by making the estimation function return the pre-specified parameters),
#' but there is limited value to this use case, as Monte Carlo is rarely necessary
#' when the parameters are known (and is certainly not necessary for the KS and AD tests).
#' It can be an useful approach for hybrid cases such as the 3-parameter Student's t distribution
#' where the number of degrees of freedom is pre-specified but the location and scale parameters are not.
#'
#' Optionally, the calculations can be parallelised over multiple cores using the `doParallel` package.
#' This is useful when the number of simulations is large and estimation of parameters is slow,
#' for example using MLE to estimate parameters from a generalised hyperbolic distribution.
#'
#' Since Monte Carlo simulation is used, the function can optionally estimate the simulation uncertainty arising from
#' a finite number of simulations, using a non-parameteric (resampling with replacement) approach from
#' the distribution of simulated test statistics produced.
#'
#' @param x The data being tested.
#' @param fn_estimate_params A function that takes the data and the extent of the overlap in the data,
#'   and returns a single object holding estimated parameters of the distribution being fitted.
#'   The method of estimation should be unbiased. Note that for many distributions, MLE
#'   only gives *asymptotically* unbiased parameters. Users should validate that their estimation
#'   functions are unbiased and if necessary adjust the threshold p-value to compensate for this.
#' @param fn_calc_test_stat A function that takes the data and the estimated parameters object,
#'   and calculates the test statistic for the distribution being tested.
#' @param fn_simulate A function takes the number of values to simulate, the estimated parameters object,
#'   and the extent of any overlap in the data,
#'   and returns that number of simulated values from the distribution being tested against.
#' @param noverlap The extent of any overlap in the data. `1` means no overlap,
#'   and \code{fn_simulate} should operate by ordinary simulation. If `noverlap > 1` then
#'   autocorrelation must be induced in the simulations that is consistent with the degree of overlap,
#'   to give unbiased test results.
#'   This is done automatically when this function is called via \code{\link{gof_test_sim}}.
#'   `fn_estimate_params` must also allow for the degree of overlap.
#' @param nreps The number of repetitions of the simulation to use.
#' @param parallelise Flag indicating whether or not to parallelise the calculations.
#' @param ncores The number of cores to use when parallelising.
#'   \code{NULL} means one fewer than the number of cores on the machine.
#' @param bs_ci The width of a confidence interval around the p-value,
#'   which will be calculated using a non-parametric bootstrap.
#'   \code{NULL} means no confidence interval will be produced.
#' @param nreps_bs_ci The number of iterations used in the bootstrapped confidence interval.
#'
#' @return A list with five components:
#' \describe{
#'   \item{ts}{The test statistic.}
#'   \item{p_value}{The p-value for the test statistic, derived by simulation.}
#'   \item{count_NA}{The number of \code{NA} values produced in the simulation of the test statistic.
#'   These generally indicate that the parameter estimation failed. These values are disregarded
#'   when calculating p-values.}
#'   \item{p_value_lower}{If \code{bs_ci} is not \code{NULL}, the lower end
#'   of the confidence interval around the p-value, calculated using
#'   a non-parametric bootstrap with \code{nreps_bs_ci} repetitions.
#'   Otherwise \code{NA}.}
#'   \item{p_value_upper}{If \code{bs_ci} is not \code{NULL}, the upper end
#'   of the confidence interval around the p-value, calculated using
#'   a non-parametric bootstrap with \code{nreps_bs_ci} repetitions.
#'   Otherwise \code{NA}.}
#' }
#' @export
#' @importFrom foreach %dopar%
#'
#' @examples
#' fn_estimate_params <- function(x, noverlap = 1) list(mean = mean(x), sd = sd(x))
#' fn_p <- function(x, params) pnorm(x, params$mean, params$sd)
#' fn_test_statistic <- function(x, est_params) calc_ks_test_stat(x, est_params, fn_p)
#' fn_simulate <- function(N, est_params) rnorm(N, est_params$mean, est_params$sd)
#' gof_test_sim_uniparam(rnorm(100), fn_estimate_params, fn_test_statistic, fn_simulate)
gof_test_sim_uniparam <- function(x, fn_estimate_params, fn_calc_test_stat, fn_simulate, noverlap = 1,
                                  nreps=999, parallelise = FALSE, ncores = NULL,
                                  bs_ci = NULL, nreps_bs_ci = 10000) {

  .validate_parallelise_ncores(parallelise, ncores)

  # Estimate the parameters from the data
  estimated_params <- tryCatch(fn_estimate_params(x, noverlap), error = function(e) NULL)

  # If it's not possible to estimate the parameters (signalled by replacing an error with an NULL value)
  # then use NA for all 3 return values.
  if (is.null(estimated_params)) {
    return(list(ts = NA, p_value = NA,
                p_value_lower = NA, p_value_upper = NA,
                count_NA = NA))
  }

  N <- length(x)

  # Calculate the test statistic on the data
  teststat_x <- fn_calc_test_stat(x, estimated_params)

  # Repeatedly simulate from the distribution we are testing against
  # In each repetition (controlled by replicate, with nsims repetitions) we simulate
  # as many values from the distribution as were in the input dataset, and then calculate the test statistic.
  # For some distributions the parameter estimation may occasionally fail - in this case
  # we catch the errors and output an NA value for the test statistic.
  if (parallelise) {
    if (is.null(ncores)) ncores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG() # To ensure reproducibility. Can't use %doRNG% directly.
    if (foreach::getDoParRegistered()) {
      message(sprintf("gof_test_epd is running in parallel with %d worker(s) using %s [%s]\n",
                      foreach::getDoParWorkers(), foreach::getDoParName(),
                      foreach::getDoParVersion()))
    }

    on.exit(parallel::stopCluster(cl))

    teststat_sims <- suppressMessages(suppressWarnings( # Suppression because estimation can produce lots of messages
      foreach::foreach(seq_len(nreps), .combine = c, .inorder = FALSE) %dopar% {
        tryCatch({
          sim_data <- fn_simulate(N, estimated_params)
          sim_params <- fn_estimate_params(sim_data, noverlap)
          fn_calc_test_stat(sim_data, sim_params)
        }, error = function(e) NA)
      }
    ))

  } else {

    teststat_sims <- suppressMessages(suppressWarnings( # Suppression because estimation can produce lots of messages
      replicate(nreps, {
        tryCatch({
          sim_data <- fn_simulate(N, estimated_params)
          sim_params <- fn_estimate_params(sim_data, noverlap)
          fn_calc_test_stat(sim_data, sim_params)
        }, error = function(e) NA)
      })
    ))
  }

  # Count the NA values and then remove them for the calculation of the p-value
  count_NA <- sum(is.na(teststat_sims))
  teststat_sims <- stats::na.omit(teststat_sims)

  # Determine the p-value from the simulated results
  # p-value is the proportion of sims which gave a higher test statistic
  # (If very few sims gave a higher test statistic, then the probability of the data being consistent with the
  # assumed distribution is very low.)
  # The +1 in numerator and denominator is the same approach used by KScorrect::LcKS.
  p_value <- .calc_prop_higher(teststat_sims, teststat_x)

  bs_p_value_ci <- .bootstrap_p_value_ci(bs_ci, nreps_bs_ci, teststat_sims, teststat_x)

  # Output the test statistic, p-value (and the lower and upper ends of the
  # confidence interval around it, if available)
  # and the number of NA values produced in the simulation of test statistics
  return(list(ts = teststat_x, p_value = p_value,
              p_value_lower = bs_p_value_ci$p_value_lower, p_value_upper = bs_p_value_ci$p_value_upper,
              count_NA = count_NA))
}

# Internal function to validate the values of the parameters parallelise and ncores
.validate_parallelise_ncores <- function(parallelise, ncores) {
  if (!parallelise & !is.null(ncores))
    stop("You specified to use multiple cores but not to run in parallel. Set parallelise = TRUE to run in parallel.")

  if (parallelise & !is.null(ncores)) {
    if (is.numeric(ncores)) {
      ncores <- as.integer(ncores)

      if (ncores < 2)
        stop("Set ncores to an integer greater than 1 to run in parallel.")
      if (ncores > parallel::detectCores())
        warning("You are attempting to run this function on more cores than your computer contains. Consider reducing ncores to improve efficiency.")
    }
    else {
      stop("Set ncores to an integer greater than 1 to run in parallel.")
    }
  }
}

# Internal function to calculate the proportion of test statistics exceeding the test statistic obtained from the data
# This uses the same approach as KScorrect in adding one to the numerator and denominator
.calc_prop_higher <- function(teststat_sims, teststat_x) {
  (sum(teststat_sims > teststat_x) + 1) / (length(teststat_sims) + 1)
}

# Internal function to give a confidence interval around the p-value using a non-parametric bootstrap
.bootstrap_p_value_ci <- function(bs_ci, nreps_bs_ci, teststat_sims, teststat_x) {
  if (is.null(bs_ci)) {
    p_value_lower = NA
    p_value_upper = NA
  } else {
    # Not worth parallelising as each iteration is so fast
    bs_p_values <- replicate(nreps_bs_ci, {
      bs_teststat_sims <- sample(teststat_sims, length(teststat_sims), replace = TRUE)
      .calc_prop_higher(bs_teststat_sims, teststat_x)
    })

    ci_lower <- (1 - bs_ci) / 2
    p_value_lower <- stats::quantile(bs_p_values, ci_lower)
    p_value_upper <- stats::quantile(bs_p_values, 1 - ci_lower)
  }

  return(list(p_value_lower = p_value_lower, p_value_upper = p_value_upper))
}

#' Perform a goodness of fit test using simulation
#'
#' Many statistical tests have null hypotheses that assume a distribution is fully specified
#' (with its parameters known in advance). It is common to estimate parameters from data,
#' and in this case a general method for adapting the statistical test is to use
#' Monte Carlo to produce a simulated distribution of the test statistic, and derive the p-value from this distribution.
#' This approach is used in the \code{\link[KScorrect]{LcKS}} function of the \code{KScorrect} package.
#' However, the implementation in \code{LcKS} only supports the KS test and a closed list of distributions,
#' because it has bespoke code for each supported distribution
#' for estimating parameters and simulating values using the estimated parameters.
#' This function generalises the approach in \code{\link[KScorrect]{LcKS}} as explained in the documentation
#' of \code{\link{gof_test_sim_uniparam}}. It is a higher-level function
#' that configures \code{\link{gof_test_sim_uniparam}} appropriately given the name of the distribution being
#' tested against. It is still necessary to provide an appropriate estimation function for \code{fn_estimate_params}
#' that is adapted to the distribution being tested including whether or not data is overlapping.
#'
#' As far as possible this function abstracts from
#' the lower-level implementation details of \code{\link{gof_test_sim_uniparam}}, by taking the name of
#' a distribution function (e.g. `"norm"`) and creating uniparameter versions of the CDF
#' ("p" function, so `pnorm` in this example) and simulation ("r" functions, here `rnorm`).
#' Where overlapping data is used (`noverlap > 1`) it will simulate from a Gaussian copula
#' to induce the same autocorrelation structure as overlapping data
#' (see \code{\link{calc_theo_sercor}}), convert this to autocorrelated random uniform values,
#' and apply a uniparameter version of the inverse CDF ("q" function) for `"dist"`.
#'
#' Notwithstanding this abstraction, it is necessary to supply an appropriate parameter estimation function
#' for `fn_estimate_params`. \code{\link{estimate_mean_sd_ol}} is used by default, to match with the
#' default value of `dist = "norm"`.
#'
#' The framework here can in principle also be used where parameters are known in advance
#' rather than estimated from the data (by making the estimation function return the pre-specified parameters),
#' but there is very little value to this use case, as Monte Carlo is rarely necessary
#' when the parameters are known (and is certainly not necessary for the KS and AD tests).
#'
#' Optionally, the calculations can be parallelised over multiple cores.
#' This is useful when the number of simulations is large
#' and estimation of parameters is slow,
#' for example using MLE to estimate parameters from a generalised hyperbolic distribution.
#'
#' Since Monte Carlo simulation is used, the function can optionally estimate the simulation uncertainty arising from
#' a finite number of simulations, using a non-parameteric (resampling with replacement) approach from
#' the distribution of simulated test statistics produced.
#'
#' @inheritParams gof_test_sim_uniparam
#' @param test_type The type of the test. Either a character string (KS and AD are supported)
#'   or a function that implements a different test statistic with the same signature as
#'   \code{\link{calc_ks_test_stat}} or \code{\link{calc_ad_test_stat}}.
#' @param dist The name of a distribution, such that it can be prepended by \code{"p"} to get
#'   a probability function, and by \code{"r"} to get a random simulation function.
#'   For example \code{"norm"} or \code{"unif"}.
#'
#' @return A list with five components:
#' \describe{
#'   \item{ts}{The test statistic.}
#'   \item{p_value}{The p-value for the test statistic, derived by simulation.}
#'   \item{count_NA}{The number of \code{NA} values produced in the simulation of the test statistic.
#'   These generally indicate that the parameter estimation failed.}
#'   \item{p_value_lower}{If \code{bs_ci} is not \code{NULL}, the lower end
#'   of the confidence interval around the p-value, calculated using
#'   a non-parametric bootstrap with \code{nreps_bs_ci} repetitions.
#'   Otherwise \code{NA}.}
#'   \item{p_value_upper}{If \code{bs_ci} is not \code{NULL}, the upper end
#'   of the confidence interval around the p-value, calculated using
#'   a non-parametric bootstrap with \code{nreps_bs_ci} repetitions.
#'   Otherwise \code{NA}.}
#' }
#' @export
#'
#' @examples
#' gof_test_sim(rnorm(100))
#' estimate_unif <- function(x, noverlap = 1) list(min = min(x), max = max(x))
#' gof_test_sim(runif(100), dist = "unif", fn_estimate_params = estimate_unif)
gof_test_sim <- function(x, test_type = c("KS", "AD"), dist = "norm", noverlap = 1, fn_estimate_params = estimate_mean_sd_ol,
                         nreps=999, parallelise = FALSE, ncores = NULL,
                         bs_ci = NULL, nreps_bs_ci = 10000) {
  if (is.character(test_type)) {
    test_type <- match.arg(test_type)
    fn_calc_test_stat = switch(test_type,
           KS = calc_ks_test_stat,
           AD = calc_ad_test_stat)

  } else if (is.function(test_type)) {
    fn_calc_test_stat <- test_type

  } else {
    stop("test_type should either be a character variable naming the kind of test or a function calculating the test statistic")
  }

  # In general params will be a list return by an estimation function, so just tack the x argument
  # on the front for do.call purposes later. However sometimes the params may be a bespoke object
  # of some kind, e.g. the result of calling ghyp() from the ghyp package, and needs to be wrapped in a list.
  # Define this function here so it exists for parallelised execution and does not need to be explicitly exported
  # as part of parallelisation.
  .build_uniparam_args <- function(x, params) {
    if (is.list(params)) {
      arg_list <- c(list(x), params)
    } else {
      arg_list <- c(list(x), list(params))
    }

    return(arg_list)
  }

  # Find the "p function" (e.g. pnorm) corresponding to dist,
  # and construct a uniparameter version of it using do.call
  fn_p <- get(paste0("p", dist))
  fn_p_uniparam <- function(q, params) {
    arg_list <- .build_uniparam_args(q, params)

    return(do.call(fn_p, arg_list))
  }

  # The "p" function is used as part of calculating the test statistic
  # - construct a uniparameter version of that function that captures the p function.
  fn_calc_test_stat_built <- .build_calc_test_stat(fn_calc_test_stat, fn_p = fn_p_uniparam)

  # For the simulation function, there are two cases to consider: no overlap (noverlap == 1),
  # when we can simulate by using the "r" function corresponding to dist and creating a uniparameter wrapper,
  # or with overlap, when we have to use a Gaussian copula approach per the Extreme Events Working Party paper.
  if (noverlap == 1) {
    fn_r <- get(paste0("r", dist))
    fn_simulate_uniparam <- function(n, params) {
      arg_list <- .build_uniparam_args(n, params)

      return(do.call(fn_r, arg_list))
    }
  } else {
    noverlap = as.integer(noverlap)
    if (noverlap > 1) {
      # Get the Cholesky decomposition of the theoretical autocorrelation matrix
      # so we can induce serial correlation in the simulated data
      # For speed, this is done outside the simulation function and captured.
      # This assumes that we will always be passed simulated data of the same length as the
      # data being tested, which is the case in practice.
      ol_chol <- chol(build_theo_sercor_mtx(length(x), noverlap))

      # Get the inverse CDF ("q" function) so we can simulate with autocorrelation
      # by applying a Gaussian copula to uniforms and inverting
      fn_q <- get(paste0("q", dist))

      fn_simulate_uniparam <- function(n, params) {
        # Get autocorrelated uniform values using a Gaussian copula approach
        p <- stats::pnorm((stats::rnorm(n) %*% ol_chol)[1,]) # The [1,] just converts from an n x 1 matrix to a vector

        arg_list <- .build_uniparam_args(p, params)

        return(do.call(fn_q, arg_list))
      }
    } else {
      stop("noverlap must be at least 1")
    }
  }

  # Now we have all our plugin functions and just use them with gof_test_sim_uniparam
  return(gof_test_sim_uniparam(x = x,
                               fn_estimate_params = fn_estimate_params,
                               fn_calc_test_stat = fn_calc_test_stat_built,
                               fn_simulate = fn_simulate_uniparam,
                               noverlap = noverlap,
                               nreps = nreps,
                               parallelise = parallelise, ncores = ncores,
                               bs_ci = bs_ci, nreps_bs_ci = nreps_bs_ci))


}

# Internal function that builds a function to calculate test statistics from data and estimated parameters
# by capturing the CDF (p) function. This gives a unified signature for use in gof_test_sim_uniparam.
.build_calc_test_stat <- function(fn_calc, fn_p) {
  captured_fn_p <- fn_p
  function(x, params) fn_calc(x, params, captured_fn_p)
}
