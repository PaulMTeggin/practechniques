#' Perform a goodness of fit test using simulation and uniparameter statistical functions
#'
#' Many statistical tests have null hypotheses that assume a distribution is fully specified
#' (with its parameters known in advance). It is common to estimate parameters from data,
#' and in this case a general method for adapting the statistical test is to use
#' simulation to derive the distribution of the test statistic, and derive the p-value from this distribution.
#'
#' Similarly in order to adapt to overlapping data when performing statistical tests,
#' it is necessary to use simulation, and to simulate in a way that induces the autocorrelation
#' consistent with overlapping data. This function can perform testing on overlapping data
#' by suitable choice of the plug-in function \code{fn_simulate}.
#' In this case the estimation function \code{fn_estimate_params} should also allow for bias
#' in parameter estimation induced by the overlap.
#'
#' In order for this function to be generic, it assumes all statistical functions involved
#' are uniparameter, i.e. all their parameters are put into a single object.
#' This entails wrapping (say) \code{\link[stats]{pnorm}} so the wrapper takes a list containing the
#' \code{mean} and \code{sd} parameters, and passes them on.
#' The function \code{\link{gof_test_sim}} sets up the wrapper functions automatically,
#' based on the un-prefixed name of the distribution (e.g. \code{norm}).
#'
#' TODO explain design further.
#'
#' @param x The data being tested.
#' @param fn_estimate_params A function that takes the data and returns an object representing
#'   the parameters of the distribution being fitted.
#' @param fn_calc_test_stat A function that takes the data and the estimated parameters,
#'   and calculates the test statistic for the distribution being tested.
#' @param fn_simulate A function takes the number of values to simulate, a set of estimated parameters,
#'   and the extent of any overlap in the data,
#'   and returns that number of simulated values from the distribution being tested.
#' @param noverlap The extent of any overlap in the data. \code{1} means no overlap,
#'   and \code{fn_simulate} should operate by ordinary simulation. If \code{noverlap > 1} then
#'   autocorrelation must be induced in the simulations that is consistent with the degree of overlap,
#'   to give unbiased test results. \code{fn_estimate_params} must also allow for the degree of overlap.
#' @param nreps The number of repetitions of the simulation to use.
#' @param parallelise Flag indicating whether or not to parallelise the calculations.
#' @param ncores The number of cores to use when parallelising.
#'   \code{NULL} means one fewer than the number of cores on the machine.
#' @param bs_ci The width of a confidence interval around the p value,
#'   which will be calculated using a non-parametric bootstrap.
#'   \code{NULL} means no confidence interval will be produced.
#' @param nreps_bs_ci The number of iterations used in the bootstrapped confidence interval.
#'
#' @return A list with five components:
#' \itemize{
#'   \item{ts}{The test statistic.}
#'   \item{p_value}{The p value for the test statistic, derived by simulation.}
#'   \item{count_NA}{The number of \code{NA} values produced in the simulation of the test statistic.
#'   These generally indicate that the parameter estimation failed.}
#'   \item{p_value_lower}{If \code{bs_ci} is not \code{NULL}, the lower end
#'   of the confidence interval around the p value, calculated using
#'   a non-parametric bootstrap with \code{nreps_bs_ci} repetitions.
#'   Otherwise \code{NA}.}
#'   \item{p_value_upper}{If \code{bs_ci} is not \code{NULL}, the upper end
#'   of the confidence interval around the p value, calculated using
#'   a non-parametric bootstrap with \code{nreps_bs_ci} repetitions.
#'   Otherwise \code{NA}.}
#' }
#' @export
#' @importFrom foreach %dopar%
#'
#' @examples
#' fn_estimate_params <- function(x) list(mean = mean(x), sd = sd(x))
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

  # Count the NA values and then remove them for the calculation of the p value
  count_NA <- sum(is.na(teststat_sims))
  teststat_sims <- stats::na.omit(teststat_sims)

  # Determine the p-value from the simulated results
  # p-value is the proportion of sims which gave a higher test statistic
  # (If very few sims gave a higher test statistic, then the probability of the data being consistent with the
  # assumed distribution is very low.)
  # The +1 in numerator and denominator is the same approach used by KScorrect::LcKS.
  p_value <- .calc_prop_higher(teststat_sims, teststat_x)

  bs_p_value_ci <- .bootstrap_p_value_ci(bs_ci, nreps_bs_ci, teststat_sims, teststat_x)

  # Output the test statistic, p value (and the lower and upper ends of the
  # confidence interval around it, if available)
  # and the number of NA values produced in the simulation of test statistics
  return(list(ts = teststat_x, p_value = p_value,
              p_value_lower = bs_p_value_ci$p_value_lower, p_value_upper = bs_p_value_ci$p_value_upper,
              count_NA = count_NA))
}

.validate_parallelise_ncores <- function(parallelise, ncores) {
  if (!parallelise & !is.null(ncores))
    stop("You specified to use multiple 'cores' but not to run in parallel. Set parallelise = TRUE to run in parallel.")

  if (parallelise & !is.null(ncores)) {
    if (is.numeric(ncores)) {
      ncores <- as.integer(ncores)

      if (ncores < 2)
        stop("Set 'ncores' to integer greater than 1 to run in parallel.")
      if (ncores > parallel::detectCores())
        warning("You are attempting to run this function on more cores than your computer contains. Consider reducing 'cores' to improve efficiency.")
    }
    else {
      stop("Set 'ncores' to integer greater than 1 to run in parallel.")
    }
  }
}

.calc_prop_higher <- function(teststat_sims, teststat_x) {
  (sum(teststat_sims > teststat_x) + 1) / (length(teststat_sims) + 1)
}

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
#' simulation to derive the distribution of the test statistic, and derive the p-value from this distribution.
#'
#' TODO explain relationship to gof_test_sim_uniparam
#'
#' @inheritParams gof_test_sim_uniparam
#' @param test_type The type of the test. Either a character string (KS and AD are supported)
#'   or a function that implements a different test statistic with the same signature as
#'   \code{\link{calc_ks_test_stat}} or \code{\link{calc_ad_test_stat}}.
#' @param dist The name of a distribution, such that it can be prepended by \code{"p"} to get
#'   a probability function, and by \code{"r"} to get a random simulation function.
#'   For example \code{"norm"} or \code{"unif"}.
#' @param noverlap The extent of any overlap in the data. \code{1} means no overlap,
#'   and the test operates by ordinary simulation. If \code{noverlap > 1} then
#'   autocorrelation will be induced in the simulations that is consistent with the degree of overlap,
#'   to give unbiased test results. \code{fn_estimate_params} must also allow for the degree of overlap.
#'
#' @return A list with five components:
#' \itemize{
#'   \item{ts}{The test statistic.}
#'   \item{p_value}{The p value for the test statistic, derived by simulation.}
#'   \item{count_NA}{The number of \code{NA} values produced in the simulation of the test statistic.
#'   These generally indicate that the parameter estimation failed.}
#'   \item{p_value_lower}{If \code{bs_ci} is not \code{NULL}, the lower end
#'   of the confidence interval around the p value, calculated using
#'   a non-parametric bootstrap with \code{nreps_bs_ci} repetitions.
#'   Otherwise \code{NA}.}
#'   \item{p_value_upper}{If \code{bs_ci} is not \code{NULL}, the upper end
#'   of the confidence interval around the p value, calculated using
#'   a non-parametric bootstrap with \code{nreps_bs_ci} repetitions.
#'   Otherwise \code{NA}.}
#' }
#' @export
#'
#' @examples
#' gof_test_sim(rnorm(100))
#' estimate_unif <- function(x) list(min = min(x), max = max(x))
#' gof_test_sim(runif(100), dist = "unif", fn_estimate_params = estimate_unif)
gof_test_sim <- function(x, test_type = c("KS", "AD"), dist = "norm", noverlap = 1, fn_estimate_params = estimate_MoM_2,
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

  # Define this function here so it exists for parallelised execution
  .build_uniparam_args <- function(x, params) {
    if (is.list(params)) {
      arg_list <- c(list(x), params)
    } else {
      arg_list <- c(list(x), list(params))
    }

    return(arg_list)
  }

  fn_p <- get(paste0("p", dist))
  fn_p_uniparam <- function(q, params) {
    arg_list <- .build_uniparam_args(q, params)

    return(do.call(fn_p, arg_list))
  }

  fn_calc_test_stat_built <- .build_calc_test_stat(fn_calc_test_stat, fn_p = fn_p_uniparam)

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

      # Get the inverse CDF so we can simulate with autocorrelation
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


  return(gof_test_sim_uniparam(x = x,
                               fn_estimate_params = fn_estimate_params,
                               fn_calc_test_stat = fn_calc_test_stat_built,
                               fn_simulate = fn_simulate_uniparam,
                               noverlap = noverlap,
                               nreps = nreps,
                               parallelise = parallelise, ncores = ncores,
                               bs_ci = bs_ci, nreps_bs_ci = nreps_bs_ci))


}

#' Calculate the one-sample KS test statistic of data against a fitted distribution
#'
#' The purpose of this function is to allow the KS test statistic to be calculated
#' for a wide range of distributions, i.e. to abstract the calculation of the test
#' statistic from the CDF (p) function.
#' This requires the estimated parameters to be specified as a single object.
#'
#' @param x The data.
#' @param fn_p The cumulative distribution function of the distribution,
#'   in the form of a function that takes the data and estimated parameters as a single object,
#'   and returns cumulative probability values.
#' @param params The parameters of the distribution, generally estimated from the data.
#'
#' @return The KS test statistic.
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
#' @param fn_p The cumulative distribution function of the distribution,
#'   in the form of a function that takes the data and estimated parameters as a single object,
#'   and returns cumulative probability values.
#' @param params The parameters of the distribution, generally estimated from the data.
#'
#' @return The AD test statistic.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' params <- list(mean = mean(x), sd = sd(x))
#' calc_ad_test_stat(x, params, function(x, params) pnorm(x, params$mean, params$sd))
#' # ADGofTest::ad.test(x, "pnorm", mean(x), sd(x))$statistic
calc_ad_test_stat <- function(x, params, fn_p) {
  N <- length(x)
  p <- fn_p(sort(x), params)

  A <- p * (1 - rev(p))
  A <- (2 * seq(p) - 1) * log(A)
  return(-mean(A) - N)
}

.build_calc_test_stat <- function(fn_calc, fn_p) {
  captured_fn_p <- fn_p
  function(x, params) fn_calc(x, params, captured_fn_p)
}

estimate_MoM_2 <- function(x, noverlap = 1) {
  list(mean = mean(x), sd = stats::sd(x) * calc_sd_ol_bias_fac(length(x), noverlap))
}

# ks_test_epd_norm_uniparam <- function(x, nreps = 999,
#                                       parallelise = FALSE,
#                                       ncores = NULL,
#                                       bs_ci = NULL, nreps_bs_ci = 10000) {
#   pnorm_uniparam <- function(x, params) {
#     stats::pnorm(x, params$mean, params$sd)
#   }
#
#   rnorm_uniparam <- function(n, params) {
#     stats::rnorm(n, params$mean, params$sd)
#   }
#
#   gof_test_epd(x,
#                fn_estimate_params = estimate_MoM_2,
#                fn_calc_test_stat = build_calc_test_stat(calc_ks_test_stat, pnorm_uniparam),
#                fn_simulate = rnorm_uniparam,
#                noverlap = 1,
#                nreps = nreps,
#                parallelise = parallelise,
#                ncores = ncores,
#                bs_ci = bs_ci, nreps_bs_ci = nreps_bs_ci)
# }

