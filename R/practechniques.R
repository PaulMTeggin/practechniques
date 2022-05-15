#' practechniques - Practical Actuarial Techniques
#'
#' A companion R package for the online book Practical Actuarial Techniques (in development).
#'
#' THIS IS A PRE-RELEASE VERSION AND NO RELIANCE WHATSOEVER MAY BE PLACED ON IT.
#'
#' At the time of writing (May 2022) the main feature is a framework
#' using Monte Carlo simulation to calculate p-values for statistical goodness-of-fit tests
#' where the parameters have been estimated from the data.
#'
#' This is a generalisation of the approach in the
#' [`KScorrect`](https://cran.r-project.org/web/packages/KScorrect/index.html) package.
#' The approach has been generalised in three ways:
#'
#' * General distributions are supported, rather than the closed list used by `KScorrect`.
#' * Multiple statistical tests are supported, not just [KS](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test).
#' * Testing can be performed against distributions fitted to overlapping data, not just IID data,
#' using the idea of a Gaussian copula to induce autocorrelation consistent with overlapping data suggested in section 4.2 of the
#' [2019 paper](https://www.cambridge.org/core/journals/british-actuarial-journal/article/calibration-of-var-models-with-overlapping-data/B20D66D81DB918AFD3BBDF9EDAC20863)
#' by the [Extreme Events Working Party](https://www.actuaries.org.uk/practice-areas/life/research-working-parties/extreme-events)
#' of the UK [Institute and Faculty of Actuaries](https://www.actuaries.org.uk/).
#'
#' The function \code{\link{gof_test_sim}} is the main entry point and abstracts some of the implementation details.
#' \code{\link{gof_test_sim_uniparam}} is the lower-level function that performs the calculations.
#'
#' @docType package
#' @name practechniques
NULL
