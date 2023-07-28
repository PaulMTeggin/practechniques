test_that("basic tests", {
  expect_identical(get_prefixed_dist_fn("p", "norm"), stats::pnorm)
  expect_identical(get_prefixed_dist_fn("q", "norm"), stats::qnorm)
  expect_identical(get_prefixed_dist_fn("r", "norm"), stats::rnorm)
  expect_identical(get_prefixed_dist_fn("d", "norm"), stats::dnorm)
  expect_identical(get_prefixed_dist_fn("d", "logis"), stats::dlogis)
})
