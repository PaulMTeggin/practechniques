test_that("pnorm by character", {
  uniparam_pnorm <- build_uniparam_fn("p", "norm")
  actual <- uniparam_pnorm(0.1, list(mean = 0, sd = 1))
  expect_equal(actual, pnorm(0.1, 0, 1))
})

test_that("pnorm by fn", {
  uniparam_pnorm <- build_uniparam_fn(fn = pnorm)
  actual <- uniparam_pnorm(0.1, list(mean = 0, sd = 1))
  expect_equal(actual, pnorm(0.1, 0, 1))
})

test_that("qnorm by character", {
  uniparam_qnorm <- build_uniparam_fn("q", "norm")
  actual <- uniparam_qnorm(0.1, list(mean = 0, sd = 1))
  expect_equal(actual, qnorm(0.1, 0, 1))
})

test_that("qnorm by fn", {
  uniparam_qnorm <- build_uniparam_fn(fn = qnorm)
  actual <- uniparam_qnorm(0.1, list(mean = 0, sd = 1))
  expect_equal(actual, qnorm(0.1, 0, 1))
})
