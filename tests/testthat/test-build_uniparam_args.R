test_that("params as list", {
  actual <- build_uniparam_args(1:10, list(mean = 0, sd = 1))
  expect_equal(actual, list(1:10, mean = 0, sd = 1))
})

test_that("params as object", {
  actual <- build_uniparam_args(1:10, matrix(0))
  expect_equal(actual, list(1:10, matrix(0)))
})
