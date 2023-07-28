test_that("pnorm", {
  x_data <- rnorm(100)
  pnorm_est <- build_estimated_fn("p", "norm", estimate_mean_sd_ol, x_data)
  x <- seq(-0.5, 0.5, length.out = 11)
  actual <- pnorm_est(x)
  expect_equal(actual, pnorm(x, mean(x_data), sd(x_data)))
})

test_that("qlogis", {
  x_data <- rlogis(100)
  qlogis_est <- build_estimated_fn("q", "logis", estimate_logis_ol, x_data)
  x <- seq(0.1, 0.9, by = 0.1)
  actual <- qlogis_est(x)
  expect_equal(actual, qlogis(x, mean(x_data), sd(x_data) * sqrt(3) / pi))
})

test_that("qlogis, multi", {
  df_data <- data.frame(PC = c(rep("PC1", 100), rep("PC2", 100), rep("PC3", 100)),
                        x_data = rlogis(300), stringsAsFactors = FALSE)

  df_est <- df_data %>%
    dplyr::group_by(PC) %>%
    dplyr::summarise(fn_p = list(build_estimated_fn("p", "logis", estimate_logis_ol, x_data)),
                     fn_q = list(build_estimated_fn("q", "logis", estimate_logis_ol, x_data)),
                     ks = list(gof_test_sim(x_data, test_type = "KS", dist = "logis", fn_estimate_params = estimate_logis_ol)),
                     ad = list(gof_test_sim(x_data, test_type = "AD", dist = "logis", fn_estimate_params = estimate_logis_ol)))

  expect_equal(df_est$fn_p[[1]](0.25), plogis(0.25, mean(df_data$x_data[1:100]), sd(df_data$x_data[1:100]) * sqrt(3) / pi))
  expect_equal(df_est$fn_p[[2]](0.5), plogis(0.5, mean(df_data$x_data[101:200]), sd(df_data$x_data[101:200]) * sqrt(3) / pi))
  expect_equal(df_est$fn_p[[3]](0.75), plogis(0.75, mean(df_data$x_data[201:300]), sd(df_data$x_data[201:300]) * sqrt(3) / pi))

  expect_equal(df_est$fn_q[[1]](0.25), qlogis(0.25, mean(df_data$x_data[1:100]), sd(df_data$x_data[1:100]) * sqrt(3) / pi))
  expect_equal(df_est$fn_q[[2]](0.5), qlogis(0.5, mean(df_data$x_data[101:200]), sd(df_data$x_data[101:200]) * sqrt(3) / pi))
  expect_equal(df_est$fn_q[[3]](0.75), qlogis(0.75, mean(df_data$x_data[201:300]), sd(df_data$x_data[201:300]) * sqrt(3) / pi))
})

