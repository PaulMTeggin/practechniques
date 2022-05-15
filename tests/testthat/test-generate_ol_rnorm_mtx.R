test_that("1m rows", {
    ## Arrange
    set.seed(1)

    ## Act
    actual <- generate_ol_rnorm_mtx(1e6, 24, 12)
    actual_mean <- apply(actual, 2, mean)
    actual_sd <- apply(actual, 2, sd)
    actual_cor <- cor(actual)

    ## Assert

    expect_true(max(abs(actual_mean)) < 0.002)
    expect_true(max(abs(actual_sd) - 1) < 0.002)

    expected_cor <- build_theo_sercor_mtx(24, 12)
    expect_true(max(abs(actual_cor - expected_cor)) < 0.002)
})
