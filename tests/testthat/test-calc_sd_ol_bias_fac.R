# Just a boundary test to exercise the code - proof this is the correct formula is by simulation
test_that("p = 1", {
    ## Arrange - nothing to do

    ## Act
    actual <- calc_sd_ol_bias_fac(c(24, 60, 120, 180, 240), 1)

    ## Assert
    expected <- rep(1, 5)
    expect_equal(actual, expected)
})
