test_that("p = 1", {
    ## Arrange - nothing to do

    ## Act
    actual <- calc_theo_sercor(1, 0:12)
    expected <- c(1, rep(0, 12)) # 1 at lag zero, otherwise 0

    ## Assert
    expect_equal(actual, expected)
})

test_that("p = 2", {
    ## Arrange - nothing to do

    ## Act
    actual <- calc_theo_sercor(2, 0:12)
    expected <- c(1, 0.5, rep(0, 11)) # 1 at lag zero, otherwise 0

    ## Assert
    expect_equal(actual, expected)
})

test_that("p = 12", {
    ## Arrange - nothing to do

    ## Act
    actual <- calc_theo_sercor(12, 0:12)
    expected <- 12:0 / 12

    ## Assert
    expect_equal(actual, expected)
})
