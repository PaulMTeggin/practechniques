test_that("noverlap = 1", {
    ## Arrange - nothing to do

    ## Act
    actual <- build_theo_sercor_mtx(12, noverlap = 1)

    ## Assert
    expected <- diag(nrow = 12)
    expect_equal(actual, expected)
})

test_that("noverlap = 2", {
    ## Arrange - nothing to do

    ## Act
    actual <- build_theo_sercor_mtx(12, noverlap = 2)
    expected <- diag(nrow = 12)
    for (i in 1:11) {
        expected[i + 1, i] <- 0.5
        expected[i, i + 1] <- 0.5
    }

    ## Assert
    expect_equal(actual, expected)
})

test_that("noverlap = 12", {
    ## Arrange - nothing to do

    ## Act
    actual <- build_theo_sercor_mtx(12, noverlap = 12)
    expected <- diag(nrow = 12)
    for (i in 2:12) {
        sercor <- (12 - (i - 1)) / 12
        for (j in seq_len(13 - i)) {
            expected[i + j - 1, j] <- sercor
            expected[j, i + j - 1] <- sercor
        }
    }

    ## Assert
    expect_equal(actual, expected)
})

