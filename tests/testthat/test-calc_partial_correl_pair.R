test_that("smoke test", {
  # Arrange
  risk_pair_df <- as.data.frame(matrix(rnorm(200), nrow = 100))
  risk_df <- as.data.frame(matrix(rnorm(1800), nrow = 200))

  # Act
  actual <- calc_all_partial_correl_pair(risk_pair_df)
  actual_all <- calc_all_partial_correl_df(risk_df)

  # Assert

  # Switch the risk names round for plotting purposes when
  flipped_idx <- which(actual_all$flip_mult == -1)
  stash <- actual_all$risk_1[flipped_idx]
  actual_all$risk_1[flipped_idx] <- actual_all$risk_2[flipped_idx]
  actual_all$risk_2[flipped_idx] <- stash


  # ggplot2::ggplot(actual_all, ggplot2::aes(x = threshold_p, y = cor, group = flip_mult, colour = flip_mult)) +
  #   ggplot2::geom_line() +
  #   ggplot2::geom_hline(yintercept = 0) +
  #   ggplot2::facet_grid(rows = ggplot2::vars(risk_1), cols = ggplot2::vars(risk_2)) +
  #   ggplot2::theme_light()
  #
  # ggplot2::ggplot(actual_all, ggplot2::aes(x = threshold_p, y = npoints, group = flip_mult, colour = flip_mult)) +
  #   ggplot2::geom_line() +
  #   ggplot2::geom_hline(yintercept = 20) +
  #   ggplot2::facet_grid(rows = ggplot2::vars(risk_1), cols = ggplot2::vars(risk_2)) +
  #   ggplot2::theme_light()
})
