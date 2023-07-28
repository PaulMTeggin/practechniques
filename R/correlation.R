calc_partial_correl_pair <- function(risk_pair_df, threshold_p, flip_sign = FALSE) {
  stopifnot(ncol(risk_pair_df) == 2)
  stopifnot(length(threshold_p) == 1)

  risk_pair_df <- risk_pair_df[stats::complete.cases(risk_pair_df),]
  if (flip_sign) flip_mult <- -1 else flip_mult <- 1
  risk_pair_df[,1] <- flip_mult * risk_pair_df[,1]

  threshold_q <- vapply(risk_pair_df, stats::quantile, numeric(1), probs = threshold_p)

  keep_idx <- which(risk_pair_df[,1] >= threshold_q[1] & risk_pair_df[,2] >= threshold_q[2])
  npoints_keep <- length(keep_idx)

  keep_df <- risk_pair_df[keep_idx,]
  cor_res <- flip_mult * stats::cor(x = keep_df[,1], y = keep_df[,2], method = "spearman") # TODO apply S-P adj

  return(list(risk_1 = colnames(risk_pair_df)[1], risk_2 = colnames(risk_pair_df)[2],
              threshold_p = threshold_p, npoints = npoints_keep, flip_mult = flip_mult, cor = cor_res))
}

calc_all_partial_correl_pair <- function(risk_pair_df, all_threshold_p = seq(0, 1, by = 0.01)) {
  res_not_flipped <- purrr::map_dfr(all_threshold_p, calc_partial_correl_pair,
                              risk_pair_df = risk_pair_df, flip_sign = FALSE)

  res_flipped <- purrr::map_dfr(all_threshold_p, calc_partial_correl_pair,
                          risk_pair_df = risk_pair_df, flip_sign = TRUE)

  res <- rbind(res_not_flipped, res_flipped)
  colnames(res) <- c("risk_1", "risk_2", "threshold_p", "npoints", "flip_mult", "cor")
  res <- as.data.frame(res)
  res$flip_mult <- factor(res$flip_mult)
  return(res)
}

calc_all_partial_correl_df <- function(risk_df, all_threshold_p = seq(0, 1, by = 0.01)) {
  idx <- tidyr::expand_grid(r1 = seq_along(risk_df), r2 = seq_along(risk_df)) # All combos
  idx <- idx[idx[,1] < idx[,2],] # Triangle

  res <- purrr::pmap_dfr(idx, function(r1, r2) {
    risk_pair_df <- risk_df[c(r1, r2)]
    calc_all_partial_correl_pair(risk_pair_df, all_threshold_p)
  })

  return(res)
}
