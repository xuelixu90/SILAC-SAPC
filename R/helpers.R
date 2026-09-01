# ============================================================
# Numerical helper functions used across the analysis
# ============================================================

#' Safe median that ignores non-finite values
#'
#' @param x numeric vector
#' @return median of finite values, or NA_real_ if none
#' @keywords internal
safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

#' Robust standard deviation based on MAD
#'
#' @param x numeric vector
#' @return MAD-scaled SD of finite values, or NA_real_ if < 2 values
#' @keywords internal
robust_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  mad(x, center = median(x), constant = 1.4826, na.rm = TRUE)
}

#' Safe coefficient of variation
#'
#' @param x numeric vector
#' @return CV of finite values, or NA_real_ if not computable
#' @keywords internal
safe_cv <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  x_mean <- mean(x)
  if (!is.finite(x_mean) || x_mean <= 1e-12) return(NA_real_)
  sd(x) / x_mean
}

#' Safe numerator / denominator ratio
#'
#' @param numerator numeric vector
#' @param denominator numeric vector
#' @return ratio where both inputs are finite and denominator > 0, else NA
#' @keywords internal
safe_fraction <- function(numerator, denominator) {
  result <- rep(NA_real_, length(numerator))
  valid <- is.finite(numerator) & is.finite(denominator) & denominator > 1e-12
  result[valid] <- numerator[valid] / denominator[valid]
  result
}

#' Replicate-support decision rule
#'
#' Requires at least `min_support_when_three_valid` supporting
#' replicates when 3 valid replicates exist, and `min_support_when_two_valid`
#' when 2 exist.
#'
#' @param n_valid_target number of valid target replicates (2 or 3)
#' @param supporting_reps number of supporting replicates
#' @param min_support_when_three_valid required supports for 3 valid reps
#' @param min_support_when_two_valid required supports for 2 valid reps
#' @return logical(1)
#' @keywords internal
get_support_rule <- function(n_valid_target, supporting_reps,
                             min_support_when_three_valid = 2,
                             min_support_when_two_valid = 2) {
  if (n_valid_target == 3) {
    return(supporting_reps >= min_support_when_three_valid)
  }
  if (n_valid_target == 2) {
    return(supporting_reps >= min_support_when_two_valid)
  }
  FALSE
}
