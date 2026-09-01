# ============================================================
# Parameter container for the serum-source analysis
# ============================================================

#' Create the parameter set for serum-source analysis
#'
#' All thresholds used by [analyze_one_group()] and
#' [run_serum_source_analysis()] are collected here.
#'
#' @param min_valid_reps minimum paired valid replicates required (default 2)
#' @param min_absolute_light_excess absolute Light-fraction excess floor for
#'   calling an excess (default 0.05)
#' @param background_multiplier_high multiplier of background SD defining the
#'   high-confidence threshold (default 3)
#' @param background_multiplier_medium multiplier of background SD defining
#'   the medium-confidence threshold (default 2)
#' @param min_support_reps_when_three_valid supporting replicates required
#'   when 3 valid replicates exist (default 2)
#' @param min_support_reps_when_two_valid supporting replicates required when
#'   2 valid replicates exist (default 2)
#' @param max_serum_fraction_cv_high_confidence maximum CV of the estimated
#'   serum fraction for high-confidence calls (default 0.50)
#' @param minimum_fraction_sd_floor lower floor applied to Light-fraction SDs
#'   (default 0.02)
#' @param min_light_dominance_ratio Light/(Light+Heavy) ratio above which a
#'   replicate is considered Light-dominant even without a Heavy channel
#'   (default 0.70)
#' @param eps numerical tolerance (default 1e-12)
#' @return A named list of parameters with class `serum_source_params`.
#' @export
#' @examples
#' p <- serum_source_params()
#' p$min_absolute_light_excess
serum_source_params <- function(
    min_valid_reps = 2,
    min_absolute_light_excess = 0.05,
    background_multiplier_high = 3,
    background_multiplier_medium = 2,
    min_support_reps_when_three_valid = 2,
    min_support_reps_when_two_valid = 2,
    max_serum_fraction_cv_high_confidence = 0.50,
    minimum_fraction_sd_floor = 0.02,
    min_light_dominance_ratio = 0.70,
    eps = 1e-12) {
  structure(
    list(
      min_valid_reps = min_valid_reps,
      min_absolute_light_excess = min_absolute_light_excess,
      background_multiplier_high = background_multiplier_high,
      background_multiplier_medium = background_multiplier_medium,
      min_support_reps_when_three_valid = min_support_reps_when_three_valid,
      min_support_reps_when_two_valid = min_support_reps_when_two_valid,
      max_serum_fraction_cv_high_confidence = max_serum_fraction_cv_high_confidence,
      minimum_fraction_sd_floor = minimum_fraction_sd_floor,
      min_light_dominance_ratio = min_light_dominance_ratio,
      eps = eps
    ),
    class = "serum_source_params"
  )
}

#' Classification categories treated as primary serum-derived evidence
#'
#' @return character vector of the four primary classification labels
#' @export
#' @examples
#' primary_serum_categories()
primary_serum_categories <- function() {
  c(
    "Acquired-only (Level 1)",
    "Acquired-only (Level 2)",
    "Source-corrected, high confidence",
    "Source-corrected, medium confidence"
  )
}
