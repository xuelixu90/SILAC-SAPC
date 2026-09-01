# ============================================================
# Core per-group serum-source classification
# ============================================================

#' Analyze one target group against the reference control group
#'
#' For every protein row, determines the Case (1-4 / Other), computes
#' replicate-level Light-fraction statistics, classifies the protein into
#' serum-source categories and estimates the source-corrected serum fraction.
#' Logic is a faithful port of the original hard-coded analysis script.
#'
#' @param lh_long long-format data from [prepare_long_data()]
#' @param annotation_data annotation table from [prepare_long_data()]
#' @param reference_group control group name (e.g. "UCEV")
#' @param target_group target group name (e.g. "UCEVPC")
#' @param params parameters from [serum_source_params()]
#' @param global_fraction_sd_floor background SD floor from
#'   [compute_ucev_background()]
#' @return tibble with one row per protein and all statistics
#' @export
analyze_one_group <- function(lh_long,
                              annotation_data,
                              reference_group,
                              target_group,
                              params = serum_source_params(),
                              global_fraction_sd_floor = params$minimum_fraction_sd_floor) {
  eps <- params$eps
  comparison_name <- paste0(target_group, "_vs_", reference_group)

  ref_data <- lh_long %>% dplyr::filter(group == reference_group)
  tgt_data <- lh_long %>% dplyr::filter(group == target_group)

  # Split once for O(1) lookup per protein (identical semantics to
  # filter-by-protein but much faster on large tables)
  ref_list <- split(ref_data, ref_data$Protein_Row_ID)
  tgt_list <- split(tgt_data, tgt_data$Protein_Row_ID)
  get_one <- function(lst, id, fallback) {
    x <- lst[[id]]
    if (is.null(x)) fallback else x
  }

  purrr::map_dfr(annotation_data$Protein_Row_ID, function(protein_id) {
    ref_one <- get_one(ref_list, protein_id, ref_data[0, ])
    tgt_one <- get_one(tgt_list, protein_id, tgt_data[0, ])

    n_valid_ref_lh <- sum(ref_one$Paired_LH_Valid_Rep, na.rm = TRUE)
    ref_control_absent <- (nrow(ref_one) > 0 &&
      sum(ref_one$Both_channels_zero_or_na, na.rm = TRUE) == nrow(ref_one))

    n_valid_tgt_lh <- sum(tgt_one$Paired_LH_Valid_Rep, na.rm = TRUE)
    n_tgt_high_pep_paired <- sum(tgt_one$High_Pep_Matched_Rep, na.rm = TRUE)
    n_tgt_single_pep_paired <- sum(tgt_one$Single_Pep_Matched_Rep, na.rm = TRUE)

    n_tgt_l_dom <- sum(tgt_one$Target_L_Dominant_Rep, na.rm = TRUE)
    n_tgt_high_pep_l_dom <- sum(tgt_one$High_Pep_L_Dominant_Rep, na.rm = TRUE)
    n_tgt_single_pep_l_dom <- sum(tgt_one$Single_Pep_L_Dominant_Rep, na.rm = TRUE)

    # ---- Case determination ----
    is_case1 <- (
      n_valid_ref_lh >= params$min_valid_reps &&
        n_valid_tgt_lh >= params$min_valid_reps &&
        ((n_valid_tgt_lh == 3 && n_tgt_high_pep_paired >= params$min_support_reps_when_three_valid) ||
           (n_valid_tgt_lh == 2 && n_tgt_high_pep_paired == params$min_support_reps_when_three_valid))
    )

    is_case2 <- (
      n_valid_ref_lh >= params$min_valid_reps &&
        n_valid_tgt_lh >= params$min_valid_reps &&
        !is_case1 &&
        ((n_valid_tgt_lh == 3 && n_tgt_single_pep_paired >= params$min_support_reps_when_three_valid) ||
           (n_valid_tgt_lh == 2 && n_tgt_single_pep_paired >= 1))
    )

    is_case3 <- (
      ref_control_absent &&
        n_tgt_l_dom >= params$min_valid_reps &&
        ((n_tgt_l_dom == 3 && n_tgt_high_pep_l_dom >= params$min_support_reps_when_three_valid) ||
           (n_tgt_l_dom == 2 && n_tgt_high_pep_l_dom == params$min_support_reps_when_three_valid))
    )

    is_case4 <- (
      ref_control_absent &&
        n_tgt_l_dom >= params$min_valid_reps &&
        !is_case3 &&
        ((n_tgt_l_dom == 3 && n_tgt_single_pep_l_dom >= params$min_support_reps_when_three_valid) ||
           (n_tgt_l_dom == 2 && n_tgt_single_pep_l_dom >= 1))
    )

    assigned_case <- dplyr::case_when(
      is_case1 ~ "Case 1",
      is_case2 ~ "Case 2",
      is_case3 ~ "Case 3",
      is_case4 ~ "Case 4",
      TRUE ~ "Other / Unassigned"
    )

    # ---- quantitative calculations ----
    L0 <- ref_one$Light
    H0 <- ref_one$Heavy
    r0 <- ref_one$Light_fraction
    Lx <- tgt_one$Light
    Hx <- tgt_one$Heavy
    rx <- tgt_one$Light_fraction

    p0 <- NA_real_
    baseline_light_fraction <- NA_real_
    target_light_fraction <- NA_real_
    light_fraction_excess <- NA_real_
    ref_sd_raw <- NA_real_
    tgt_sd_raw <- NA_real_
    sigma_used <- NA_real_
    high_threshold <- NA_real_
    medium_threshold <- NA_real_
    z_score <- NA_real_
    prob_excess <- NA_real_
    ci_sep <- NA
    supporting_high_reps <- NA_integer_
    supporting_med_reps <- NA_integer_
    est_serum_frac <- NA_real_
    serum_frac_cv <- NA_real_

    if (is_case3 || is_case4 || ref_control_absent) {
      # ---- control absent: acquired-only branches ----
      target_light_fraction <- safe_median(rx[tgt_one$Target_L_Dominant_Rep])
      tgt_sd_raw <- robust_sd(rx[tgt_one$Target_L_Dominant_Rep])
      if (!is.finite(tgt_sd_raw)) tgt_sd_raw <- global_fraction_sd_floor
      sigma_used <- global_fraction_sd_floor
      high_threshold <- max(
        params$min_absolute_light_excess,
        params$background_multiplier_high * sigma_used
      )
      medium_threshold <- max(
        params$min_absolute_light_excess,
        params$background_multiplier_medium * sigma_used
      )

      classification <- dplyr::case_when(
        is_case3 ~ "Acquired-only (Level 1)",
        is_case4 ~ "Acquired-only (Level 2)",
        TRUE ~ "Control_absent_other"
      )
    } else if (n_valid_ref_lh >= params$min_valid_reps &&
      n_valid_tgt_lh >= params$min_valid_reps) {
      # ---- both control and target present: source-corrected branches ----
      valid_ref_mask <- ref_one$Paired_LH_Valid_Rep
      valid_tgt_mask <- tgt_one$Paired_LH_Valid_Rep

      ref_heavy_frac_rep <- safe_fraction(H0, L0 + H0)
      p0 <- safe_median(ref_heavy_frac_rep[valid_ref_mask])
      baseline_light_fraction <- safe_median(r0[valid_ref_mask])
      target_light_fraction <- safe_median(rx[valid_tgt_mask])
      light_fraction_excess <- target_light_fraction - baseline_light_fraction

      ref_sd_raw <- robust_sd(r0[valid_ref_mask])
      tgt_sd_raw <- robust_sd(rx[valid_tgt_mask])
      if (!is.finite(ref_sd_raw)) ref_sd_raw <- global_fraction_sd_floor
      if (!is.finite(tgt_sd_raw)) tgt_sd_raw <- global_fraction_sd_floor

      sigma_used <- max(ref_sd_raw, global_fraction_sd_floor)
      high_threshold <- max(
        params$min_absolute_light_excess,
        params$background_multiplier_high * sigma_used
      )
      medium_threshold <- max(
        params$min_absolute_light_excess,
        params$background_multiplier_medium * sigma_used
      )

      pooled_se <- sqrt(
        (ref_sd_raw^2 / n_valid_ref_lh) +
          (tgt_sd_raw^2 / n_valid_tgt_lh) +
          global_fraction_sd_floor^2
      )
      if (is.finite(light_fraction_excess) && pooled_se > eps) {
        z_score <- light_fraction_excess / pooled_se
        prob_excess <- stats::pnorm(z_score)
      }

      ref_upper_95 <- baseline_light_fraction +
        1.96 * ref_sd_raw / sqrt(n_valid_ref_lh)
      tgt_lower_95 <- target_light_fraction -
        1.96 * tgt_sd_raw / sqrt(n_valid_tgt_lh)
      ci_sep <- is.finite(tgt_lower_95) && is.finite(ref_upper_95) &&
        tgt_lower_95 > ref_upper_95

      tgt_excess_rep <- rx - baseline_light_fraction
      supporting_high_reps <- sum(is.finite(tgt_excess_rep) &
        tgt_excess_rep >= high_threshold)
      supporting_med_reps <- sum(is.finite(tgt_excess_rep) &
        tgt_excess_rep >= medium_threshold)

      high_supp_rule <- get_support_rule(
        n_valid_tgt_lh, supporting_high_reps,
        params$min_support_reps_when_three_valid,
        params$min_support_reps_when_two_valid
      )
      med_supp_rule <- get_support_rule(
        n_valid_tgt_lh, supporting_med_reps,
        params$min_support_reps_when_three_valid,
        params$min_support_reps_when_two_valid
      )

      if (!is.finite(p0) || p0 <= eps || p0 >= 1) {
        serum_frac_rep <- rep(NA_real_, length(Hx))
      } else {
        A <- Hx / p0
        expected_residual_L <- A * (1 - p0)
        estimated_serum_L <- pmax(Lx - expected_residual_L, 0)
        serum_frac_rep <- safe_fraction(estimated_serum_L, A + estimated_serum_L)
      }

      est_serum_frac <- safe_median(serum_frac_rep[valid_tgt_mask])
      serum_frac_cv <- safe_cv(serum_frac_rep[valid_tgt_mask])

      strong_quant <- (
        is.finite(light_fraction_excess) &&
          light_fraction_excess >= high_threshold &&
          high_supp_rule &&
          is.finite(serum_frac_cv) &&
          serum_frac_cv <= params$max_serum_fraction_cv_high_confidence
      )

      med_quant <- (
        is.finite(light_fraction_excess) &&
          light_fraction_excess >= medium_threshold &&
          med_supp_rule
      )

      classification <- dplyr::case_when(
        is_case1 && strong_quant && ci_sep ~ "Source-corrected, high confidence",
        is_case1 && med_quant ~ "Source-corrected, medium confidence",
        is_case2 && strong_quant && ci_sep ~
          "Source-corrected, single-peptide strong candidate",
        is_case2 && med_quant ~
          "Source-corrected, single-peptide candidate",
        is.finite(light_fraction_excess) && light_fraction_excess <= 0 ~
          "Likely_residual",
        TRUE ~ "Ambiguous"
      )
    } else {
      classification <- ifelse(
        n_valid_ref_lh < params$min_valid_reps,
        "Control_missing",
        "Target_group_missing"
      )
    }

    tibble::tibble(
      Protein_Row_ID = protein_id,
      comparison = comparison_name,
      assigned_case = assigned_case,
      classification = classification,

      Target_Both_LH_Valid_Reps = n_valid_tgt_lh,
      Control_Both_LH_Valid_Reps = n_valid_ref_lh,

      Target_L_Dominant_Reps = n_tgt_l_dom,
      Target_High_Pep_L_Dom_Reps = n_tgt_high_pep_l_dom,
      Target_Single_Pep_L_Dom_Reps = n_tgt_single_pep_l_dom,

      Control_Is_Absent = ref_control_absent,

      UCEV_heavy_labeling_p0 = p0,
      UCEV_baseline_Light_fraction = baseline_light_fraction,
      target_Light_fraction = target_light_fraction,
      Light_fraction_excess = light_fraction_excess,

      sigma_used = sigma_used,
      high_threshold = high_threshold,
      medium_threshold = medium_threshold,

      z_score = z_score,
      probability_excess_Light = prob_excess,
      CI_separated = ci_sep,

      estimated_serum_fraction = est_serum_frac,
      serum_fraction_CV = serum_frac_cv
    )
  }) %>%
    dplyr::left_join(annotation_data, by = "Protein_Row_ID")
}

#' Calibrate serum Light intensity for each target replicate
#'
#' Subtracts the UCEV baseline Light fraction from each target replicate and
#' converts the excess back to intensity units, adding one
#' `Serum_L_Calibrated_<sample>` column per target replicate.
#'
#' @param cmp_results results tibble from [analyze_one_group()]
#' @param data_converted converted wide data from [prepare_long_data()]
#' @param sample_info sample mapping tibble
#' @param target_group target group name
#' @param reference_group reference group name
#' @return cmp_results joined with raw intensity columns and the calibrated
#'   serum Light columns
#' @export
calibrate_serum_light <- function(cmp_results,
                                  data_converted,
                                  sample_info,
                                  target_group,
                                  reference_group) {
  eps <- 1e-12

  tgt_sample_table <- sample_info %>% dplyr::filter(group == target_group)
  ref_sample_table <- sample_info %>% dplyr::filter(group == reference_group)

  output_raw_columns <- unique(c(
    tgt_sample_table$light_column,
    tgt_sample_table$heavy_column,
    tgt_sample_table$peptide_column,
    ref_sample_table$light_column,
    ref_sample_table$heavy_column,
    ref_sample_table$peptide_column
  ))

  pair_raw_data <- data_converted %>%
    dplyr::select(Protein_Row_ID, dplyr::all_of(output_raw_columns))

  cmp_full <- cmp_results %>%
    dplyr::left_join(pair_raw_data, by = "Protein_Row_ID")

  for (i in seq_len(nrow(tgt_sample_table))) {
    sample_id <- tgt_sample_table$sample_id[i]
    light_col <- tgt_sample_table$light_column[i]
    heavy_col <- tgt_sample_table$heavy_column[i]

    output_column <- paste("Serum_L_Calibrated", sample_id)

    light_value <- cmp_full[[light_col]]
    heavy_value <- cmp_full[[heavy_col]]
    heavy_value_fixed <- ifelse(is.na(heavy_value), 0, heavy_value)
    total_value <- light_value + heavy_value_fixed

    light_fraction_rep <- ifelse(
      is.finite(light_value) & light_value > 0 &
        is.finite(total_value) & total_value > 0,
      light_value / total_value,
      NA_real_
    )

    baseline_fraction <- cmp_full$UCEV_baseline_Light_fraction

    calibrated_value <- rep(NA_real_, nrow(cmp_full))

    control_absent_mask <- (
      cmp_full$Control_Is_Absent &
        is.finite(light_value) &
        light_value > 0
    )
    calibrated_value[control_absent_mask] <- light_value[control_absent_mask]

    calibration_mask <- (
      !control_absent_mask &
        is.finite(light_fraction_rep) &
        light_fraction_rep > eps &
        is.finite(baseline_fraction)
    )

    delta_fraction <- light_fraction_rep - baseline_fraction
    calibrated_value[calibration_mask] <- pmax(
      (delta_fraction[calibration_mask] / light_fraction_rep[calibration_mask]) *
        light_value[calibration_mask],
      0
    )

    cmp_full[[output_column]] <- calibrated_value
  }

  cmp_full
}
