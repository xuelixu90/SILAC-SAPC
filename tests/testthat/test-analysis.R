# Mini synthetic dataset mirroring the real workbook structure
make_mini_data <- function() {
  tibble::tibble(
    `Protein IDs` = c("P_ALB", "P_RES", "P_ACQ"),
    `Representive ID` = c("P_ALB", "P_RES", "P_ACQ"),
    `Protein names` = c("Albumin", "Residual", "Acquired"),
    `Gene names` = c("ALB", "RES", "ACQ"),
    In_house_db = c("in_house_serum", NA, NA),
    Plasma_Altas_db = c("PA001", NA, "PA003"),
    # control group CTR, 2 reps
    "Intensity L C1" = c(1000, 1000, NA),
    "Intensity H C1" = c(500, 500, NA),
    "Intensity L C2" = c(1100, 950, NA),
    "Intensity H C2" = c(550, 480, NA),
    # target group TGT, 2 reps
    "Intensity L T1" = c(2000, 1000, 1500),
    "Intensity H T1" = c(400, 480, NA),
    "Intensity L T2" = c(2100, 980, 1400),
    "Intensity H T2" = c(420, 500, NA),
    "Unique peptides C1" = c(5, 5, NA),
    "Unique peptides C2" = c(5, 5, NA),
    "Unique peptides T1" = c(5, 5, 3),
    "Unique peptides T2" = c(5, 5, 3)
  )
}

make_mini_sample_info <- function() {
  tibble::tibble(
    sample_id = c("C1", "C2", "T1", "T2"),
    group = c("CTR", "CTR", "TGT", "TGT"),
    light_column = paste("Intensity L", c("C1", "C2", "T1", "T2")),
    heavy_column = paste("Intensity H", c("C1", "C2", "T1", "T2")),
    peptide_column = paste("Unique peptides", c("C1", "C2", "T1", "T2"))
  )
}

test_that("numeric helpers behave as expected", {
  expect_equal(safe_median(c(1, NA, 3)), 2)
  expect_true(is.na(safe_median(c(NA, NA))))
  expect_equal(robust_sd(c(5, 5, 5)), 0)
  expect_true(is.na(robust_sd(5)))
  expect_equal(safe_fraction(1, 4), 0.25)
  expect_true(is.na(safe_fraction(1, 0)))
  expect_true(is.na(safe_fraction(NA, 4)))
  expect_false(get_support_rule(3, 1))
  expect_true(get_support_rule(3, 2))
  expect_true(get_support_rule(2, 2))
  expect_false(get_support_rule(1, 5))
})

test_that("params and categories factory", {
  p <- serum_source_params()
  expect_s3_class(p, "serum_source_params")
  expect_equal(p$min_absolute_light_excess, 0.05)
  expect_length(primary_serum_categories(), 4)
})

test_that("sample mapping contains en dash columns", {
  si <- build_sample_info()
  expect_true(any(grepl("\u2013", si$light_column)))
  expect_equal(nrow(si), 15)
})

test_that("prepare_long_data builds valid flags and fractions", {
  prep <- prepare_long_data(make_mini_data(), make_mini_sample_info())
  lh <- prep$lh_long

  expect_setequal(unique(lh$group), c("CTR", "TGT"))
  alb <- lh %>% dplyr::filter(Protein_Row_ID == "P_ALB", group == "CTR")
  expect_true(all(alb$Paired_LH_Valid_Rep))
  expect_true(all(abs(alb$Light_fraction - c(2 / 3, 2 / 3)) < 1e-6))

  acq <- lh %>% dplyr::filter(Protein_Row_ID == "P_ACQ", group == "CTR")
  expect_true(all(acq$Both_channels_zero_or_na))

  acq_t <- lh %>% dplyr::filter(Protein_Row_ID == "P_ACQ", group == "TGT")
  expect_true(all(acq_t$Target_L_Dominant_Rep))

  expect_true(all(c("Database_evidence", "Gene_names") %in% names(lh)))
})

test_that("compute_ucev_background returns floor", {
  prep <- prepare_long_data(make_mini_data(), make_mini_sample_info())
  bg <- compute_ucev_background(prep$lh_long, "CTR")
  expect_equal(bg$global_fraction_sd_floor, 0.02)
  # P_RES has two slightly different control fractions -> sd > 0;
  # P_ALB has identical fractions (sd = 0) and P_ACQ has no valid reps,
  # so exactly one protein passes the filter
  expect_equal(nrow(bg$sd_table), 1)
})

test_that("analyze_one_group classifies mini proteins correctly", {
  prep <- prepare_long_data(make_mini_data(), make_mini_sample_info())
  bg <- compute_ucev_background(prep$lh_long, "CTR")
  res <- analyze_one_group(
    prep$lh_long, prep$annotation, "CTR", "TGT",
    global_fraction_sd_floor = bg$global_fraction_sd_floor
  )

  expect_setequal(res$Protein_Row_ID, c("P_ALB", "P_RES", "P_ACQ"))

  alb <- res %>% dplyr::filter(Protein_Row_ID == "P_ALB")
  expect_equal(alb$classification, "Source-corrected, high confidence")
  expect_true(alb$CI_separated)
  expect_equal(round(alb$UCEV_heavy_labeling_p0, 3), 0.333)
  expect_equal(round(alb$estimated_serum_fraction, 3), 0.5)

  res_alb_num <- res %>% dplyr::filter(Protein_Row_ID == "P_ALB")
  expect_true(abs(res_alb_num$Light_fraction_excess - 1 / 6) < 1e-3)

  acq <- res %>% dplyr::filter(Protein_Row_ID == "P_ACQ")
  expect_equal(acq$classification, "Acquired-only (Level 1)")
  expect_true(acq$Control_Is_Absent)

  resi <- res %>% dplyr::filter(Protein_Row_ID == "P_RES")
  expect_true(resi$classification %in% c("Ambiguous", "Likely_residual"))
})

test_that("calibrate_serum_light adds calibrated columns", {
  prep <- prepare_long_data(make_mini_data(), make_mini_sample_info())
  bg <- compute_ucev_background(prep$lh_long, "CTR")
  res <- analyze_one_group(
    prep$lh_long, prep$annotation, "CTR", "TGT",
    global_fraction_sd_floor = bg$global_fraction_sd_floor
  )
  full <- calibrate_serum_light(res, prep$data_converted,
                                make_mini_sample_info(), "TGT", "CTR")
  expect_true("Serum_L_Calibrated T1" %in% names(full))
  expect_true("Serum_L_Calibrated T2" %in% names(full))

  # ACQ (control absent) keeps raw light intensity
  acq <- full %>% dplyr::filter(Protein_Row_ID == "P_ACQ")
  expect_equal(acq[["Serum_L_Calibrated T1"]], 1500)
})

test_that("plot functions return ggplot objects", {
  prep <- prepare_long_data(make_mini_data(), make_mini_sample_info())
  bg <- compute_ucev_background(prep$lh_long, "CTR")
  res <- analyze_one_group(
    prep$lh_long, prep$annotation, "CTR", "TGT",
    global_fraction_sd_floor = bg$global_fraction_sd_floor
  )
  full <- calibrate_serum_light(res, prep$data_converted,
                                make_mini_sample_info(), "TGT", "CTR")
  cs <- res %>% dplyr::count(comparison, assigned_case, classification,
                             name = "protein_count")

  expect_s3_class(plot_classification_summary(cs), "ggplot")
  expect_s3_class(plot_excess_distribution(res), "ggplot")
  expect_s3_class(plot_serum_fraction(full), "ggplot")
  expect_s3_class(plot_database_evidence(full), "ggplot")
  expect_s3_class(plot_excess_vs_probability(res), "ggplot")
  expect_s3_class(plot_excess_heatmap(res), "ggplot")
  expect_s3_class(plot_ucev_background(
    tibble::tibble(UCEV_Light_fraction_sd = c(0.01, 0.03, 0.05)), 0.02
  ), "ggplot")
})
