# ============================================================
# Core computation pipeline (no downstream I/O)
# ============================================================

#' Run the core serum-source computation pipeline
#'
#' Reads the annotated MaxQuant workbook, builds long-format Heavy/Light
#' data, estimates UCEV background variation, classifies proteins for every
#' target group and performs source correction.
#'
#' This function performs **no** downstream I/O: it does not generate plots,
#' write CSV/XLSX files, produce a Markdown report, or render a notebook.
#' Use the accompanying downstream analysis script for those steps.
#'
#' @param input_file path to the annotated Excel workbook
#' @param sheet_name sheet index/name (default 1)
#' @param sample_info sample mapping (default [build_sample_info()])
#' @param reference_group control group name (default "UCEV"). The target
#'   groups are taken from `sample_info` excluding this group when
#'   `target_groups` is `NULL`.
#' @param target_groups character vector of target group names, or `NULL`
#'   (default) to use all groups in `sample_info` except `reference_group`.
#' @param params parameters from [serum_source_params()]
#' @param annotation_cols named character vector passed to
#'   [prepare_long_data()] mapping logical annotation keys to workbook column
#'   names (defaults match the NP_1 MaxQuant workbook).
#' @param na character vector of values to treat as missing when reading the
#'   workbook (passed to [read_input_data()]).
#' @return a list with elements `lh_long`, `annotation`, `data_converted`,
#'   `ucev_sd`, `global_fraction_sd_floor`, `results` (all comparisons,
#'   without calibrated columns), `all_merged` (with calibrated columns),
#'   `primary_serum`, `case_summary`, `database_summary`.
#' @export
run_core_analysis <- function(input_file,
                              sheet_name = 1,
                              sample_info = build_sample_info(),
                              reference_group = "UCEV",
                              target_groups = NULL,
                              params = serum_source_params(),
                              annotation_cols = c(
                                protein_ids = "Protein IDs",
                                representative_id = "Representive ID",
                                protein_names = "Protein names",
                                gene_names = "Gene names",
                                in_house_db = "In_house_db",
                                plasma_atlas_db = "Plasma_Altas_db"
                              ),
                              na = c("", " ", "NA", "NaN", "#N/A", "#DIV/0!")) {
  # Auto-infer target groups from the sample mapping if not supplied
  if (is.null(target_groups)) {
    all_groups <- unique(sample_info$group)
    target_groups <- setdiff(all_groups, reference_group)
    if (length(target_groups) == 0) {
      stop("No target groups found in sample_info after excluding '",
           reference_group, "'.", call. = FALSE)
    }
  }

  message("Reading input file...")
  raw_data <- read_input_data(input_file, sheet_name, na = na)
  message("Input dimensions: ", nrow(raw_data), " rows x ", ncol(raw_data), " columns")

  validate_sample_info(sample_info, raw_data)

  prepared <- prepare_long_data(
    raw_data,
    sample_info,
    min_light_dominance_ratio = params$min_light_dominance_ratio,
    annotation_cols = annotation_cols
  )

  ucev <- compute_ucev_background(
    prepared$lh_long,
    reference_group,
    params$min_valid_reps,
    params$minimum_fraction_sd_floor
  )

  message("Starting serum-source analysis for all groups...")
  message("Reference group : ", reference_group)
  message("Target groups   : ", paste(target_groups, collapse = ", "))
  results_list <- list()
  full_list <- list()

  for (tgt_grp in target_groups) {
    comparison_name <- paste0(tgt_grp, "_vs_", reference_group)
    message("--> Processing: ", comparison_name)

    cmp_results <- analyze_one_group(
      lh_long = prepared$lh_long,
      annotation_data = prepared$annotation,
      reference_group = reference_group,
      target_group = tgt_grp,
      params = params,
      global_fraction_sd_floor = ucev$global_fraction_sd_floor
    )

    cmp_full <- calibrate_serum_light(
      cmp_results = cmp_results,
      data_converted = prepared$data_converted,
      sample_info = sample_info,
      target_group = tgt_grp,
      reference_group = reference_group
    )

    results_list[[comparison_name]] <- cmp_results
    full_list[[comparison_name]] <- cmp_full
  }

  all_merged <- dplyr::bind_rows(full_list)
  all_results <- dplyr::bind_rows(results_list)
  primary <- primary_serum_categories()
  all_primary <- all_merged %>% dplyr::filter(classification %in% primary)

  case_summary <- all_results %>%
    dplyr::count(comparison, assigned_case, classification, name = "protein_count") %>%
    dplyr::arrange(comparison, assigned_case, classification)

  database_summary <- all_merged %>%
    dplyr::count(comparison, Database_evidence, classification, name = "protein_count") %>%
    dplyr::arrange(comparison, classification, dplyr::desc(protein_count))

  invisible(list(
    lh_long = prepared$lh_long,
    annotation = prepared$annotation,
    data_converted = prepared$data_converted,
    ucev_sd = ucev$sd_table,
    global_fraction_sd_floor = ucev$global_fraction_sd_floor,
    results = all_results,
    all_merged = all_merged,
    primary_serum = all_primary,
    case_summary = case_summary,
    database_summary = database_summary
  ))
}
