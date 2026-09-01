# ============================================================
# End-to-end pipeline wrapper
# ============================================================

#' Run the complete serum-source analysis pipeline
#'
#' Reproduces the full workflow of the original hard-coded script:
#' reads the annotated MaxQuant workbook, builds long-format Heavy/Light
#' data, estimates UCEV background variation, classifies proteins for every
#' target group, performs source correction, writes all CSV/XLSX outputs,
#' generates key visualizations and writes a Markdown report.
#'
#' @param input_file path to the annotated Excel workbook
#' @param sheet_name sheet index/name (default 1)
#' @param output_dir output directory (created if needed)
#' @param sample_info sample mapping (default [build_sample_info()])
#' @param reference_group control group (default "UCEV")
#' @param target_groups target groups (default the four NP_1 groups)
#' @param params parameters from [serum_source_params()]
#' @param save_outputs write CSV/XLSX outputs (default TRUE)
#' @param make_plots write PNG figures into `output_dir/plots` (default TRUE)
#' @param write_report write `Analysis_Report.md` into `output_dir`
#'   (default TRUE)
#' @param render_notebook logical. If TRUE, render the package's
#'   `analysis_notebook.Rmd` into `output_dir/Analysis_Notebook.nb.html`
#'   and `Analysis_Notebook.md` so the full execution log (code, messages,
#'   intermediate tables, every plot) can be opened directly in a browser
#'   or Markdown viewer (default TRUE).
#' @return (invisibly) a list with elements `lh_long`, `annotation`,
#'   `data_converted`, `ucev_sd`, `global_fraction_sd_floor`, `results`
#'   (all comparisons), `all_merged`, `primary_serum`, `case_summary`,
#'   `database_summary`, `plots` (named list of ggplot objects), and
#'   (if rendered) `notebook_output` with paths to the generated HTML and
#'   Markdown notebook files.
#' @export
run_serum_source_analysis <- function(input_file,
                                      sheet_name = 1,
                                      output_dir,
                                      sample_info = build_sample_info(),
                                      reference_group = "UCEV",
                                      target_groups = c(
                                        "DGEV_PC_1-3",
                                        "DGEV_PC_10-12",
                                        "DGEV_PC_4-9",
                                        "UCEVPC"
                                      ),
                                      params = serum_source_params(),
                                      save_outputs = TRUE,
                                      make_plots = TRUE,
                                      write_report = TRUE,
                                      render_notebook = TRUE) {
  message("Reading input file...")
  raw_data <- read_input_data(input_file, sheet_name)
  message("Input dimensions: ", nrow(raw_data), " rows x ", ncol(raw_data), " columns")

  validate_sample_info(sample_info, raw_data)

  prepared <- prepare_long_data(
    raw_data,
    sample_info,
    min_light_dominance_ratio = params$min_light_dominance_ratio
  )

  ucev <- compute_ucev_background(
    prepared$lh_long,
    reference_group,
    params$min_valid_reps,
    params$minimum_fraction_sd_floor
  )

  message("Starting serum-source analysis for all groups...")
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

  plots <- list()
  if (make_plots) {
    plots <- list(
      classification_summary = plot_classification_summary(case_summary),
      excess_distribution = plot_excess_distribution(all_results),
      serum_fraction = plot_serum_fraction(all_primary),
      database_evidence = plot_database_evidence(all_primary),
      ucev_background = plot_ucev_background(
        ucev$sd_table, ucev$global_fraction_sd_floor
      ),
      excess_vs_probability = plot_excess_vs_probability(all_results),
      excess_heatmap = plot_excess_heatmap(all_results)
    )
  }

  if (save_outputs) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    save_analysis_outputs(
      output_dir = output_dir,
      sample_info = sample_info,
      annotation = prepared$annotation,
      ucev_sd = ucev$sd_table,
      all_merged = all_merged,
      all_primary = all_primary,
      case_summary = case_summary,
      database_summary = database_summary,
      target_groups = target_groups,
      reference_group = reference_group,
      primary = primary
    )
  }

  if (make_plots && save_outputs) {
    save_analysis_plots(plots, file.path(output_dir, "plots"))
  }

  if (write_report && save_outputs) {
    write_analysis_report(
      output_dir = output_dir,
      input_file = input_file,
      params = params,
      case_summary = case_summary,
      database_summary = database_summary,
      all_primary = all_primary,
      all_results = all_results,
      ucev_sd = ucev$sd_table,
      global_fraction_sd_floor = ucev$global_fraction_sd_floor
    )
  }

  notebook_output <- NULL
  if (render_notebook && save_outputs) {
    notebook_output <- render_analysis_notebook(
      input_file  = input_file,
      output_dir  = output_dir,
      sample_info = sample_info,
      params      = params
    )
  }

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
    database_summary = database_summary,
    plots = plots,
    notebook_output = notebook_output
  ))
}
