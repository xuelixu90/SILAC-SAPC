# ============================================================
# Markdown analysis report generator
# ============================================================

# Convert a data.frame to a pipe-table string
md_table <- function(df, digits = 3) {
  if (nrow(df) == 0) return("_No rows._")
  num_cols <- vapply(df, is.numeric, logical(1))
  df[num_cols] <- lapply(df[num_cols], function(x) {
    ifelse(is.na(x), "NA",
      ifelse(x == round(x) & abs(x) < 1e15,
        formatC(x, format = "d", big.mark = ","),
        formatC(x, format = "g", digits = digits)
      )
    )
  })
  lines <- c(
    paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  )
  for (i in seq_len(nrow(df))) {
    lines <- c(lines, paste0("| ", paste(as.character(unlist(df[i, ])),
      collapse = " | "
    ), " |"))
  }
  paste(lines, collapse = "\n")
}

#' Write the Markdown analysis report
#'
#' Produces `Analysis_Report.md` summarizing methods, parameters,
#' classification results, database evidence, top candidate proteins and
#' embedded figures from the `plots/` directory.
#'
#' @param output_dir output directory containing `plots/`
#' @param input_file input workbook path (for provenance)
#' @param params parameters from [serum_source_params()]
#' @param case_summary,all_primary,database_summary analysis tables
#' @param all_results merged per-protein results (all comparisons)
#' @param ucev_sd UCEV background table
#' @param global_fraction_sd_floor global SD floor
#' @return invisibly, the path of the written report
#' @export
write_analysis_report <- function(output_dir,
                                  input_file,
                                  params,
                                  case_summary,
                                  all_primary,
                                  database_summary,
                                  all_results,
                                  ucev_sd,
                                  global_fraction_sd_floor) {
  primary <- primary_serum_categories()

  # wide classification count table
  cls_wide <- case_summary %>%
    dplyr::group_by(comparison, classification) %>%
    dplyr::summarise(n = sum(protein_count), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = classification, values_from = n, values_fill = 0)

  primary_counts <- all_primary %>%
    dplyr::count(comparison, name = "Primary_serum_proteins")

  top_proteins <- all_primary %>%
    dplyr::group_by(comparison) %>%
    dplyr::arrange(dplyr::desc(Light_fraction_excess), .by_group = TRUE) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      Group = comparison,
      Protein = substr(Protein_IDs, 1, 40),
      Gene = Gene_names,
      Classification = classification,
      `Light excess` = round(Light_fraction_excess, 3),
      `Serum fraction` = round(estimated_serum_fraction, 3),
      `CV` = round(serum_fraction_CV, 3),
      `Database evidence` = Database_evidence
    )

  ucev_stats <- c(
    n = nrow(ucev_sd),
    median_sd = stats::median(ucev_sd$UCEV_Light_fraction_sd, na.rm = TRUE),
    p90_sd = as.numeric(stats::quantile(ucev_sd$UCEV_Light_fraction_sd,
                                        0.9, na.rm = TRUE))
  )

  par_tbl <- tibble::tibble(
    Parameter = c(
      "min_valid_reps", "min_absolute_light_excess",
      "background_multiplier_high", "background_multiplier_medium",
      "max_serum_fraction_cv_high_confidence", "minimum_fraction_sd_floor",
      "min_light_dominance_ratio", "global_fraction_sd_floor (computed)"
    ),
    Value = c(
      params$min_valid_reps, params$min_absolute_light_excess,
      params$background_multiplier_high, params$background_multiplier_medium,
      params$max_serum_fraction_cv_high_confidence, params$minimum_fraction_sd_floor,
      params$min_light_dominance_ratio, round(global_fraction_sd_floor, 4)
    )
  )

  fig <- function(name) paste0("![", name, "](plots/", name, ".png)")

  report <- paste0(
    "# Heavy/Light Isotope-Based Serum Source Analysis Report\n\n",
    "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M"), "  \n",
    "Input file: `", input_file, "`  \n",
    "Tool: R package `SerumSourceAnalyzer` v",
    utils::packageVersion("SerumSourceAnalyzer"), "\n\n",
    "## 1. Method Overview\n\n",
    "Using UCEV as the control group, each target group is compared ",
    "protein-by-protein on the Light channel fraction ",
    "(Light / (Light + Heavy)). Proteins whose target-group Light ",
    "fraction significantly exceeds the UCEV background are classified ",
    "as serum-derived, and p0 (the UCEV heavy-labeling rate) is used to ",
    "source-correct and estimate the serum fraction. ",
    "Classification confidence is determined jointly by the Light-fraction ",
    "excess, replicate support count, CV, and CI separation.\n",
    "When both channels are completely absent in the control, ",
    "Acquired-only classification (Level 1/2) is applied instead.\n\n",
    "## 2. Analysis Parameters\n\n",
    md_table(par_tbl), "\n\n",
    "## 3. UCEV Background Variation\n\n",
    sprintf(
      "Proteins used for background estimation: %d; median robust SD of Light fraction = %.4f, P90 = %.4f, global SD floor = %.4f.\n\n",
      ucev_stats[["n"]], ucev_stats[["median_sd"]], ucev_stats[["p90_sd"]],
      global_fraction_sd_floor
    ),
    fig("ucev_background"), "\n\n",
    "## 4. Classification Summary\n\n",
    md_table(as.data.frame(cls_wide)), "\n\n",
    "### Primary serum-derived protein counts per group\n\n",
    md_table(as.data.frame(primary_counts)), "\n\n",
    fig("classification_summary"), "\n\n",
    "## 5. Light-Fraction Excess Distribution\n\n",
    fig("excess_distribution"), "\n\n",
    "### Excess vs Statistical Probability (volcano-style)\n\n",
    fig("excess_vs_probability"), "\n\n",
    "### Top protein excess heatmap\n\n",
    fig("excess_heatmap"), "\n\n",
    "## 6. Serum Fraction Estimates (primary categories)\n\n",
    fig("serum_fraction"), "\n\n",
    "## 7. Database Evidence Composition\n\n",
    md_table(as.data.frame(
      database_summary %>%
        dplyr::filter(classification %in% primary) %>%
        dplyr::count(comparison, Database_evidence, wt = protein_count,
                     name = "protein_count")
    )), "\n\n",
    fig("database_evidence"), "\n\n",
    "## 8. Top 10 Candidate Proteins per Group\n\n",
    md_table(as.data.frame(top_proteins)), "\n\n",
    "## 9. Output Files\n\n",
    "- `00_All_Comparisons_Merged_Results.csv` - per-protein results for all comparisons\n",
    "- `00_All_Comparisons_Primary_Serum_Derived_Merged.csv` - merged primary serum-derived proteins\n",
    "- `00_Classification_Summary_by_Group.csv` - classification summary\n",
    "- `<comparison>/` subdirectories - full result tables by category for each group\n",
    "- `Serum_Source_Analysis_All_Groups_Summary.xlsx` - key-table summary workbook\n",
    "- `plots/` - all figures as PNG\n"
  )

  report_path <- file.path(output_dir, "Analysis_Report.md")
  writeLines(enc2utf8(report), report_path, useBytes = TRUE)
  invisible(report_path)
}
