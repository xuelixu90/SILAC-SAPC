# ============================================================================
# downstream_analysis.R
# ----------------------------------------------------------------------------
# What is this script?
#   It is NOT part of the serumsource package. The package does the pure
#   computation; this script does everything around it:
#     1. draws the 7 key ggplot2 figures
#     2. writes all CSV / XLSX result tables
#     3. builds the Markdown analysis report
#     4. renders an HTML / Markdown execution notebook
#
# How to use it (two options):
#   A. Two-step: run the core computation yourself, then hand the result over
#        library(serumsource)
#        core <- run_core_analysis(input_file = "data.xlsx")
#        source("downstream/downstream_analysis.R")
#        run_downstream_analysis(core_result = core,
#                                input_file  = "data.xlsx",
#                                output_dir  = "out/")
#   B. One call does everything (core + downstream):
#        run_full_analysis(input_file = "data.xlsx", output_dir = "out/")
#
# Requirements: the packages below are checked on load; install any that
# are missing with install.packages().
# ============================================================================

# ---- Dependencies (not part of the serumsource package) -------------------
required_pkgs <- c("ggplot2", "dplyr", "tibble", "stringr", "tidyr", "writexl")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "The following packages are required for downstream analysis but are ",
    "not installed:\n  ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them with:\n  install.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(tidyr)
})

# ============================================================================
# 1. Visualisations
# ============================================================================

# Internal color mapping for classification labels
.classification_palette <- c(
  "Source-corrected, high confidence" = "#C0392B",
  "Source-corrected, medium confidence" = "#E67E22",
  "Source-corrected, single-peptide strong candidate" = "#B0A8B9",
  "Source-corrected, single-peptide candidate" = "#C9C3D0",
  "Acquired-only (Level 1)" = "#2E86AB",
  "Acquired-only (Level 2)" = "#7FB3D5",
  "Likely_residual" = "#95A5A6",
  "Ambiguous" = "#D5D8DC",
  "Control_absent_other" = "#F2F3F4",
  "Control_missing" = "#F2F3F4",
  "Target_group_missing" = "#F2F3F4"
)

#' Stacked bar chart of classification counts per comparison
plot_classification_summary <- function(case_summary) {
  primary <- serumsource::primary_serum_categories()
  ggplot(case_summary,
         aes(x = comparison, y = protein_count, fill = classification)) +
    geom_col(width = 0.7, color = "grey30", linewidth = 0.2) +
    scale_fill_manual(values = .classification_palette) +
    coord_flip() +
    labs(
      title = "Protein classification by group",
      subtitle = "Primary serum-derived categories in strong colors",
      x = NULL, y = "Number of proteins", fill = "Classification"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")
}

#' Distribution of Light-fraction excess per comparison
plot_excess_distribution <- function(results_df) {
  primary <- serumsource::primary_serum_categories()
  results_df <- results_df %>%
    mutate(is_primary = classification %in% primary)
  ggplot(results_df,
         aes(x = comparison, y = Light_fraction_excess)) +
    geom_violin(fill = "grey90", color = NA, scale = "width") +
    geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white") +
    geom_point(
      data = ~ filter(.x, is_primary),
      aes(color = classification),
      position = position_jitter(width = 0.15),
      size = 1.4, alpha = 0.9
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = .classification_palette, na.translate = FALSE) +
    labs(
      title = "Light-fraction excess over UCEV background",
      x = NULL, y = "Light fraction excess (target - UCEV baseline)",
      color = "Primary classification"
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

#' Estimated serum fractions for primary serum-derived proteins
plot_serum_fraction <- function(results_df,
                                categories = serumsource::primary_serum_categories()) {
  d <- results_df %>%
    filter(
      classification %in% categories,
      is.finite(estimated_serum_fraction)
    )
  ggplot(d, aes(x = comparison, y = estimated_serum_fraction)) +
    geom_hline(yintercept = c(0, 1), linetype = "dashed", color = "grey50") +
    geom_boxplot(aes(fill = classification), width = 0.4,
                 outlier.shape = NA, alpha = 0.6) +
    geom_point(aes(color = classification),
               position = position_jitter(width = 0.12),
               size = 1.5, alpha = 0.8) +
    scale_fill_manual(values = .classification_palette) +
    scale_color_manual(values = .classification_palette) +
    coord_cartesian(ylim = c(0, 1.05)) +
    labs(
      title = "Source-corrected serum fraction estimates",
      subtitle = "Median per-protein serum fraction within primary categories",
      x = NULL, y = "Estimated serum fraction", fill = "Classification",
      color = "Classification"
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

# optional ggrepel support (falls back to geom_text when unavailable)
ggrepel_if_available <- function(label_df, label_col) {
  if (nrow(label_df) == 0) return(NULL)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    return(
      ggrepel::geom_text_repel(
        data = label_df,
        aes(x = Light_fraction_excess, y = probability_excess_Light,
            label = .data[[label_col]]),
        size = 2.6, max.overlaps = 20, segment.color = "grey60"
      )
    )
  }
  geom_text(
    data = label_df,
    aes(x = Light_fraction_excess, y = probability_excess_Light,
        label = .data[[label_col]]),
    size = 2.4, vjust = -0.7
  )
}

#' Volcano-style plot: excess vs statistical significance
plot_excess_vs_probability <- function(results_df,
                                       comparison = NULL,
                                       label_top_n = 10) {
  d <- results_df
  if (!is.null(comparison)) {
    d <- d %>% filter(comparison == !!comparison)
  }
  d <- d %>%
    filter(
      is.finite(Light_fraction_excess),
      is.finite(probability_excess_Light)
    ) %>%
    mutate(
      gene_short = ifelse(
        is.na(Gene_names) | trimws(Gene_names) == "",
        Protein_Row_ID,
        trimws(str_split(Gene_names, ";", simplify = TRUE)[, 1])
      )
    )

  top_idx <- order(
    d$Light_fraction_excess * d$probability_excess_Light,
    decreasing = TRUE
  )[seq_len(min(label_top_n, nrow(d)))]

  ggplot(d, aes(x = Light_fraction_excess, y = probability_excess_Light)) +
    geom_point(aes(color = classification), size = 1.6, alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggrepel_if_available(d[top_idx, ], "gene_short") +
    scale_color_manual(values = .classification_palette) +
    labs(
      title = ifelse(is.null(comparison),
        "Light-fraction excess vs probability (all groups)",
        paste0("Light-fraction excess vs probability: ", comparison)
      ),
      x = "Light fraction excess", y = "P(excess is real)", color = "Classification"
    ) +
    theme_bw(base_size = 12)
}

#' Database-evidence composition of primary serum-derived proteins
plot_database_evidence <- function(primary_df) {
  d <- primary_df %>%
    count(comparison, Database_evidence, name = "protein_count")
  ggplot(d, aes(x = comparison, y = protein_count, fill = Database_evidence)) +
    geom_col(width = 0.6, color = "grey30", linewidth = 0.2) +
    coord_flip() +
    labs(
      title = "Database evidence for primary serum-derived proteins",
      x = NULL, y = "Number of proteins", fill = "Database evidence"
    ) +
    theme_bw(base_size = 12)
}

#' Distribution of UCEV background Light-fraction SD
plot_ucev_background <- function(ucev_sd, global_floor = NULL) {
  p <- ggplot(ucev_sd, aes(x = UCEV_Light_fraction_sd)) +
    geom_histogram(bins = 40, fill = "#4E79A7", color = "white") +
    labs(
      title = "UCEV background variation of Light fraction",
      x = "Robust SD of Light fraction (per protein)",
      y = "Number of proteins"
    ) +
    theme_bw(base_size = 12)
  if (!is.null(global_floor)) {
    p <- p + geom_vline(
      xintercept = global_floor, linetype = "dashed", color = "#C0392B"
    ) +
      annotate(
        "text", x = global_floor, y = Inf, vjust = 2, hjust = -0.05,
        label = sprintf("Global SD floor = %.3f", global_floor),
        color = "#C0392B", size = 3.2
      )
  }
  p
}

#' Heatmap of Light-fraction excess for top proteins across comparisons
plot_excess_heatmap <- function(results_df, top_n = 40) {
  d <- results_df %>%
    filter(is.finite(Light_fraction_excess)) %>%
    mutate(
      gene_short = trimws(str_split(Gene_names, ";", simplify = TRUE)[, 1]),
      gene_short = ifelse(
        is.na(gene_short) | gene_short == "",
        Protein_Row_ID, gene_short
      )
    )

  top_proteins <- d %>%
    group_by(Protein_Row_ID) %>%
    summarise(max_excess = max(Light_fraction_excess), .groups = "drop") %>%
    arrange(desc(max_excess)) %>%
    slice_head(n = top_n) %>%
    pull(Protein_Row_ID)

  d <- d %>%
    filter(Protein_Row_ID %in% top_proteins) %>%
    mutate(
      gene_short = substr(gene_short, 1, 15),
      label = make.unique(gene_short)
    )

  ggplot(d, aes(x = comparison, y = label, fill = Light_fraction_excess)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, na.value = "grey90"
    ) +
    labs(
      title = paste0("Top ", top_n, " proteins by Light-fraction excess"),
      x = NULL, y = NULL, fill = "Light fraction\nexcess"
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 7))
}

# ============================================================================
# 2. Output saving (CSV / XLSX / PNG)
# ============================================================================

save_analysis_outputs <- function(output_dir,
                                  sample_info,
                                  annotation,
                                  ucev_sd,
                                  all_merged,
                                  all_primary,
                                  case_summary,
                                  database_summary,
                                  target_groups,
                                  reference_group,
                                  primary) {
  write.csv(sample_info, file.path(output_dir, "01_sample_information.csv"),
    row.names = FALSE, na = ""
  )
  write.csv(annotation, file.path(output_dir, "02_protein_database_annotation.csv"),
    row.names = FALSE, na = ""
  )
  write.csv(ucev_sd, file.path(output_dir, "04_UCEV_background_variation.csv"),
    row.names = FALSE, na = ""
  )
  write.csv(all_merged, file.path(output_dir, "00_All_Comparisons_Merged_Results.csv"),
    row.names = FALSE, na = ""
  )
  write.csv(
    all_primary,
    file.path(output_dir, "00_All_Comparisons_Primary_Serum_Derived_Merged.csv"),
    row.names = FALSE, na = ""
  )
  write.csv(case_summary, file.path(output_dir, "00_Classification_Summary_by_Group.csv"),
    row.names = FALSE, na = ""
  )
  write.csv(
    database_summary,
    file.path(output_dir, "00_Database_Annotation_Summary_All_Groups.csv"),
    row.names = FALSE, na = ""
  )

  for (tgt_grp in target_groups) {
    comparison_name <- paste0(tgt_grp, "_vs_", reference_group)
    comparison_dir <- file.path(output_dir, comparison_name)
    dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)

    cmp_full <- all_merged %>% filter(comparison == comparison_name)

    write.csv(
      cmp_full,
      file.path(comparison_dir, paste0(comparison_name, "_00_All_Results.csv")),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% filter(classification %in% primary),
      file.path(
        comparison_dir,
        paste0(comparison_name, "_01_Merged_Primary_Serum_Derived_Proteins.csv")
      ),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% filter(classification == "Acquired-only (Level 1)"),
      file.path(comparison_dir, paste0(comparison_name, "_02_Acquired_only_Level1.csv")),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% filter(classification == "Acquired-only (Level 2)"),
      file.path(comparison_dir, paste0(comparison_name, "_03_Acquired_only_Level2.csv")),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% filter(classification == "Source-corrected, high confidence"),
      file.path(
        comparison_dir,
        paste0(comparison_name, "_04_Source_corrected_High_Confidence.csv")
      ),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% filter(classification == "Source-corrected, medium confidence"),
      file.path(
        comparison_dir,
        paste0(comparison_name, "_05_Source_corrected_Medium_Confidence.csv")
      ),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% filter(grepl("single-peptide", classification, ignore.case = TRUE)),
      file.path(
        comparison_dir,
        paste0(comparison_name, "_06_Source_corrected_Single_Peptide_Candidates.csv")
      ),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% filter(classification == "Likely_residual"),
      file.path(comparison_dir, paste0(comparison_name, "_07_Likely_Residual.csv")),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% filter(classification == "Ambiguous"),
      file.path(comparison_dir, paste0(comparison_name, "_08_Ambiguous.csv")),
      row.names = FALSE, na = ""
    )

    primary_df <- cmp_full %>%
      filter(classification %in% primary) %>%
      select(
        Protein_Row_ID, Original_Excel_Row,
        Protein_IDs, Representative_ID,
        Protein_names, Gene_names,
        In_house_db, Plasma_Altas_db,
        In_house_serum_annotation,
        Human_plasma_PeptideAtlas_annotation,
        Database_evidence,
        assigned_case, classification,
        estimated_serum_fraction, serum_fraction_CV,
        Light_fraction_excess,
        everything()
      )
    write.csv(
      primary_df,
      file.path(
        output_dir,
        paste0(comparison_name, "_Primary_Serum_Proteins_with_Database_Annotation.csv")
      ),
      row.names = FALSE, na = ""
    )
  }

  excel_output <- list(
    Sample_information = sample_info,
    Protein_annotation = annotation,
    UCEV_background = ucev_sd,
    Classification_summary = case_summary,
    Database_summary = database_summary,
    All_primary_serum_proteins = all_primary
  )
  writexl::write_xlsx(
    excel_output,
    path = file.path(output_dir, "Serum_Source_Analysis_All_Groups_Summary.xlsx")
  )

  invisible(output_dir)
}

save_analysis_plots <- function(plots,
                                plots_dir,
                                width = 9,
                                height = 6,
                                dpi = 200) {
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(plots)) {
    p <- plots[[nm]]
    if (is.null(p) || !inherits(p, "ggplot")) next
    h <- height
    if (nm %in% c("classification_summary", "database_evidence")) {
      h <- max(height, 7)
    }
    if (nm == "excess_heatmap") {
      h <- max(height, 10)
    }
    ggsave(
      filename = file.path(plots_dir, paste0(nm, ".png")),
      plot = p, width = width, height = h, dpi = dpi, bg = "white"
    )
  }
  invisible(plots_dir)
}

# ============================================================================
# 3. Markdown analysis report
# ============================================================================

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

write_analysis_report <- function(output_dir,
                                  input_file,
                                  params,
                                  case_summary,
                                  all_primary,
                                  database_summary,
                                  all_results,
                                  ucev_sd,
                                  global_fraction_sd_floor) {
  primary <- serumsource::primary_serum_categories()

  cls_wide <- case_summary %>%
    group_by(comparison, classification) %>%
    summarise(n = sum(protein_count), .groups = "drop") %>%
    pivot_wider(names_from = classification, values_from = n, values_fill = 0)

  primary_counts <- all_primary %>%
    count(comparison, name = "Primary_serum_proteins")

  top_proteins <- all_primary %>%
    group_by(comparison) %>%
    arrange(desc(Light_fraction_excess), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    transmute(
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

  par_tbl <- tibble(
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
    "Tool: R package `serumsource` v",
    utils::packageVersion("serumsource"), "\n\n",
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
        filter(classification %in% primary) %>%
        count(comparison, Database_evidence, wt = protein_count,
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

# ============================================================================
# 4. Notebook rendering
# ============================================================================

# Locate the analysis notebook template. It lives next to this script
# (downstream/analysis_notebook.Rmd).
find_notebook_rmd <- function() {
  here <- tryCatch(
    normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
    error = function(e) NA_character_
  )
  candidates <- c(
    if (!is.na(here)) file.path(dirname(here), "analysis_notebook.Rmd"),
    file.path(getwd(), "downstream", "analysis_notebook.Rmd"),
    file.path(dirname(getwd()), "downstream", "analysis_notebook.Rmd"),
    file.path(getwd(), "analysis_notebook.Rmd")
  )
  found <- candidates[file.exists(candidates) & nchar(candidates) > 0]
  if (length(found) == 0) {
    stop("Cannot find analysis_notebook.Rmd template. Looked in:\n  ",
         paste(candidates, collapse = "\n  "))
  }
  normalizePath(found[1], winslash = "/", mustWork = TRUE)
}

render_analysis_notebook <- function(input_file,
                                     output_dir,
                                     sample_info = serumsource::build_sample_info(),
                                     params = serumsource::serum_source_params()) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    warning("Package 'rmarkdown' is not installed. Skipping notebook render. ",
            "Install with:\n  install.packages(c('rmarkdown','knitr'))")
    return(NULL)
  }
  if (!requireNamespace("knitr", quietly = TRUE)) {
    warning("Package 'knitr' is not installed. Skipping notebook render.")
    return(NULL)
  }

  probe_pandoc <- function() {
    if (rmarkdown::pandoc_available()) return(TRUE)
    candidates <- c(
      "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/pandoc.exe",
      "C:/Program Files/RStudio/resources/app/bin/quarto/bin/pandoc.exe",
      "C:/Program Files (x86)/RStudio/resources/app/bin/quarto/bin/tools/pandoc.exe",
      "C:/Program Files/Pandoc/pandoc.exe",
      file.path(Sys.getenv("LOCALAPPDATA"), "Pandoc", "pandoc.exe", fsep = "/"),
      file.path(Sys.getenv("ProgramFiles"),     "Pandoc", "pandoc.exe", fsep = "/")
    )
    found <- candidates[file.exists(candidates)]
    if (length(found) == 0) return(FALSE)
    pandoc_dir <- dirname(normalizePath(found[1], winslash = "/", mustWork = TRUE))
    Sys.setenv(RSTUDIO_PANDOC = pandoc_dir)
    tryCatch({
      unlink(Sys.getenv("RMARKDOWN_PANDOC"), recursive = TRUE)
    }, error = function(e) NULL)
    rmarkdown::pandoc_available()
  }
  have_pandoc <- probe_pandoc()
  if (!have_pandoc) {
    warning("pandoc >= 1.12.3 is required but not found. HTML notebook ",
            "cannot be rendered. Install RStudio or Pandoc (https://pandoc.org).")
    return(NULL)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  rmd_src  <- find_notebook_rmd()
  rmd_copy <- file.path(output_dir, "Analysis_Notebook.Rmd")
  file.copy(rmd_src, rmd_copy, overwrite = TRUE)

  lib_dir <- if (requireNamespace("serumsource", quietly = TRUE)) {
    dirname(system.file(package = "serumsource"))
  } else {
    .libPaths()[1]
  }

  params_list <- list(
    INPUT_FILE  = normalizePath(input_file, winslash = "/", mustWork = TRUE),
    OUTPUT_DIR  = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
    PARAMS      = params,
    SAMPLE_INFO = sample_info,
    LIB_DIR     = lib_dir
  )

  message("Rendering HTML notebook ...")
  nb_html <- tryCatch({
    rmarkdown::render(
      input          = rmd_copy,
      output_format  = "html_notebook",
      output_file    = "Analysis_Notebook.nb.html",
      output_dir     = output_dir,
      params         = params_list,
      envir          = new.env(parent = globalenv()),
      quiet          = FALSE,
      knit_root_dir   = output_dir
    )
  }, error = function(e) {
    warning("html_notebook render failed: ", conditionMessage(e),
            "\nFalling back to html_document.")
    NULL
  })

  message("Rendering html_document (fallback + print-friendly) ...")
  reg_html <- tryCatch({
    rmarkdown::render(
      input          = rmd_copy,
      output_format  = "html_document",
      output_file    = "Analysis_Notebook.html",
      output_dir     = output_dir,
      params         = params_list,
      envir          = new.env(parent = globalenv()),
      quiet          = FALSE,
      knit_root_dir   = output_dir
    )
  }, error = function(e) {
    warning("html_document render failed: ", conditionMessage(e))
    NULL
  })

  message("Rendering Markdown notebook ...")
  md_file <- tryCatch({
    rmarkdown::render(
      input          = rmd_copy,
      output_format  = "md_document",
      output_file    = "Analysis_Notebook.md",
      output_dir     = output_dir,
      params         = params_list,
      envir          = new.env(parent = globalenv()),
      quiet          = FALSE,
      knit_root_dir   = output_dir
    )
  }, error = function(e) {
    warning("md_document render failed: ", conditionMessage(e))
    NULL
  })

  html_out <- if (!is.null(nb_html) && file.exists(nb_html)) {
    nb_html
  } else if (!is.null(reg_html) && file.exists(reg_html)) {
    reg_html
  } else {
    NA_character_
  }

  out <- list(
    rmd      = rmd_copy,
    html     = if (is.na(html_out)) NULL else html_out,
    nb_html  = if (!is.null(nb_html) && file.exists(nb_html)) nb_html else NULL,
    html_doc = if (!is.null(reg_html) && file.exists(reg_html)) reg_html else NULL,
    markdown = if (!is.null(md_file)  && file.exists(md_file))  md_file  else NULL
  )
  class(out) <- "notebook_output"
  message("Notebook render finished.")
  message("  HTML (open in browser): ", if (!is.null(out$html)) out$html else "<failed>")
  message("  Markdown:              ", if (!is.null(out$markdown)) out$markdown else "<failed>")
  invisible(out)
}

# ============================================================================
# 5. Orchestration: run the full downstream pipeline
# ============================================================================

#' Run the downstream analysis on top of a core result
#'
#' Generates all plots, CSV/XLSX outputs, the Markdown report and (optionally)
#' the execution notebook from the result of `serumsource::run_core_analysis()`.
#'
#' @param core_result list returned by `serumsource::run_core_analysis()`
#' @param input_file path to the input workbook (for report provenance)
#' @param output_dir destination directory (created if needed)
#' @param sample_info sample mapping table (default `build_sample_info()`)
#' @param params a `serumsource::serum_source_params()` object
#' @param target_groups target group names, or `NULL` (default) to use all
#'   groups in `sample_info` except `reference_group`
#' @param reference_group reference group name (default "UCEV")
#' @param make_plots logical, generate and save plots (default TRUE)
#' @param save_outputs logical, write CSV/XLSX outputs (default TRUE)
#' @param write_report logical, write the Markdown report (default TRUE)
#' @param render_notebook logical, render the HTML/MD notebook (default TRUE)
#' @return invisibly, a list with `plots` and `notebook_output`
run_downstream_analysis <- function(core_result,
                                    input_file,
                                    output_dir,
                                    sample_info = serumsource::build_sample_info(),
                                    params = serumsource::serum_source_params(),
                                    target_groups = NULL,
                                    reference_group = "UCEV",
                                    make_plots = TRUE,
                                    save_outputs = TRUE,
                                    write_report = TRUE,
                                    render_notebook = TRUE) {
  primary <- serumsource::primary_serum_categories()

  # Auto-infer target groups from the sample mapping when not supplied
  if (is.null(target_groups)) {
    target_groups <- setdiff(unique(sample_info$group), reference_group)
    if (length(target_groups) == 0) {
      stop("No target groups found in sample_info after excluding '",
           reference_group, "'.", call. = FALSE)
    }
  }

  plots <- list()
  if (make_plots) {
    plots <- list(
      classification_summary = plot_classification_summary(core_result$case_summary),
      excess_distribution = plot_excess_distribution(core_result$results),
      serum_fraction = plot_serum_fraction(core_result$primary_serum),
      database_evidence = plot_database_evidence(core_result$primary_serum),
      ucev_background = plot_ucev_background(
        core_result$ucev_sd, core_result$global_fraction_sd_floor
      ),
      excess_vs_probability = plot_excess_vs_probability(core_result$results),
      excess_heatmap = plot_excess_heatmap(core_result$results)
    )
  }

  if (save_outputs) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    save_analysis_outputs(
      output_dir = output_dir,
      sample_info = sample_info,
      annotation = core_result$annotation,
      ucev_sd = core_result$ucev_sd,
      all_merged = core_result$all_merged,
      all_primary = core_result$primary_serum,
      case_summary = core_result$case_summary,
      database_summary = core_result$database_summary,
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
      case_summary = core_result$case_summary,
      all_primary = core_result$primary_serum,
      database_summary = core_result$database_summary,
      all_results = core_result$results,
      ucev_sd = core_result$ucev_sd,
      global_fraction_sd_floor = core_result$global_fraction_sd_floor
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
    plots = plots,
    notebook_output = notebook_output
  ))
}

#' One-call full analysis (core computation + downstream)
#'
#' Convenience wrapper that runs `serumsource::run_core_analysis()` followed
#' by `run_downstream_analysis()`.
#'
#' @param input_file path to the annotated Excel workbook
#' @param output_dir destination directory
#' @param sheet_name sheet index/name (default 1)
#' @param sample_info sample mapping (default `build_sample_info()`)
#' @param reference_group control group name (default "UCEV")
#' @param target_groups character vector of target group names, or `NULL`
#'   (default) to use all groups in `sample_info` except `reference_group`.
#' @param params parameters from `serumsource::serum_source_params()`
#' @param annotation_cols named character vector mapping logical annotation
#'   keys to workbook column names (passed to `run_core_analysis()`).
#' @param na character vector of values to treat as missing (passed to
#'   `run_core_analysis()`).
#' @param make_plots,save_outputs,write_report,render_notebook downstream
#'   toggles passed to [run_downstream_analysis()]
#' @return invisibly, a list with `core` and `downstream` elements
run_full_analysis <- function(input_file,
                              output_dir,
                              sheet_name = 1,
                              sample_info = serumsource::build_sample_info(),
                              reference_group = "UCEV",
                              target_groups = NULL,
                              params = serumsource::serum_source_params(),
                              annotation_cols = c(
                                protein_ids = "Protein IDs",
                                representative_id = "Representive ID",
                                protein_names = "Protein names",
                                gene_names = "Gene names",
                                in_house_db = "In_house_db",
                                plasma_atlas_db = "Plasma_Altas_db"
                              ),
                              na = c("", " ", "NA", "NaN", "#N/A", "#DIV/0!"),
                              make_plots = TRUE,
                              save_outputs = TRUE,
                              write_report = TRUE,
                              render_notebook = TRUE) {
  core <- serumsource::run_core_analysis(
    input_file = input_file,
    sheet_name = sheet_name,
    sample_info = sample_info,
    reference_group = reference_group,
    target_groups = target_groups,
    params = params,
    annotation_cols = annotation_cols,
    na = na
  )

  # Determine target groups for downstream output splitting
  if (is.null(target_groups)) {
    target_groups <- setdiff(unique(sample_info$group), reference_group)
  }

  downstream <- run_downstream_analysis(
    core_result = core,
    input_file = input_file,
    output_dir = output_dir,
    sample_info = sample_info,
    params = params,
    target_groups = target_groups,
    reference_group = reference_group,
    make_plots = make_plots,
    save_outputs = save_outputs,
    write_report = write_report,
    render_notebook = render_notebook
  )

  invisible(list(core = core, downstream = downstream))
}
