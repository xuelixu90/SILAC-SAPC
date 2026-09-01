# ============================================================
# ggplot2-based visualizations of key analysis results
# ============================================================

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
#'
#' @param case_summary summary tibble with columns `comparison`,
#'   `classification`, `protein_count` (see [run_serum_source_analysis()])
#' @return ggplot object
#' @export
plot_classification_summary <- function(case_summary) {
  primary <- primary_serum_categories()
  ggplot2::ggplot(
    case_summary,
    ggplot2::aes(x = comparison, y = protein_count, fill = classification)
  ) +
    ggplot2::geom_col(width = 0.7, color = "grey30", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = .classification_palette) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Protein classification by group",
      subtitle = "Primary serum-derived categories in strong colors",
      x = NULL, y = "Number of proteins", fill = "Classification"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "right")
}

#' Distribution of Light-fraction excess per comparison
#'
#' Violin + boxplot of `Light_fraction_excess`; values above 0 indicate more
#' Light channel than the UCEV background baseline.
#'
#' @param results_df merged results tibble (all comparisons)
#' @return ggplot object
#' @export
plot_excess_distribution <- function(results_df) {
  primary <- primary_serum_categories()
  results_df <- results_df %>%
    dplyr::mutate(is_primary = classification %in% primary)
  ggplot2::ggplot(
    results_df,
    ggplot2::aes(x = comparison, y = Light_fraction_excess)
  ) +
    ggplot2::geom_violin(fill = "grey90", color = NA, scale = "width") +
    ggplot2::geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white") +
    ggplot2::geom_point(
      data = ~ dplyr::filter(.x, is_primary),
      ggplot2::aes(color = classification),
      position = ggplot2::position_jitter(width = 0.15),
      size = 1.4, alpha = 0.9
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::scale_color_manual(values = .classification_palette, na.translate = FALSE) +
    ggplot2::labs(
      title = "Light-fraction excess over UCEV background",
      x = NULL, y = "Light fraction excess (target - UCEV baseline)",
      color = "Primary classification"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
}

#' Estimated serum fractions for primary serum-derived proteins
#'
#' @param results_df merged results tibble
#' @param categories classification categories to show (default: primary)
#' @return ggplot object
#' @export
plot_serum_fraction <- function(results_df,
                                categories = primary_serum_categories()) {
  d <- results_df %>%
    dplyr::filter(
      classification %in% categories,
      is.finite(estimated_serum_fraction)
    )
  ggplot2::ggplot(
    d,
    ggplot2::aes(x = comparison, y = estimated_serum_fraction)
  ) +
    ggplot2::geom_hline(yintercept = c(0, 1), linetype = "dashed", color = "grey50") +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = classification),
      width = 0.4, outlier.shape = NA, alpha = 0.6
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = classification),
      position = ggplot2::position_jitter(width = 0.12),
      size = 1.5, alpha = 0.8
    ) +
    ggplot2::scale_fill_manual(values = .classification_palette) +
    ggplot2::scale_color_manual(values = .classification_palette) +
    ggplot2::coord_cartesian(ylim = c(0, 1.05)) +
    ggplot2::labs(
      title = "Source-corrected serum fraction estimates",
      subtitle = "Median per-protein serum fraction within primary categories",
      x = NULL, y = "Estimated serum fraction", fill = "Classification",
      color = "Classification"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
}

#' Volcano-style plot: excess vs statistical significance
#'
#' x = Light-fraction excess, y = probability that the excess is real
#' (`probability_excess_Light`, one-sided z-test vs UCEV background).
#'
#' @param results_df merged results tibble
#' @param comparison optional comparison name to subset to
#' @param label_top_n label the top N proteins by excess x probability
#' @return ggplot object
#' @export
plot_excess_vs_probability <- function(results_df,
                                       comparison = NULL,
                                       label_top_n = 10) {
  d <- results_df
  if (!is.null(comparison)) {
    d <- d %>% dplyr::filter(comparison == !!comparison)
  }
  d <- d %>%
    dplyr::filter(
      is.finite(Light_fraction_excess),
      is.finite(probability_excess_Light)
    ) %>%
    dplyr::mutate(
      gene_short = ifelse(
        is.na(Gene_names) | trimws(Gene_names) == "",
        Protein_Row_ID,
        trimws(stringr::str_split(Gene_names, ";", simplify = TRUE)[, 1])
      )
    )

  top_idx <- order(
    d$Light_fraction_excess * d$probability_excess_Light,
    decreasing = TRUE
  )[seq_len(min(label_top_n, nrow(d)))]

  ggplot2::ggplot(d, ggplot2::aes(x = Light_fraction_excess, y = probability_excess_Light)) +
    ggplot2::geom_point(
      ggplot2::aes(color = classification), size = 1.6, alpha = 0.85
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggrepel_if_available(d[top_idx, ], "gene_short") +
    ggplot2::scale_color_manual(values = .classification_palette) +
    ggplot2::labs(
      title = ifelse(is.null(comparison),
        "Light-fraction excess vs probability (all groups)",
        paste0("Light-fraction excess vs probability: ", comparison)
      ),
      x = "Light fraction excess", y = "P(excess is real)", color = "Classification"
    ) +
    ggplot2::theme_bw(base_size = 12)
}

# optional ggrepel support (falls back to geom_text when unavailable)
ggrepel_if_available <- function(label_df, label_col) {
  if (nrow(label_df) == 0) return(NULL)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    return(
      ggrepel::geom_text_repel(
        data = label_df,
        ggplot2::aes(x = Light_fraction_excess, y = probability_excess_Light,
                     label = .data[[label_col]]),
        size = 2.6, max.overlaps = 20, segment.color = "grey60"
      )
    )
  }
  ggplot2::geom_text(
    data = label_df,
    ggplot2::aes(x = Light_fraction_excess, y = probability_excess_Light,
                 label = .data[[label_col]]),
    size = 2.4, vjust = -0.7
  )
}

#' Database-evidence composition of primary serum-derived proteins
#'
#' @param primary_df primary serum-derived proteins (merged across groups)
#' @return ggplot object
#' @export
plot_database_evidence <- function(primary_df) {
  d <- primary_df %>%
    dplyr::count(comparison, Database_evidence, name = "protein_count")
  ggplot2::ggplot(
    d,
    ggplot2::aes(x = comparison, y = protein_count, fill = Database_evidence)
  ) +
    ggplot2::geom_col(width = 0.6, color = "grey30", linewidth = 0.2) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Database evidence for primary serum-derived proteins",
      x = NULL, y = "Number of proteins", fill = "Database evidence"
    ) +
    ggplot2::theme_bw(base_size = 12)
}

#' Distribution of UCEV background Light-fraction SD
#'
#' @param ucev_sd per-protein UCEV background table from
#'   [compute_ucev_background()]
#' @param global_floor the global SD floor applied downstream
#' @return ggplot object
#' @export
plot_ucev_background <- function(ucev_sd, global_floor = NULL) {
  p <- ggplot2::ggplot(ucev_sd, ggplot2::aes(x = UCEV_Light_fraction_sd)) +
    ggplot2::geom_histogram(bins = 40, fill = "#4E79A7", color = "white") +
    ggplot2::labs(
      title = "UCEV background variation of Light fraction",
      x = "Robust SD of Light fraction (per protein)",
      y = "Number of proteins"
    ) +
    ggplot2::theme_bw(base_size = 12)
  if (!is.null(global_floor)) {
    p <- p + ggplot2::geom_vline(
      xintercept = global_floor, linetype = "dashed", color = "#C0392B"
    ) +
      ggplot2::annotate(
        "text", x = global_floor, y = Inf, vjust = 2, hjust = -0.05,
        label = sprintf("Global SD floor = %.3f", global_floor),
        color = "#C0392B", size = 3.2
      )
  }
  p
}

#' Heatmap of Light-fraction excess for top proteins across comparisons
#'
#' Proteins are ranked by the maximum excess across comparisons and the
#' top `top_n` are drawn as a heatmap.
#'
#' @param results_df merged results tibble
#' @param top_n number of proteins to display (default 40)
#' @return ggplot object
#' @export
plot_excess_heatmap <- function(results_df, top_n = 40) {
  d <- results_df %>%
    dplyr::filter(is.finite(Light_fraction_excess)) %>%
    dplyr::mutate(
      gene_short = trimws(stringr::str_split(Gene_names, ";", simplify = TRUE)[, 1]),
      gene_short = dplyr::if_else(
        is.na(gene_short) | gene_short == "",
        Protein_Row_ID, gene_short
      )
    )

  top_proteins <- d %>%
    dplyr::group_by(Protein_Row_ID) %>%
    dplyr::summarise(max_excess = max(Light_fraction_excess), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(max_excess)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(Protein_Row_ID)

  d <- d %>%
    dplyr::filter(Protein_Row_ID %in% top_proteins) %>%
    dplyr::mutate(
      gene_short = substr(gene_short, 1, 15),
      label = make.unique(gene_short)
    )

  ggplot2::ggplot(
    d,
    ggplot2::aes(x = comparison, y = label, fill = Light_fraction_excess)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, na.value = "grey90"
    ) +
    ggplot2::labs(
      title = paste0("Top ", top_n, " proteins by Light-fraction excess"),
      x = NULL, y = NULL, fill = "Light fraction\nexcess"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))
}
