# ============================================================
# Save analysis outputs (CSV / XLSX / PNG), mirroring the
# original script's output structure
# ============================================================

#' Save all CSV and XLSX outputs of the analysis
#'
#' @param output_dir output directory
#' @param sample_info,annotation,ucev_sd,all_merged,all_primary,case_summary,
#'   database_summary analysis tables
#' @param target_groups target group names
#' @param reference_group reference group name
#' @param primary primary classification categories
#' @return invisibly, the output directory
#' @keywords internal
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

    cmp_full <- all_merged %>% dplyr::filter(comparison == comparison_name)

    write.csv(
      cmp_full,
      file.path(comparison_dir, paste0(comparison_name, "_00_All_Results.csv")),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% dplyr::filter(classification %in% primary),
      file.path(
        comparison_dir,
        paste0(comparison_name, "_01_Merged_Primary_Serum_Derived_Proteins.csv")
      ),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% dplyr::filter(classification == "Acquired-only (Level 1)"),
      file.path(comparison_dir, paste0(comparison_name, "_02_Acquired_only_Level1.csv")),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% dplyr::filter(classification == "Acquired-only (Level 2)"),
      file.path(comparison_dir, paste0(comparison_name, "_03_Acquired_only_Level2.csv")),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% dplyr::filter(classification == "Source-corrected, high confidence"),
      file.path(
        comparison_dir,
        paste0(comparison_name, "_04_Source_corrected_High_Confidence.csv")
      ),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% dplyr::filter(classification == "Source-corrected, medium confidence"),
      file.path(
        comparison_dir,
        paste0(comparison_name, "_05_Source_corrected_Medium_Confidence.csv")
      ),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% dplyr::filter(grepl("single-peptide", classification, ignore.case = TRUE)),
      file.path(
        comparison_dir,
        paste0(comparison_name, "_06_Source_corrected_Single_Peptide_Candidates.csv")
      ),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% dplyr::filter(classification == "Likely_residual"),
      file.path(comparison_dir, paste0(comparison_name, "_07_Likely_Residual.csv")),
      row.names = FALSE, na = ""
    )
    write.csv(
      cmp_full %>% dplyr::filter(classification == "Ambiguous"),
      file.path(comparison_dir, paste0(comparison_name, "_08_Ambiguous.csv")),
      row.names = FALSE, na = ""
    )

    # per-group primary proteins with database annotation first
    primary_df <- cmp_full %>%
      dplyr::filter(classification %in% primary) %>%
      dplyr::select(
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
        dplyr::everything()
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

#' Save analysis plots as PNG files
#'
#' @param plots named list of ggplot objects
#' @param plots_dir directory to write PNGs into (created)
#' @param width,height,dpi passed to ggplot2::ggsave
#' @return invisibly, the plots directory
#' @export
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
    ggplot2::ggsave(
      filename = file.path(plots_dir, paste0(nm, ".png")),
      plot = p, width = width, height = h, dpi = dpi, bg = "white"
    )
  }
  invisible(plots_dir)
}
