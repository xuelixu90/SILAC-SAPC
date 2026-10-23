# ============================================================
# Data input, annotation preparation and long-format assembly
# ============================================================

#' Read the annotated MaxQuant Excel workbook
#'
#' @param input_file path to the .xlsx workbook
#' @param sheet_name sheet index or name (default 1)
#' @param na character vector of values to treat as missing
#'   (default MaxQuant-style NA strings)
#' @return tibble with raw data
#' @export
read_input_data <- function(input_file,
                            sheet_name = 1,
                            na = c("", " ", "NA", "NaN", "#N/A", "#DIV/0!")) {
  readxl::read_excel(
    path = input_file,
    sheet = sheet_name,
    na = na
  )
}

#' Build the protein annotation table
#'
#' Adds normalized database-annotation columns and a combined
#' `Database_evidence` label.
#'
#' @param data converted data frame (must contain `Protein_Row_ID`,
#'   `Original_Excel_Row`, `Protein_IDs_original`,
#'   `Representative_ID_original` and the annotation columns)
#' @param annotation_cols named character vector mapping logical annotation
#'   keys to the actual column names in the data. The required keys are
#'   `protein_names`, `gene_names`, `in_house_db`, `plasma_atlas_db`.
#'   Defaults match the NP_1 MaxQuant annotation workbook.
#' @return tibble annotation table
#' @export
prepare_annotation <- function(data, annotation_cols = c(
                                 protein_names = "Protein names",
                                 gene_names = "Gene names",
                                 in_house_db = "In_house_db",
                                 plasma_atlas_db = "Plasma_Altas_db"
                               )) {
  ann <- as.list(annotation_cols)
  annotation <- data %>%
    dplyr::transmute(
      Protein_Row_ID,
      Original_Excel_Row,
      Protein_IDs = Protein_IDs_original,
      Representative_ID = Representative_ID_original,
      Protein_names = as.character(.data[[ann$protein_names]]),
      Gene_names = as.character(.data[[ann$gene_names]]),
      In_house_db = as.character(.data[[ann$in_house_db]]),
      Plasma_Altas_db = as.character(.data[[ann$plasma_atlas_db]])
    ) %>%
    dplyr::mutate(
      In_house_serum_annotation = dplyr::case_when(
        is.na(In_house_db) | trimws(In_house_db) == "" ~ "Not annotated",
        TRUE ~ In_house_db
      ),
      Human_plasma_PeptideAtlas_annotation = dplyr::case_when(
        is.na(Plasma_Altas_db) | trimws(Plasma_Altas_db) == "" ~ "Not annotated",
        TRUE ~ Plasma_Altas_db
      ),
      Database_evidence = dplyr::case_when(
        !is.na(In_house_db) & trimws(In_house_db) != "" &
          !is.na(Plasma_Altas_db) & trimws(Plasma_Altas_db) != "" ~
          "In-house serum database + Human Plasma PeptideAtlas",
        !is.na(In_house_db) & trimws(In_house_db) != "" ~
          "In-house serum database only",
        !is.na(Plasma_Altas_db) & trimws(Plasma_Altas_db) != "" ~
          "Human Plasma PeptideAtlas only",
        TRUE ~ "Not annotated in either database"
      )
    )
  annotation
}

#' Prepare long-format Heavy/Light data and annotation
#'
#' Performs numeric conversion, builds unique protein row IDs, constructs the
#' per-replicate long table with replicate validity / dominance flags and the
#' Light fraction, and attaches protein annotation.
#'
#' @param raw_data raw tibble from [read_input_data()]
#' @param sample_info sample mapping from [build_sample_info()]
#' @param min_light_dominance_ratio Light dominance ratio (default 0.70)
#' @param annotation_cols named character vector mapping logical keys to the
#'   actual column names in the workbook. Required keys: `protein_ids`,
#'   `representative_id`, `protein_names`, `gene_names`, `in_house_db`,
#'   `plasma_atlas_db`. Defaults match the NP_1 MaxQuant annotation workbook.
#' @return list with elements `data_converted`, `annotation`, `lh_long`
#' @export
prepare_long_data <- function(raw_data,
                              sample_info = build_sample_info(),
                              min_light_dominance_ratio = 0.70,
                              annotation_cols = c(
                                protein_ids = "Protein IDs",
                                representative_id = "Representive ID",
                                protein_names = "Protein names",
                                gene_names = "Gene names",
                                in_house_db = "In_house_db",
                                plasma_atlas_db = "Plasma_Altas_db"
                              )) {
  eps <- 1e-12
  ann <- as.list(annotation_cols)

  required_anno <- unname(unlist(ann))
  missing_anno <- setdiff(required_anno, names(raw_data))
  if (length(missing_anno) > 0) {
    stop(
      "Missing annotation columns: ",
      paste(missing_anno, collapse = ", "),
      call. = FALSE
    )
  }

  data <- raw_data %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c(
          sample_info$light_column,
          sample_info$heavy_column,
          sample_info$peptide_column
        )),
        ~ suppressWarnings(as.numeric(.x))
      ),
      Original_Excel_Row = dplyr::row_number(),
      Protein_IDs_original = as.character(.data[[ann$protein_ids]]),
      Representative_ID_original = as.character(.data[[ann$representative_id]]),
      Protein_Row_ID = make.unique(
        ifelse(
          is.na(.data[[ann$protein_ids]]) |
            trimws(as.character(.data[[ann$protein_ids]])) == "",
          paste0("Protein_row_", dplyr::row_number()),
          trimws(as.character(.data[[ann$protein_ids]]))
        )
      )
    )

  annotation <- prepare_annotation(data, annotation_cols = annotation_cols)

  lh_long <- purrr::map_dfr(
    seq_len(nrow(sample_info)),
    function(i) {
      tibble::tibble(
        Protein_Row_ID = data$Protein_Row_ID,
        sample_id = sample_info$sample_id[i],
        group = sample_info$group[i],
        Light_raw = data[[sample_info$light_column[i]]],
        Heavy_raw = data[[sample_info$heavy_column[i]]],
        Peptides_raw = data[[sample_info$peptide_column[i]]]
      )
    }
  ) %>%
    dplyr::mutate(
      Paired_LH_Valid_Rep =
        is.finite(Light_raw) & Light_raw > 0 &
          is.finite(Heavy_raw) & Heavy_raw > 0,

      Target_L_Dominant_Rep =
        is.finite(Light_raw) & Light_raw > 0 &
          (is.na(Heavy_raw) | Heavy_raw == 0 |
             (is.finite(Heavy_raw) &
                Light_raw / (Light_raw + Heavy_raw) >= min_light_dominance_ratio)),

      Both_channels_zero_or_na =
        (is.na(Light_raw) | Light_raw == 0) &
          (is.na(Heavy_raw) | Heavy_raw == 0),

      High_Pep_Matched_Rep =
        Paired_LH_Valid_Rep & is.finite(Peptides_raw) & Peptides_raw >= 2,

      Single_Pep_Matched_Rep =
        Paired_LH_Valid_Rep & is.finite(Peptides_raw) & Peptides_raw == 1,

      High_Pep_L_Dominant_Rep =
        Target_L_Dominant_Rep & is.finite(Peptides_raw) & Peptides_raw >= 2,

      Single_Pep_L_Dominant_Rep =
        Target_L_Dominant_Rep & is.finite(Peptides_raw) & Peptides_raw == 1,

      Light = ifelse(Paired_LH_Valid_Rep | Target_L_Dominant_Rep, Light_raw, NA_real_),
      Heavy = ifelse(
        Paired_LH_Valid_Rep | Target_L_Dominant_Rep,
        ifelse(is.na(Heavy_raw), 0, Heavy_raw),
        NA_real_
      ),
      Total_Light_Heavy = Light + Heavy,
      Light_fraction = safe_fraction(Light, Total_Light_Heavy)
    ) %>%
    dplyr::left_join(annotation, by = "Protein_Row_ID")

  list(
    data_converted = data,
    annotation = annotation,
    lh_long = lh_long
  )
}

#' Estimate UCEV background variation of the Light fraction
#'
#' @param lh_long long-format table from [prepare_long_data()]
#' @param reference_group name of the control group (default "UCEV")
#' @param min_valid_reps minimum valid replicates (default 2)
#' @param minimum_fraction_sd_floor minimum SD floor (default 0.02)
#' @return list with `sd_table` (per-protein background SD) and
#'   `global_fraction_sd_floor` (median-based floor used downstream)
#' @export
compute_ucev_background <- function(lh_long,
                                    reference_group = "UCEV",
                                    min_valid_reps = 2,
                                    minimum_fraction_sd_floor = 0.02) {
  ucev_protein_sd <- lh_long %>%
    dplyr::filter(group == reference_group) %>%
    dplyr::group_by(Protein_Row_ID) %>%
    dplyr::summarise(
      n_valid_reps = sum(Paired_LH_Valid_Rep, na.rm = TRUE),
      UCEV_Light_fraction_sd = robust_sd(Light_fraction),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      n_valid_reps >= min_valid_reps,
      is.finite(UCEV_Light_fraction_sd),
      UCEV_Light_fraction_sd > 0
    )

  global_fraction_sd_floor <- if (nrow(ucev_protein_sd) == 0) {
    minimum_fraction_sd_floor
  } else {
    max(
      minimum_fraction_sd_floor,
      stats::median(ucev_protein_sd$UCEV_Light_fraction_sd, na.rm = TRUE)
    )
  }

  list(
    sd_table = ucev_protein_sd,
    global_fraction_sd_floor = global_fraction_sd_floor
  )
}
