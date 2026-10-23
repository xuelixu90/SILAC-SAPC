# ============================================================
# Sample mapping for the NP_1 heavy/light isotope dataset
# ============================================================

#' Build the default sample mapping table
#'
#' The column names in the MaxQuant annotation workbook contain en dashes
#' ("–") rather than hyphens ("-"); this factory returns the mapping between
#' logical samples/groups and the exact Excel column names.
#'
#' @return A tibble with columns `sample_id`, `group`, `light_column`,
#'   `heavy_column`, `peptide_column`.
#' @export
#' @examples
#' si <- build_sample_info()
#' head(si$sample_id)
build_sample_info <- function() {
  tibble::tibble(
    sample_id = c(
      # UCEV controls (3 reps)
      "UCEV_1", "UCEV_2", "UCEV_3",
      # DGEV_PC_1-3 (3 reps)
      "DGEV\u2013PC_1-3_1", "DGEV\u2013PC_1-3_2", "DGEV\u2013PC_1-3_3",
      # DGEV_PC_10-12 (3 reps)
      "DGEV\u2013PC_10-12_1", "DGEV\u2013PC_10-12_2", "DGEV\u2013PC_10-12_3",
      # DGEV_PC_4-9 (3 reps)
      "DGEV\u2013PC_4-9_1", "DGEV\u2013PC_4-9_2", "DGEV\u2013PC_4-9_3",
      # UCEVPC (3 reps)
      "UCEV\u2013PC_1", "UCEV\u2013PC_2", "UCEV\u2013PC_3"
    ),
    group = c(
      rep("UCEV", 3),
      rep("DGEV_PC_1-3", 3),
      rep("DGEV_PC_10-12", 3),
      rep("DGEV_PC_4-9", 3),
      rep("UCEVPC", 3)
    ),
    light_column = c(
      "Intensity L UCEV_1",
      "Intensity L UCEV_2",
      "Intensity L UCEV_3",
      "Intensity L DGEV\u2013PC_1-3_1",
      "Intensity L DGEV\u2013PC_1-3_2",
      "Intensity L DGEV\u2013PC_1-3_3",
      "Intensity L DGEV\u2013PC_10-12_1",
      "Intensity L DGEV\u2013PC_10-12_2",
      "Intensity L DGEV\u2013PC_10-12_3",
      "Intensity L DGEV\u2013PC_4-9_1",
      "Intensity L DGEV\u2013PC_4-9_2",
      "Intensity L DGEV\u2013PC_4-9_3",
      "Intensity L UCEV\u2013PC_1",
      "Intensity L UCEV\u2013PC_2",
      "Intensity L UCEV\u2013PC_3"
    ),
    heavy_column = c(
      "Intensity H UCEV_1",
      "Intensity H UCEV_2",
      "Intensity H UCEV_3",
      "Intensity H DGEV\u2013PC_1-3_1",
      "Intensity H DGEV\u2013PC_1-3_2",
      "Intensity H DGEV\u2013PC_1-3_3",
      "Intensity H DGEV\u2013PC_10-12_1",
      "Intensity H DGEV\u2013PC_10-12_2",
      "Intensity H DGEV\u2013PC_10-12_3",
      "Intensity H DGEV\u2013PC_4-9_1",
      "Intensity H DGEV\u2013PC_4-9_2",
      "Intensity H DGEV\u2013PC_4-9_3",
      "Intensity H UCEV\u2013PC_1",
      "Intensity H UCEV\u2013PC_2",
      # NOTE: kept exactly as in the original analysis (no UCEV prefix)
      "Intensity H EV\u2013PC_3"
    ),
    peptide_column = c(
      "Unique peptides UCEV_1",
      "Unique peptides UCEV_2",
      "Unique peptides UCEV_3",
      "Unique peptides DGEV\u2013PC_1-3_1",
      "Unique peptides DGEV\u2013PC_1-3_2",
      "Unique peptides DGEV\u2013PC_1-3_3",
      "Unique peptides DGEV\u2013PC_10-12_1",
      "Unique peptides DGEV\u2013PC_10-12_2",
      "Unique peptides DGEV\u2013PC_10-12_3",
      "Unique peptides DGEV\u2013PC_4-9_1",
      "Unique peptides DGEV\u2013PC_4-9_2",
      "Unique peptides DGEV\u2013PC_4-9_3",
      "Unique peptides UCEV\u2013PC_1",
      "Unique peptides UCEV\u2013PC_2",
      "Unique peptides UCEV\u2013PC_3"
    )
  )
}

#' Validate a sample mapping against the raw data columns
#'
#' @param sample_info sample mapping tibble (see [build_sample_info()])
#' @param raw_data raw data as read from the Excel workbook
#' @return invisibly returns `sample_info`; errors when columns are missing
#' @keywords internal
validate_sample_info <- function(sample_info, raw_data) {
  required_cols <- c(
    sample_info$light_column,
    sample_info$heavy_column,
    sample_info$peptide_column
  )
  missing_cols <- setdiff(required_cols, names(raw_data))
  if (length(missing_cols) > 0) {
    stop(
      "The following columns are missing from the input file:\n",
      paste(missing_cols, collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(sample_info)
}
