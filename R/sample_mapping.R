# ============================================================
# Sample mapping for the NP_1 heavy/light isotope dataset
# ============================================================

#' Build the sample mapping table
#'
#' Maps logical samples/groups to the exact Excel column names. The column
#' names in the default MaxQuant annotation workbook contain en dashes
#' ("–") rather than hyphens ("-").
#'
#' Called without arguments this returns the NP_1 default mapping. For other
#' datasets, pass `sample_names` together with `groups`; the workbook column
#' names are then built as `<prefix> <sample_name>` (MaxQuant-style prefixes
#' by default). If your workbook uses completely different column names,
#' build the mapping tibble by hand with the columns `sample_id`, `group`,
#' `light_column`, `heavy_column`, `peptide_column` and pass it through the
#' `sample_info` argument of [run_core_analysis()].
#'
#' @param sample_names character vector of sample names exactly as they
#'   appear after the prefix in the workbook column names; `NULL` (default)
#'   uses the NP_1 samples
#' @param groups character vector, same length as `sample_names`, giving the
#'   group of each sample; required when `sample_names` is supplied
#' @param light_prefix,heavy_prefix,peptide_prefix prefixes used to build the
#'   workbook column names; defaults are the standard MaxQuant prefixes
#' @return A tibble with columns `sample_id`, `group`, `light_column`,
#'   `heavy_column`, `peptide_column`.
#' @export
#' @examples
#' si <- build_sample_info()
#' head(si$sample_id)
#' \dontrun{
#' # A different dataset with column names such as "Intensity H WT_rep1":
#' si2 <- build_sample_info(
#'   sample_names = c("WT_rep1", "WT_rep2", "WT_rep3",
#'                    "KO_rep1", "KO_rep2", "KO_rep3"),
#'   groups = rep(c("WT", "KO"), each = 3)
#' )
#' }
build_sample_info <- function(sample_names = NULL,
                              groups = NULL,
                              light_prefix = "Intensity L",
                              heavy_prefix = "Intensity H",
                              peptide_prefix = "Unique peptides") {
  default_ids <- c(
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
  )
  default_groups <- c(
    rep("UCEV", 3),
    rep("DGEV_PC_1-3", 3),
    rep("DGEV_PC_10-12", 3),
    rep("DGEV_PC_4-9", 3),
    rep("UCEVPC", 3)
  )

  if (is.null(sample_names) && is.null(groups)) {
    sample_names <- default_ids
    groups <- default_groups
  } else if (is.null(sample_names) || is.null(groups)) {
    stop(
      "'sample_names' and 'groups' must be provided together ",
      "(or omit both for the NP_1 defaults).",
      call. = FALSE
    )
  }
  if (length(sample_names) != length(groups)) {
    stop("'sample_names' and 'groups' must have the same length.", call. = FALSE)
  }

  tibble::tibble(
    sample_id = as.character(sample_names),
    group = as.character(groups),
    light_column = paste(light_prefix, sample_names),
    heavy_column = paste(heavy_prefix, sample_names),
    peptide_column = paste(peptide_prefix, sample_names)
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
