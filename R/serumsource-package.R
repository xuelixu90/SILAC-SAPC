#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom magrittr %>%
#' @importFrom stats median sd mad pnorm
## usethis namespace: end
NULL

# Suppress R CMD check NOTEs for dplyr non-standard evaluation
utils::globalVariables(c(
  "Protein_Row_ID", "Original_Excel_Row", "Protein_IDs_original",
  "Representative_ID_original", "Protein names", "Gene names",
  "In_house_db", "Plasma_Altas_db", "Protein IDs", "Representative ID",
  "Light_raw", "Heavy_raw", "Paired_LH_Valid_Rep", "Peptides_raw",
  "Target_L_Dominant_Rep", "Light", "Heavy", "Total_Light_Heavy",
  "group", "n_valid_reps", "UCEV_Light_fraction_sd", "Light_fraction",
  "classification", "comparison", "assigned_case", "Database_evidence",
  "Peptide_Evidence_Level", "Target_High_Pep_Paired_Reps",
  "Target_Single_Pep_Paired_Reps",
  "protein_count", ".data"
))

