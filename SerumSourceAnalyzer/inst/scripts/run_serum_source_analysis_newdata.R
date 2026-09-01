#!/usr/bin/env Rscript
# Reproduce the full analysis of the original script (using the SerumSourceAnalyzer package)
# Input:  F:/WYF/NP_1/MaxQuant_annotated_with_total.xlsx
# Output: F:/WYF/NP_1/Serum_source_analysis_all_groups_1

library(SerumSourceAnalyzer)

result <- run_serum_source_analysis(
  input_file = "F:/WYF/NP_1/MaxQuant_annotated_with_total.xlsx",
  sheet_name = 1,
  output_dir = "F:/WYF/NP_1/Serum_source_analysis_all_groups_1",
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
  write_report = TRUE
)

# Classification summary preview
print(result$case_summary, n = 100)
