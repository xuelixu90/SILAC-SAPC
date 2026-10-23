#!/usr/bin/env Rscript
# ============================================================================
# Reproduce the full NP_1 serum-source analysis with the serumsource package
# ----------------------------------------------------------------------------
# What this script does (2 steps):
#   1. Core computation      -> run_core_analysis() from the serumsource package
#                               (pure math: read, classify, source-correct)
#   2. Downstream analysis   -> plots, CSV/XLSX tables, Markdown report, notebook
#                               (run_downstream_analysis() from downstream_analysis.R)
#
# Before running, edit these three paths:
#   - input_file          : your MaxQuant annotated workbook (.xlsx)
#   - output_dir          : where results should be written
#   - downstream_script   : location of downstream/downstream_analysis.R
#
# Run with:  Rscript run_serum_source_analysis_newdata.R
# ============================================================================

library(serumsource)

# --- 1. Core computation (no files are written here) -----------------------
# target_groups = NULL means: analyse every group in the sample mapping
# except the reference group. Give a character vector to analyse only some
# groups, e.g. target_groups = c("DGEV_PC_1-3", "UCEVPC").
core <- run_core_analysis(
  input_file      = "path/to/MaxQuant_annotated_with_total.xlsx",
  sheet_name      = 1,
  sample_info     = build_sample_info(),
  reference_group = "UCEV",
  target_groups   = NULL,
  params          = serum_source_params()
)

# Classification summary preview
print(core$case_summary, n = 100)

# --- 2. Downstream analysis (plots, tables, report, notebook) --------------
# Point downstream_script to wherever downstream/downstream_analysis.R lives.
downstream_script <- "downstream/downstream_analysis.R"
if (!file.exists(downstream_script)) {
  stop("Cannot find downstream script: ", downstream_script,
       "\nPlease point downstream_script to downstream/downstream_analysis.R")
}
source(downstream_script)

# IMPORTANT: sample_info / params must be the same objects as in step 1,
# otherwise the two steps would describe different experiments.
run_downstream_analysis(
  core_result = core,
  input_file  = "path/to/MaxQuant_annotated_with_total.xlsx",
  output_dir  = "path/to/output",
  sample_info = build_sample_info(),
  params      = serum_source_params()
)