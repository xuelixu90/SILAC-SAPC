# SerumSourceAnalyzer

**Heavy/Light Isotope-Based Serum Source Analysis R Package**

Analyzes SILAC-style Heavy/Light dual-channel labeled EV (extracellular vesicle) proteomics data. Using UCEV as the control group, the package classifies each protein as serum-derived versus residual/background, performs p0-based source correction to estimate serum fractions, annotates candidates against an in-house serum database and the Human Plasma PeptideAtlas, and outputs structured result tables, figures, and a Markdown analysis report.

This package is a modular refactor of the `Heavy_Light-Isotope-Based_Serum_Source_Analysis_newdata.R` hard-coded script. The analysis logic is fully aligned with the original, with added visualization and reporting capabilities.

> **Consistency validation**: On `MaxQuant_annotated_with_total.xlsx` (5298 proteins x 116 columns),
> the package output matches the original script value-by-value -- classification summary 54 rows with
> zero differences, 535 primary serum-derived proteins with identical sets, and maximum numeric
> difference of 0 across all numeric columns (excess / serum fraction / p0 / CV, etc.).

---

## 1. Workflow Overview

```
MaxQuant annotated table (xlsx)
        |
        v
+---------------------------------------------+
| read_input_data()      read + NA handling   |
| prepare_long_data()    long table + rep     |
|   validity flags                           |
|   +-- Paired_LH_Valid_Rep  both > 0         |
|   +-- Target_L_Dominant_Rep L/(L+H) >= 0.70 |
+---------------------------------------------+
        |
        v
+---------------------------------------------+
| compute_ucev_background()                   |
|   per-protein robust SD of Light fraction   |
|   (MAD-based)                               |
|   -> global SD floor (max(0.02, median))    |
+---------------------------------------------+
        |
        v
+---------------------------------------------+
| analyze_one_group()   for each target group |
|   Case 1/2: control + target both have      |
|     paired data                             |
|     -> Source-corrected high/medium          |
|       confidence, single-peptide candidates,|
|       Likely_residual, Ambiguous            |
|   Case 3/4: control both channels missing   |
|     -> Acquired-only Level 1 / Level 2      |
| calibrate_serum_light()                     |
|     -> Serum_L_Calibrated_<sample> columns  |
+---------------------------------------------+
        |
        +--> CSV / XLSX result tables (save_analysis_outputs)
        +--> 7 key figures (visualization.R, PNG)
        +--> Analysis_Report.md (write_analysis_report)
```

## 2. Installation

```r
# Dependencies: readxl, dplyr, magrittr, purrr, tibble, stringr, tidyr, writexl, ggplot2
install.packages(c("readxl", "dplyr", "magrittr", "purrr", "tibble",
                   "stringr", "tidyr", "writexl", "ggplot2"))

# Install from local source (optional --build to produce a zip)
install.packages("./code/SerumSourceAnalyzer", repos = NULL, type = "source")
```

## 3. Quick Start (one-call full analysis)

```r
library(SerumSourceAnalyzer)

res <- run_serum_source_analysis(
  input_file  = "./MaxQuant_annotated_with_total.xlsx",
  output_dir  = "./Serum_source_analysis_all_groups_1",
  make_plots  = TRUE,
  write_report = TRUE
)
```

After completion, `output_dir` contains:

| File | Content |
|---|---|
| `Analysis_Report.md` | **Markdown analysis report (with all figures)** |
| `Analysis_Notebook.nb.html` | **Interactive HTML notebook — open in any browser to review the full execution log (code + output + figures)** |
| `Analysis_Notebook.html` / `Analysis_Notebook.md` | Print-friendly HTML and plain Markdown versions of the notebook |
| `plots/*.png` | 7 key visualizations |
| `00_All_Comparisons_Merged_Results.csv` | Per-protein results for all comparisons |
| `00_All_Comparisons_Primary_Serum_Derived_Merged.csv` | Merged primary serum-derived proteins |
| `00_Classification_Summary_by_Group.csv` | Classification summary |
| `00_Database_Annotation_Summary_All_Groups.csv` | Database evidence summary |
| `<group>_vs_UCEV/` subdirectories | Full result tables by category for each group |
| `<group>_vs_UCEV_Primary_Serum_Proteins_with_Database_Annotation.csv` | Per-group primary serum proteins + annotation |
| `Serum_Source_Analysis_All_Groups_Summary.xlsx` | Key-table summary workbook |

The return value `res` is a list containing `results` (per-protein statistics), `all_merged`, `primary_serum`, `case_summary`, `plots`, etc., which can be used for further analysis.

### 3.1 Viewing the HTML Notebook Interface

The package ships an execution notebook (`inst/rmarkdown/analysis_notebook.Rmd`) that re-runs the whole pipeline and captures code, console messages, intermediate tables, and every figure into one **self-contained HTML file**. Three ways to view it:

**(a) Automatic (default)** — `run_serum_source_analysis()` renders the notebook as part of the run (`render_notebook = TRUE` is the default). Double-click the generated file to open it in any browser:

```
output_dir/
├── Analysis_Notebook.nb.html   <- interactive notebook (code folding, TOC)
├── Analysis_Notebook.html      <- print-friendly HTML
└── Analysis_Notebook.md        <- plain Markdown
```

**(b) Standalone function** — render the notebook for an existing analysis without re-running the full pipeline wrapper:

```r
SerumSourceAnalyzer::render_analysis_notebook(
  input_file = "F:/WYF/NP_1/MaxQuant_annotated_with_total.xlsx",
  output_dir = "F:/WYF/NP_1/Serum_source_analysis_all_groups_1"
)
```

**(c) Manual `rmarkdown::render()`** — full control over parameters:

```r
rmarkdown::render(
  "inst/rmarkdown/analysis_notebook.Rmd",
  output_dir = "output_dir",
  params = list(
    INPUT_FILE  = "./MaxQuant_annotated_with_total.xlsx",
    OUTPUT_DIR  = "output_dir",
    PARAMS      = list(),
    SAMPLE_INFO = NULL,
    LIB_DIR     = NULL
  )
)
```

**Requirements**: `rmarkdown` + `knitr` (`install.packages(c("rmarkdown", "knitr"))`) and pandoc (bundled with RStudio; the package auto-detects RStudio's pandoc in non-interactive sessions). If rendering is skipped due to missing dependencies, a warning is emitted and the rest of the analysis is unaffected.

## 4. Analysis Logic Details

### 4.1 Replicate validity flags (long table `lh_long`)

| Flag | Definition |
|---|---|
| `Paired_LH_Valid_Rep` | Both Light and Heavy are finite and > 0 |
| `Target_L_Dominant_Rep` | Light > 0 and (Heavy missing/zero, or L/(L+H) >= 0.70) |
| `High_Pep_Matched_Rep` | Valid pair and unique peptides >= 2 |
| `Single_Pep_Matched_Rep` | Valid pair and unique peptides = 1 |

Light fraction: `Light_fraction = Light / (Light + Heavy)`.

### 4.2 UCEV background

For each protein in the UCEV group, the robust SD (MAD x 1.4826) of the Light fraction is computed. The global floor `global_fraction_sd_floor = max(0.02, median(per-protein SD))` stabilizes thresholds for low-variability proteins.

### 4.3 Case determination and classification

**Case 1**: Control and target both have >= 2 valid paired replicates, and >= 2 target replicates have >= 2 unique peptides.
Within Case 1/2, further quantitative grading applies (threshold = max(0.05, k x sigma), k=3 for high confidence / k=2 for medium):

| Classification | Condition |
|---|---|
| Source-corrected, high confidence | Excess >= high threshold and >= 2 reps support and CV <= 0.50 and 95% CI separated |
| Source-corrected, medium confidence | Excess >= medium threshold and >= 2 reps support |
| Source-corrected, single-peptide (strong) candidate | Case 2 version of the above rules |
| Likely_residual | Excess <= 0 |
| Ambiguous | All others |

**Case 3/4**: Both channels completely absent in UCEV (target has Light-dominant signal) -> `Acquired-only (Level 1)` (peptides >= 2 support) / `Acquired-only (Level 2)` (single-peptide support).

`Primary serum categories` = the 4 strong-evidence categories above (Acquired L1/L2 + corrected high/medium confidence).

### 4.4 p0 source correction and serum fraction

`p0` = median of Heavy/(Heavy+Light) in UCEV (the heavy-labeling rate). For each target replicate:

```
A = Heavy / p0                     # total EV estimate
expected_residual_L = A * (1 - p0)
serum_L = max(Light - expected_residual_L, 0)
serum_fraction = serum_L / (A + serum_L)
```

`calibrate_serum_light()` similarly converts each replicate's excess back to intensity units, producing `Serum_L_Calibrated_<sample>` columns (raw Light intensity is used directly when the control is absent).

## 5. Step-by-Step Usage (flexible API)

```r
raw  <- read_input_data("./MaxQuant_annotated_with_total.xlsx")
si   <- build_sample_info()                # sample mapping (en dash column names)
prep <- prepare_long_data(raw, si)         # $data_converted / $annotation / $lh_long
bg   <- compute_ucev_background(prep$lh_long)

r1 <- analyze_one_group(prep$lh_long, prep$annotation,
                        "UCEV", "UCEVPC",
                        global_fraction_sd_floor = bg$global_fraction_sd_floor)
```

## 6. Visualization Functions

| Function | Figure | Key information |
|---|---|---|
| `plot_classification_summary()` | Stacked bar chart | Protein counts per classification per group; primary categories highlighted |
| `plot_excess_distribution()` | Violin + boxplot | Light-fraction excess distribution per group; primary categories color-coded |
| `plot_excess_vs_probability()` | Volcano-style scatter | Excess vs P(excess is real); top genes labeled |
| `plot_serum_fraction()` | Boxplot + scatter | Serum fraction estimates for primary categories |
| `plot_database_evidence()` | Stacked bar chart | Database evidence composition for primary-category proteins |
| `plot_ucev_background()` | Histogram | UCEV background SD distribution with global floor |
| `plot_excess_heatmap()` | Heatmap | Top N proteins' excess across comparisons |

```r
p <- plot_classification_summary(res$case_summary)
p
ggplot2::ggsave("cls.png", p, width = 9, height = 6, dpi = 300)
```

## 7. Custom Parameters

```r
p <- serum_source_params(
  min_absolute_light_excess = 0.05,   # absolute excess floor
  background_multiplier_high = 3,     # high-confidence sigma multiplier
  background_multiplier_medium = 2,    # medium-confidence sigma multiplier
  min_light_dominance_ratio = 0.70    # Light dominance ratio
)
res <- run_serum_source_analysis(..., params = p)
```

Adding new experimental groups: simply extend the mapping table returned by `build_sample_info()` (keep en dash column names), and pass the group names to `target_groups`.

## 8. Testing

```r
testthat::test_local("./code/SerumSourceAnalyzer")
```
