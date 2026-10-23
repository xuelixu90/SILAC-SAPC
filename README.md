# serumsource

**Find out which EV proteins come from serum — a core algorithm R package**

`serumsource` analyses SILAC Heavy/Light labelled extracellular-vesicle (EV)
proteomics data and answers one question per protein: **is this protein
really carried by serum, or is it just background / residual EV material?**

## How it works in one minute

1. Your input is a MaxQuant annotation table (xlsx). Each sample has a
   **Light** and a **Heavy** intensity channel.
2. For every protein and replicate we compute the **Light fraction** =
   Light / (Light + Heavy). Serum proteins are enriched in the Light
   channel, so a high Light fraction points towards serum.
3. We measure how much the Light fraction **wobbles in the control group
   (UCEV)** — that wobble is the background noise floor.
4. A target group is called **serum-derived** only when its Light fraction
   clearly exceeds that noise floor, with enough replicate support.
5. Finally, a **p0-based source correction** converts the excess into an
   estimated **serum fraction** (0-1) per protein.

> **Design principle**: this package does the **math only** — no plots, no
> file output. Plotting, CSV/XLSX export, the Markdown report and the HTML
> notebook live in a separate downstream script:
> [`downstream/downstream_analysis.R`](downstream/downstream_analysis.R).

> **Consistency check**: on the NP_1 input (5298 proteins x 116 columns)
> the package reproduces the original hard-coded script value-by-value
> (identical classification summary, identical 535 primary serum-derived
> proteins, maximum numeric difference 0).

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
        v
   run_core_analysis()  ->  list of results (no I/O)
        |
        v  (hand off to the downstream script)
+---------------------------------------------+
| downstream/downstream_analysis.R            |
|   +--> CSV / XLSX result tables             |
|   +--> 7 key figures (PNG)                  |
|   +--> Analysis_Report.md                   |
|   +--> Analysis_Notebook.nb.html            |
+---------------------------------------------+
```

## 2. Installation

```r
# Core dependencies
install.packages(c("readxl", "dplyr", "magrittr", "purrr", "tibble"))

# Option A: install from GitHub
# install.packages("remotes")
remotes::install_github("xuelixu90/SILAC-SAPC")

# Option B: install from local source (point pkg_dir to the folder
# containing DESCRIPTION)
pkg_dir <- "path/to/serumsource"
install.packages(pkg_dir, repos = NULL, type = "source")
```

> Downstream dependencies (only needed for the script, **not** the package):
> `ggplot2`, `writexl`, `stringr`, `tidyr`, `rmarkdown`, `knitr`.

## 3. Quick Start

### 3.1 Core computation (package)

```r
library(serumsource)

# Only input_file is required -- every other argument has a default matching
# the NP_1 dataset. Shown here with explicit arguments for clarity:
core <- run_core_analysis(
  input_file      = "path/to/MaxQuant_annotated_with_total.xlsx",  # required
  sample_info     = build_sample_info(),    # group/column mapping (NP_1 defaults)
  reference_group = "UCEV",                  # control group name
  target_groups   = NULL,                    # NULL = auto: all groups except reference
  params          = serum_source_params()    # thresholds; see Section 7 to customise
)

# New dataset with different groups or column names? Override what differs, e.g.:
#   run_core_analysis(input_file, sample_info = my_mapping,
#                     reference_group = "CTRL", target_groups = c("A", "B"))

# core$results        : per-protein classification (all comparisons)
# core$all_merged     : results with calibrated Serum_L columns
# core$primary_serum  : primary serum-derived proteins
# core$case_summary   : classification counts per group
print(core$case_summary, n = 100)
```

**All arguments of `run_core_analysis()` at a glance:**

| Argument | Default | When to change it |
|---|---|---|
| `input_file` | (required) | - |
| `sheet_name` | `1` | Your data is not on the first sheet |
| `sample_info` | `build_sample_info()` | Different samples/groups than NP_1 |
| `reference_group` | `"UCEV"` | A different control group name |
| `target_groups` | `NULL` (auto) | Analyse only a subset of groups |
| `params` | `serum_source_params()` | Loosen/tighten thresholds (Section 7) |
| `annotation_cols` | NP_1 column names | Your workbook uses different column names |
| `na` | MaxQuant NA strings | Other missing-value codes in your table |

`annotation_cols` is the one to watch for new datasets - it maps logical
keys to the actual column names in **your** workbook:

```r
run_core_analysis(input_file, annotation_cols = c(
  protein_ids       = "Protein IDs",
  representative_id = "Representive ID",
  protein_names     = "Protein names",
  gene_names        = "Gene names",
  in_house_db       = "In_house_db",
  plasma_atlas_db   = "Plasma_Altas_db"
))
```

### 3.2 Full analysis (core + downstream script)

```r
library(serumsource)

# 1. Core computation (all parameters have NP_1 defaults -- see 3.1)
core <- run_core_analysis(
  input_file = "path/to/MaxQuant_annotated_with_total.xlsx"
)

# 2. Downstream: plots, CSV/XLSX, report, notebook
source("downstream/downstream_analysis.R")

# Only core_result / input_file / output_dir are required.
# sample_info / params / group names MUST match the values used in step 1.
run_downstream_analysis(
  core_result     = core,                   # result list from run_core_analysis()
  input_file      = "path/to/MaxQuant_annotated_with_total.xlsx",
  output_dir      = "path/to/output",       # created if it does not exist
  sample_info     = build_sample_info(),    # group mapping (must match step 1)
  reference_group = "UCEV",                 # control group name
  target_groups   = NULL,                   # NULL = auto: all groups except
                                            # reference_group; or set explicitly,
                                            # e.g. c("DGEV_PC_1-3", "UCEVPC")
  params          = serum_source_params(),  # thresholds (must match step 1)
  make_plots      = TRUE,                   # generate & save the 7 figures
  save_outputs    = TRUE,                   # write CSV / XLSX tables
  write_report    = TRUE,                   # build Analysis_Report.md
  render_notebook = TRUE                    # render Analysis_Notebook.nb.html
)
```

After completion, `output_dir` contains:

| File | Content |
|---|---|
| `Analysis_Report.md` | **Start here** — Markdown report with all figures |
| `Analysis_Notebook.nb.html` | Interactive HTML notebook — full execution log |
| `Analysis_Notebook.html` / `.md` | Print-friendly HTML and plain Markdown |
| `plots/*.png` | 7 key visualizations |
| `00_All_Comparisons_Merged_Results.csv` | Per-protein results for all comparisons |
| `00_All_Comparisons_Primary_Serum_Derived_Merged.csv` | Merged primary serum-derived proteins |
| `00_Classification_Summary_by_Group.csv` | Classification summary |
| `00_Database_Annotation_Summary_All_Groups.csv` | Database evidence summary |
| `<group>_vs_UCEV/` subdirectories | Full result tables by category for each group |
| `Serum_Source_Analysis_All_Groups_Summary.xlsx` | Key-table summary workbook |

## 4. How the classification works

### 4.1 Per-replicate flags (long table `lh_long`)

For every protein x replicate we first decide whether the numbers can be
trusted:

| Flag | Plain-language meaning |
|---|---|
| `Paired_LH_Valid_Rep` | Both Light and Heavy intensities are present and positive |
| `Target_L_Dominant_Rep` | Light is clearly dominant (Light > 0 and L/(L+H) >= 0.70, or Heavy is absent) |
| `High_Pep_Matched_Rep` | Valid pair measured with >= 2 unique peptides |
| `Single_Pep_Matched_Rep` | Valid pair measured with only 1 unique peptide |

Light fraction: `Light_fraction = Light / (Light + Heavy)`.

### 4.2 The background noise floor (UCEV)

The control group UCEV should contain **no** serum-derived signal, so its
Light-fraction variation is pure background. For each protein we compute a
robust SD (MAD x 1.4826) across UCEV replicates; a global floor
`global_fraction_sd_floor = max(0.02, median(per-protein SD))` keeps
thresholds stable for very quiet proteins.

### 4.3 Classification rules

**Case 1/2** (control and target both have usable paired data): a protein is
called serum-derived when its Light-fraction **excess** over the background
passes `threshold = max(0.05, k x sigma)` (k = 3 high confidence, k = 2
medium):

| Classification | In plain words |
|---|---|
| Source-corrected, high confidence | Strong excess, >= 2 replicates agree, low CV, CI separated |
| Source-corrected, medium confidence | Clear excess, >= 2 replicates agree |
| Source-corrected, single-peptide (strong) candidate | Same rules but only 1 unique peptide available |
| Likely_residual | Excess <= 0 — looks like background |
| Ambiguous | Everything else |

**Case 3/4** (control has no signal at all, target is Light-dominant): the
protein was **acquired** by the EVs, not corrected against a background ->
`Acquired-only (Level 1)` (>= 2 peptides) / `Acquired-only (Level 2)`
(1 peptide).

`Primary serum categories` = the 4 strong-evidence labels above
(Acquired L1/L2 + corrected high/medium confidence).

### 4.4 p0 source correction: from "excess" to "serum fraction"

`p0` = median of Heavy/(Heavy+Light) in UCEV — the fraction of EV material
that stayed on the Heavy channel. For each target replicate:

```
A = Heavy / p0                     # total EV estimate
expected_residual_L = A * (1 - p0) # Light we would expect from residual EVs
serum_L = max(Light - expected_residual_L, 0)   # Light above that = serum
serum_fraction = serum_L / (A + serum_L)
```

`calibrate_serum_light()` does the same conversion in intensity units,
producing `Serum_L_Calibrated_<sample>` columns (raw Light intensity is used
directly when the control is absent).

## 5. Step-by-Step Usage (flexible API)

Prefer to run the pipeline piece by piece? Each exported function is
independent:

```r
raw  <- read_input_data("path/to/MaxQuant_annotated_with_total.xlsx")  # 1. read
si   <- build_sample_info()                # 2. sample->group mapping (en dash names)
prep <- prepare_long_data(raw, si)         # 3. long table + annotation
bg   <- compute_ucev_background(prep$lh_long)  # 4. background noise floor

r1 <- analyze_one_group(prep$lh_long, prep$annotation,   # 5. classify one group
                        "UCEV", "UCEVPC",
                        global_fraction_sd_floor = bg$global_fraction_sd_floor)
```

## 6. Exported Core Functions

| Function | What it does |
|---|---|
| `read_input_data()` | Read the annotated MaxQuant Excel workbook |
| `build_sample_info()` | Default sample-to-column mapping for the NP_1 dataset |
| `prepare_long_data()` | Build the long-format Heavy/Light table + annotation |
| `compute_ucev_background()` | Estimate the background noise floor from UCEV |
| `analyze_one_group()` | Classify one target group against one control |
| `calibrate_serum_light()` | Convert excess into calibrated serum Light intensity |
| `run_core_analysis()` | Run the whole computation in one call (no I/O) |
| `serum_source_params()` | All threshold parameters in one object |
| `primary_serum_categories()` | The four strong-evidence category labels |

## 7. Custom Parameters

All thresholds live in one object. The ones you are most likely to touch:

```r
p <- serum_source_params(
  min_absolute_light_excess = 0.05,    # smallest excess that counts as "real"
  background_multiplier_high = 3,      # sigma multiplier for high confidence
  background_multiplier_medium = 2,    # sigma multiplier for medium confidence
  min_light_dominance_ratio = 0.70     # L/(L+H) cutoff for Light dominance
)
core <- run_core_analysis(..., params = p)
```

<details>
<summary>All 10 thresholds (defaults shown)</summary>

| Parameter | Default | Meaning |
|---|---|---|
| `min_valid_reps` | 2 | Minimum paired valid replicates per protein |
| `min_absolute_light_excess` | 0.05 | Smallest excess that counts as "real" |
| `background_multiplier_high` | 3 | Sigma multiplier for high confidence |
| `background_multiplier_medium` | 2 | Sigma multiplier for medium confidence |
| `min_support_reps_when_three_valid` | 2 | Replicates needed when 3 valid reps exist |
| `min_support_reps_when_two_valid` | 2 | Replicates needed when 2 valid reps exist |
| `max_serum_fraction_cv_high_confidence` | 0.50 | Max CV for high-confidence calls |
| `minimum_fraction_sd_floor` | 0.02 | Lower floor on the Light-fraction SD |
| `min_light_dominance_ratio` | 0.70 | L/(L+H) cutoff for Light dominance |
| `eps` | 1e-12 | Numerical tolerance |

</details>

Adding new experimental groups: simply extend the mapping table returned by
`build_sample_info()` (keep en dash column names), and pass the group names
to `target_groups`.

## 8. Testing

```r
testthat::test_local("path/to/serumsource")
```