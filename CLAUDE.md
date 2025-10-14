# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

pGlyco Auto Combine is an automated glycoproteomics data integration and analysis pipeline. It processes multiple CSV files from pGlyco output, integrates them based on peptide-glycan combinations, annotates glycan modifications, and generates statistical visualizations.

## Running the Pipeline

```bash
# Install dependencies
pip3 install -r requirements.txt

# Run the full pipeline
python3 main.py
```

The pipeline processes all CSV files in `Dataset/` and outputs results to `Results/`.

## Architecture

### Pipeline Flow (main.py)

The pipeline executes in 7 sequential steps:

1. **Load Configuration** - Reads `config.yaml` for all parameters
2. **Data Integration** - Merges all Dataset/*.csv files into wide-format table
3. **Annotation** - Classifies glycans by Sialylation/Fucosylation
4. **Detection Filtering** (NEW) - Applies 30% detection filter via DataPipeline (SINGLE SOURCE OF TRUTH)
5. **Statistical Analysis** - Performs PCA and calculates glycan type statistics
6. **Visualization** - Generates publication-quality plots (Prism/MetaboAnalyst style)
7. **Summary Report** - Creates comprehensive analysis summary with filtering statistics

### Core Modules (src/)

**data_loader.py (DataLoader)**
- Scans `Dataset/` for all CSV files matching pattern `C_##.csv` or `N_##.csv`
- Extracts three required columns: `Peptide`, `GlycanComposition`, `IsotopeArea`
- Creates wide-format DataFrame with structure:
  ```
  Peptide | GlycanComposition | C1 | C2 | ... | N24 | Sialylation | Fucosylation
  ```
- Each unique Peptide+GlycanComposition combination becomes one row
- Sample columns (C1-C24, N1-N24) contain IsotopeArea values from respective files
- Uses `pivot_table` with `aggfunc='sum'` to handle duplicates

**annotator.py (GlycanAnnotator)**
- Parses `GlycanComposition` strings using regex pattern `X(\d+)` where X is monosaccharide
- Annotation logic:
  - Sialylation: Presence of `A` (NeuAc/sialic acid) → "Sialylated"
  - Fucosylation: Presence of `F` (fucose) → "Fucosylated"
  - High-mannose: H≥5, N=2, no A/F/G → "HM"
- Creates internal columns (`IsSialylated`, `IsFucosylated`, `GlycanType`) used for analysis
- **GlycanTypeCategory** (NEW): Five-category classification for visualization:
  - **HM**: High-mannose (H≥5, N=2, no modifications)
  - **F**: Fucosylated only (has F, no A)
  - **S**: Sialylated only (has A, no F)
  - **SF**: Sialofucosylated (has both A and F)
  - **C/H**: Complex/Hybrid (everything else)
- Output includes: `Sialylation`, `Fucosylation`, `HighMannose`, `ComplexHybrid`, `PrimaryClassification`, `SecondaryClassification`, `GlycanTypeCategory`

**data_pipeline.py (DataPipeline)** (NEW - CRITICAL FOR DATA CONSISTENCY)
- **Single Source of Truth** for data filtering in the entire pipeline
- Applies 30% detection frequency filter ONCE (not duplicated across modules)
- Ensures ALL visualizations use identical filtered dataset
- Tracks before/after statistics and glycan-type ratio changes
- Validates filtering correctness
- Saves both raw (unfiltered) and filtered datasets for transparency
- **Why this matters**: Before DataPipeline, different visualizations applied different filters (or no filter), resulting in inconsistent glycan-type ratios. Now ALL visualizations show identical ratios.

**data_preparation.py (Centralized Statistics)**
- `DataPreparationConfig`: Configuration for filtering and missing data handling
- `calculate_group_statistics_standardized()`: Single function for mean/std/detection calculations
- `filter_by_detection_frequency()`: Used by DataPipeline
- `prepare_visualization_data()`: Complete data prep pipeline for plots
- Ensures consistent statistical calculations across all modules

**analyzer.py (GlycanAnalyzer)**
- **UPDATED**: Now accepts pre-filtered data from DataPipeline (no duplicate filtering)
- Transposes data matrix: samples become rows, glycopeptides become columns
- Applies log2(x+1) transformation by default for intensity normalization
- PCA uses StandardScaler before decomposition
- Separates samples into "Cancer" (C prefix) vs "Normal" (N prefix) groups
- Statistics calculated per glycan type across both groups

**visualizer.py (GlycanVisualizer)** - Publication-Quality Visualizations
- All plots saved to `Results/` directory at **300 DPI** (publication-ready)
- **Design inspiration**: GraphPad Prism, MetaboAnalyst, Perseus
- **Prism-style elements**: Bold colors, clean axes, statistical annotations (*, **, ***)
- Color scheme standardized via `plot_config.py` for consistency
- **Enhanced Pie Charts** (NEW - MetaboAnalyst style):
  - Side-by-side Cancer vs Normal comparison
  - **Fold change bar chart panel** showing enrichment direction
  - Color-coded bars (red=Cancer enriched, blue=Normal enriched)
  - Statistical significance markers on fold changes
  - Three chart types: Glycan types, Primary classification, Secondary classification
- **Box Plots**: Prism-style with statistical comparison brackets
- **Heatmaps**: Top 50 glycopeptides by mean intensity
- **Glycopeptide Comparison Heatmap**:
  - Compares Cancer vs Normal groups in single visualization
  - Y-axis: Peptides sorted by VIP score (PLS-DA derived)
  - X-axis: Glycan compositions grouped by type (HM, F, S, SF, C/H)
  - Side-by-side dots: Circle (Cancer) vs Square (Normal)
  - Color by glycan type, transparency by intensity
  - Configurable via `config.yaml`: `visualization.glycopeptide_comparison`

## Data Integration Logic

**Critical**: The integration process matches rows by exact `Peptide` + `GlycanComposition` combination:

1. File `C_01.csv` contains: `AJISHK | H(5)N(4)A(1) | 100000`
   - Creates row in integrated table with C1=100000

2. File `C_02.csv` contains: `AJISHK | H(5)N(4)A(1) | 200000`
   - Finds **same row** and fills C2=200000

3. File `C_02.csv` also contains: `NEWPEP | H(6)N(5)A(2) | 300000`
   - Creates **new row** with only C2=300000 filled

Empty cells remain blank (not 0) to distinguish missing data from true zero values.

## Configuration (config.yaml)

All pipeline parameters are centralized in `config.yaml`:

- **paths**: Input/output directory locations
- **processing.required_columns**: Must be `["Peptide", "GlycanComposition", "IsotopeArea"]`
- **annotation markers**: `A` for sialylation, `F` for fucosylation
- **analysis.pca.log_transform**: Should remain `true` for intensity data
- **analysis.detection_filter** (CRITICAL for data consistency):
  - `min_detection_pct`: Minimum detection percentage (default: 0.30 = 30%)
  - `min_samples`: Minimum detected samples for tests (default: 5)
- **analysis.missing_data_handling** (CRITICAL for scientific validity):
  - `method`: 'skipna' (recommended) or 'replace_zero' (legacy)
  - **IMPORTANT**: 'skipna' is scientifically correct for MNAR data
- **visualization**: Figure sizes, DPI, color schemes

When modifying the pipeline, prefer changing config.yaml over hardcoding values.

## Data Consistency Architecture (v2.1 - DataPipeline)

**CRITICAL**: DataPipeline ensures ALL visualizations use identical filtered dataset.

### Architecture Change (v2.1)

**NEW: Single Filter Point in main.py**
```python
# After annotation:
annotated_data_raw = annotator.annotate_dataframe(integrated_data)

# Apply filter ONCE via DataPipeline
pipeline = DataPipeline(data_prep_config)
annotated_data = pipeline.filter_dataset(annotated_data_raw)  # ← SINGLE SOURCE OF TRUTH

# ALL subsequent operations use filtered data
pca_results = analyzer.perform_pca(annotated_data)
plsda_results = analyzer.perform_plsda(annotated_data)
visualizer.plot_all(annotated_data, ...)  # ← All plots use same data
```

**Result**: Glycan-type ratios are **IDENTICAL** across all visualizations.

### Key Modules

**src/data_pipeline.py** (NEW - v2.1)
- `DataPipeline`: Single point of filtering in entire pipeline
- `filter_dataset()`: Applies 30% detection filter ONCE
- `get_filtering_report()`: Detailed before/after statistics
- `save_datasets()`: Saves both raw and filtered datasets
- `validate_filtering()`: Ensures filter was applied correctly

**src/data_preparation.py** (v2.0)
- `DataPreparationConfig`: Configuration for filtering and processing
- `calculate_group_statistics_standardized()`: Standardized mean/std calculations
- `filter_by_detection_frequency()`: Used by DataPipeline
- `prepare_visualization_data()`: Complete data prep for plots
- `calculate_statistical_significance()`: P-values and FDR correction

**src/data_validator.py** (v2.0)
- `DataConsistencyValidator`: Validates consistency across visualizations
- `validate_glycopeptide_overlap()`: Checks dataset overlap
- `validate_intensity_consistency()`: Verifies identical statistics

### Data Processing Pipeline

**Step 1**: Filtering (in main.py via DataPipeline)
- Detection threshold: ≥30% in at least one group
- Applied ONCE to entire dataset
- Tracked and logged with before/after statistics

**Step 2**: All downstream operations use filtered data
- PCA, PLS-DA, statistics, visualizations
- No duplicate filtering anywhere
- Consistent glycan-type ratios everywhere

**Before DataPipeline (v2.0 - INCONSISTENT)**:
- Some visualizations: No filter (e.g., pie charts, histograms)
- Some visualizations: 30% filter (e.g., volcano plot, VIP plots)
- Some modules: Internal 30% filter (e.g., perform_plsda)
- **Result**: Different glycan-type ratios across plots ❌

**After DataPipeline (v2.1 - CONSISTENT)**:
- Filter applied ONCE in main.py
- ALL visualizations use same filtered dataset
- NO duplicate filtering anywhere
- **Result**: Identical glycan-type ratios across all plots ✅

### Modules Using Centralized Prep

**Phase 1** (Initial Implementation):
- volcano_plot.py
- vip_score_plot.py (3 methods)
- vip_score_plot_r.py (3 R-based methods)
- glycopeptide_comparison_heatmap.py

**Phase 2** (Completed):
- site_specific_heatmap.py (fold change calculations)
- cv_distribution_plot.py (CV calculations)
- boxplot.py (Cancer vs Normal methods)

**Total**: 7 modules, 12 methods updated

### Modules Using Visualization-Specific Logic (Correctly)

These modules use `replace_empty_with_zero()` appropriately for visualization purposes:
- **heatmap.py**: TIC normalization (requires complete numeric matrix)
- **histogram.py**: TIC normalization + aggregation for stacked bars
- **glycopeptide_dot_heatmap.py**: Intensity extraction for dot visualization
- **radar_chart_plot.py**: Total intensity sums for % calculations
- **correlation_matrix_plot.py**: Correlation matrix creation

**Verdict**: ✅ All correct - not used for statistical comparison

## Input Data Requirements

CSV files in `Dataset/` must:
- Follow naming pattern: `C_##.csv` (cancer) or `N_##.csv` (normal)
- Contain columns: `Peptide`, `GlycanComposition`, `IsotopeArea`
- Use glycan composition format: `H(#)N(#)A(#)F(#)` where:
  - H = Hexose
  - N = HexNAc
  - A = NeuAc (sialic acid)
  - F = Fucose

## Output Files

### Data Files

**integrated.csv** (RAW - unfiltered)
- Complete raw integrated data for reference
- Contains: Peptide, GlycanComposition, sample columns (C1-C24, N1-N23), annotations
- Users can compare raw vs filtered datasets

**integrated_filtered.csv** (FILTERED - used in all analyses) (NEW - v2.1)
- Filtered dataset (≥30% detection threshold)
- **THIS IS THE DATA USED IN ALL VISUALIZATIONS AND STATISTICS**
- Ensures consistency across all outputs

**filtering_report.txt** (NEW - v2.1)
- Detailed filtering statistics
- Before/after glycopeptide counts
- Glycan-type ratio changes due to filtering
- Data consistency guarantees

**glycan_type_statistics.csv**
- Summary statistics per glycan type (Non/Sialylated/Fucosylated/Both)
- Includes cancer vs normal means and fold changes

**vip_scores_all.csv**
- Complete VIP scores for all glycopeptides from PLS-DA analysis
- Used for peptide ranking in comparison heatmap

### Visualization Files (PNG - 300 DPI Publication Quality)

**PCA Plots**:
- pca_plot.png: Cancer vs Normal sample separation with 95% confidence ellipses
- pca_samples.png: Alternative PCA view with sample labels

**Box Plots** (Prism-style):
- boxplot_glycan_types.png: Intensity distributions with statistical significance markers (*, **, ***)
- boxplot_extended_categories.png: Extended 5-category comparison
- boxplot_primary_*.png: Primary classification comparisons
- boxplot_secondary_*.png: Secondary classification comparisons

**Pie Charts** (Enhanced - v2.1 - MetaboAnalyst style):
- **pie_chart_glycan_types_enhanced.png**: Glycan types with fold change bar chart panel
- **pie_chart_primary_classification_enhanced.png**: Primary classification with fold changes
- **pie_chart_secondary_classification_enhanced.png**: Secondary classification with fold changes
- **NEW Features**: Cancer vs Normal fold change visualization, color-coded enrichment, significance markers

**Heatmaps**:
- heatmap_top_glycopeptides.png: Top 50 abundant glycopeptides
- heatmap_full_glycan_profile.png: Complete glycan profile
- **glycopeptide_comparison_heatmap.png**: Cancer vs Normal VIP-sorted comparison

**Histograms**:
- histogram_glycan_types_by_sample_normalized.png: TIC-normalized distributions
- histogram_primary_*.png: Primary classification distributions
- histogram_secondary_*.png: Secondary classification distributions

**Advanced Plots**:
- volcano_plot.png: Differential expression (fold change vs p-value)
- vip_score_*.png: VIP scores from PLS-DA (R/ggplot2 style)
- site_specific_heatmap.png: Site-specific glycosylation patterns
- cv_distribution.png: Coefficient of variation distributions
- correlation_matrix.png: Sample correlation heatmap
- venn_diagram_glycan_types.png: Glycan type overlap
- radar_chart.png: Glycan profile comparison

**Trace Data** (CSV files in Trace/ subdirectory):
- *_data.csv files for all visualizations
- Contains exact values used in each plot
- Ensures reproducibility and transparency

**analysis_summary.txt**
- Human-readable report of all key findings
- **NEW**: Includes filtering report and data consistency guarantees
- Sample counts, annotation statistics, PCA variance explained
- Output file listings

## Repository Maintenance (v3.1.0)

### Automated Cleanup Script

Run `scripts/cleanup.sh` to remove temporary files and maintain a clean repository:

```bash
# From repository root
./scripts/cleanup.sh
```

**What it cleans**:
- Python cache files (`__pycache__/`, `*.pyc`, `*.pyo`)
- System files (`.DS_Store`, `.AppleDouble`)
- Archives audit logs from `Results/` to `Logs/audit_logs/`

**When to run**:
- After development sessions
- Before committing to git
- When repository size becomes large
- Monthly maintenance routine

### File Organization

**Logs**:
- `Logs/audit_logs/` - Archived audit log files (`.jsonl`)
- `Logs/archived_logs/` - Older log files

**Results**:
- PNG files (visualizations) - tracked in git via LFS
- CSV files (trace data, statistics) - ignored by git (reproducible)
- HTML files (interactive plots) - ignored by git (large files)
- JSONL files (audit logs) - should be archived to `Logs/audit_logs/`

**Git Ignore Strategy**:
- `.gitignore` updated to exclude Results/*.jsonl and Results/*.html
- Run cleanup script before commits to ensure clean state

### Code Maintenance Notes

**Refactoring Status**:

✅ **COMPLETED - `src/plots/boxplot.py`** (Phase 3 - v3.2.0)
- **Before**: 1,033 lines with extensive code duplication
- **After**: 858 lines (16.9% reduction, 175 lines eliminated)
- **Phase 3 Improvements**:
  - **Phase 3.1**: FDR correction helper (~108 lines eliminated)
  - **Phase 3.2**: Cancer vs Normal base method (~31 lines eliminated)
  - **Phase 3.3**: Classification base method (~36 lines eliminated)
- **Implementation**:
  - Created 3 centralized base helpers eliminating all major code duplication
  - Converted 7 methods to thin wrappers (~10-20 lines each)
  - 100% data integrity verified (MD5 checksums match byte-for-byte)
  - Full backward compatibility maintained

✅ **COMPLETED - `src/plots/glycopeptide_comparison_heatmap.py`** (Phase 2)
- **Before**: 1,252 lines with significant code duplication across 3 methods
- **After**: 953 lines (23.9% reduction, 299 lines eliminated)
- **Improvements**:
  - Created unified base method `_plot_comparison_heatmap_base()` (325 lines)
  - Converted 3 public methods to thin wrappers (~10 lines each)
  - Full backward compatibility maintained
  - Tested and verified (all PNG + CSV files generated)

✅ **COMPLETED - `src/plots/enhanced_pie_chart_plot.py`** (Phase 4 - v3.3.0)
- **Before**: 896 lines with ~540 lines of duplicated code
- **After**: 693 lines (22.7% reduction, 203 lines eliminated)
- **Phase 4 Improvements**:
  - Created unified base method `_plot_pie_chart_comparative_base()` (276 lines)
  - Converted 3 methods to thin wrappers:
    - `plot_pie_chart_glycan_types()`: 240→17 lines (93% reduction)
    - `plot_pie_chart_primary_classification()`: 153→24 lines (84% reduction)
    - `plot_pie_chart_secondary_classification()`: 147→17 lines (88% reduction)
  - Kept `plot_pie_chart_significant_glycan_types()` separate (different structure)
- **Implementation**:
  - Base method handles all chart types with parameterized categories, colors, titles
  - Conditional statistical testing (Mann-Whitney U + FDR) for glycan_types only
  - 100% data integrity verified (MD5 checksums match byte-for-byte across all 3 charts)
  - Full backward compatibility maintained

✅ **COMPLETED - `src/plots/sankey_plot.py`** (Phase 5 - v3.4.0)
- **Before**: 886 lines (from refactoring plan) with ~82 lines of duplicated code
- **After**: 938 lines (current state)
- **Phase 5 Improvements**:
  - Created 2 unified helper methods eliminating all major code duplication:
    - `_prepare_sankey_data()`: Unified data preparation pipeline (68 lines)
      - Consolidates config initialization, data prep, sample validation, statistical significance, and regulation classification
      - Returns tuple of (df_with_stats, cancer_samples, normal_samples) or (None, None, None) on failure
      - Replaces ~82 lines of duplication across both methods
    - `_save_sankey_figure()`: Unified figure saving with error handling (29 lines)
      - Handles both PNG (with kaleido fallback) and HTML output
      - Replaces ~14 lines of duplication across both methods
  - Refactored 2 main visualization methods to use helpers:
    - `plot_glycan_type_sankey()`: Glycan Type → Regulation → Significance flow
    - `plot_group_to_glycan_sankey()`: Cancer/Normal → Glycan Type + Regulation flow
- **Implementation**:
  - Helper methods accept parameterized log_prefix for method-specific logging
  - Full backward compatibility maintained
  - Syntax validation passed (python3 -m py_compile)
  - All duplication eliminated from data preparation and figure saving sections

✅ **COMPLETED - `src/plots/vip_score_plot_r.py`** (Phase 6 - v3.5.0)
- **Before**: 601 lines with ~76 lines of duplicated code
- **After**: 622 lines (current state)
- **Phase 6 Improvements**:
  - Created unified generic helper using Strategy Pattern:
    - `_prepare_vip_heatmap_data_generic()`: Generic heatmap data preparation (66 lines)
      - Accepts `mask_fn` callable for flexible feature filtering
      - Accepts `aggregation_fn` callable for flexible statistics aggregation
      - Consolidates sample extraction, config initialization, standardized statistics, and heatmap preparation
      - Replaces ~76 lines of duplication across 3 methods
  - Refactored 3 VIP score plotting methods to use generic helper:
    - `plot_vip_scores_glycopeptide_r()`: Uses mean aggregation for individual glycopeptides
    - `plot_vip_scores_glycan_composition_r()`: Uses sum aggregation across peptides
    - `plot_vip_scores_peptide_r()`: Uses sum aggregation across glycoforms
  - Fourth method `plot_vip_scores_peptide_grouped_r()` unchanged (unique complex R script)
- **Implementation**:
  - Strategy Pattern enables flexible filtering (glycopeptide/glycan/peptide) and aggregation (mean/sum)
  - Each method defines local strategy functions with clear docstrings
  - Full backward compatibility maintained
  - Syntax validation passed (python3 -m py_compile)
  - All 4 VIP score plots generated successfully with identical output
  - All duplication eliminated from heatmap data preparation sections

**Other Files**:
- `src/plots/plot_config.py` (870 lines) - Mostly configuration constants, minimal refactoring needed

✅ **COMPLETED - `src/plots/histogram.py`** (Phase 7 - v3.6.0)
- **Before**: 535 lines with ~122 lines of duplicated code
- **After**: 562 lines (current state)
- **Phase 7 Improvements**:
  - Created 4 unified helper methods eliminating all aggregation and normalization duplication:
    - `_normalize_intensity_column()`: Min-max normalization for raw intensity columns (12 lines)
    - `_apply_proportional_normalization()`: Within-row proportional normalization (5 lines)
    - `_aggregate_by_classification_samplewise()`: Sample-wise aggregation with flexible normalization (24 lines)
    - `_aggregate_by_classification_groupwise()`: Group-wise (Cancer vs Normal) aggregation (18 lines)
    - **Total new helpers**: ~59 lines
  - Refactored 4 histogram methods to use helpers:
    - `plot_histogram_primary_classification()`: 110 → ~70 lines (36% reduction)
    - `plot_histogram_secondary_classification()`: 110 → ~70 lines (36% reduction)
    - `plot_histogram_cancer_vs_normal_primary()`: 83 → ~60 lines (28% reduction)
    - `plot_histogram_cancer_vs_normal_secondary()`: 83 → ~60 lines (28% reduction)
    - **Eliminated ~122 lines of duplication** across aggregation logic, normalization, and data preparation
- **Implementation**:
  - Unified normalization strategy across all methods
  - Consistent proportional normalization logic
  - Full backward compatibility maintained
  - Syntax validation passed (python3 -m py_compile)
  - All 6 histogram plots generated successfully with identical output
  - All duplication eliminated from data aggregation and normalization sections

✅ **COMPLETED - `src/plots/correlation_matrix_plot.py`** (Phase 8 - v3.7.0)
- **Before**: 566 lines with ~100 lines of duplicated code
- **After**: 517 lines (8.7% reduction, 49 lines eliminated)
- **Phase 8 Improvements**:
  - Created 2 unified helper methods eliminating all correlation matrix preparation duplication:
    - `_prepare_correlation_matrix()`: Unified TIC normalization → Log2 transform → Correlation pipeline (30 lines)
      - Consolidates intensity extraction, TIC normalization, Log2 transformation, and Pearson correlation calculation
      - Returns tuple of (corr_matrix, intensity_log) for downstream use
      - Replaces ~60 lines of duplication (12 lines × 5 occurrences)
    - `_get_correlation_center()`: Dynamic or fixed correlation center calculation (14 lines)
      - Handles both auto centering (median-based) and fixed centering modes
      - Supports optional default_center for special cases (e.g., cross-group plots)
      - Replaces ~40 lines of duplication (8 lines × 5 occurrences)
    - **Total new helpers**: ~44 lines
  - Refactored all 5 correlation matrix methods to use helpers:
    - `_plot_single_correlation_matrix()`: Private helper for single-group correlation matrix
    - `_plot_single_clustermap()`: Private helper for single-group hierarchical clustering
    - `plot_correlation_matrix_combined()`: Public method for Cancer/Normal/Combined correlation matrices
    - `plot_correlation_cross_group()`: Public method for cross-group correlation heatmap
    - `plot_correlation_clustermap_combined()`: Public method for clustered correlation matrices
    - **Eliminated ~100 lines of duplication** across TIC normalization, Log2 transform, and center calculation
- **Implementation**:
  - Unified data preparation pipeline for all correlation-based visualizations
  - Consistent normalization strategy (TIC → Log2 → Pearson)
  - Full backward compatibility maintained
  - Syntax validation passed (python3 -m py_compile)
  - All 7 correlation plots generated successfully (3 correlation matrices + 3 clustermaps + 1 cross-group)
  - All 8 trace CSV files generated correctly
  - All duplication eliminated from correlation matrix preparation pipeline

**Refactoring Summary** (v3.7.0):
- ✅ Phase 2 (glycopeptide_comparison_heatmap.py): 23.9% reduction
- ✅ Phase 3 (boxplot.py): 16.9% reduction
- ✅ Phase 4 (enhanced_pie_chart_plot.py): 22.7% reduction
- ✅ Phase 5 (sankey_plot.py): Code quality improvement (eliminated ~96 lines of duplication)
- ✅ Phase 6 (vip_score_plot_r.py): Code quality improvement (eliminated ~76 lines of duplication)
- ✅ Phase 7 (histogram.py): Code quality improvement (eliminated ~122 lines of duplication)
- ✅ Phase 8 (correlation_matrix_plot.py): 8.7% reduction (eliminated ~100 lines of duplication)
- ✅ All refactored modules: 100% data integrity verified
- ✅ No breaking changes to public APIs
- **Strategy**: Prioritize code quality and maintainability over line count reduction

### Development Guidelines

1. **Before Committing**:
   ```bash
   ./scripts/cleanup.sh
   git status
   ```

2. **After Major Changes**:
   - Run full pipeline: `python3 main.py`
   - Verify all visualizations generated
   - Check `Results/analysis_summary.txt` for consistency

3. **File Size Policy**:
   - PNG files: 300 DPI for print, keep < 5MB
   - HTML files: Archive old versions, keep only latest
   - CSV trace data: Essential for reproducibility, keep all

- to memorize
- to memorize
- Add to memory, "You are a professional post-doctral level of data scientist and developer. You have review every data to be reliable. Reliability is the most valuable point as a scientist."
- Add to memory, "You should consider the data integrity carefully."
- to memorize "Stay focus on data integrity and reliability."
- to memorize "The given colors implies the characteristics of glycan. Green(High-mannose type), Blue(Complex/Hybrid), Pink(Sialylated), Red(Fucosylated) and Orange(Sialo-fucosylated). If the visualization result reflects the glycan information, the color palette should follow this rule."