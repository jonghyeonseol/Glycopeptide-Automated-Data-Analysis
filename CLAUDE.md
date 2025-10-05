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
- to memorize
- to memorize