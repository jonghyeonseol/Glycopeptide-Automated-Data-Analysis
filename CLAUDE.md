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

The pipeline executes in 6 sequential steps:

1. **Load Configuration** - Reads `config.yaml` for all parameters
2. **Data Integration** - Merges all Dataset/*.csv files into wide-format table
3. **Annotation** - Classifies glycans by Sialylation/Fucosylation
4. **Statistical Analysis** - Performs PCA and calculates glycan type statistics
5. **Visualization** - Generates PCA plots, boxplots, heatmaps
6. **Summary Report** - Creates comprehensive analysis summary

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

**analyzer.py (GlycanAnalyzer)**
- Transposes data matrix: samples become rows, glycopeptides become columns
- Applies log2(x+1) transformation by default for intensity normalization
- PCA uses StandardScaler before decomposition
- Separates samples into "Cancer" (C prefix) vs "Normal" (N prefix) groups
- Statistics calculated per glycan type across both groups

**visualizer.py (GlycanVisualizer)**
- All plots saved to `Results/` directory at 300 DPI
- Uses seaborn "whitegrid" style
- Color scheme for glycan types defined in config.yaml
- Heatmap shows top 50 glycopeptides by mean intensity
- **Glycopeptide Comparison Heatmap** (NEW):
  - Compares Cancer vs Normal groups in single visualization
  - Y-axis: Peptides sorted by VIP score (PLS-DA derived)
  - X-axis: Glycan compositions grouped by type (HM, F, S, SF, C/H)
  - Side-by-side dots: Circle (Cancer) vs Square (Normal)
  - Color by glycan type, transparency by intensity
  - Bottom panel: Gradient colored bar showing glycan type regions
  - Top panel: Aggregated intensity comparison between groups
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

## Centralized Data Preparation (NEW - v2.0)

**CRITICAL**: All visualizations now use standardized data preparation to ensure consistency.

### Key Modules

**src/data_preparation.py**
- `DataPreparationConfig`: Configuration for filtering and processing
- `prepare_visualization_data()`: Single source of truth for data prep
- `calculate_group_statistics_standardized()`: Standardized mean/std calculations
- `filter_by_detection_frequency()`: Uniform detection filtering
- `calculate_statistical_significance()`: P-values and FDR correction

**src/data_validator.py**
- `DataConsistencyValidator`: Validates consistency across visualizations
- `validate_glycopeptide_overlap()`: Checks dataset overlap
- `validate_intensity_consistency()`: Verifies identical statistics

### Data Processing Pipeline

All visualizations follow this standardized pipeline:

1. **Detection Filtering**: 30% detection OR 5 samples in at least one group
2. **Mean Calculation**: Uses `skipna=True` to exclude missing values
3. **Statistics**: Same method for all Cancer_Mean, Normal_Mean, etc.
4. **Fold Change**: Log2(Cancer+1 / Normal+1) for symmetry

**Before (INCONSISTENT)**:
- Volcano plot: 20% filter, non-zero mean
- VIP plots: No filter, include-zero mean
- Comparison heatmap: 30% filter (inner join) + 50% filter = double filtering

**After (CONSISTENT)**:
- All visualizations: 30% filter, skipna mean
- Same glycopeptides appear in all plots
- Same intensity values across all visualizations

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

**integrated.csv**
- Primary output with all integrated sample data
- Contains: Peptide, GlycanComposition, sample columns (C1-C24, N1-N23), Sialylation, Fucosylation, HighMannose, ComplexHybrid, PrimaryClassification, SecondaryClassification, GlycanTypeCategory, Proteins
- Users should be able to open and analyze this file directly

**glycan_type_statistics.csv**
- Summary statistics per glycan type (Non/Sialylated/Fucosylated/Both)
- Includes cancer vs normal means and fold changes

**vip_scores_all.csv**
- Complete VIP scores for all glycopeptides from PLS-DA analysis
- Used for peptide ranking in comparison heatmap

**Visualization files** (PNG):
- pca_plot.png: Cancer vs Normal sample separation
- pca_samples.png: Alternative PCA view
- boxplot_glycan_types.png: Intensity distributions by glycan type
- glycan_type_distribution.png: Count of each glycan type
- heatmap_top_glycopeptides.png: Top 50 abundant glycopeptides
- **pie_chart_glycan_types.png** (v2.0): Glycan type distribution pie charts (Cancer vs Normal)
- **pie_chart_primary_classification.png** (v2.0): Primary classification pie charts
- **pie_chart_secondary_classification.png** (v2.0): Secondary classification pie charts
- **glycopeptide_comparison_heatmap.png**: Cancer vs Normal comparison with VIP-sorted peptides

**analysis_summary.txt**
- Human-readable report of all key findings
- Includes sample counts, annotation statistics, PCA variance explained
