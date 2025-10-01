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
- Creates internal columns (`IsSialylated`, `IsFucosylated`, `GlycanType`) used for analysis
- Only `Sialylation` and `Fucosylation` columns are included in final output

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
- **visualization**: Figure sizes, DPI, color schemes

When modifying the pipeline, prefer changing config.yaml over hardcoding values.

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

**integrated_example.csv**
- Primary output with all integrated sample data
- Contains only: Peptide, GlycanComposition, sample columns, Sialylation, Fucosylation
- Users should be able to open and analyze this file directly

**glycan_type_statistics.csv**
- Summary statistics per glycan type (Non/Sialylated/Fucosylated/Both)
- Includes cancer vs normal means and fold changes

**Visualization files** (PNG):
- pca_plot.png: Cancer vs Normal sample separation
- pca_samples.png: Alternative PCA view
- boxplot_glycan_types.png: Intensity distributions by glycan type
- glycan_type_distribution.png: Count of each glycan type
- heatmap_top_glycopeptides.png: Top 50 abundant glycopeptides

**analysis_summary.txt**
- Human-readable report of all key findings
- Includes sample counts, annotation statistics, PCA variance explained
