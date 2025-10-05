# pGlyco Auto Combine - Architecture

## Overview

This document describes the architecture, design patterns, and code organization of the pGlyco Auto Combine pipeline.

---

## Directory Structure

```
pGlyco_auto_combine/
├── Root Configuration
│   ├── README.md              # User-facing documentation
│   ├── CLAUDE.md              # AI assistant instructions
│   ├── ARCHITECTURE.md        # This file
│   ├── config.yaml            # Pipeline configuration
│   ├── main.py                # Entry point
│   ├── requirements.txt       # Python dependencies
│   └── .gitignore             # Git ignore patterns
│
├── docs/                      # Documentation
│   ├── README.md              # Documentation index
│   ├── CHANGELOG.md           # Version history
│   ├── glycan-sorting-guide.md
│   ├── normalization.md
│   ├── trace-data-reference.md
│   ├── verification-guide.md
│   ├── visualization-enhancements.md
│   ├── visualization-guide.md
│   └── visualization-updates.md
│
├── scripts/                   # Utility scripts
│   ├── README.md
│   └── verify_trace_data.py  # Data verification
│
├── src/                       # Source code
│   ├── Core Modules
│   │   ├── data_loader.py     # CSV integration
│   │   ├── annotator.py       # Glycan classification
│   │   ├── analyzer.py        # Statistical analysis
│   │   ├── visualizer.py      # Visualization coordinator
│   │   └── utils.py           # Shared utilities
│   │
│   └── plots/                 # Visualization modules
│       ├── boxplot.py
│       ├── correlation_matrix_plot.py
│       ├── cv_distribution_plot.py
│       ├── distribution_plot.py
│       ├── glycopeptide_comparison_heatmap.py ⭐ NEW
│       ├── glycopeptide_dot_heatmap.py
│       ├── heatmap.py
│       ├── histogram.py
│       ├── pca_plot.py
│       ├── radar_chart_plot.py
│       ├── site_specific_heatmap.py
│       ├── venn_diagram_plot.py
│       ├── vip_score_plot.py
│       ├── vip_score_plot_r.py
│       └── volcano_plot.py
│
├── Dataset/                   # Input data (user-provided)
│   ├── C_01.csv ... C_24.csv (Cancer samples)
│   └── N_01.csv ... N_24.csv (Normal samples)
│
└── Results/                   # Generated outputs
    ├── integrated.csv
    ├── vip_scores_all.csv
    ├── analysis_summary.txt
    ├── glycan_type_statistics.csv
    ├── *.png (visualizations)
    └── Trace/                 # Verification data
        ├── *_data.csv
        └── *_summary.csv
```

---

## Pipeline Architecture

### Execution Flow

```
main.py
  │
  ├─[1]─> Load Configuration (config.yaml)
  │
  ├─[2]─> Data Integration (DataLoader)
  │       ├─ Scan Dataset/ for CSV files
  │       ├─ Extract Peptide, GlycanComposition, IsotopeArea
  │       ├─ Create wide-format table
  │       └─ Output: integrated.csv
  │
  ├─[3]─> Annotation (GlycanAnnotator)
  │       ├─ Parse glycan compositions
  │       ├─ Classify into 5 categories (HM, F, S, SF, C/H)
  │       ├─ Add GlycanTypeCategory column
  │       └─ Output: annotated data
  │
  ├─[4]─> Statistical Analysis (GlycanAnalyzer)
  │       ├─ TIC normalization
  │       ├─ Log2 transformation
  │       ├─ PCA analysis
  │       ├─ PLS-DA → VIP scores
  │       └─ Output: vip_scores_all.csv, statistics
  │
  ├─[5]─> Visualization (GlycanVisualizer)
  │       ├─ PCA plots
  │       ├─ Boxplots
  │       ├─ Heatmaps
  │       ├─ Volcano plot
  │       ├─ Glycopeptide comparison heatmap ⭐
  │       └─ Output: *.png files + trace data
  │
  └─[6]─> Summary Report
          └─ Output: analysis_summary.txt
```

---

## Core Modules

### 1. data_loader.py (DataLoader)

**Purpose**: Integrate multiple CSV files into single wide-format table

**Key Methods**:
- `load_csv_files()` - Scan and load all CSV files
- `integrate_data()` - Pivot to wide format
- `save_integrated_data()` - Export integrated.csv

**Data Flow**:
```
C_01.csv   PEPTIDE | GLYCAN | AREA
C_02.csv   PEPTIDE | GLYCAN | AREA   →  PEPTIDE | GLYCAN | C1 | C2 | ... | N24
...        ...                           ...
N_24.csv   PEPTIDE | GLYCAN | AREA
```

### 2. annotator.py (GlycanAnnotator)

**Purpose**: Classify glycan compositions into functional categories

**Classification System**:
```
HM (High-Mannose):      H≥5, N=2, no A/F/G
F (Fucosylated):        Has F, no A
S (Sialylated):         Has A, no F
SF (Sialofucosylated):  Has A and F
C/H (Complex/Hybrid):   Everything else
```

**Key Methods**:
- `parse_glycan_composition()` - Extract H, N, A, F, G counts
- `is_high_mannose()` - HM detection logic
- `get_glycan_type_category()` - 5-category classification

### 3. analyzer.py (GlycanAnalyzer)

**Purpose**: Statistical analysis and dimensionality reduction

**Key Methods**:
- `normalize_tic()` - Total Ion Current normalization
- `perform_pca()` - Principal Component Analysis
- `perform_plsda()` - PLS-DA with VIP scores
- `calculate_statistics()` - Group statistics

**Normalization Pipeline**:
```
Raw intensities
  → TIC normalization (median-based)
  → Log2(x+1) transformation
  → StandardScaler
  → PCA/PLS-DA
```

### 4. visualizer.py (GlycanVisualizer)

**Purpose**: Coordinate all visualization generation

**Design Pattern**: Mixin-based architecture
```python
class GlycanVisualizer(
    PCAPlotMixin,
    BoxplotMixin,
    HeatmapMixin,
    VolcanoPlotMixin,
    GlycopeptideComparisonHeatmapMixin,  # NEW
    ...
):
    pass
```

**Benefits**:
- Modular: Each visualization in separate file
- Maintainable: Easy to add/remove plots
- Clean: No monolithic class

---

## Visualization Architecture

### Mixin Pattern

Each visualization is a separate mixin class:

```python
# src/plots/pca_plot.py
class PCAPlotMixin:
    def plot_pca(self, ...):
        # PCA visualization logic
        pass

# src/visualizer.py
class GlycanVisualizer(PCAPlotMixin, BoxplotMixin, ...):
    def __init__(self, config, output_dir):
        self.config = config
        self.output_dir = output_dir
        self.dpi = config['visualization']['dpi']
```

### Key Visualizations

**1. PCA Plot** (`pca_plot.py`)
- Cancer vs Normal separation
- Variance explained
- Sample clustering

**2. Glycopeptide Comparison Heatmap** (`glycopeptide_comparison_heatmap.py`) ⭐
- **NEW**: Advanced comparison visualization
- × symbol = Cancer (red)
- \+ symbol = Normal (blue)
- Glycan type grouping (HM, F, S, SF, C/H)
- VIP score sorting
- Complete trace data export

**3. Volcano Plot** (`volcano_plot.py`)
- Differential expression
- Top 3 increases/decreases
- Glycan type coloring
- Significance thresholds

**4. Heatmaps** (`heatmap.py`)
- Clustered heatmap (top 50)
- Full glycan profile
- Hierarchical clustering

---

## Data Flow

### Input → Processing → Output

```
[Dataset/]
  C_01.csv, ..., N_24.csv
        ↓
[DataLoader] → integrated.csv
        ↓
[GlycanAnnotator] → annotated data (in memory)
        ↓
[GlycanAnalyzer] → vip_scores_all.csv, statistics.csv
        ↓
[GlycanVisualizer] → *.png files
        ↓              ↓
        ↓         [Trace Data]
        ↓              ↓
        ↓         Results/Trace/*.csv
        ↓
[Summary Report] → analysis_summary.txt
```

---

## Configuration System

### config.yaml Structure

```yaml
paths:
  dataset_dir: "Dataset"
  output_dir: "Results"

processing:
  required_columns: ["Peptide", "GlycanComposition", "IsotopeArea"]

annotation:
  sialylation_marker: "A"
  fucosylation_marker: "F"

analysis:
  pca:
    n_components: 2
    log_transform: true
  normalization:
    method: "tic"

visualization:
  figsize: [12, 8]
  dpi: 300
  glycan_type_colors:
    HM: "#00CC00"
    F: "#FF0000"
    S: "#FF69B4"
    SF: "#FFA500"
    "C/H": "#0000FF"

glycopeptide_comparison:
  enabled: true
  max_peptides: 20
  max_glycans_per_type: 15
```

---

## Design Patterns

### 1. Mixin Pattern (Visualizations)
- **Why**: Separate concerns, modular design
- **Where**: `src/visualizer.py` + `src/plots/*.py`

### 2. Factory Pattern (Plot Creation)
- **Why**: Consistent interface for all plots
- **Where**: Each plot mixin has standard methods

### 3. Pipeline Pattern (Main Flow)
- **Why**: Sequential processing steps
- **Where**: `main.py` execution flow

### 4. Configuration Pattern
- **Why**: Centralized settings, no hardcoding
- **Where**: `config.yaml` loaded by all modules

---

## Trace Data System

### Purpose
Enable manual verification of all visualizations

### Structure
```
Results/Trace/
├── pca_plot_data.csv
├── heatmap_top_glycopeptides_data.csv
├── volcano_plot_data.csv
├── glycopeptide_comparison_heatmap_summary.csv  ⭐
└── glycopeptide_comparison_heatmap_data.csv     ⭐
```

### Comparison Heatmap Trace Data

**Summary file** (22 columns):
- Position info: Plot_X_Position, Plot_Y_Position
- Statistics: Means, StdDev, Fold_Change, etc.
- Plot properties: Cancer_Alpha, Normal_Alpha

**Full data file** (69 columns):
- All summary columns
- Individual samples: C1-C24, N1-N24
- Raw intensity values

---

## Testing Strategy

### Verification Scripts

**scripts/verify_trace_data.py**:
- Checks all calculations
- Validates sorting
- Verifies grouping
- 8 comprehensive checks

**Run**:
```bash
python3 scripts/verify_trace_data.py
```

---

## Extension Points

### Adding New Visualizations

1. Create new file: `src/plots/my_plot.py`
2. Define mixin class:
   ```python
   class MyPlotMixin:
       def plot_my_visualization(self, df, ...):
           # Implementation
           pass
   ```
3. Add to `src/visualizer.py`:
   ```python
   class GlycanVisualizer(MyPlotMixin, ...):
       pass
   ```
4. Call in `main.py`:
   ```python
   visualizer.plot_my_visualization(data)
   ```

### Adding New Glycan Types

1. Update `src/annotator.py`:
   - Add detection logic in `get_glycan_type_category()`
2. Update `config.yaml`:
   - Add color for new type
3. Update visualizations as needed

---

## Dependencies

### Python Packages

**Core**:
- pandas - Data manipulation
- numpy - Numerical operations
- matplotlib - Plotting
- seaborn - Statistical visualizations

**Analysis**:
- scikit-learn - PCA, PLS-DA, scaling
- scipy - Statistical tests
- statsmodels - FDR correction

**R Integration**:
- rpy2 - R interface for ggplot2-based VIP plots

See `requirements.txt` for complete list.

---

## Performance Considerations

### Large Datasets
- TIC normalization: O(n*m) where n=samples, m=features
- PLS-DA: O(n*m*k) where k=components
- Trace data: Complete data export for verification

### Memory
- Wide-format table kept in memory
- Visualizations generated sequentially
- Results saved incrementally

---

## Error Handling

### Input Validation
- CSV file existence check
- Required column verification
- Data type validation

### Graceful Degradation
- Skip missing samples (e.g., N19)
- Handle empty cells (fillna with 0)
- Continue pipeline if individual plots fail

---

## Future Enhancements

### Potential Additions
1. Interactive visualizations (Plotly)
2. Web dashboard (Streamlit/Dash)
3. Batch processing multiple datasets
4. Machine learning classification
5. Database backend for large datasets

---

**Version**: 2.0
**Last Updated**: 2025-10-05
**Status**: Production-ready
