# pGlyco Auto Combine

Automated glycoproteomics data integration and analysis pipeline for cancer vs normal tissue comparison.

## Features

1. **Data Integration**: Automatically integrates all CSV files from the Dataset folder
2. **Glycan Annotation**: 5-category glycan type classification (HM, F, S, SF, C/H)
3. **Statistical Analysis**: PLS-DA, PCA, and statistical comparison between groups
4. **Advanced Visualizations**:
   - PCA plots with group separation
   - Glycopeptide comparison heatmaps (Cancer vs Normal)
   - Volcano plots with significance annotations
   - Boxplots and distribution analysis
5. **Complete Traceability**: All visualizations include trace data for manual verification

## 📚 Documentation

- **[Architecture](ARCHITECTURE.md)** - System design and code organization
- **[User Guides](docs/)** - Comprehensive documentation index
- **[Visualization Guide](docs/visualization-guide.md)** - All plot types explained
- **[Verification Guide](docs/verification-guide.md)** - Manual data verification
- **[Change Log](docs/CHANGELOG.md)** - Recent updates and improvements

## Project Structure

```
pGlyco_auto_combine/
├── README.md                # This file
├── CLAUDE.md                # AI assistant instructions
├── config.yaml              # Configuration file
├── main.py                  # Main pipeline execution script
├── requirements.txt         # Python package dependencies
├── docs/                    # Documentation
│   ├── README.md           # Documentation index
│   ├── visualization-guide.md
│   ├── verification-guide.md
│   ├── trace-data-reference.md
│   └── normalization.md
├── src/                     # Source code
│   ├── data_loader.py
│   ├── annotator.py
│   ├── analyzer.py
│   ├── visualizer.py
│   ├── plots/              # Visualization modules
│   └── utils.py
├── Dataset/                 # Input CSV files (C_01.csv ~ N_24.csv)
└── Results/                 # Output results
    ├── integrated.csv
    ├── glycan_type_statistics.csv
    ├── vip_scores_all.csv
    ├── *.png               # All visualizations
    └── Trace/              # Trace data for verification
        ├── *_data.csv
        └── *_summary.csv
```

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### 1. Data Preparation
- Place CSV files in the `Dataset/` folder (filename format: `C_01.csv`, `C_02.csv`, ..., `N_01.csv`, ...)
- Each CSV file must contain the following columns:
  - `Peptide`: Peptide sequence
  - `GlycanComposition`: Glycan composition (e.g., H(5)N(4)A(2)F(1))
  - `IsotopeArea`: Quantification value

### 2. Run Pipeline

```bash
python3 main.py
```

### 3. Check Results
- All output files will be generated in the `Results/` folder
- Trace data for verification: `Results/Trace/`
- See [Verification Guide](docs/verification-guide.md) for manual verification steps

## Annotation Rules

### Sialylation
- If `A` is present in GlycanComposition → **Sialylated**
- Example: `H(5)N(4)A(2)` → Sialylated (2 sialic acids)

### Fucosylation
- If `F` is present in GlycanComposition → **Fucosylated**
- Example: `H(5)N(4)F(1)` → Fucosylated (1 fucose)

### Glycan Type
- `Non`: No Sialylation, No Fucosylation
- `Sialylated`: Sialylation only
- `Fucosylated`: Fucosylation only
- `Both`: Both Sialylation and Fucosylation present

## Output File Descriptions

### integrated_example.csv
- File integrating data from all samples
- Structure:
  ```
  Peptide | GlycanComposition | C1 | C2 | ... | N1 | N2 | ... | Sialylation | Fucosylation | GlycanType
  ```

### analysis_summary.txt
- Summary of overall analysis results
- Number of samples, glycan type distribution, PCA results, statistics, etc.

### Visualization Files
- **pca_plot.png**: PCA distribution of Cancer vs Normal samples
- **boxplot_glycan_types.png**: Intensity distribution by glycan type
- **glycan_type_distribution.png**: Glycan type count distribution
- **heatmap_top_glycopeptides.png**: Heatmap of top 50 glycopeptides

## Configuration

You can adjust various parameters by modifying the `config.yaml` file:

- Path settings
- QC filtering criteria
- Annotation markers
- PCA parameters
- Visualization options

## Developer

Development Date: 2025-10-01

## License

MIT License
