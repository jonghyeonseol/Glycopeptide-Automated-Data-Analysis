# Data Normalization Pipeline

## Overview

This document describes the standardized data normalization pipeline applied consistently across all visualizations in the pGlyco Auto Combine project.

## Standard Pipeline

All sample-to-sample comparison visualizations follow this pipeline:

```
Raw Intensity → TIC Normalization → Log2 Transform → Visualization
```

### Step 1: TIC (Total Ion Current) Normalization

**Purpose**: Remove sample loading differences and total protein amount variations

**Method**:
```python
sample_sums = intensity_matrix.sum(axis=0)
median_sum = sample_sums.median()
intensity_normalized = intensity_matrix / sample_sums * median_sum
```

**Why**: Sample C1 may have 2x more total protein than Sample C2. Without TIC normalization, all glycopeptides in C1 would appear 2x higher, creating false biological differences.

### Step 2: Log2 Transformation

**Purpose**:
- Stabilize variance across intensity ranges
- Make fold changes symmetric (2-fold up = -2-fold down)
- Reduce impact of extreme outliers

**Method**:
```python
intensity_log = np.log2(intensity_normalized + 1)
```

**Pseudocount (+1)**: Prevents log(0) = -∞

### Step 3 (PCA only): RobustScaler

**Purpose**: Scale features to comparable ranges using median/IQR (robust to outliers)

**Method**:
```python
from sklearn.preprocessing import RobustScaler
scaler = RobustScaler()
intensity_scaled = scaler.fit_transform(intensity_log)
```

## Visualization-Specific Pipelines

### Sample Comparison (TIC normalization required)

| Visualization | Pipeline | Notes |
|--------------|----------|-------|
| **PCA** | TIC → Log2 → RobustScaler → PCA | Feature scaling essential |
| **Heatmap** | TIC → Log2 → Hierarchical Clustering | Top 50 or full profile |
| **Correlation Matrix** | TIC → Log2 → Pearson Correlation | Sample-to-sample |
| **Histogram** | TIC → Sum by Category | Stacked bar plot |
| **Boxplot (Cancer vs Normal)** | TIC → Non-zero mean → Log2 | QC: detection ≥10% |

### Statistical Testing (Raw intensity preferred)

| Visualization | Pipeline | Rationale |
|--------------|----------|-----------|
| **Volcano Plot** | Raw → Non-zero only → Mann-Whitney U | Preserve biological magnitude |
| **CV Distribution** | Raw → Non-zero only → CV calculation | Measure technical variation |
| **Old Boxplot** | Raw → Log2 → Remove zeros | Legacy glycan type comparison |

### Key Differences

**Missing Value (0) Handling**:
- **Sample comparison**: Include zeros (0 = not detected)
- **Statistical tests**: Exclude zeros (detection-based analysis)

**Normalization**:
- **Sample comparison**: TIC normalization (removes total amount differences)
- **Statistical tests**: Raw intensity (preserves biological differences)

## Code Implementation

### Heatmap (src/plots/heatmap.py)

```python
def plot_heatmap(self, df: pd.DataFrame, figsize: tuple = (16, 12), top_n: int = 50):
    """
    Pipeline: TIC Normalization → Log2 Transform → Hierarchical Clustering
    """
    # Step 1: TIC Normalization
    intensity_matrix = replace_empty_with_zero(df[sample_cols])
    sample_sums = intensity_matrix.sum(axis=0)
    median_sum = sample_sums.median()
    intensity_normalized = intensity_matrix / sample_sums.replace(0, 1) * median_sum

    # Step 2: Log2 transform
    heatmap_data = np.log2(intensity_normalized + 1)

    # Step 3: Hierarchical clustering
    sns.clustermap(heatmap_data, ...)
```

### Correlation Matrix (src/plots/correlation_matrix_plot.py)

```python
def _plot_single_correlation_matrix(self, df: pd.DataFrame, samples: list, ...):
    """
    Pipeline: TIC Normalization → Log2 Transform → Pearson Correlation
    """
    # Step 1: TIC Normalization
    intensity_data = replace_empty_with_zero(df[samples])
    sample_sums = intensity_data.sum(axis=0)
    intensity_normalized = intensity_data / sample_sums.replace(0, 1) * median_sum

    # Step 2: Log2 transform
    intensity_log = np.log2(intensity_normalized + 1)

    # Step 3: Pearson correlation
    corr_matrix = intensity_log.corr(method='pearson')
```

### Boxplot Cancer vs Normal (src/plots/boxplot.py)

```python
def plot_boxplot_cancer_vs_normal_primary(self, df: pd.DataFrame, ...):
    """
    Pipeline: TIC Normalization → Non-zero mean → Log2 Transform
    """
    # Step 1: TIC Normalization (entire dataset)
    intensity_matrix = replace_empty_with_zero(df[sample_cols])
    sample_sums = intensity_matrix.sum(axis=0)
    intensity_normalized = intensity_matrix / sample_sums.replace(0, 1) * median_sum

    # Step 2: Calculate non-zero mean per sample
    for sample in sample_cols:
        values = subset_df[sample].values
        nonzero_values = values[values > 0]
        mean_intensity = nonzero_values.mean()

        # Step 3: Log2 transform the mean
        log_intensity = np.log2(mean_intensity + 1)
```

### PCA (src/analyzer.py)

```python
def perform_pca(self, df: pd.DataFrame) -> dict:
    """
    Pipeline: TIC Normalization → Log2 Transform → RobustScaler → PCA
    """
    intensity_matrix_t, sample_names, _ = self.prepare_intensity_matrix(df)

    # prepare_intensity_matrix already applies:
    # Step 1: TIC Normalization
    # Step 2: Log2 Transform

    # Step 3: RobustScaler
    self.scaler = RobustScaler()
    intensity_scaled = self.scaler.fit_transform(intensity_matrix_t)

    # Step 4: PCA
    self.pca = PCA(n_components=self.n_components)
    pca_coords = self.pca.fit_transform(intensity_scaled)
```

## Rationale

### Why TIC Normalization?

**Problem without TIC normalization**:
- Sample C1: 100 µg total protein → IsotopeArea = 1,000,000
- Sample C2: 50 µg total protein → IsotopeArea = 500,000
- **Same glycopeptide, different abundance due to loading only**

**Solution**:
- Normalize to median total intensity
- Now both samples show same relative glycopeptide levels

### Why Log2 Transform?

1. **Variance stabilization**: High-abundance glycopeptides have higher variance in raw space
2. **Symmetric fold changes**:
   - 2-fold increase: log2(2) = +1
   - 2-fold decrease: log2(0.5) = -1
3. **Outlier reduction**: Extreme values compressed

### Why Non-zero Only for Boxplots?

**Dataset characteristics**:
- ComplexHybrid: 79% zeros (non-detection)
- High Mannose: 74.7% zeros

**Problem with including zeros**:
- Q1 = 0, Median = 0, Q3 = 0 → Box height = 0 (invisible)
- All detected values appear as outliers

**Solution**:
- Calculate mean of detected values only
- Represents "typical detected intensity"
- Boxes are now visible and interpretable

## Quality Control (QC) in Boxplots

**Detection Rate Filter**: Exclude samples with <10% detection rate

**Example**:
- High Mannose in sample C22: 388 glycopeptides, 15 detected → 3.9% detection rate
- **Decision**: Exclude (unreliable mean due to sparse data)

## Consistency Verification

All sample comparison visualizations now use:
✅ Same TIC normalization method
✅ Same Log2 transformation
✅ Same median_sum target

This ensures:
- Heatmap clustering reflects true biological patterns
- Correlation matrix shows sample similarities correctly
- Boxplots compare Cancer vs Normal fairly
- PCA separates groups based on biology, not loading differences

## References

- TIC normalization: Removes technical variation in total sample amount
- RobustScaler: More resistant to outliers than StandardScaler
- Log2 transform: Standard in omics data analysis
- Non-zero filtering: Common in proteomics when dealing with missing values
