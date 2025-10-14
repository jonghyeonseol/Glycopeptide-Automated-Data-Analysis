# Enhanced Pie Chart Refactoring Plan

## Current State
- **File**: `src/plots/enhanced_pie_chart_plot.py`
- **Lines**: 896
- **Methods**: 4 total
  1. `plot_pie_chart_glycan_types()` (lines 69-317, 248 lines)
  2. `plot_pie_chart_primary_classification()` (lines 319-479, 160 lines)
  3. `plot_pie_chart_secondary_classification()` (lines 481-635, 154 lines)
  4. `plot_pie_chart_significant_glycan_types()` (lines 637-896, 259 lines) ← KEEP SEPARATE

## Duplication Analysis

### Methods 1-3 Share Identical Structure:
1. **Data Aggregation** (~30 lines each)
   - Get sample columns
   - Loop through categories
   - Calculate cancer/normal totals

2. **Figure Setup** (~10 lines each)
   - Create 2x2 GridSpec
   - 3 subplots (Cancer pie, Normal pie, Fold change bars)

3. **Pie Charts** (~30 lines each)
   - Cancer pie chart with percentages
   - Normal pie chart with percentages

4. **Fold Change Analysis** (~60-110 lines)
   - Calculate fold changes (Cancer/Normal)
   - Color by enrichment direction
   - Statistical testing (only glycan_types)

5. **Bar Chart** (~30 lines each)
   - Plot fold change bars
   - Add value labels
   - Reference line at FC=1.0

6. **Styling & Saving** (~30 lines each)
   - Titles, labels, legend
   - Publication theme
   - Save figure and trace data

**Total Duplicated**: ~190-240 lines × 3 = 570-720 lines

## Base Method Design

```python
def _plot_pie_chart_comparative_base(
    self,
    df: pd.DataFrame,
    categories: list,
    column_name: str,
    colors_dict: dict,
    chart_type: str,
    title_main: str,
    title_cancer: str,
    title_normal: str,
    title_fc: str,
    output_filename: str,
    perform_statistical_test: bool = False,
    rotation: int = 0,
    figsize: tuple = (16, 12)
) -> None:
    """
    Unified base method for comparative pie charts with fold change

    Creates:
    - Cancer pie chart (ax1, top-left)
    - Normal pie chart (ax2, top-right)
    - Fold change bar chart (ax3, bottom span)

    Args:
        df: Annotated DataFrame with intensity data
        categories: List of category names
        column_name: DataFrame column for filtering
        colors_dict: Category -> color mapping
        chart_type: Type identifier for logging
        title_main: Overall figure title
        title_cancer: Cancer pie title
        title_normal: Normal pie title
        title_fc: Fold change title
        output_filename: Output file name (without .png)
        perform_statistical_test: Enable Mann-Whitney U + FDR correction
        rotation: X-axis label rotation (degrees)
        figsize: Figure size
    """
```

## Wrapper Implementations

### Method 1: glycan_types
```python
def plot_pie_chart_glycan_types(self, df: pd.DataFrame, figsize: tuple = (16, 12)):
    """Create enhanced glycan type pie charts with fold change"""
    logger.info("Creating enhanced glycan type pie charts with fold change...")

    return self._plot_pie_chart_comparative_base(
        df=df,
        categories=['Non', 'Sialylated', 'Fucosylated', 'Both'],
        column_name='GlycanType',
        colors_dict=LEGACY_GLYCAN_COLORS,
        chart_type='glycan_types',
        title_main='Glycan Type Distribution & Comparative Analysis',
        title_cancer='Cancer Group\nGlycan Type Distribution',
        title_normal='Normal Group\nGlycan Type Distribution',
        title_fc='Fold Change Analysis (Cancer vs Normal)',
        output_filename='pie_chart_glycan_types_enhanced',
        perform_statistical_test=True,  # ← Only this method has stats
        rotation=0,
        figsize=figsize
    )
```

### Method 2: primary_classification
```python
def plot_pie_chart_primary_classification(self, df: pd.DataFrame, figsize: tuple = (16, 12)):
    """Create enhanced primary classification pie charts"""
    logger.info("Creating enhanced primary classification pie charts with fold change...")

    colors_primary = {
        'Truncated': '#CCCCCC',
        'High Mannose': '#2ECC71',
        'ComplexHybrid': '#3498DB',
        'Outlier': '#95A5A6'
    }

    return self._plot_pie_chart_comparative_base(
        df=df,
        categories=['Truncated', 'High Mannose', 'ComplexHybrid', 'Outlier'],
        column_name='PrimaryClassification',
        colors_dict=colors_primary,
        chart_type='primary_classification',
        title_main='Primary Classification Distribution & Comparative Analysis',
        title_cancer='Cancer Group\nPrimary Classification',
        title_normal='Normal Group\nPrimary Classification',
        title_fc='Fold Change Analysis (Cancer vs Normal)',
        output_filename='pie_chart_primary_classification_enhanced',
        perform_statistical_test=False,
        rotation=15,
        figsize=figsize
    )
```

### Method 3: secondary_classification
```python
def plot_pie_chart_secondary_classification(self, df: pd.DataFrame, figsize: tuple = (18, 12)):
    """Create enhanced secondary classification pie charts"""
    logger.info("Creating enhanced secondary classification pie charts with fold change...")

    return self._plot_pie_chart_comparative_base(
        df=df,
        categories=['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated'],
        column_name='SecondaryClassification',
        colors_dict=EXTENDED_CATEGORY_COLORS,
        chart_type='secondary_classification',
        title_main='Secondary Classification Distribution & Comparative Analysis',
        title_cancer='Cancer Group\nSecondary Classification',
        title_normal='Normal Group\nSecondary Classification',
        title_fc='Fold Change Analysis (Cancer vs Normal)',
        output_filename='pie_chart_secondary_classification_enhanced',
        perform_statistical_test=False,
        rotation=20,
        figsize=figsize
    )
```

## Expected Results

**Before Refactoring**:
- 896 lines total
- 562 lines for methods 1-3 (with duplication)

**After Refactoring**:
- Base method: ~250 lines
- 3 thin wrappers: ~20 lines each = 60 lines
- Method 4 unchanged: 259 lines
- **Total**: ~569 lines (36% reduction from 896)
- **Lines Eliminated**: ~327 lines

## Implementation Steps

1. ✅ Analyze structure
2. ✅ Design base method
3. ⏳ Create `_plot_pie_chart_comparative_base()` at line 69
4. ⏳ Replace method 1 with thin wrapper
5. ⏳ Replace method 2 with thin wrapper
6. ⏳ Replace method 3 with thin wrapper
7. ⏳ Test all three outputs (save baselines first)
8. ⏳ Validate data integrity (MD5 checksums)
9. ⏳ Update CLAUDE.md with new file sizes
