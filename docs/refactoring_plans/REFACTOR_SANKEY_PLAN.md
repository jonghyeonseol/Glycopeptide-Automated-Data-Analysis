# Phase 5: sankey_plot.py Refactoring Plan

## Current Status
- **File**: `src/plots/sankey_plot.py`
- **Lines**: 886 lines
- **Methods**: 2 main visualization methods with significant duplication

## Duplication Analysis

### Method 1: `plot_glycan_type_sankey()` (lines 189-498, ~310 lines)
### Method 2: `plot_group_to_glycan_sankey()` (lines 500-886, ~386 lines)

### Identified Duplications

#### 1. Config Initialization (6 lines × 2 = 12 lines)
**Lines 228-233 and 542-547**: Identical
```python
if config is None:
    config = DataPreparationConfig(
        min_detection_pct=0.30,
        min_samples=5,
        missing_data_method='skipna'
    )
```

#### 2. Data Preparation Pipeline (12 lines × 2 = 24 lines)
**Lines 236-247 and 550-561**: Nearly identical (only log_prefix differs)
```python
df_prepared = prepare_visualization_data(
    df=df, config=config, vip_scores=vip_scores,
    merge_method='left', apply_detection_filter=False,
    log_prefix="[Sankey]" / "[Group Sankey]"  # Only difference
)
if len(df_prepared) == 0:
    logger.error("No glycopeptides available!")
    return
```

#### 3. Sample Validation (3 lines × 2 = 6 lines)
**Lines 250-252 and 564-566**: Identical
```python
cancer_samples, normal_samples = self._get_sample_lists(df_prepared)
if not self._validate_samples(cancer_samples, normal_samples):
    return None
```

#### 4. Statistical Significance Calculation (8 lines × 2 = 16 lines)
**Lines 256-263 and 570-577**: Identical
```python
logger.info("Calculating statistical significance...")
df_with_stats = calculate_statistical_significance(
    df_prep=df_prepared,
    cancer_samples=cancer_samples,
    normal_samples=normal_samples,
    method='mannwhitneyu',
    fdr_correction=True
)
```

#### 5. Regulation Classification (5 lines × 2 = 10 lines)
**Lines 265-269 and 579-583**: Identical
```python
df_with_stats = self._calculate_regulation_status(
    df_with_stats,
    log2fc_threshold=log2fc_threshold,
    fdr_threshold=fdr_threshold
)
```

#### 6. Figure Saving Pattern (7 lines × 2 = 14 lines)
**Lines 471-477 and 848-854**: Similar pattern
```python
try:
    fig.write_image(str(output_png), width=W, height=H, scale=2)
    logger.info(f"Saved ... to {output_png}")
except Exception as e:
    logger.warning(f"Could not save PNG (requires kaleido): {e}")
fig.write_html(str(output_html))
logger.info(f"Saved ... to {output_html}")
```

### Total Duplication
**~82 lines of direct duplication**

## Refactoring Strategy

### Helper Method 1: `_prepare_sankey_data()`
**Purpose**: Unified data preparation pipeline for Sankey diagrams

**Parameters**:
- `df`: Input DataFrame
- `vip_scores`: VIP scores DataFrame
- `config`: DataPreparationConfig (optional)
- `log2fc_threshold`: Fold change threshold
- `fdr_threshold`: FDR significance threshold
- `log_prefix`: Prefix for logging messages

**Returns**: Tuple of (df_with_stats, cancer_samples, normal_samples) or (None, None, None) on failure

**Lines**: ~50 lines

**Replaces**: Lines 228-269 and 542-583 (~82 lines total)

### Helper Method 2: `_save_sankey_figure()`
**Purpose**: Save Sankey figure as PNG and HTML with error handling

**Parameters**:
- `fig`: Plotly Figure object
- `output_png`: Path object for PNG file
- `output_html`: Path object for HTML file
- `width`: Figure width in pixels
- `height`: Figure height in pixels
- `diagram_type`: Description for logging ("Sankey diagram", "Group → Glycan Type Sankey")

**Returns**: None

**Lines**: ~15 lines

**Replaces**: Lines 471-477 and 848-854 (~14 lines total)

## Expected Results

### Before Refactoring
- **Total lines**: 886
- **Duplication**: ~82 lines
- **Method 1**: ~310 lines
- **Method 2**: ~386 lines

### After Refactoring
- **New helper methods**:
  - `_prepare_sankey_data()`: ~50 lines
  - `_save_sankey_figure()`: ~15 lines
  - **Total new code**: ~65 lines

- **Method 1 reduction**: 310 → ~270 lines (40 lines eliminated)
- **Method 2 reduction**: 386 → ~346 lines (40 lines eliminated)

- **Expected total**: 886 - 82 + 65 = **~869 lines**
- **Net reduction**: **~17 lines (1.9%)**
- **Code quality improvement**: Significant (eliminates all major duplication)

## Implementation Steps

1. ✅ Create refactoring plan document
2. ⏳ Implement `_prepare_sankey_data()` helper method
3. ⏳ Implement `_save_sankey_figure()` helper method
4. ⏳ Refactor `plot_glycan_type_sankey()` to use helpers
5. ⏳ Refactor `plot_group_to_glycan_sankey()` to use helpers
6. ⏳ Run syntax check
7. ⏳ Run full pipeline test
8. ⏳ Validate outputs (check HTML/PNG files exist and are valid)
9. ⏳ Update CLAUDE.md with Phase 5 completion

## Validation Criteria

✅ **Syntax Check**: Python compile successful
✅ **Pipeline Test**: Executes without errors
✅ **Output Files**: Both HTML and PNG files generated
✅ **Data Integrity**: Same number of nodes/links in diagrams
✅ **No Breaking Changes**: Public API unchanged
