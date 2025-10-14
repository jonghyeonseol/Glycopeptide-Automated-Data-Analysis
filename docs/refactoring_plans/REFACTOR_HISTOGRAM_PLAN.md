# Phase 7: histogram.py Refactoring Plan

## Current Status
- **File**: `src/plots/histogram.py`
- **Lines**: 535 lines
- **Methods**: 5 visualization methods with significant duplication

## Duplication Analysis

### Identified Duplications

#### 1. Sample Intensity Normalization (13 lines × 2 = 26 lines)
**Lines 172-184 and 281-293**: IDENTICAL
```python
if normalization == 'raw':
    intensity_col = replace_empty_with_zero(df[sample])
    min_val = intensity_col.min()
    max_val = intensity_col.max()
    if max_val > min_val:
        intensity_col = (intensity_col - min_val) / (max_val - min_val)
    else:
        intensity_col = intensity_col * 0
else:
    intensity_col = replace_empty_with_zero(df[sample])
```

#### 2. Sample-wise Aggregation Loop (~21 lines × 2 = ~42 lines)
**Lines 169-190 and 278-299**: Nearly identical (only category lists differ)
```python
for sample in sample_cols:
    sample_data = {'Sample': sample}
    [normalization logic]
    for category in categories:
        mask = df[classification_col] == category
        sample_data[category] = intensity_col[mask].sum()
    data_for_plot.append(sample_data)
```

#### 3. Proportional Normalization (~8 lines × 4 = ~32 lines)
**Lines 197-204, 306-313, 404-410, 488-494**: Nearly identical
```python
if normalization == 'aggregated':  # or unconditional in cancer_vs_normal methods
    row_sums = plot_df[categories].sum(axis=1)
    for col in categories:
        plot_df[col] = plot_df[col] / row_sums
    plot_df = plot_df.fillna(0)
```

#### 4. Group Aggregation Loop (~11 lines × 2 = ~22 lines)
**Lines 388-398 and 472-482**: Nearly identical
```python
for group_name, group_samples in [('Cancer', cancer_samples), ('Normal', normal_samples)]:
    group_data = {'Group': group_name}
    intensity_matrix = replace_empty_with_zero(df[group_samples])
    for category in categories:
        mask = df[classification_col] == category
        group_data[category] = intensity_matrix[mask].sum().sum()
    data_for_plot.append(group_data)
```

### Total Duplication
**~122 lines of duplication**

## Refactoring Strategy

### Helper Method 1: `_normalize_intensity_column()`
**Purpose**: Normalize a single intensity column
**Parameters**: `intensity_col`, `normalization` ('raw' or other)
**Returns**: Normalized Series
**Lines**: ~12 lines
**Replaces**: Lines 172-184, 281-293 (~26 lines)

### Helper Method 2: `_aggregate_by_classification_samplewise()`
**Purpose**: Aggregate intensities by classification across samples
**Parameters**: `df`, `sample_cols`, `categories`, `classification_col`, `normalization`
**Returns**: DataFrame with aggregated data
**Lines**: ~25 lines
**Replaces**: Lines 169-204, 278-313 (~90 lines)

### Helper Method 3: `_apply_proportional_normalization()`
**Purpose**: Apply within-row proportional normalization
**Parameters**: `plot_df`, `categories`
**Returns**: Normalized DataFrame
**Lines**: ~5 lines
**Replaces**: Lines 197-204, 306-313, 404-410, 488-494 (~32 lines)

### Helper Method 4: `_aggregate_by_classification_groupwise()`
**Purpose**: Aggregate intensities by classification for cancer vs normal groups
**Parameters**: `df`, `cancer_samples`, `normal_samples`, `categories`, `classification_col`
**Returns**: DataFrame with aggregated data
**Lines**: ~15 lines
**Replaces**: Lines 388-410, 472-494 (~44 lines)

## Expected Results

### Before Refactoring
- **Total lines**: 535
- **Duplication**: ~122 lines
- **Method 2**: ~110 lines
- **Method 3**: ~110 lines
- **Method 4**: ~83 lines
- **Method 5**: ~83 lines

### After Refactoring
- **New helper methods**: ~57 lines
- **Method 2 reduction**: 110 → ~70 lines (40 lines eliminated)
- **Method 3 reduction**: 110 → ~70 lines (40 lines eliminated)
- **Method 4 reduction**: 83 → ~60 lines (23 lines eliminated)
- **Method 5 reduction**: 83 → ~60 lines (23 lines eliminated)

- **Expected total**: 535 - 122 + 57 = **~470 lines**
- **Net reduction**: **~65 lines (12.1%)**
- **Code quality improvement**: Significant (eliminates all major duplication)

## Implementation Steps

1. ✅ Create refactoring plan document
2. ⏳ Implement helper methods
3. ⏳ Refactor classification methods (2 methods)
4. ⏳ Refactor cancer_vs_normal methods (2 methods)
5. ⏳ Run syntax check
6. ⏳ Run full pipeline test
7. ⏳ Validate outputs
8. ⏳ Update CLAUDE.md

## Validation Criteria

✅ **Syntax Check**: Python compile successful
✅ **Pipeline Test**: Executes without errors
✅ **Output Files**: All histogram PNG files generated
✅ **Data Integrity**: Identical data in trace CSV files
✅ **No Breaking Changes**: Public API unchanged
