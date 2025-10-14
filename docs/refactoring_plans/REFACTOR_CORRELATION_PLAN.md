# Phase 8: correlation_matrix_plot.py Refactoring Plan

## Current Status
- **File**: `src/plots/correlation_matrix_plot.py`
- **Lines**: 566 lines
- **Methods**: 7 methods (3 public wrappers + 2 private helpers + 2 combined methods)

## Duplication Analysis

### Identified Duplications

#### 1. TIC Normalization + Log2 Transform (11 lines × 5 = 55 lines)
**Lines 77-88, 200-210, 280-290, 390-400, 487-496**: IDENTICAL
```python
# Step 1: TIC (Total Ion Current) Normalization
sample_sums = intensity_data.sum(axis=0)
median_sum = sample_sums.median()
sample_sums_safe = sample_sums.replace(0, 1)
intensity_normalized = intensity_data / sample_sums_safe * median_sum

# Step 2: Log2 transform
intensity_log = np.log2(intensity_normalized + 1)
```

#### 2. Correlation Matrix Calculation (1 line × 5 = 5 lines)
**Lines 91, 213, 293, 402, 499**: IDENTICAL
```python
corr_matrix = intensity_log.corr(method='pearson')
```

#### 3. Dynamic Center Calculation (8 lines × 5 = 40 lines)
**Lines 96-106, 216-223, 295-303, 407-417, 501-509**: Nearly identical
```python
if CORR_CENTER_AUTO:
    corr_center = np.median(corr_matrix.values)
    logger.debug(f"  Using dynamic center: {corr_center:.3f}")
else:
    corr_center = CORR_CENTER_FIXED
    logger.debug(f"  Using fixed center: {corr_center:.3f}")
```

### Total Duplication
**~100 lines of duplication**

## Refactoring Strategy

### Helper Method: `_prepare_correlation_matrix()`
**Purpose**: Unified correlation matrix preparation pipeline
**Parameters**:
- `df`: Annotated DataFrame
- `samples`: List of sample columns
- `method`: Correlation method (default='pearson')

**Returns**: Tuple of (corr_matrix, intensity_log)

**Pipeline**:
1. Extract intensity data with `replace_empty_with_zero()`
2. Apply TIC normalization
3. Apply Log2 transform
4. Calculate correlation matrix
5. Return both correlation matrix and log-transformed data

**Lines**: ~20 lines

**Replaces**: Lines 77-91, 200-213, 280-293, 390-402, 487-499 (~60 lines total)

### Helper Method 2: `_get_correlation_center()`
**Purpose**: Calculate correlation center (dynamic or fixed)
**Parameters**:
- `corr_matrix`: Correlation matrix
- `default_center`: Optional default center for special cases

**Returns**: float (center value)

**Lines**: ~8 lines

**Replaces**: Lines 96-106, 216-223, 295-303, 407-417, 501-509 (~40 lines total)

## Expected Results

### Before Refactoring
- **Total lines**: 566
- **Duplication**: ~100 lines
- **Helper methods**: 2 existing (_plot_single_correlation_matrix, _plot_single_clustermap)

### After Refactoring
- **New helper methods**: 2 helpers (~28 lines total)
- **Net reduction**: ~72 lines eliminated from main methods
- **Expected total**: 566 - 72 = **~494 lines**
- **Net reduction**: **~72 lines (12.7%)**

## Implementation Steps

1. ✅ Create refactoring plan document
2. ⏳ Implement `_prepare_correlation_matrix()` helper
3. ⏳ Implement `_get_correlation_center()` helper
4. ⏳ Refactor `_plot_single_correlation_matrix()` to use helpers
5. ⏳ Refactor `_plot_single_clustermap()` to use helpers
6. ⏳ Refactor `plot_correlation_matrix_combined()` to use helpers
7. ⏳ Refactor `plot_correlation_cross_group()` to use helpers
8. ⏳ Refactor `plot_correlation_clustermap_combined()` to use helpers
9. ⏳ Run syntax check
10. ⏳ Run full pipeline test
11. ⏳ Validate outputs
12. ⏳ Update CLAUDE.md

## Validation Criteria

✅ **Syntax Check**: Python compile successful
✅ **Pipeline Test**: Executes without errors
✅ **Output Files**: All correlation PNG files generated
✅ **Data Integrity**: Identical correlation matrices and trace CSV files
✅ **No Breaking Changes**: Public API unchanged
