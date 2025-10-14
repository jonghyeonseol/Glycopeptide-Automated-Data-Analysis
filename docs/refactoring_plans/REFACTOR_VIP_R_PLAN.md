# Phase 6: vip_score_plot_r.py Refactoring Plan

## Current Status
- **File**: `src/plots/vip_score_plot_r.py`
- **Lines**: 601 lines
- **Methods**: 4 visualization methods (3 with significant duplication + 1 unique complex method)

## Duplication Analysis

### Method 1: `plot_vip_scores_glycopeptide_r()` (lines 232-287, ~56 lines)
### Method 2: `plot_vip_scores_glycan_composition_r()` (lines 289-339, ~51 lines)
### Method 3: `plot_vip_scores_peptide_r()` (lines 341-391, ~51 lines)
### Method 4: `plot_vip_scores_peptide_grouped_r()` (lines 393-601, ~209 lines) - UNIQUE, NO DUPLICATION

### Identified Duplications (Methods 1-3 Only)

#### 1. Sample Column Extraction (1 line × 3 = 3 lines)
**Lines 255, 304, 356**: Identical
```python
cancer_samples, normal_samples = get_sample_columns(df)
```

#### 2. Config Initialization (1 line × 3 = 3 lines)
**Lines 258, 307, 359**: Identical
```python
config = DataPreparationConfig(missing_data_method='skipna')
```

#### 3. Heatmap Data Preparation Loop (~20 lines × 3 = ~60 lines)
**Lines 259-280, 308-331, 360-383**: Similar structure with minor differences

**Method 1** (Glycopeptide - uses MEAN):
```python
heatmap_data = []
for _, row in top_n_data.iterrows():
    mask = (df['Peptide'] == row['Peptide']) & (df['GlycanComposition'] == row['GlycanComposition'])
    if mask.sum() > 0:
        glycopeptide_row = df[mask]
        cancer_stats = calculate_group_statistics_standardized(
            glycopeptide_row, cancer_samples, method=config.missing_data_method
        )
        normal_stats = calculate_group_statistics_standardized(
            glycopeptide_row, normal_samples, method=config.missing_data_method
        )
        cancer_mean = cancer_stats['mean'].iloc[0] if not cancer_stats['mean'].isna().all() else 0
        normal_mean = normal_stats['mean'].iloc[0] if not normal_stats['mean'].isna().all() else 0
        heatmap_data.append({'Feature': row['Feature'], 'Cancer': cancer_mean, 'Normal': normal_mean})
    else:
        heatmap_data.append({'Feature': row['Feature'], 'Cancer': 0, 'Normal': 0})
```

**Method 2** (Glycan Composition - uses SUM):
```python
heatmap_data = []
for _, row in glycan_vip.iterrows():
    glycan_comp = row['GlycanComposition']
    mask = df['GlycanComposition'] == glycan_comp
    if mask.sum() > 0:
        glycan_rows = df[mask]
        cancer_stats = calculate_group_statistics_standardized(
            glycan_rows, cancer_samples, method=config.missing_data_method
        )
        normal_stats = calculate_group_statistics_standardized(
            glycan_rows, normal_samples, method=config.missing_data_method
        )
        cancer_total = cancer_stats['sum'].sum()  # SUM across peptides
        normal_total = normal_stats['sum'].sum()
        heatmap_data.append({'Feature': row['Feature'], 'Cancer': cancer_total, 'Normal': normal_total})
    else:
        heatmap_data.append({'Feature': row['Feature'], 'Cancer': 0, 'Normal': 0})
```

**Method 3** (Peptide - uses SUM):
```python
heatmap_data = []
for _, row in peptide_vip.iterrows():
    peptide = row['Peptide']
    mask = df['Peptide'] == peptide
    if mask.sum() > 0:
        peptide_rows = df[mask]
        cancer_stats = calculate_group_statistics_standardized(
            peptide_rows, cancer_samples, method=config.missing_data_method
        )
        normal_stats = calculate_group_statistics_standardized(
            peptide_rows, normal_samples, method=config.missing_data_method
        )
        cancer_total = cancer_stats['sum'].sum()  # SUM across glycoforms
        normal_total = normal_stats['sum'].sum()
        heatmap_data.append({'Feature': row['Feature'], 'Cancer': cancer_total, 'Normal': normal_total})
    else:
        heatmap_data.append({'Feature': row['Feature'], 'Cancer': 0, 'Normal': 0})
```

**Differences**:
- Mask creation logic (filter by peptide+glycan vs glycan vs peptide)
- Aggregation method (mean for single glycopeptide, sum for grouped features)

#### 4. Output File and Plot Creation (3-4 lines × 3 = ~10 lines)
**Lines 284-287, 336-339, 388-391**: Similar pattern
```python
output_file = self.output_dir / 'vip_score_*.png'
self._create_vip_plot_r(vip_plot_data, heatmap_df,
                        f'Top {top_n} ... by VIP Score',
                        'Label', str(output_file))
```

### Total Duplication
**~76 lines of duplication** (excluding the unique grouped method)

## Refactoring Strategy

### Helper Method 1: `_prepare_vip_heatmap_data_generic()`
**Purpose**: Unified heatmap data preparation for VIP score plots with flexible aggregation

**Parameters**:
- `df`: Annotated DataFrame
- `features_df`: DataFrame with features and VIP scores
- `feature_col`: Column name for feature identification
- `cancer_samples`: List of cancer sample columns
- `normal_samples`: List of normal sample columns
- `config`: DataPreparationConfig
- `mask_fn`: Callable that takes (df, row) and returns boolean mask
- `aggregation_fn`: Callable that takes (stats_dict, row) and returns aggregated value

**Returns**: Tuple of (heatmap_df, vip_plot_data)

**Lines**: ~30 lines

**Replaces**: Lines 255-282, 304-334, 356-386 (~85 lines total)

### Design Pattern: Strategy Pattern
The helper method uses two callable strategies:
1. **mask_fn**: Defines how to filter data for each feature
   - Glycopeptide: `lambda df, row: (df['Peptide'] == row['Peptide']) & (df['GlycanComposition'] == row['GlycanComposition'])`
   - Glycan: `lambda df, row: df['GlycanComposition'] == row['GlycanComposition']`
   - Peptide: `lambda df, row: df['Peptide'] == row['Peptide']`

2. **aggregation_fn**: Defines how to aggregate statistics
   - Mean: `lambda stats, row: stats['mean'].iloc[0] if not stats['mean'].isna().all() else 0`
   - Sum: `lambda stats, row: stats['sum'].sum()`

## Expected Results

### Before Refactoring
- **Total lines**: 601
- **Duplication**: ~76 lines across methods 1-3
- **Method 1**: ~56 lines
- **Method 2**: ~51 lines
- **Method 3**: ~51 lines
- **Method 4**: ~209 lines (no changes)

### After Refactoring
- **New helper method**:
  - `_prepare_vip_heatmap_data_generic()`: ~30 lines
  - **Total new code**: ~30 lines

- **Method 1 reduction**: 56 → ~30 lines (26 lines eliminated)
- **Method 2 reduction**: 51 → ~25 lines (26 lines eliminated)
- **Method 3 reduction**: 51 → ~25 lines (26 lines eliminated)
- **Method 4**: ~209 lines (unchanged)

- **Expected total**: 601 - 76 + 30 = **~555 lines**
- **Net reduction**: **~46 lines (7.7%)**
- **Code quality improvement**: Significant (eliminates all heatmap prep duplication)

## Implementation Steps

1. ✅ Create refactoring plan document
2. ⏳ Implement `_prepare_vip_heatmap_data_generic()` helper method
3. ⏳ Refactor `plot_vip_scores_glycopeptide_r()` to use helper
4. ⏳ Refactor `plot_vip_scores_glycan_composition_r()` to use helper
5. ⏳ Refactor `plot_vip_scores_peptide_r()` to use helper
6. ⏳ Run syntax check
7. ⏳ Run full pipeline test
8. ⏳ Validate outputs (check PNG files exist and are valid)
9. ⏳ Update CLAUDE.md with Phase 6 completion

## Validation Criteria

✅ **Syntax Check**: Python compile successful
✅ **Pipeline Test**: Executes without errors
✅ **Output Files**: All VIP score PNG files generated
✅ **Data Integrity**: Same VIP scores and heatmap values in plots
✅ **No Breaking Changes**: Public API unchanged
