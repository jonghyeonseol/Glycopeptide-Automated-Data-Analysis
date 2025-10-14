# Phase 2 Refactoring Plan - glycopeptide_comparison_heatmap.py

**Goal**: Reduce file from 1,272 lines to ~600 lines by eliminating code duplication

## Current Structure Analysis

### Three Methods with Significant Duplication:
1. **plot_glycopeptide_comparison_heatmap** (lines 50-513, 464 lines)
   - Top N peptides/glycans view
   - 3-panel layout: line plot + colorbar + heatmap

2. **plot_glycopeptide_comparison_heatmap_full** (lines 514-915, 402 lines)
   - Full-scale view (ALL glycopeptides)
   - 3-panel layout: line plot + colorbar + heatmap
   - Dynamic figure sizing

3. **plot_glycopeptide_comparison_heatmap_by_type** (lines 916-1272, 357 lines)
   - Glycan-type-specific view
   - 2-panel layout: line plot + heatmap (no colorbar)
   - Single glycan type color

### Identified Duplication Patterns:

#### A. Helper Functions (Defined 3 Times - IDENTICAL)
- **glycan_sort_key(glycan_comp)** - Lines 167, 586, 990
  - Natural/numeric sorting of glycan compositions
  - Extracts monosaccharide counts via regex
  - ~14 lines each = 42 lines total → Can reduce to 14 lines (1 definition)

- **get_symbol_info(row)** - Lines 442, 853, 1214
  - Determines symbols (×, +) and linewidth based on group presence
  - Checks for qualitative differences (Cancer-only, Normal-only)
  - ~21 lines each = 63 lines total → Can reduce to 21 lines (1 definition)

#### B. Data Preparation Logic (Similar ~50 lines each = 150 lines)
- Config initialization with defaults
- Glycan type order/colors definition
- Data preparation via `prepare_visualization_data()`
- Empty data checks
- VIP score sorting

#### C. Ordering Logic (Similar ~30 lines each = 90 lines)
- Peptide ordering by VIP score (max per peptide, descending)
- Glycan ordering (grouped by type, naturally sorted via glycan_sort_key)
- Index dictionaries (peptide_to_idx, glycan_to_idx)
- Glycan type position tracking

#### D. Top Panel Plotting (Similar ~30 lines each = 90 lines)
- Average intensity line plots (Cancer vs Normal)
- X-axis limits and ticks
- Y-axis label
- Grid configuration
- Legend

#### E. Symbol Plotting Logic (Similar ~50 lines each = 150 lines)
- Iterate through df_plot
- Check group presence (has_cancer, has_normal)
- Determine linewidth (qualitative difference → bold)
- Plot × (Cancer) and + (Normal) symbols
- Apply glycan type colors

#### F. Trace Data Preparation (Similar ~90 lines each = 270 lines)
- Plot position mapping
- Calculate group statistics via `calculate_group_statistics_standardized()`
- Fold change calculations
- Symbol info extraction via `get_symbol_info()`
- Column reordering
- Save comprehensive and summary CSV files

#### G. Legend Creation (Similar ~30 lines each = 90 lines)
- Glycan type color patches
- Symbol indicators with linewidth variations
- Legend positioning and styling

**Total Duplication**: ~945 lines of duplicated/similar code

---

## Refactoring Strategy

### Phase 2.1: Extract Module-Level Helper Functions
**File**: Create `src/plots/heatmap_helpers.py` (new file)

```python
def glycan_sort_key(glycan_comp: str) -> tuple:
    """Natural/numeric sorting for glycan compositions"""
    # Lines 167-180 (extract once, reuse 3 times)

def get_symbol_info(row, linewidth_bold: float = 3.0, linewidth_normal: float = 2.5) -> tuple:
    """Determine symbol markers and linewidth based on group presence"""
    # Lines 442-462 (extract once, reuse 3 times)

def create_peptide_order(df: pd.DataFrame, vip_col: str = 'VIP_Score') -> list:
    """Create peptide order sorted by VIP score"""
    # Extract common peptide ordering logic

def create_glycan_order(df: pd.DataFrame, glycan_types: list, include_positions: bool = True) -> tuple:
    """Create glycan order grouped by type with natural sorting"""
    # Extract common glycan ordering logic
    # Returns (glycan_order, glycan_type_positions)
```

**Impact**: Reduces ~150 lines of duplication

---

### Phase 2.2: Extract Class Helper Methods
**File**: Modify `src/plots/glycopeptide_comparison_heatmap.py`

Add private helper methods to `GlycopeptideComparisonHeatmapMixin`:

```python
def _prepare_heatmap_data(self, df, vip_scores, config, glycan_types, filter_config) -> pd.DataFrame:
    """Standardized data preparation for all heatmap methods"""
    # Common: config defaults, data prep, ordering
    # Returns df_plot ready for visualization

def _plot_top_panel(self, ax, glycan_order, df, marker_size='large'):
    """Plot average intensity comparison panel (Cancer vs Normal)"""
    # Lines 208-239, 632-663, 1024-1055 (extract once, reuse 3 times)

def _plot_colorbar_panel(self, ax, glycan_order, glycan_type_positions, glycan_type_colors):
    """Plot glycan type color bar panel"""
    # Lines 241-266, 665-688 (extract once, reuse 2 times)

def _plot_symbol_heatmap(self, ax, df, peptide_to_idx, glycan_to_idx,
                         glycan_type_colors, marker_size='large', single_color=None):
    """Plot symbol-based heatmap (× Cancer, + Normal)"""
    # Lines 268-350, 690-771, 1057-1131 (extract once, reuse 3 times)
    # Parameter single_color: If provided, use single color instead of glycan_type_colors

def _create_legend(self, glycan_types, glycan_type_positions, glycan_type_colors,
                   marker_size='large', single_type=None):
    """Create legend for heatmap"""
    # Lines 352-385, 773-806, 1133-1159 (extract once, reuse 3 times)

def _prepare_and_save_trace_data(self, df, peptide_to_idx, glycan_to_idx,
                                  filename_base: str, include_type_info: bool = True):
    """Prepare and save comprehensive trace data"""
    # Lines 397-502, 818-903, 1179-1265 (extract once, reuse 3 times)
```

**Impact**: Reduces ~600 lines of duplication

---

### Phase 2.3: Create Unified Base Method
**File**: Modify `src/plots/glycopeptide_comparison_heatmap.py`

Add private base method:

```python
def _plot_comparison_heatmap_base(self, df, vip_scores, config,
                                  filter_mode: str,  # 'top_n', 'full', or 'by_type'
                                  max_peptides: int = None,
                                  max_glycans_per_type: int = None,
                                  glycan_type_filter: str = None,
                                  figsize: tuple = None):
    """
    Unified base method for all comparison heatmap variants

    filter_mode options:
    - 'top_n': Top N peptides/glycans (original method)
    - 'full': All glycopeptides (full-scale method)
    - 'by_type': Single glycan type (by-type method)
    """
    # Use helper methods to build heatmap
    # Conditional logic for layout (2-panel vs 3-panel)
    # Conditional logic for figure sizing (fixed vs dynamic)
```

**Impact**: Reduces ~200 lines by unifying control flow

---

### Phase 2.4: Refactor Public Methods
**File**: Modify `src/plots/glycopeptide_comparison_heatmap.py`

Simplify public methods to thin wrappers:

```python
def plot_glycopeptide_comparison_heatmap(self, df, vip_scores, config=None,
                                         figsize=(24, 16), max_peptides=50,
                                         max_glycans_per_type=15):
    """Top N peptides/glycans view"""
    return self._plot_comparison_heatmap_base(
        df, vip_scores, config,
        filter_mode='top_n',
        max_peptides=max_peptides,
        max_glycans_per_type=max_glycans_per_type,
        figsize=figsize
    )

def plot_glycopeptide_comparison_heatmap_full(self, df, vip_scores, config=None):
    """Full-scale view (ALL glycopeptides)"""
    return self._plot_comparison_heatmap_base(
        df, vip_scores, config,
        filter_mode='full',
        figsize=None  # Dynamic sizing
    )

def plot_glycopeptide_comparison_heatmap_by_type(self, df, vip_scores, glycan_type, config=None):
    """Glycan-type-specific view"""
    return self._plot_comparison_heatmap_base(
        df, vip_scores, config,
        filter_mode='by_type',
        glycan_type_filter=glycan_type,
        figsize=None  # Dynamic sizing
    )
```

**Impact**: Each public method becomes ~10 lines (down from 400+ lines)

---

## Expected Results

### Before Refactoring:
- **Total lines**: 1,272
- **Public methods**: 3 × ~400 lines = 1,200 lines
- **Code duplication**: ~945 lines (~74%)

### After Refactoring:
- **New file**: `src/plots/heatmap_helpers.py` (~100 lines)
- **Private helpers**: ~400 lines (extracted from 3× duplication)
- **Base method**: ~150 lines (unified control flow)
- **Public methods**: 3 × ~15 lines = 45 lines (thin wrappers)
- **Total lines**: ~695 lines (45% reduction)

**Benefits**:
1. ✅ Single source of truth for helper functions
2. ✅ Easier to maintain (fix once, apply everywhere)
3. ✅ Improved testability (can test helpers independently)
4. ✅ Better code readability (clear separation of concerns)
5. ✅ Reduced risk of bugs from copy-paste errors

---

## Implementation Steps

1. **Step 1**: Create `heatmap_helpers.py` with module-level functions
2. **Step 2**: Extract class helper methods (`_plot_top_panel`, etc.)
3. **Step 3**: Create unified `_plot_comparison_heatmap_base` method
4. **Step 4**: Refactor public methods to thin wrappers
5. **Step 5**: Run pipeline test to verify correctness
6. **Step 6**: Visual comparison of outputs (before vs after)

---

## Testing Strategy

**Critical**: Must ensure NO changes to visual output

1. **Baseline Generation**: Save current pipeline outputs (3 PNG files)
2. **Refactored Execution**: Run pipeline with refactored code
3. **Visual Comparison**: Pixel-by-pixel comparison of PNG files
4. **Trace Data Validation**: CSV files must be identical
5. **Log Verification**: Check for warnings/errors

**Acceptance Criteria**:
- ✅ All 3 heatmap PNG files visually identical
- ✅ All 6 trace CSV files (data + summary) byte-identical
- ✅ No new warnings or errors in logs
- ✅ File size reduction: 1,272 → ~600 lines

---

## Risk Assessment

**Low Risk**:
- Module-level helper functions (pure functions, no side effects)
- Private class methods (not exposed to external callers)

**Medium Risk**:
- Base method unification (complex conditional logic)
- Need careful testing of all 3 execution paths

**Mitigation**:
- Keep original file as backup
- Implement incrementally (test after each step)
- Use git for version control

---

**Status**: Plan approved, ready for implementation
**Next**: Phase 2.1 - Create heatmap_helpers.py
