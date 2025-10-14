# Phase 2 Refactoring Progress Report

## Phase 2.1: Module-Level Helper Functions ✅ COMPLETED

**Created**: `src/plots/heatmap_helpers.py` (214 lines)

**Extracted Functions**:
1. `glycan_sort_key()` - Natural sorting for glycan compositions
2. `get_symbol_info()` - Determine symbols and linewidth based on group presence
3. `create_peptide_order()` - Sort peptides by VIP score
4. `create_glycan_order()` - Group and sort glycans by type
5. `create_index_mappings()` - Create position dictionaries

**Impact**: Created single source of truth for helper functions used across all 3 methods

---

## Phase 2.2a: Import Helper Functions and Remove Duplicates ✅ COMPLETED

**File**: `src/plots/glycopeptide_comparison_heatmap.py`

### Changes Made

1. **Added imports** (lines 19-25):
```python
from .heatmap_helpers import (
    glycan_sort_key,
    get_symbol_info,
    create_peptide_order,
    create_glycan_order,
    create_index_mappings
)
```

2. **Replaced duplicate code in all 3 methods**:
   - Method 1 (`plot_glycopeptide_comparison_heatmap`):
     - Lines 157-198 (42 lines) → Lines 164-173 (10 lines) = **32 lines removed**
     - Lines 442-463 (22 lines) → Lines 417-422 (6 lines) = **16 lines removed**

   - Method 2 (`plot_glycopeptide_comparison_heatmap_full`):
     - Lines 574-611 (38 lines) → Lines 534-543 (10 lines) = **28 lines removed**
     - Lines 853-871 (19 lines) → Lines 783-788 (6 lines) = **13 lines removed**

   - Method 3 (`plot_glycopeptide_comparison_heatmap_by_type`):
     - Lines 983-1007 (25 lines) → Lines 901-909 (9 lines) = **16 lines removed**
     - Lines 1214-1232 (19 lines) → Lines 1116-1121 (6 lines) = **13 lines removed**

### Results

**File Size Reduction**:
- **Before**: 1,273 lines
- **After**: 1,160 lines
- **Removed**: 113 lines (8.9% reduction)

**Testing**:
- ✅ Syntax check passed
- ✅ Full pipeline test passed (71 seconds)
- ✅ All 7 glycopeptide comparison heatmap PNG files generated correctly:
  - glycopeptide_comparison_heatmap.png (445K)
  - glycopeptide_comparison_heatmap_full.png (4.0M)
  - glycopeptide_comparison_heatmap_HM.png (551K)
  - glycopeptide_comparison_heatmap_F.png (918K)
  - glycopeptide_comparison_heatmap_S.png (2.1M)
  - glycopeptide_comparison_heatmap_SF.png (1.4M)
  - glycopeptide_comparison_heatmap_C_H.png (505K)

**Code Quality Improvements**:
- ✅ Eliminated duplicate function definitions (5 functions × 3 methods = 15 definitions → 5 definitions)
- ✅ Single source of truth for all helper functions
- ✅ Improved maintainability (changes in one place apply to all methods)
- ✅ Better testability (helper functions can be tested independently)
- ✅ Cleaner, more readable code

---

## Next Steps

### Phase 2.2b: Extract Class Helper Methods (Pending)
Extract common plotting logic into private class methods:
- `_plot_top_panel()` - Average intensity line plots
- `_plot_colorbar_panel()` - Glycan type color bar
- `_plot_symbol_heatmap()` - Symbol-based heatmap
- `_create_legend()` - Legend creation
- `_prepare_and_save_trace_data()` - Trace data preparation

**Expected Impact**: ~400-500 lines reduction

### Phase 2.3: Create Unified Base Method (Pending)
Create `_plot_comparison_heatmap_base()` with conditional logic for all 3 variants

**Expected Impact**: ~200 lines reduction

### Phase 2.4: Refactor Public Methods (Pending)
Convert public methods to thin wrappers calling base method

**Expected Impact**: Each method reduced to ~10-15 lines

### Phase 2.5: Final Testing and Verification (Pending)
- Comprehensive pipeline testing
- Visual comparison of outputs
- CSV trace data validation
- Performance benchmarking

---

## Total Progress

**Current State**:
- Lines reduced: 113 (8.9%)
- Current size: 1,160 lines
- Target size: ~600 lines

**Progress**: 16.6% towards target (113 out of 673 lines to remove)

**Files Modified**:
- ✅ Created: `src/plots/heatmap_helpers.py`
- ✅ Modified: `src/plots/glycopeptide_comparison_heatmap.py`

**Status**: Phase 2.1 and 2.2a completed successfully. Ready to proceed with Phase 2.2b.
