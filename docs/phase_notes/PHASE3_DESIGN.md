# Phase 3: Boxplot.py Refactoring Design

## Executive Summary

**File**: `src/plots/boxplot.py`
**Current Size**: 951 lines
**Target Size**: ~650 lines (31% reduction)
**Estimated Reduction**: ~300 lines

## Current Structure Analysis

### Methods Overview

| Method | Lines | Description | Type |
|--------|-------|-------------|------|
| `calculate_cohens_d()` | 40 | Helper function | Utility |
| `plot_boxplot()` | 185 | Glycan types boxplot | Main |
| `plot_boxplot_extended()` | 193 | Extended categories boxplot | Main |
| `plot_boxplot_primary_classification()` | 106 | Primary classification | Secondary |
| `plot_boxplot_secondary_classification()` | 106 | Secondary classification | Secondary |
| `plot_boxplot_cancer_vs_normal_primary()` | 130 | Cancer vs Normal (primary) | Cancer/Normal |
| `plot_boxplot_cancer_vs_normal_secondary()` | 140 | Cancer vs Normal (secondary) | Cancer/Normal |

**Total**: 7 methods, 900 lines (excluding class definition)

---

## Code Duplication Analysis

### Pattern 1: FDR Correction Logic (HIGH PRIORITY)

**Location**: Methods 1-2 (`plot_boxplot()` and `plot_boxplot_extended()`)

**Duplicated Code**:
- Lines 134-242 (Method 1): ~108 lines
- Lines 319-434 (Method 2): ~116 lines
- **Total Duplication**: ~110 lines × 2 = ~220 lines

**Similarities**:
1. P-value collection loop (identical structure)
2. FDR correction (Benjamini-Hochberg)
3. Significance bracket plotting
4. Effect size calculation (Cohen's d)
5. Logging results

**Differences**:
- Column name: `'GlycanType'` vs `'ExtendedCategory'`
- Category list: `existing_types` vs `existing_categories`
- Log prefix: `'FDR correction applied'` vs `'FDR correction applied (extended)'`

**Refactoring Potential**: **~200 lines** (can be reduced to single method)

---

### Pattern 2: Cancer vs Normal Comparison Logic (MEDIUM PRIORITY)

**Location**: Methods 5-6 (`plot_boxplot_cancer_vs_normal_primary()` and `plot_boxplot_cancer_vs_normal_secondary()`)

**Duplicated Code**:
- Lines 700-810 (Method 5): ~111 lines
- Lines 830-951 (Method 6): ~122 lines
- **Total Duplication**: ~115 lines × 2 = ~230 lines

**Similarities**:
1. TIC normalization (lines 700-709 vs 830-839) - IDENTICAL
2. Data preparation loop (lines 717-749 vs 847-879) - 95% similar
3. Plotting logic (lines 756-810 vs 884-951) - 90% similar
4. QC filtering loop (lines 717-749) - IDENTICAL

**Differences**:
- Classification categories: `primary_categories` vs `secondary_categories`
- Column name: `'PrimaryClassification'` vs `'SecondaryClassification'`
- Figure size: `(10, 6)` vs `(12, 6)`

**Refactoring Potential**: **~180 lines** (can be reduced to single base method + 2 wrappers)

---

### Pattern 3: Classification Boxplot Logic (LOW PRIORITY)

**Location**: Methods 3-4 (`plot_boxplot_primary_classification()` and `plot_boxplot_secondary_classification()`)

**Duplicated Code**:
- Lines 491-574 (Method 3): ~84 lines
- Lines 596-680 (Method 4): ~85 lines
- **Total Duplication**: ~85 lines × 2 = ~170 lines

**Similarities**:
1. Data preparation loop (lines 491-528 vs 596-633) - 95% similar
2. Plotting logic (lines 529-574 vs 635-680) - 98% similar
3. Normalization options (identical)

**Differences**:
- Classification categories: `primary_categories` vs `secondary_categories`
- Column name: `'PrimaryClassification'` vs `'SecondaryClassification'`
- Figure size: `(12, 8)` vs `(14, 8)`

**Refactoring Potential**: **~150 lines** (can be reduced to single base method + 2 wrappers)

---

## Refactoring Strategy

### Phase 3.1: Extract FDR Correction Helper (~200 lines reduction)

**Create**: `_perform_fdr_correction_and_plot_brackets()`

**Purpose**: Consolidate FDR correction logic for Methods 1-2

**Parameters**:
- `ax`: Matplotlib axes
- `boxplot_data`: DataFrame
- `category_column`: Column name ('GlycanType' or 'ExtendedCategory')
- `existing_categories`: List of categories
- `log_prefix`: Logging prefix string

**Impact**:
- Method 1: 185 → ~90 lines
- Method 2: 193 → ~90 lines
- New helper: ~120 lines
- **Net reduction**: (185 + 193) - (90 + 90 + 120) = **78 lines**

---

### Phase 3.2: Extract Cancer vs Normal Base Method (~180 lines reduction)

**Create**: `_plot_boxplot_cancer_vs_normal_base()`

**Purpose**: Unified base method for Methods 5-6

**Parameters**:
- `df`: Annotated DataFrame
- `classification_type`: 'primary' or 'secondary'
- `figsize`: Figure size tuple
- `categories`: List of categories
- `column_name`: Classification column name

**Impact**:
- Method 5: 130 → ~15 lines (thin wrapper)
- Method 6: 140 → ~15 lines (thin wrapper)
- New base method: ~140 lines
- **Net reduction**: (130 + 140) - (15 + 15 + 140) = **100 lines**

---

### Phase 3.3: Extract Classification Base Method (~150 lines reduction)

**Create**: `_plot_boxplot_classification_base()`

**Purpose**: Unified base method for Methods 3-4

**Parameters**:
- `df`: Annotated DataFrame
- `classification_type`: 'primary' or 'secondary'
- `normalization`: 'raw' or 'aggregated'
- `figsize`: Figure size tuple
- `categories`: List of categories
- `column_name`: Classification column name

**Impact**:
- Method 3: 106 → ~10 lines (thin wrapper)
- Method 4: 106 → ~10 lines (thin wrapper)
- New base method: ~95 lines
- **Net reduction**: (106 + 106) - (10 + 10 + 95) = **97 lines**

---

## Expected Final Structure

### Helper Methods (All Private)

1. `calculate_cohens_d()` - Existing utility function (40 lines)
2. `_perform_fdr_correction_and_plot_brackets()` - NEW (Phase 3.1, ~120 lines)
3. `_plot_boxplot_cancer_vs_normal_base()` - NEW (Phase 3.2, ~140 lines)
4. `_plot_boxplot_classification_base()` - NEW (Phase 3.3, ~95 lines)

**Total Helpers**: ~395 lines

### Public Methods (Thin Wrappers or Simplified)

1. `plot_boxplot()` - Simplified (~90 lines after Phase 3.1)
2. `plot_boxplot_extended()` - Simplified (~90 lines after Phase 3.1)
3. `plot_boxplot_primary_classification()` - Thin wrapper (~10 lines after Phase 3.3)
4. `plot_boxplot_secondary_classification()` - Thin wrapper (~10 lines after Phase 3.3)
5. `plot_boxplot_cancer_vs_normal_primary()` - Thin wrapper (~15 lines after Phase 3.2)
6. `plot_boxplot_cancer_vs_normal_secondary()` - Thin wrapper (~15 lines after Phase 3.2)

**Total Public Methods**: ~230 lines

---

## Total Impact Summary

**Before**: 951 lines
**After**:
- Helpers: 395 lines
- Public methods: 230 lines
- Class definition: ~50 lines (unchanged)
- **Total**: ~675 lines

**Reduction**: 951 → 675 = **276 lines (29% reduction)**

---

## Implementation Order (Recommended)

### Option A: Bottom-Up Approach (Safest)
1. Phase 3.3: Classification base method (smallest, least complex)
2. Phase 3.2: Cancer vs Normal base method (medium complexity)
3. Phase 3.1: FDR correction helper (most complex, high impact)

**Advantages**: Build confidence, test incrementally, easier rollback

### Option B: Top-Down Approach (Fastest Impact)
1. Phase 3.1: FDR correction helper (highest impact first)
2. Phase 3.2: Cancer vs Normal base method
3. Phase 3.3: Classification base method

**Advantages**: Maximum line reduction early, immediate value

---

## Testing Strategy

After each phase:
1. ✅ Syntax check: `python3 -m py_compile src/plots/boxplot.py`
2. ✅ Full pipeline test: `python3 main.py`
3. ✅ Verify PNG outputs:
   - `boxplot_glycan_types.png`
   - `boxplot_extended_categories.png`
   - `boxplot_primary_raw_normalized.png`
   - `boxplot_primary_aggregated_normalized.png`
   - `boxplot_secondary_raw_normalized.png`
   - `boxplot_secondary_aggregated_normalized.png`
   - `boxplot_primary_cancer_vs_normal.png`
   - `boxplot_primary_cancer_vs_normal_qc.png`
   - `boxplot_secondary_cancer_vs_normal.png`
   - `boxplot_secondary_cancer_vs_normal_qc.png`
4. ✅ Verify trace data CSVs (10 files)
5. ✅ Visual comparison (before/after screenshots)

---

## Risk Assessment

**Low Risk**:
- Phase 3.3 (classification base): Highly similar code, clear parameters
- All refactorings maintain exact same logic

**Medium Risk**:
- Phase 3.2 (Cancer vs Normal base): TIC normalization must be preserved exactly
- QC filtering loop must work identically

**High Risk**:
- Phase 3.1 (FDR correction): Statistical calculations must be identical
- Bracket positioning must be pixel-perfect

**Mitigation**: Test after each phase, compare outputs numerically and visually

---

## Next Steps

**Decision Required**: Which implementation order?
- **A) Bottom-Up** (3.3 → 3.2 → 3.1): Safest, incremental
- **B) Top-Down** (3.1 → 3.2 → 3.3): Fastest impact

**Recommendation**: **Option A (Bottom-Up)** for first refactoring of this file type, then Option B for subsequent files once pattern is established.

---

## Success Criteria

✅ All 10 PNG files generate identically
✅ All 10 trace data CSV files match
✅ Statistical values unchanged (p-values, FDR, Cohen's d)
✅ File size reduced by 250-300 lines (26-31%)
✅ No breaking changes to public APIs
✅ Code is more maintainable and DRY

---

*Phase 3 Design Document*
*Created: 2025-10-14*
*Target: boxplot.py refactoring*
