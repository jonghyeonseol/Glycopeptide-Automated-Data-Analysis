# Data Reliability Analysis Report
## pGlyco Auto Combine Pipeline

**Date**: 2025-10-05
**Status**: CRITICAL ISSUES IDENTIFIED

---

## Executive Summary

Analysis of the pGlyco Auto Combine pipeline has revealed **critical inconsistencies** in data filtering and processing across different visualizations. These inconsistencies lead to:

1. Different glycopeptides appearing in different plots
2. Same glycopeptides showing different intensity values across visualizations
3. Non-reproducible scientific conclusions

---

## Issue 1: Inconsistent Detection Frequency Filters

### Current Implementation

| Visualization | Source File | Line | Detection Threshold |
|--------------|-------------|------|-------------------|
| **PLS-DA / VIP Scores** | `src/analyzer.py` | 312-323 | ≥30% in at least one group |
| **Volcano Plot** | `src/plots/volcano_plot.py` | 64-79 | ≥20% in at least one group **AND** ≥5 samples |
| **Glycopeptide Comparison Heatmap** | `src/plots/glycopeptide_comparison_heatmap.py` | 79-97 | ≥50% in at least one group |

### Impact Example

Consider a glycopeptide detected in:
- 10/24 Cancer samples (41.7%)
- 5/23 Normal samples (21.7%)
- Max detection: 41.7%

**Where does it appear?**
- ✅ Volcano plot (>20%)
- ✅ VIP score visualizations (>30%)
- ❌ Comparison heatmap (<50%)

**Result**: User sees this glycopeptide as "significant" in volcano plot but it's missing from the comparison heatmap, creating confusion.

---

## Issue 2: Inconsistent Mean Intensity Calculations

### Three Different Methods Used

#### Method A: VIP Score Plots
**File**: `src/plots/vip_score_plot.py:48-49`

```python
cancer_mean = replace_empty_with_zero(glycopeptide_row[cancer_samples]).mean()
```

- Converts empty strings to 0
- **Includes zeros in mean calculation**
- Denominator = total number of samples (including missing)

#### Method B: Glycopeptide Comparison Heatmap
**File**: `src/plots/glycopeptide_comparison_heatmap.py:73-77`

```python
cancer_stats = calculate_group_statistics(df_with_vip, cancer_samples)
# Uses skipna=True internally
```

- **Excludes missing values** from mean calculation (skipna=True)
- Denominator = number of detected samples only

#### Method C: Volcano Plot
**File**: `src/plots/volcano_plot.py:57-83`

```python
cancer_values = replace_empty_with_zero(row[cancer_samples]).values.astype(float)
cancer_nonzero = cancer_values[cancer_values > 0]
cancer_mean = np.mean(cancer_nonzero)
```

- Explicitly filters out ALL zeros
- **Calculates mean of only non-zero values**
- Denominator = number of non-zero samples

### Impact Example

**Sample intensities**: [100, 200, 0, 0, 300] (3 detected, 2 missing)

| Visualization | Mean Calculation | Result |
|--------------|-----------------|--------|
| VIP Score Plots | (100+200+0+0+300) / 5 | **120** |
| Comparison Heatmap | (100+200+300) / 3 | **200** |
| Volcano Plot | (100+200+300) / 3 | **200** |

**Result**: Same glycopeptide appears with **40% lower intensity** in VIP plots compared to other visualizations.

---

## Issue 3: Different Data Sources and Pre-filtering

### Data Flow Inconsistency

#### VIP Score Plots (Direct Access)
```python
# main.py line 161
visualizer.plot_vip_scores_glycopeptide_r(annotated_data, plsda_results['vip_scores'])
# Uses full annotated_data (no pre-filtering)
```

#### Volcano Plot (Inline Filtering)
```python
# volcano_plot.py recalculates everything
for idx, row in df.iterrows():
    # Applies 20% detection filter inline
```

#### Glycopeptide Comparison Heatmap (Double Filtering)
```python
# glycopeptide_comparison_heatmap.py:69
df_with_vip = df.merge(vip_scores, on=['Peptide', 'GlycanComposition'], how='inner')
# INNER JOIN = only glycopeptides that passed PLS-DA filter (30%)

# Line 91
df_with_vip = df_with_vip[df_with_vip['Max_Detection_Pct'] >= 0.5]
# Additional 50% filter applied
```

**Result**: Comparison heatmap uses a **highly selective double-filtered subset** not used anywhere else in the pipeline.

---

## Scientific Validity Concerns

### 1. Zero-Inflation Bias (Method A)

**Problem**: Including zeros/missing values in mean calculation (VIP plots) violates assumptions for Missing Not At Random (MNAR) data in proteomics.

**Scientific Issue**: Glycopeptides with low detection rates get artificially low mean intensities, biasing downstream analysis.

**Standard Practice**: Use `skipna=True` to exclude missing values from statistical calculations (Methods B and C are correct).

### 2. Non-Reproducible Comparisons

**Problem**: Users cannot validate results across visualizations:
- Volcano plot shows glycopeptide X as "significantly up-regulated"
- VIP plot shows same glycopeptide with low intensity
- Comparison heatmap doesn't show it at all

**Scientific Issue**: Contradictory evidence in the same analysis undermines credibility of findings.

### 3. Selection Bias (Double Filtering)

**Problem**: Glycopeptide comparison heatmap uses cascading filters:
1. PLS-DA filter: 30% detection (from `inner join` with vip_scores)
2. Heatmap filter: 50% detection (explicit filter)

**Scientific Issue**: Only shows glycopeptides that are:
- Detected in ≥30% samples (passed PLS-DA)
- **AND** detected in ≥50% samples (passed heatmap filter)

This creates a **highly selective subset** (approximately top 20-30% of data) not representative of the full dataset.

---

## Recommendations

### Option 1: Quick Fix (Standardize Parameters)

**Time**: 1-2 hours
**Risk**: Low

1. **Standardize detection threshold** to 30% across all visualizations
2. **Standardize mean calculation** to use `skipna=True` (exclude missing values)
3. **Remove double-filtering** in glycopeptide comparison heatmap

**Files to modify**:
- `src/plots/volcano_plot.py` (change 20% → 30%)
- `src/plots/vip_score_plot.py` (change mean calculation to skipna=True)
- `src/plots/glycopeptide_comparison_heatmap.py` (use `left join` instead of `inner join`)

---

### Option 2: Comprehensive Fix (Recommended)

**Time**: 4-6 hours
**Risk**: Medium
**Long-term benefit**: High

#### Step 1: Create Centralized Data Preparation Module

**New file**: `src/data_preparation.py`

```python
def filter_by_detection(df, sample_cols, min_detection_pct=0.30):
    """
    Standardized detection frequency filter

    Args:
        df: DataFrame with sample columns
        sample_cols: List of sample column names
        min_detection_pct: Minimum detection % in at least one group

    Returns:
        Filtered DataFrame
    """
    # Implementation with detailed logging
    pass

def calculate_group_means(df, sample_cols, method='skipna'):
    """
    Standardized mean calculation with consistent missing data handling

    Args:
        df: DataFrame
        sample_cols: Sample column names
        method: 'skipna' (default) or 'replace_zero'

    Returns:
        Series of mean values
    """
    # Implementation
    pass

def prepare_visualization_data(df, vip_scores=None, detection_threshold=0.30):
    """
    Centralized data preparation for all visualizations

    Args:
        df: Annotated DataFrame
        vip_scores: Optional VIP scores for merging
        detection_threshold: Uniform detection threshold

    Returns:
        Prepared DataFrame with consistent filtering and statistics
    """
    # Implementation
    pass
```

#### Step 2: Update All Visualization Modules

**Files to modify**:
- `src/plots/volcano_plot.py`
- `src/plots/vip_score_plot.py`
- `src/plots/vip_score_plot_r.py`
- `src/plots/glycopeptide_comparison_heatmap.py`

Replace inline filtering/calculation logic with calls to centralized functions.

#### Step 3: Add Data Validation

**New file**: `src/data_validator.py`

```python
def validate_visualization_consistency(df_volcano, df_vip, df_heatmap):
    """
    Validate that visualizations use consistent data

    Checks:
    1. Same glycopeptides included (within tolerance)
    2. Same mean intensities (within numerical precision)
    3. Same detection statistics

    Raises:
        ValidationError if inconsistencies detected
    """
    pass
```

#### Step 4: Add Configuration Parameter

**Update**: `config.yaml`

```yaml
analysis:
  detection_filter:
    min_detection_pct: 0.30  # Uniform across all visualizations
    min_samples: 5           # Minimum samples for statistical tests

  missing_data_handling:
    method: 'skipna'  # 'skipna' or 'replace_zero'
```

---

### Option 3: Detailed Technical Specification

Create comprehensive specification document before implementation, including:
1. Detailed data flow diagrams
2. Unit test requirements
3. Validation protocols
4. Backwards compatibility plan

---

## Immediate Actions Required

### Priority 1 (Critical)
- [ ] Decide on fix approach (Option 1 vs 2)
- [ ] Document current behavior for users (warning in documentation)
- [ ] Add disclaimer to analysis_summary.txt about potential inconsistencies

### Priority 2 (High)
- [ ] Implement chosen fix
- [ ] Create unit tests for data consistency
- [ ] Re-run full pipeline to validate fixes

### Priority 3 (Medium)
- [ ] Update CLAUDE.md with standardized data handling guidelines
- [ ] Create developer documentation for data processing standards

---

## Questions for User

1. **Detection threshold**: Should we use 20%, 30%, or 50% uniformly?
   - 20%: Most inclusive, may include noisy data
   - 30%: Balanced (current PLS-DA threshold)
   - 50%: Most stringent, reduces false positives but loses data

2. **Missing data handling**: Confirm `skipna=True` is scientifically appropriate for your data?

3. **Implementation timeline**: Prefer quick fix (Option 1) or comprehensive solution (Option 2)?

4. **Backwards compatibility**: Need to preserve exact numerical results from previous runs?

---

## Appendix: Code Locations

### Detection Filter Implementations
- **PLS-DA**: `src/analyzer.py:312-345` (30% threshold)
- **Volcano**: `src/plots/volcano_plot.py:64-79` (20% threshold + 5 samples)
- **Comparison Heatmap**: `src/plots/glycopeptide_comparison_heatmap.py:79-97` (50% threshold)

### Mean Calculation Implementations
- **VIP Plots**: `src/plots/vip_score_plot.py:48-49` (includes zeros)
- **VIP Plots (R)**: `src/plots/vip_score_plot_r.py:221` (includes zeros)
- **Comparison Heatmap**: `src/plots/glycopeptide_comparison_heatmap.py:73-77` (skipna)
- **Volcano Plot**: `src/plots/volcano_plot.py:57-83` (non-zero only)

### Data Merge Operations
- **Comparison Heatmap**: `src/plots/glycopeptide_comparison_heatmap.py:69` (inner join creates pre-filtering)

---

**Report Generated**: 2025-10-05
**Analysis Tool**: Claude Code (Sonnet 4.5)
