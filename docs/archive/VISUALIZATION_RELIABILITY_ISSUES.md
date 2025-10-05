# Visualization Reliability Issues - Scientific Analysis

**Date**: 2025-10-05
**Issue**: User reports "visualized results are so suspicious"
**Root Cause**: **Sparse feature detection** - many glycopeptides detected in too few samples

---

## Executive Summary

**CRITICAL FINDING**: The visualizations include comparisons of glycopeptides that are detected in very few samples, making the results **statistically unreliable** and potentially misleading.

### Key Statistics

**Glycopeptide Comparison Heatmap** (Top 75 glycopeptides by VIP):
- **Median Cancer detection**: 4/24 samples (17%)
- **Median Normal detection**: 9/23 samples (39%)
- **79% of glycopeptides**: < 50% detection in Cancer
- **47% of glycopeptides**: < 50% detection in BOTH groups

**Volcano Plot** (1971 glycopeptides tested):
- Minimum detection threshold: 3 samples per group (12.5%)
- **90 significant features** (FDR < 0.05) - detection rates unknown
- No detection count columns exported for verification

---

## Problem 1: Sparse Detection in Heatmap

### Current Behavior

The `glycopeptide_comparison_heatmap.py` shows glycopeptides with **no detection frequency filtering**:

```python
# Example from actual data (Results/Trace/glycopeptide_comparison_heatmap_data.csv):
# Peptide: TVDKSTGKPTLYJVSLVMSDTAGTCY, Glycan: H(8)N(2)
# Cancer_SampleCount: 21/24 (88%)  ← Good
# Normal_SampleCount: 2/23 (9%)    ← UNRELIABLE!
# Fold_Change: 2.50
# This comparison is scientifically questionable
```

### Why This Is A Problem

**Scientific Issue**: Comparing a feature detected in 21 Cancer samples vs. 2 Normal samples:
1. **Biased estimate**: Normal mean is based on only 2 measurements (high variance)
2. **Detection bias**: The feature might be truly absent in most Normal samples OR below detection limit
3. **No statistical power**: Cannot distinguish biological difference from technical noise

**Mass Spectrometry Standard Practice**:
- **Require 50-70% detection in at least one group** before including in comparisons
- Or require **minimum absolute count** (e.g., ≥5 samples AND ≥30% in at least one group)

### Impact on Results

Looking at the heatmap data:
```
Cancer_SampleCount distribution:
  25th percentile: 1 sample
  Median: 4 samples
  75th percentile: 10 samples

Normal_SampleCount distribution:
  25th percentile: 1 sample
  Median: 9 samples
  75th percentile: 21 samples
```

**Interpretation**: Most glycopeptides in the heatmap are **under-detected** in at least one group.

---

## Problem 2: Minimal Detection Threshold in Volcano Plot

### Current Behavior

The `volcano_plot.py` filters at line 56-57:

```python
# Skip if insufficient data
if len(cancer_nonzero) < 3 or len(normal_nonzero) < 3:
    continue
```

**Threshold**: 3 non-zero values per group
- 3/24 Cancer samples = **12.5% detection**
- 3/23 Normal samples = **13% detection**

### Why This Is Too Permissive

**Statistical Issue**: Mann-Whitney U test with n=3 per group:
- **Minimum p-value possible**: ~0.05 (even if groups are perfectly separated)
- **Power**: Extremely low to detect true differences
- **High Type II error**: Will miss many real differences

**FDR Correction Issue**: Running 1971 tests with many under-powered tests inflates false negatives.

### Recommendation

**Minimum detection threshold should be**:
- **Option 1 (Conservative)**: ≥50% detection in at least one group
  - Cancer: ≥12/24 samples OR Normal: ≥11/23 samples
- **Option 2 (Moderate)**: ≥30% detection in at least one group
  - Cancer: ≥7/24 samples OR Normal: ≥7/23 samples
- **Option 3 (Minimal)**: ≥5 samples AND ≥20% in at least one group

Current threshold (n≥3, 12.5%) is **below acceptable standards**.

---

## Problem 3: PLS-DA and VIP Scores on Sparse Data

### Current Behavior

The `analyzer.py` performs PLS-DA on the **full dataset** (6434 glycopeptides) without filtering for detection frequency.

```python
# From analyzer.py line 86-88
intensity_matrix = df[sample_cols].copy()
intensity_matrix = replace_empty_with_zero(intensity_matrix)
# → Zeros replace ALL missing values before PLS-DA
```

### Why This Is A Problem

**PLS-DA Assumption**: Features should be measured consistently across samples
- **Violates assumption**: Many glycopeptides have 0 in 50-90% of samples
- **Model bias**: PLS-DA weights are influenced by sparse detection patterns
- **VIP score reliability**: VIP scores for sparsely detected features are unstable

### Example Impact

Looking at heatmap glycopeptides:
```
VIP_Score: 3.43 (appears important)
Cancer detection: 21/24 (good)
Normal detection: 2/23 (sparse!)

→ High VIP might reflect detection bias, not biological importance
```

### Recommendation

**Option 1**: Filter data before PLS-DA
```python
# Before PLS-DA, filter to glycopeptides with:
# - Detected in ≥50% of samples in at least one group
# - This reduces feature count but increases reliability
```

**Option 2**: Add detection rate as covariate
```python
# Include detection frequency as additional information
# Helps PLS-DA account for detection bias
```

---

## Problem 4: Inconsistent Detection Metrics

### Current State

**What's MISSING from trace data**:
1. **Detection percentage** (e.g., "detected in 75% of Cancer samples")
2. **Detection difference** (e.g., "50% more detected in Cancer vs Normal")
3. **Minimum reliable sample count** warnings

**What EXISTS**:
- `Cancer_SampleCount`, `Normal_SampleCount` (raw counts)
- But no interpretation or filtering based on these

### Why This Matters

Users (and reviewers) need to know:
- "This fold change is based on 20/24 Cancer vs 18/23 Normal" ← Reliable
- "This fold change is based on 3/24 Cancer vs 2/23 Normal" ← **Unreliable, should not be shown**

---

## Recommended Fixes (Priority Order)

### HIGH PRIORITY (Affects Scientific Validity)

#### 1. Add Detection Frequency Filter to Heatmap

**File**: `src/plots/glycopeptide_comparison_heatmap.py`

**Implementation**:
```python
def plot_glycopeptide_comparison_heatmap(self, df, vip_scores, output_dir,
                                         max_peptides=20, max_glycans_per_type=15,
                                         min_detection_pct=0.5):  # NEW PARAMETER
    """
    Args:
        min_detection_pct: Minimum detection percentage in at least one group (default: 0.5 = 50%)
    """

    # After merging VIP scores
    cancer_samples, normal_samples = get_sample_columns(df)

    # Calculate detection statistics
    cancer_stats = calculate_group_statistics(df_with_vip, cancer_samples)
    normal_stats = calculate_group_statistics(df_with_vip, normal_samples)

    df_with_vip['Cancer_Detection_Pct'] = cancer_stats['count'] / len(cancer_samples)
    df_with_vip['Normal_Detection_Pct'] = normal_stats['count'] / len(normal_samples)

    # Filter: require min_detection_pct in at least one group
    df_with_vip['Max_Detection_Pct'] = df_with_vip[['Cancer_Detection_Pct', 'Normal_Detection_Pct']].max(axis=1)
    df_filtered = df_with_vip[df_with_vip['Max_Detection_Pct'] >= min_detection_pct].copy()

    logger.info(f"Detection filtering: {len(df_filtered)}/{len(df_with_vip)} glycopeptides pass ≥{min_detection_pct*100}% detection")

    # Continue with existing logic on df_filtered
```

#### 2. Increase Volcano Plot Detection Threshold

**File**: `src/plots/volcano_plot.py` line 56-57

**Current**:
```python
if len(cancer_nonzero) < 3 or len(normal_nonzero) < 3:
    continue
```

**Recommended**:
```python
# Require minimum 5 samples AND 20% detection in at least one group
min_samples = 5
min_detection_pct = 0.20

cancer_detection_pct = len(cancer_nonzero) / len(cancer_samples)
normal_detection_pct = len(normal_nonzero) / len(normal_samples)
max_detection_pct = max(cancer_detection_pct, normal_detection_pct)

if (len(cancer_nonzero) < min_samples and len(normal_nonzero) < min_samples) or \
   (max_detection_pct < min_detection_pct):
    continue
```

**Also add detection counts to trace data**:
```python
volcano_data.append({
    # ... existing fields ...
    'Cancer_Count': len(cancer_nonzero),
    'Normal_Count': len(normal_nonzero),
    'Cancer_Detection_Pct': len(cancer_nonzero) / len(cancer_samples),
    'Normal_Detection_Pct': len(normal_nonzero) / len(normal_samples)
})
```

#### 3. Filter Data Before PLS-DA

**File**: `src/analyzer.py` in `perform_plsda()` method

**Add before line 86**:
```python
def perform_plsda(self, df: pd.DataFrame, n_components: int = DEFAULT_PLSDA_COMPONENTS,
                  min_detection_pct: float = 0.3) -> Dict:
    """
    Perform PLS-DA with detection frequency filtering

    Args:
        min_detection_pct: Minimum detection % in at least one group (default: 0.3 = 30%)
    """

    # Get sample columns
    cancer_samples, normal_samples = get_sample_columns(df)

    # Calculate detection frequency
    cancer_detection = (df[cancer_samples] != '').sum(axis=1) / len(cancer_samples)
    normal_detection = (df[normal_samples] != '').sum(axis=1) / len(normal_samples)
    max_detection = pd.concat([cancer_detection, normal_detection], axis=1).max(axis=1)

    # Filter
    df_filtered = df[max_detection >= min_detection_pct].copy()

    logger.info(f"PLS-DA detection filter: {len(df_filtered)}/{len(df)} glycopeptides (≥{min_detection_pct*100}% detection)")

    # Continue with existing PLS-DA on df_filtered
    intensity_matrix, sample_names, feature_info = self.prepare_intensity_matrix(df_filtered)
```

### MEDIUM PRIORITY (Improves Transparency)

#### 4. Add Detection Metrics to Summary Report

**File**: `main.py` in summary generation

**Add**:
```python
# Detection frequency analysis
cancer_samples, normal_samples = get_sample_columns(annotated_data)
cancer_detection = (annotated_data[cancer_samples] != '').sum(axis=1) / len(cancer_samples)
normal_detection = (annotated_data[normal_samples] != '').sum(axis=1) / len(normal_samples)

summary_lines.append(f"\nDetection Frequency Analysis:")
summary_lines.append(f"  Glycopeptides detected in ≥50% of Cancer samples: {(cancer_detection >= 0.5).sum()} ({(cancer_detection >= 0.5).sum()/len(annotated_data)*100:.1f}%)")
summary_lines.append(f"  Glycopeptides detected in ≥50% of Normal samples: {(normal_detection >= 0.5).sum()} ({(normal_detection >= 0.5).sum()/len(annotated_data)*100:.1f}%)")
summary_lines.append(f"  Glycopeptides detected in ≥50% of at least one group: {((cancer_detection >= 0.5) | (normal_detection >= 0.5)).sum()}")
```

#### 5. Add Detection Warnings to Visualizations

**All visualization files**: Add text annotation when showing sparse features

```python
# Example for heatmap
if sparse_features_shown:
    ax.text(0.02, 0.98,
           f"⚠️ {sparse_count} glycopeptides shown have <50% detection",
           transform=ax.transAxes, fontsize=10,
           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
```

---

## Configuration Changes Needed

**File**: `config.yaml`

**Add new section**:
```yaml
visualization:
  # ... existing ...

  detection_filtering:
    # Minimum detection frequency for including features in visualizations
    # Set to 0.0 to disable filtering (not recommended)
    min_detection_pct: 0.5  # 50% in at least one group

    # Minimum absolute sample count
    min_sample_count: 5

    # Apply filtering to specific plots
    apply_to:
      - glycopeptide_comparison_heatmap
      - volcano_plot
      - pls_da  # Filters data before PLS-DA

statistical_tests:
  # Minimum requirements for volcano plot statistical tests
  min_samples_per_group: 5
  min_detection_pct: 0.2  # 20% in at least one group
```

---

## Testing After Fixes

**Expected Results**:

1. **Heatmap**: Should show ~20-40 glycopeptides instead of 75 (with 50% detection filter)
2. **Volcano Plot**: Should test ~500-800 glycopeptides instead of 1971
3. **Significant Features**: Should have higher reliability (fewer false positives from sparse data)
4. **VIP Scores**: Should be more stable (calculated on consistently detected features)

**Verification**:
```python
# Check detection stats in filtered heatmap
df = pd.read_csv('Results/Trace/glycopeptide_comparison_heatmap_data.csv')
print(f"Cancer detection median: {df['Cancer_Detection_Pct'].median()}")
print(f"Normal detection median: {df['Normal_Detection_Pct'].median()}")
# Should both be >0.5 with recommended filters
```

---

## Summary

**Current State**: Visualizations include many glycopeptides with **sparse detection** (median 17-39%), violating mass spectrometry best practices.

**Why Suspicious**: Fold changes and statistical tests on sparsely detected features are **unreliable** and can mislead interpretation.

**Recommended Action**: Implement HIGH PRIORITY fixes (detection filtering) before using for publication or biological interpretation.

**Impact**: After filtering, you'll see:
- Fewer glycopeptides in visualizations (but more reliable)
- Higher confidence in fold changes (based on adequate sample sizes)
- More interpretable results (detection patterns visible in data)

---

**Report Generated**: 2025-10-05
**Severity**: HIGH - Affects scientific validity of all comparative visualizations
