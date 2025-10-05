# Trace Data Guide - Manual Verification

## 📊 Overview

The visualization generates **comprehensive trace data** to ensure complete traceability from the plot back to the original data. This allows you to manually verify every dot, calculation, and statistic shown in the heatmap.

---

## 📁 Output Files

## ⚠️ Important Note About Data Values

**The trace data contains the EXACT values used to create the visualization:**

- `Cancer_Mean` = Aggregated mean intensity (averaged across C1-C24 samples)
- `Normal_Mean` = Aggregated mean intensity (averaged across N1-N23 samples)
- Individual sample columns (`C1`, `C2`, ..., `N24`) = **Original raw intensities** from `integrated.csv`

**These means are calculated directly from the raw sample values** (after replacing empty cells with 0).
They represent the actual data points plotted as dots in the visualization.

---

### 1. Summary File (Recommended for Quick Review)

**File**: `Results/Trace/glycopeptide_comparison_heatmap_summary.csv`
**Size**: ~19 KB
**Rows**: Same as number of dots in the plot

**Columns** (22 total):

| Column | Description | Use for Verification |
|--------|-------------|---------------------|
| `Peptide` | Peptide sequence | Identify the feature |
| `GlycanComposition` | Glycan structure (e.g., H(5)N(4)A(2)) | Identify the feature |
| `GlycanTypeCategory` | Type: HM, F, S, SF, C/H | Verify color grouping |
| `Plot_X_Position` | X coordinate in plot (0-based) | Find exact column position |
| `Plot_Y_Position` | Y coordinate in plot (0-based) | Find exact row position |
| `PeptideVIP` | VIP score from PLS-DA | Verify Y-axis sorting |
| **Cancer Group** | | |
| `Cancer_Mean` | Mean intensity across C1-C24 | Verify left dot (○) intensity |
| `Cancer_StdDev` | Standard deviation | Assess variability |
| `Cancer_SampleCount` | # samples with non-zero value | Check detection frequency |
| `Cancer_Min` | Minimum value across samples | See range |
| `Cancer_Max` | Maximum value across samples | See range |
| **Normal Group** | | |
| `Normal_Mean` | Mean intensity across N1-N23 | Verify right dot (□) intensity |
| `Normal_StdDev` | Standard deviation | Assess variability |
| `Normal_SampleCount` | # samples with non-zero value | Check detection frequency |
| `Normal_Min` | Minimum value across samples | See range |
| `Normal_Max` | Maximum value across samples | See range |
| **Comparison** | | |
| `Fold_Change` | Cancer_Mean / Normal_Mean | Ratio comparison |
| `Log2_Fold_Change` | Log2 fold change | Standard metric |
| **Plot Properties** | | |
| `Cancer_Alpha` | Transparency value (0-1) | Verify dot darkness |
| `Normal_Alpha` | Transparency value (0-1) | Verify dot darkness |
| `Cancer_Dot_Plotted` | TRUE/FALSE | Was dot shown? |
| `Normal_Dot_Plotted` | TRUE/FALSE | Was dot shown? |

### 2. Full Data File (Complete Individual Samples)

**File**: `Results/Trace/glycopeptide_comparison_heatmap_data.csv`
**Size**: ~55 KB
**Rows**: Same as number of dots in the plot

**Contains**: All 22 summary columns PLUS all individual sample intensities:
- `C1`, `C2`, ..., `C24` (24 cancer samples)
- `N1`, `N2`, ..., `N24` (23 normal samples, N19 missing in your data)

**Total Columns**: 22 + 47 = **69 columns**

---

## 🔍 Manual Verification Examples

### Example 1: Verify a Specific Dot

**Goal**: Check if the dot at position (X=5, Y=2) is correctly plotted.

#### Step 1: Find the data row
```bash
# Using the summary file
grep "Plot_X_Position.*5" Results/Trace/glycopeptide_comparison_heatmap_summary.csv | \
grep "Plot_Y_Position.*2"
```

Or open in Excel/spreadsheet and filter:
- `Plot_X_Position = 5`
- `Plot_Y_Position = 2`

#### Step 2: Verify the values

Check the row shows:
- ✓ Peptide sequence
- ✓ Glycan composition
- ✓ Cancer_Mean (should match left dot darkness)
- ✓ Normal_Mean (should match right dot darkness)
- ✓ Cancer_Dot_Plotted = TRUE (if left dot is visible)
- ✓ Normal_Dot_Plotted = TRUE (if right dot is visible)

#### Step 3: Verify transparency calculation

**Formula used in code**:
```
Alpha = min(0.3 + (Intensity / Max_Intensity) * 0.7, 1.0)
```

Where:
- Min alpha = 0.3 (lightest)
- Max alpha = 1.0 (darkest)
- Scales linearly with intensity

**Manual check**:
1. Find `max_intensity` = highest value in Cancer_Mean and Normal_Mean columns
2. Calculate: `expected_alpha = min(0.3 + (Cancer_Mean / max_intensity) * 0.7, 1.0)`
3. Compare with `Cancer_Alpha` in trace data

### Example 2: Verify Group Averages

**Goal**: Check if Cancer_Mean is correctly calculated.

#### Using the full data file:

```python
import pandas as pd

# Load full data
df = pd.read_csv('Results/Trace/glycopeptide_comparison_heatmap_data.csv')

# Select a row
row = df.iloc[0]

# Cancer sample columns
cancer_cols = [f'C{i}' for i in range(1, 25)]

# Calculate mean manually
cancer_values = row[cancer_cols].astype(float)
manual_mean = cancer_values.mean()

# Compare with saved value
saved_mean = row['Cancer_Mean']

print(f"Manual: {manual_mean}")
print(f"Saved:  {saved_mean}")
print(f"Match:  {abs(manual_mean - saved_mean) < 0.01}")
```

### Example 3: Verify VIP Score Sorting

**Goal**: Check if peptides are sorted by VIP score (descending).

```python
import pandas as pd

df = pd.read_csv('Results/Trace/glycopeptide_comparison_heatmap_summary.csv')

# Group by Y position and get VIP scores
vip_by_position = df.groupby('Plot_Y_Position')['PeptideVIP'].mean().sort_index()

# Check if descending (higher VIP at lower Y positions)
is_descending = all(vip_by_position.iloc[i] >= vip_by_position.iloc[i+1]
                   for i in range(len(vip_by_position)-1))

print(f"VIP scores properly sorted: {is_descending}")
print(vip_by_position)
```

### Example 4: Verify Glycan Type Grouping

**Goal**: Check if glycans are grouped correctly by type.

```python
import pandas as pd

df = pd.read_csv('Results/Trace/glycopeptide_comparison_heatmap_summary.csv')

# Check grouping by X position
glycan_order = df.sort_values('Plot_X_Position')[['Plot_X_Position', 'GlycanComposition', 'GlycanTypeCategory']]

print(glycan_order)

# Verify that glycans of same type are together
for gtype in ['HM', 'F', 'S', 'SF', 'C/H']:
    positions = df[df['GlycanTypeCategory'] == gtype]['Plot_X_Position'].values
    if len(positions) > 1:
        is_contiguous = all(positions[i+1] == positions[i]+1 for i in range(len(positions)-1))
        print(f"{gtype}: Contiguous = {is_contiguous}")
```

### Example 5: Verify Statistical Calculations

**Goal**: Manually recalculate fold change and verify.

```python
import pandas as pd
import numpy as np

df = pd.read_csv('Results/Trace/glycopeptide_comparison_heatmap_summary.csv')

# Select a row
row = df.iloc[0]

# Manual fold change
if row['Normal_Mean'] > 0:
    manual_fc = row['Cancer_Mean'] / row['Normal_Mean']
    manual_log2fc = np.log2(manual_fc) if manual_fc > 0 else np.nan
else:
    manual_fc = np.inf
    manual_log2fc = np.nan

print(f"Saved FC:     {row['Fold_Change']}")
print(f"Manual FC:    {manual_fc}")
print(f"Saved Log2FC: {row['Log2_Fold_Change']}")
print(f"Manual Log2FC:{manual_log2fc}")
```

---

## 📋 Verification Checklist

Use this checklist to validate the visualization:

### Data Integrity
- [ ] All rows have valid Peptide and GlycanComposition
- [ ] Plot positions (X, Y) are within expected range
- [ ] No duplicate positions (each X,Y combination is unique)
- [ ] All VIP scores are positive numbers
- [ ] Sample counts ≤ total sample number (24 for Cancer, 23 for Normal)

### Statistical Calculations
- [ ] Cancer_Mean matches manual calculation from individual samples
- [ ] Normal_Mean matches manual calculation from individual samples
- [ ] Standard deviations are reasonable (< mean for most cases)
- [ ] Fold changes are correctly calculated
- [ ] Log2 fold changes are correctly calculated

### Plot Properties
- [ ] Alpha values are between 0.3 and 1.0
- [ ] Dots marked as plotted (TRUE) have Mean > 0
- [ ] Dots not plotted (FALSE) have Mean = 0 or very low
- [ ] Transparency increases with intensity

### Sorting and Grouping
- [ ] Y-axis (peptides) sorted by VIP score (descending)
- [ ] X-axis (glycans) grouped by type (HM, F, S, SF, C/H)
- [ ] Within each type, glycans are ordered by total intensity

---

## 🔧 Excel/Spreadsheet Verification

### Quick Steps in Excel:

1. **Open summary file**: `glycopeptide_comparison_heatmap_summary.csv`

2. **Enable filtering**: Select all → Data → Filter

3. **Find a specific dot**:
   - Filter `Plot_X_Position` = column number
   - Filter `Plot_Y_Position` = row number

4. **Verify calculations**:
   ```
   # Add a new column to verify fold change
   = Cancer_Mean / Normal_Mean

   # Add a new column to verify log2 fold change
   = LOG(Cancer_Mean / Normal_Mean, 2)
   ```

5. **Check sorting**:
   - Sort by `Plot_Y_Position` (ascending)
   - Check `PeptideVIP` column decreases

6. **Check grouping**:
   - Sort by `Plot_X_Position` (ascending)
   - Check `GlycanTypeCategory` groups together

---

## 🐍 Python Verification Script

Create `verify_trace_data.py`:

```python
import pandas as pd
import numpy as np

def verify_heatmap_data():
    """Comprehensive verification of trace data"""

    print("="*80)
    print("Verification of Glycopeptide Comparison Heatmap Trace Data")
    print("="*80)

    # Load data
    summary = pd.read_csv('Results/Trace/glycopeptide_comparison_heatmap_summary.csv')
    full_data = pd.read_csv('Results/Trace/glycopeptide_comparison_heatmap_data.csv')

    print(f"\n1. Data loaded successfully")
    print(f"   - Summary rows: {len(summary)}")
    print(f"   - Full data rows: {len(full_data)}")

    # Check 1: Verify means match individual samples
    print(f"\n2. Verifying Cancer_Mean calculations...")
    cancer_cols = [f'C{i}' for i in range(1, 25)]

    mismatches = 0
    for idx, row in full_data.iterrows():
        cancer_values = row[cancer_cols].astype(float)
        manual_mean = cancer_values.mean()
        saved_mean = row['Cancer_Mean']

        if abs(manual_mean - saved_mean) > 0.01:
            mismatches += 1

    print(f"   - Mismatches: {mismatches}")
    print(f"   - Status: {'PASS' if mismatches == 0 else 'FAIL'}")

    # Check 2: Verify VIP sorting
    print(f"\n3. Verifying VIP score sorting...")
    vip_by_pos = summary.groupby('Plot_Y_Position')['PeptideVIP'].mean().sort_index()
    is_sorted = all(vip_by_pos.iloc[i] >= vip_by_pos.iloc[i+1]
                   for i in range(len(vip_by_pos)-1))
    print(f"   - Properly sorted (descending): {is_sorted}")
    print(f"   - Status: {'PASS' if is_sorted else 'FAIL'}")

    # Check 3: Verify fold changes
    print(f"\n4. Verifying fold change calculations...")
    fc_mismatches = 0
    for idx, row in summary.iterrows():
        if row['Normal_Mean'] > 0:
            manual_fc = row['Cancer_Mean'] / row['Normal_Mean']
            if abs(manual_fc - row['Fold_Change']) > 0.01:
                fc_mismatches += 1

    print(f"   - Mismatches: {fc_mismatches}")
    print(f"   - Status: {'PASS' if fc_mismatches == 0 else 'FAIL'}")

    # Check 4: Verify plot flags
    print(f"\n5. Verifying plot flags...")
    cancer_flag_errors = sum((summary['Cancer_Mean'] > 0) != summary['Cancer_Dot_Plotted'])
    normal_flag_errors = sum((summary['Normal_Mean'] > 0) != summary['Normal_Dot_Plotted'])

    print(f"   - Cancer flag errors: {cancer_flag_errors}")
    print(f"   - Normal flag errors: {normal_flag_errors}")
    print(f"   - Status: {'PASS' if (cancer_flag_errors + normal_flag_errors) == 0 else 'FAIL'}")

    # Summary
    print(f"\n" + "="*80)
    print(f"VERIFICATION COMPLETE")
    print(f"="*80)

    all_pass = (mismatches == 0 and is_sorted and
                fc_mismatches == 0 and
                (cancer_flag_errors + normal_flag_errors) == 0)

    if all_pass:
        print(f"✓ All checks PASSED - Data is reliable")
    else:
        print(f"✗ Some checks FAILED - Review above for details")

    return all_pass

if __name__ == "__main__":
    verify_heatmap_data()
```

**Run the script**:
```bash
python3 verify_trace_data.py
```

---

## 📊 Understanding the Data Flow

### From Raw Data to Visualization:

```
1. Raw CSV files (C_01.csv, ..., N_24.csv)
   ↓
2. Integrated data (integrated.csv)
   ↓
3. Annotated data (with GlycanTypeCategory)
   ↓
4. PLS-DA analysis → VIP scores
   ↓
5. Filter top peptides by VIP
   ↓
6. Filter top glycans by type and intensity
   ↓
7. Calculate group averages (Cancer_Mean, Normal_Mean)
   ↓
8. Plot dots with transparency based on intensity
   ↓
9. Save comprehensive trace data ← YOU ARE HERE
```

Every step is traceable through the saved files.

---

## 🎯 Key Verification Points

### Critical Values to Check:

1. **VIP Scores**: Should match `vip_scores_all.csv`
2. **Individual Samples**: Should match `integrated.csv`
3. **Group Means**: Should be average of individual samples
4. **Fold Changes**: Should be Cancer_Mean / Normal_Mean
5. **Plot Positions**: Should reflect VIP sorting and type grouping

### Red Flags:

❌ Fold change = infinity but Normal_Mean > 0
❌ Sample count > total samples
❌ Standard deviation > mean (for most cases)
❌ Alpha values outside [0.3, 1.0]
❌ Dots plotted but Mean = 0

---

## 💾 Additional Trace Files

All visualizations in the pipeline save trace data:

- `pca_plot_data.csv`
- `heatmap_top_glycopeptides_data.csv`
- `boxplot_glycan_types_data.csv`
- `glycopeptide_comparison_heatmap_data.csv` ← NEW
- `glycopeptide_comparison_heatmap_summary.csv` ← NEW

Located in: `Results/Trace/`

---

## 📞 Support

If you find discrepancies:

1. Check this guide's verification examples
2. Run the Python verification script
3. Review the calculation formulas in the code
4. Open an issue with specific examples

**All data is fully traceable - if something looks wrong, you can verify it!**

---

## ✅ Summary

**Two trace files ensure complete traceability**:

1. **Summary file** (19 KB): Quick review with statistics
2. **Full data file** (55 KB): All individual sample values

**You can verify**:
- Every dot position
- Every intensity calculation
- Every statistical measure
- Every sorting decision
- Every grouping choice

**The visualization is fully reproducible from the trace data.**
