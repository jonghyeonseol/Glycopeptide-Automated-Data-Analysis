# Manual Verification Guide - Quick & Simple

## 🎯 Purpose

This guide helps you **manually verify** that each dot in the visualization correctly represents your data.

---

## 📊 What Data is Available

### Two CSV Files for Verification:

1. **Summary File** (Recommended): `Results/Trace/glycopeptide_comparison_heatmap_summary.csv`
   - Contains all key statistics
   - 22 columns, easy to review in Excel
   - **Use this for quick spot-checks**

2. **Full Data File**: `Results/Trace/glycopeptide_comparison_heatmap_data.csv`
   - Contains all individual sample values (C1-C24, N1-N23)
   - 69 columns total
   - **Use this for detailed verification**

---

## ✅ Simple Verification Steps

### Step 1: Open in Excel/Spreadsheet

```
Open: Results/Trace/glycopeptide_comparison_heatmap_summary.csv
```

### Step 2: Find the Dot You Want to Verify

**In the visualization:**
- Count from left to right to get X position (starts at 0)
- Count from top to bottom to get Y position (starts at 0)

**In Excel:**
- Enable filters (Data → Filter)
- Filter `Plot_X_Position` = your X value
- Filter `Plot_Y_Position` = your Y value

**Example:**
- Third column from left, second row from top
- Filter: `Plot_X_Position = 2`, `Plot_Y_Position = 1`

### Step 3: Verify the Values

For the filtered row, check:

| Column | What to Verify |
|--------|----------------|
| `Peptide` | Correct peptide sequence? |
| `GlycanComposition` | Correct glycan (e.g., H(5)N(4)A(2))? |
| `GlycanTypeCategory` | Correct type (HM/F/S/SF/C/H)? |
| `Cancer_Mean` | Average of cancer samples |
| `Normal_Mean` | Average of normal samples |
| `Cancer_Dot_Plotted` | TRUE if you see left dot (○) |
| `Normal_Dot_Plotted` | TRUE if you see right dot (□) |

### Step 4: Check Dot Darkness

**Darker dots should have higher means:**

If visualizing shows:
- **Dark Cancer dot (○)** → `Cancer_Mean` should be high
- **Light Cancer dot (○)** → `Cancer_Mean` should be low
- **Dark Normal dot (□)** → `Normal_Mean` should be high
- **Light Normal dot (□)** → `Normal_Mean` should be low

**Transparency formula:**
```
Alpha = 0.3 + (Intensity / Max_Intensity) × 0.7

Where:
- Min darkness (most transparent) = 0.3
- Max darkness (fully opaque) = 1.0
```

The `Cancer_Alpha` and `Normal_Alpha` columns show these calculated values.

---

## 🔍 Quick Verification Examples

### Example 1: Verify Top-Left Dot

**Visual**: Top-left corner of heatmap
**Position**: X=0, Y=0

**Excel Steps:**
1. Filter `Plot_X_Position = 0`
2. Filter `Plot_Y_Position = 0`
3. Check `Peptide`, `GlycanComposition`
4. Check `VIP_Score` (should be highest, as it's top row)
5. Check `Cancer_Mean` and `Normal_Mean`

### Example 2: Verify Cancer-Specific Dots

**Look for**: Rows where only Cancer dot is shown

**Excel Steps:**
1. Filter `Cancer_Dot_Plotted = TRUE`
2. Filter `Normal_Dot_Plotted = FALSE`
3. Verify `Cancer_Mean > 0` and `Normal_Mean = 0` (or very low)

### Example 3: Verify Fold Changes

**Check**: Do high fold-change rows show darker Cancer dots?

**Excel Steps:**
1. Sort by `Log2_Fold_Change` (descending)
2. Top rows should have much higher `Cancer_Mean` than `Normal_Mean`
3. Visual should show darker ○ (Cancer) and lighter □ (Normal)

### Example 4: Verify Glycan Type Colors

**Check**: Do all glycans in the HM region have correct type?

**Excel Steps:**
1. Filter `GlycanTypeCategory = HM`
2. Check all rows have consecutive `Plot_X_Position` values
3. In visualization, these positions should be in the green-colored region

---

## 📋 Spot-Check Checklist

Perform these quick checks to validate data integrity:

### Data Consistency
- [ ] All rows have valid `Peptide` and `GlycanComposition`
- [ ] `Plot_X_Position` ranges from 0 to (total glycans - 1)
- [ ] `Plot_Y_Position` ranges from 0 to (total peptides - 1)
- [ ] No duplicate (X, Y) positions

### Intensity Values
- [ ] `Cancer_Mean >= 0` for all rows
- [ ] `Normal_Mean >= 0` for all rows
- [ ] If `Cancer_Dot_Plotted = TRUE`, then `Cancer_Mean > 0`
- [ ] If `Normal_Dot_Plotted = TRUE`, then `Normal_Mean > 0`

### Statistics
- [ ] `Fold_Change > 0` (or inf) for all rows
- [ ] `Cancer_Alpha` between 0 and 1
- [ ] `Normal_Alpha` between 0 and 1
- [ ] `Cancer_SampleCount <= 24`
- [ ] `Normal_SampleCount <= 23`

### Sorting and Grouping
- [ ] Top row (Y=0) has highest `PeptideVIP`
- [ ] `PeptideVIP` decreases as Y increases
- [ ] Glycans of same type have consecutive X positions

---

## 🔢 Manual Calculation Example

Want to verify the mean calculation? Use the full data file:

### Using Excel

1. **Open**: `Results/Trace/glycopeptide_comparison_heatmap_data.csv`

2. **Find your row** (e.g., row 2)

3. **Calculate Cancer Mean manually**:
   ```
   =AVERAGE(C1:Z1)  # Adjust columns for C1-C24
   ```

4. **Compare** with `Cancer_Mean` column
   - Should match (small rounding differences OK)

### Using Python

```python
import pandas as pd

# Load data
df = pd.read_csv('Results/Trace/glycopeptide_comparison_heatmap_data.csv')

# Select a row (e.g., first row)
row = df.iloc[0]

# Cancer samples
cancer_cols = [f'C{i}' for i in range(1, 25)]
cancer_values = row[cancer_cols].replace('', 0).astype(float)

# Calculate mean
manual_mean = cancer_values.mean()
saved_mean = row['Cancer_Mean']

print(f"Manual calculation: {manual_mean:,.0f}")
print(f"Saved in file:      {saved_mean:,.0f}")
print(f"Match: {abs(manual_mean - saved_mean) < 1}")
```

---

## 💡 Tips for Efficient Verification

### 1. Spot-Check Strategy

Don't verify every single dot - use statistical sampling:

- **Verify 5-10 random dots**: Should all be correct
- **Verify highest VIP scores**: Top 3 rows
- **Verify each glycan type**: 1-2 dots per type (HM, F, S, SF, C/H)
- **Verify extreme cases**: Highest/lowest fold changes

### 2. Visual Pattern Verification

Instead of checking individual numbers:

**Pattern**: Do all dots in HM region (green) show `GlycanTypeCategory = HM`?
- Filter by `GlycanTypeCategory = HM`
- Check `Plot_X_Position` values form contiguous block
- Visually confirm these positions are in green region

**Pattern**: Are high VIP peptides at the top?
- Sort by `Plot_Y_Position` (ascending)
- Check `PeptideVIP` decreases as you go down
- Visually confirm top rows in heatmap match top peptides in data

### 3. Cross-Reference with Source Data

**Verify end-to-end**:
1. Open `Results/integrated.csv` (source data)
2. Find a specific `Peptide` + `GlycanComposition` combination
3. Check C1-C24 and N1-N23 values
4. Calculate means manually
5. Find same combination in trace data
6. Compare values

---

## 🎓 Understanding the Columns

### Position Columns

| Column | Range | Meaning |
|--------|-------|---------|
| `Plot_X_Position` | 0 to 34 | Column in heatmap (left to right) |
| `Plot_Y_Position` | 0 to 19 | Row in heatmap (top to bottom) |

Position (0, 0) = **Top-left corner**
Position (34, 19) = **Bottom-right corner** (for 35 glycans × 20 peptides)

### Statistical Columns

| Column | Meaning |
|--------|---------|
| `Cancer_StdDev` | How much cancer samples vary |
| `Cancer_SampleCount` | How many cancer samples detected this glycopeptide |
| `Fold_Change` | Ratio of Cancer / Normal (>1 = higher in cancer) |
| `Log2_Fold_Change` | Log2 of fold change (standard in biology) |

---

## ❓ Common Questions

**Q: Why don't the individual sample values add up to the mean?**
A: The mean is calculated from C1-C24 (or N1-N23). Make sure you're averaging the correct columns.

**Q: Why is Normal_SampleCount sometimes 23, not 24?**
A: Sample N19 is missing from your dataset (common in real data).

**Q: Can I verify the VIP scores?**
A: Yes! Compare with `Results/vip_scores_all.csv` - should match for the same peptide-glycan combination.

**Q: What if I find a mismatch?**
A: Double-check your calculation methodology. If still mismatched, note the specific row and report it.

---

## ✅ Verification Complete!

Once you've spot-checked several dots and they all match:

**✓ Your visualization is reliable and accurate**
**✓ Every dot traces back to your original data**
**✓ You can confidently use this for analysis and publication**

---

## 📞 Need Help?

1. Check column descriptions above
2. Try the manual calculation example
3. Review the full `TRACE_DATA_GUIDE.md` for detailed information
4. Contact support with specific row numbers if issues persist

**Remember**: Complete traceability means you can verify ANY value in the visualization against the source data!
