# Phase 3.1 Dry Run Analysis - FDR Correction Helper

## Executive Summary

**Target**: Extract FDR correction logic from `plot_boxplot()` and `plot_boxplot_extended()`
**Method**: Create `_perform_fdr_correction_and_plot_brackets()` helper
**Expected Impact**: ~78 lines net reduction (220 duplicated → 120 helper + 20 calls)

---

## 1. Current Code Analysis - Side-by-Side Comparison

### Method 1: `plot_boxplot()` (Lines 134-242)

```python
# Step 1: Collect all p-values and effect sizes
p_values_dict = {}  # {glycan_type: p_value}
cohens_d_dict = {}  # {glycan_type: effect_size}

for glycan_type in existing_types:                                    # ← DIFFERENCE: variable name
    # Get data for each group
    cancer_data = boxplot_data[
        (boxplot_data['GlycanType'] == glycan_type) &                # ← DIFFERENCE: column name
        (boxplot_data['Group'] == 'Cancer')
    ]['Intensity'].values

    normal_data = boxplot_data[
        (boxplot_data['GlycanType'] == glycan_type) &                # ← DIFFERENCE: column name
        (boxplot_data['Group'] == 'Normal')
    ]['Intensity'].values

    # Skip if either group has insufficient data
    if len(cancer_data) < 3 or len(normal_data) < 3:
        p_values_dict[glycan_type] = np.nan                           # ← DIFFERENCE: key name
        cohens_d_dict[glycan_type] = np.nan                           # ← DIFFERENCE: key name
        logger.warning(f"{glycan_type}: Insufficient data ...")       # ← DIFFERENCE: log text
        continue

    # Perform Mann-Whitney U test (non-parametric)
    try:
        statistic, p_value = stats.mannwhitneyu(cancer_data, normal_data, alternative='two-sided')
        cohens_d = calculate_cohens_d(cancer_data, normal_data)

        p_values_dict[glycan_type] = p_value                          # ← DIFFERENCE: key name
        cohens_d_dict[glycan_type] = cohens_d                         # ← DIFFERENCE: key name

    except Exception as e:
        logger.warning(f"Statistical test failed for {glycan_type}: {str(e)}")
        p_values_dict[glycan_type] = np.nan
        cohens_d_dict[glycan_type] = np.nan

# Step 2: Apply FDR correction (Benjamini-Hochberg)
valid_glycan_types = [gt for gt, p in p_values_dict.items() if not np.isnan(p)]
valid_p_values = [p_values_dict[gt] for gt in valid_glycan_types]

fdr_dict = {}
if len(valid_p_values) > 0:
    # Apply Benjamini-Hochberg FDR correction
    reject, fdr_values, _, _ = multipletests(valid_p_values, method='fdr_bh')
    fdr_dict = dict(zip(valid_glycan_types, fdr_values))

    # Log FDR correction results
    n_significant = sum(fdr < 0.05 for fdr in fdr_values)
    logger.info(f"FDR correction applied: {len(valid_p_values)} tests, "  # ← DIFFERENCE: log text
               f"{n_significant} significant (FDR < 0.05)")

    for gt, raw_p, fdr in zip(valid_glycan_types, valid_p_values, fdr_values):
        logger.info(f"  {gt}: p={raw_p:.4f} → FDR={fdr:.4f}")
else:
    logger.warning("No valid p-values for FDR correction")           # ← DIFFERENCE: log text

# Step 3: Plot with FDR-corrected significance markers
y_max = boxplot_data['Intensity'].max()
y_range = boxplot_data['Intensity'].max() - boxplot_data['Intensity'].min()

for i, glycan_type in enumerate(existing_types):                     # ← DIFFERENCE: variable name
    # Get FDR-corrected p-value
    fdr = fdr_dict.get(glycan_type, np.nan)                          # ← DIFFERENCE: key name
    cohens_d = cohens_d_dict.get(glycan_type, np.nan)                # ← DIFFERENCE: key name

    # Skip if no valid FDR value
    if np.isnan(fdr):
        continue

    # Determine significance level using FDR (not raw p-value)
    if fdr < 0.001:
        sig_marker = '***'
    elif fdr < 0.01:
        sig_marker = '**'
    elif fdr < 0.05:
        sig_marker = '*'
    else:
        sig_marker = 'ns'

    # ✨ ENHANCED: Add significance bracket if significant
    if sig_marker != 'ns':
        # Calculate x positions for the connecting line
        n_types = len(existing_types)                                 # ← DIFFERENCE: variable name
        x_offset = (i - (n_types - 1) / 2) * (BOXPLOT_WIDTH / n_types)

        x1 = 0 + x_offset  # Cancer position
        x2 = 1 + x_offset  # Normal position

        y_position = y_max + y_range * 0.05 * (1 + i * 0.3)

        # Format annotation text with effect size
        if not np.isnan(cohens_d):
            annotation_text = f"{sig_marker}\n(d={cohens_d:.2f})"
        else:
            annotation_text = sig_marker

        # ✨ Use enhanced statistical bracket (rounded ends, fancy box)
        enhance_statistical_bracket(
            ax, x1, x2, y_position,
            text=annotation_text,
            color='black',
            fontsize=ANNOTATION_SIZE
        )

        logger.info(
            f"{glycan_type}: Cancer vs Normal FDR={fdr:.4f} "        # ← DIFFERENCE: log text
            f"({sig_marker}), Cohen's d={cohens_d:.3f}"
        )
```

### Method 2: `plot_boxplot_extended()` (Lines 319-434)

**IDENTICAL STRUCTURE** with only these differences:
- `glycan_type` → `category`
- `existing_types` → `existing_categories`
- `'GlycanType'` → `'ExtendedCategory'`
- Log prefix: `"FDR correction applied"` → `"FDR correction applied (extended)"`

---

## 2. Proposed Helper Method Design

### Signature

```python
def _perform_fdr_correction_and_plot_brackets(
    self,
    ax: plt.Axes,
    boxplot_data: pd.DataFrame,
    category_column: str,
    existing_categories: list,
    log_context: str = ""
) -> None:
    """
    Perform FDR-corrected statistical testing and plot significance brackets

    This method consolidates the FDR correction workflow:
    1. Collect p-values and effect sizes for each category
    2. Apply Benjamini-Hochberg FDR correction
    3. Plot significance brackets with effect sizes
    4. Log results

    Args:
        ax: Matplotlib axes for plotting brackets
        boxplot_data: Long-format DataFrame with columns:
            - {category_column}: Category identifier (e.g., 'GlycanType', 'ExtendedCategory')
            - 'Group': 'Cancer' or 'Normal'
            - 'Intensity': Numerical values for comparison
        category_column: Name of column containing categories
        existing_categories: Ordered list of categories to test
        log_context: Prefix for logging (e.g., "extended" for extended boxplot)

    Returns:
        None (modifies ax in-place, logs results)

    Notes:
        - Uses Mann-Whitney U test (non-parametric)
        - Calculates Cohen's d effect size
        - Applies Benjamini-Hochberg FDR correction for multiple testing
        - Only plots brackets for FDR < 0.05
        - Requires ≥3 samples per group for testing

    Scientific Integrity:
        - All statistical calculations are deterministic
        - FDR correction is standard Benjamini-Hochberg
        - Effect sizes use pooled standard deviation
        - No data transformation occurs
    """
```

### Implementation

```python
def _perform_fdr_correction_and_plot_brackets(
    self,
    ax: plt.Axes,
    boxplot_data: pd.DataFrame,
    category_column: str,
    existing_categories: list,
    log_context: str = ""
) -> None:
    """[Full docstring as above]"""

    # Step 1: Collect all p-values and effect sizes
    p_values_dict = {}  # {category: p_value}
    cohens_d_dict = {}  # {category: effect_size}

    for category in existing_categories:
        # Get data for each group
        cancer_data = boxplot_data[
            (boxplot_data[category_column] == category) &
            (boxplot_data['Group'] == 'Cancer')
        ]['Intensity'].values

        normal_data = boxplot_data[
            (boxplot_data[category_column] == category) &
            (boxplot_data['Group'] == 'Normal')
        ]['Intensity'].values

        # Skip if either group has insufficient data
        if len(cancer_data) < 3 or len(normal_data) < 3:
            p_values_dict[category] = np.nan
            cohens_d_dict[category] = np.nan
            logger.warning(
                f"{category}: Insufficient data (Cancer: {len(cancer_data)}, "
                f"Normal: {len(normal_data)})"
            )
            continue

        # Perform Mann-Whitney U test (non-parametric)
        try:
            statistic, p_value = stats.mannwhitneyu(
                cancer_data, normal_data, alternative='two-sided'
            )
            cohens_d = calculate_cohens_d(cancer_data, normal_data)

            p_values_dict[category] = p_value
            cohens_d_dict[category] = cohens_d

        except Exception as e:
            logger.warning(f"Statistical test failed for {category}: {str(e)}")
            p_values_dict[category] = np.nan
            cohens_d_dict[category] = np.nan

    # Step 2: Apply FDR correction (Benjamini-Hochberg)
    valid_categories = [cat for cat, p in p_values_dict.items() if not np.isnan(p)]
    valid_p_values = [p_values_dict[cat] for cat in valid_categories]

    fdr_dict = {}
    if len(valid_p_values) > 0:
        # Apply Benjamini-Hochberg FDR correction
        reject, fdr_values, _, _ = multipletests(valid_p_values, method='fdr_bh')
        fdr_dict = dict(zip(valid_categories, fdr_values))

        # Log FDR correction results
        n_significant = sum(fdr < 0.05 for fdr in fdr_values)
        log_prefix = f"FDR correction applied{' (' + log_context + ')' if log_context else ''}"
        logger.info(
            f"{log_prefix}: {len(valid_p_values)} tests, "
            f"{n_significant} significant (FDR < 0.05)"
        )

        for cat, raw_p, fdr in zip(valid_categories, valid_p_values, fdr_values):
            logger.info(f"  {cat}: p={raw_p:.4f} → FDR={fdr:.4f}")
    else:
        log_suffix = f" ({log_context})" if log_context else ""
        logger.warning(f"No valid p-values for FDR correction{log_suffix}")

    # Step 3: Plot with FDR-corrected significance markers
    y_max = boxplot_data['Intensity'].max()
    y_range = boxplot_data['Intensity'].max() - boxplot_data['Intensity'].min()

    for i, category in enumerate(existing_categories):
        # Get FDR-corrected p-value
        fdr = fdr_dict.get(category, np.nan)
        cohens_d = cohens_d_dict.get(category, np.nan)

        # Skip if no valid FDR value
        if np.isnan(fdr):
            continue

        # Determine significance level using FDR (not raw p-value)
        if fdr < 0.001:
            sig_marker = '***'
        elif fdr < 0.01:
            sig_marker = '**'
        elif fdr < 0.05:
            sig_marker = '*'
        else:
            sig_marker = 'ns'

        # ✨ ENHANCED: Add significance bracket if significant
        if sig_marker != 'ns':
            # Calculate x positions for the connecting line
            n_categories = len(existing_categories)
            x_offset = (i - (n_categories - 1) / 2) * (BOXPLOT_WIDTH / n_categories)

            x1 = 0 + x_offset  # Cancer position
            x2 = 1 + x_offset  # Normal position

            y_position = y_max + y_range * 0.05 * (1 + i * 0.3)

            # Format annotation text with effect size
            if not np.isnan(cohens_d):
                annotation_text = f"{sig_marker}\n(d={cohens_d:.2f})"
            else:
                annotation_text = sig_marker

            # ✨ Use enhanced statistical bracket (rounded ends, fancy box)
            enhance_statistical_bracket(
                ax, x1, x2, y_position,
                text=annotation_text,
                color='black',
                fontsize=ANNOTATION_SIZE
            )

            logger.info(
                f"{category}: Cancer vs Normal FDR={fdr:.4f} "
                f"({sig_marker}), Cohen's d={cohens_d:.3f}"
            )
```

---

## 3. Modified Method Calls

### Method 1: `plot_boxplot()` (After Refactoring)

**BEFORE** (Lines 134-242, ~108 lines):
```python
# [All the FDR correction code shown above]
```

**AFTER** (~5 lines):
```python
# Apply FDR correction and plot significance brackets
self._perform_fdr_correction_and_plot_brackets(
    ax=ax,
    boxplot_data=boxplot_data,
    category_column='GlycanType',
    existing_categories=existing_types,
    log_context=""  # No prefix for standard boxplot
)
```

### Method 2: `plot_boxplot_extended()` (After Refactoring)

**BEFORE** (Lines 319-434, ~116 lines):
```python
# [All the FDR correction code, nearly identical]
```

**AFTER** (~5 lines):
```python
# Apply FDR correction and plot significance brackets
self._perform_fdr_correction_and_plot_brackets(
    ax=ax,
    boxplot_data=boxplot_data,
    category_column='ExtendedCategory',
    existing_categories=existing_categories,
    log_context="extended"  # Add "extended" prefix to logs
)
```

---

## 4. Data Integrity Guarantees

### Statistical Correctness

✅ **Mann-Whitney U Test**:
- Input: `cancer_data`, `normal_data` (unchanged)
- Parameters: `alternative='two-sided'` (unchanged)
- Output: `p_value` (deterministic for same input)

✅ **Cohen's d Calculation**:
- Formula: `(mean1 - mean2) / pooled_std` (unchanged)
- Pooled std calculation: Standard formula (unchanged)
- Output: Deterministic for same input

✅ **FDR Correction**:
- Method: Benjamini-Hochberg (`method='fdr_bh'`)
- Library: `statsmodels.stats.multitest.multipletests`
- Output: Deterministic for same input p-values

✅ **Bracket Positioning**:
- Formula: `y_position = y_max + y_range * 0.05 * (1 + i * 0.3)` (unchanged)
- X offset: `(i - (n - 1) / 2) * (BOXPLOT_WIDTH / n)` (unchanged)
- Output: Pixel-perfect identical

### Parameterization Safety

✅ **No Hard-Coding**:
- Column names passed as parameters (not hard-coded)
- Category lists passed as parameters
- All constants from `plot_config.py` (unchanged)

✅ **No Over-Fitting**:
- Works with ANY categorical column
- Works with ANY list of categories
- No assumptions about data structure beyond:
  - Column exists in DataFrame
  - Contains 'Cancer'/'Normal' groups
  - Has 'Intensity' values

✅ **Backward Compatible**:
- Public method signatures unchanged
- Output PNG files identical
- Trace CSV files identical

---

## 5. Testing Strategy

### Pre-Refactoring Baseline

1. **Run current pipeline**:
```bash
python3 main.py
```

2. **Save baseline outputs**:
```bash
cp Results/boxplot_glycan_types.png Results/baseline_boxplot_glycan_types.png
cp Results/boxplot_extended_categories.png Results/baseline_boxplot_extended_categories.png
cp Results/Trace/boxplot_glycan_types_data.csv Results/Trace/baseline_boxplot_glycan_types_data.csv
cp Results/Trace/boxplot_extended_categories_data.csv Results/Trace/baseline_boxplot_extended_categories_data.csv
```

3. **Extract statistical values**:
```bash
grep "FDR=" logs/latest.log > baseline_fdr_values.txt
```

### Post-Refactoring Validation

1. **Syntax check**:
```bash
python3 -m py_compile src/plots/boxplot.py
```

2. **Run refactored pipeline**:
```bash
python3 main.py
```

3. **Compare PNG files** (visual inspection):
```bash
# Open side-by-side comparison
open Results/baseline_boxplot_glycan_types.png Results/boxplot_glycan_types.png
open Results/baseline_boxplot_extended_categories.png Results/boxplot_extended_categories.png
```

4. **Compare CSV files** (numerical validation):
```bash
diff Results/Trace/baseline_boxplot_glycan_types_data.csv Results/Trace/boxplot_glycan_types_data.csv
diff Results/Trace/baseline_boxplot_extended_categories_data.csv Results/Trace/boxplot_extended_categories_data.csv
# Expected output: No differences (files identical)
```

5. **Compare statistical values**:
```bash
grep "FDR=" logs/latest.log > refactored_fdr_values.txt
diff baseline_fdr_values.txt refactored_fdr_values.txt
# Expected output: Identical FDR values
```

### Success Criteria

✅ **PNG files**: Visually identical (pixel-perfect if possible)
✅ **CSV files**: `diff` shows no differences
✅ **FDR values**: Identical p-values and FDR corrections in logs
✅ **Bracket positions**: Same coordinates in both versions
✅ **Effect sizes**: Identical Cohen's d values
✅ **File size**: ~78 lines reduction (951 → 873)

### Failure Rollback Plan

If ANY test fails:
1. ❌ **STOP** - Do not proceed
2. 🔍 **Analyze** - Identify what changed
3. 🔄 **Rollback** - `git checkout src/plots/boxplot.py`
4. 📝 **Document** - Record what went wrong
5. 🔧 **Fix** - Adjust helper method design
6. ♻️ **Retry** - Test again with corrected approach

---

## 6. Over-Fitting Prevention

### What We're Avoiding

❌ **Hard-coded assumptions**:
```python
# BAD: Assumes specific column names
if 'GlycanType' in df.columns:
    # ...

# GOOD: Accepts column name as parameter
def helper(category_column: str):
    # Works with any column
```

❌ **Data-specific logic**:
```python
# BAD: Assumes specific categories
if category in ['Non', 'Sialylated', 'Fucosylated', 'Both']:
    # ...

# GOOD: Works with any category list
for category in existing_categories:
    # Generic iteration
```

❌ **Magic numbers**:
```python
# BAD: Hard-coded thresholds
if len(cancer_data) < 3:  # Why 3?

# GOOD: Documented scientific rationale
MIN_SAMPLES_FOR_TEST = 3  # Mann-Whitney U requires ≥3 samples
if len(cancer_data) < MIN_SAMPLES_FOR_TEST:
```

### What We're Preserving

✅ **Generic algorithms**: Mann-Whitney U, FDR correction (standard methods)
✅ **Configurable parameters**: Column names, category lists (passed in)
✅ **Scientific constants**: Significance thresholds (0.001, 0.01, 0.05 - standard)
✅ **Layout constants**: From `plot_config.py` (centralized)

---

## 7. Expected Changes Summary

### File Structure

**BEFORE**:
```
calculate_cohens_d()                          40 lines
plot_boxplot()                               185 lines (including 108 lines FDR)
plot_boxplot_extended()                      193 lines (including 116 lines FDR)
[other methods...]                           533 lines
Total:                                       951 lines
```

**AFTER**:
```
calculate_cohens_d()                          40 lines
_perform_fdr_correction_and_plot_brackets()  120 lines (NEW)
plot_boxplot()                                82 lines (108 → 5 lines FDR)
plot_boxplot_extended()                       82 lines (116 → 5 lines FDR)
[other methods...]                           533 lines
Total:                                       873 lines
```

**Net Change**: 951 → 873 = **-78 lines (8.2% reduction)**

### Code Quality Improvements

✅ **DRY Principle**: Eliminated 220 lines of duplication
✅ **Single Source of Truth**: FDR logic in one place
✅ **Maintainability**: Bug fixes apply to both methods
✅ **Testability**: Can unit test FDR helper independently
✅ **Readability**: Method intent clearer (delegate to helper)

---

## 8. Risk Mitigation

### High-Risk Areas

🔴 **Statistical calculations**:
- Mitigation: Use exact same functions (`mannwhitneyu`, `multipletests`, `calculate_cohens_d`)
- Validation: Compare log outputs line-by-line

🔴 **Bracket positioning**:
- Mitigation: Preserve exact formulas for x_offset, y_position
- Validation: Visual comparison of PNG files

🟡 **Logging changes**:
- Mitigation: Parameterize log context, preserve format
- Validation: Compare log outputs (except whitespace)

### Low-Risk Areas

🟢 **Method signatures**: Public APIs unchanged
🟢 **Data flow**: Same input → same output
🟢 **Dependencies**: No new imports required

---

## 9. Next Steps Decision

**OPTION 1: Proceed with Phase 3.1**
- ✅ You've reviewed the dry run analysis
- ✅ You understand the changes
- ✅ You agree with testing strategy
- ✅ You're comfortable with risk level

**OPTION 2: Request modifications**
- ❓ Any concerns about helper design?
- ❓ Want different parameterization?
- ❓ Need more detailed testing?
- ❓ Prefer different approach?

**OPTION 3: Abort and try Option A**
- ❌ Too risky to start with FDR correction
- ❌ Prefer safer classification pattern first
- ❌ Want more incremental approach

---

## 10. Your Decision Required

Please advise:

1. ✅ **Approve Phase 3.1 implementation?** (Yes/No/Modifications needed)
2. 📋 **Any concerns about the helper method design?**
3. 🧪 **Is the testing strategy sufficient?**
4. 🔍 **Do you want to see actual code diffs before implementation?**
5. 📊 **Should I extract baseline outputs first, then implement?**

---

*Phase 3.1 Dry Run Analysis*
*Created: 2025-10-14*
*Purpose: Pre-implementation review for data integrity assurance*
