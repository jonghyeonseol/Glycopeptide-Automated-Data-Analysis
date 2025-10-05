# Scientific Validity Review - pGlyco Auto Combine

**Date**: 2025-10-05
**Reviewer**: Code Analysis
**Principle**: "Every value used for visualization should meet high standards driven by scientific validity"

---

## Executive Summary

**Overall Grade: B+ (Good, with room for improvement)**

The pipeline employs scientifically sound core methods (PCA, PLS-DA, FDR correction), with some excellent choices (RobustScaler). However, there are important issues with missing data handling that could introduce systematic bias.

**Critical Finding**: Missing values are treated as zeros throughout, which systematically underestimates intensities and biases statistical comparisons.

---

## 1. Data Integration

### Pivot Table Aggregation
```python
df.pivot_table(values='IsotopeArea', index=['Peptide', 'GlycanComposition'],
               columns='Sample', aggfunc='sum')
```

**Validity**: ✅ **VALID**

**Rationale**:
- Sum aggregation is correct for technical replicates (multiple MS scans of same sample)
- Would use mean for biological replicates

**Recommendation**:
- Add logging to report if/when duplicates are summed
- Document whether data contains technical or biological replicates

---

## 2. Normalization Methods

### TIC (Total Ion Current) Normalization

```python
sample_sum = data_matrix.sum(axis=0)
target_median = np.median(sample_sum)
scaling_factors = target_median / sample_sum
normalized = data_matrix * scaling_factors
```

**Validity**: ✅ **VALID** (with important caveats)

**Scientific Basis**:
- Standard in proteomics/glycoproteomics
- Corrects for technical variation (loading differences, instrument variation)
- Assumes global intensity differences are technical, not biological

**Critical Assumption**: Total glycoprotein abundance is similar across groups
- **Valid for**: Comparing glycosylation patterns between groups
- **Invalid for**: Diseases with massive changes in total glycoprotein production

**For your cancer vs normal study**: ✅ Appropriate

### Log2 Transformation

```python
log_transformed = np.log2(normalized + 1)
```

**Validity**: ✅ **VALID**

**Scientific Basis**:
- Stabilizes variance (high-intensity values have high variance)
- Makes data more normally distributed (many tests assume normality)
- +1 pseudocount correctly handles zeros

**Standard Practice**: Yes (universal in omics)

---

## 3. Statistical Analysis

### PCA with RobustScaler

```python
scaler = RobustScaler()  # Uses median and IQR, not mean and SD
X_scaled = scaler.fit_transform(intensity_matrix)
pca = PCA(n_components=2)
```

**Validity**: ✅ **EXCELLENT**

**Why This is Better Than Standard**:
- **RobustScaler** uses median (50th percentile) and IQR (25th-75th percentile)
- **StandardScaler** uses mean and standard deviation (sensitive to outliers)
- Omics data often has outliers → RobustScaler is superior choice

**Recommendation**: This is actually BETTER than what many published papers use!

### PLS-DA + VIP Scores

```python
plsda = PLSRegression(n_components=2)
plsda.fit(X_scaled, y)

# VIP calculation (Wold et al., 2001)
vip[i] = sqrt(p * sum(w[i,j]^2 * s[j]) / total_s)
```

**Validity**: ✅ **VALID**

**Scientific Basis**:
- PLS-DA = Partial Least Squares Discriminant Analysis
- Supervised method (uses group labels to find discriminative features)
- VIP = Variable Importance in Projection (standard biomarker metric)

**Concern**: ⚠️ No cross-validation
- **Risk**: May overfit to current samples
- **Impact**: Selected biomarkers might not generalize to new samples
- **Recommendation**: Add k-fold cross-validation or permutation testing

**For your current analysis**: Acceptable for exploratory research, but would strengthen conclusions with validation

---

## 4. Critical Issue: Missing Data Handling

### Current Approach

```python
replace_empty_with_zero(df[sample_cols])
```

**Validity**: ⚠️ **QUESTIONABLE - SYSTEMATIC BIAS**

**Problem**:
Missing values in mass spectrometry can mean:
1. **True biological absence** (glycopeptide not present)
2. **Below detection limit** (present but too low to detect)
3. **Technical failure** (sample not measured)

Treating all as zero conflates these meanings.

**Impact on Results**:
```python
# Example:
Cancer samples: [100, 200, 150, 0, 0, 0]  # 3 real values, 3 missing → treated as zero
Normal samples: [120, 180, 140, 110, 130, 125]  # all measured

# Current calculation (missing = 0):
Cancer mean = (100+200+150+0+0+0)/6 = 75
Normal mean = (120+180+140+110+130+125)/6 = 134
Fold change = 75/134 = 0.56  # Appears DOWN in cancer

# Correct calculation (skip missing):
Cancer mean = (100+200+150)/3 = 150
Normal mean = (120+180+140+110+130+125)/6 = 134
Fold change = 150/134 = 1.12  # Actually UP in cancer!
```

**Systematic Bias**: Always underestimates group with more missing values

**Recommended Fix**:
```python
# Option 1: Use only non-missing values
cancer_values = df[cancer_cols].replace('', np.nan).astype(float)
cancer_mean = cancer_values.mean(axis=1, skipna=True)

# Option 2: Impute using minimum detected value
min_value = df[df > 0].min().min()
imputed = df.replace('', min_value)
```

---

## 5. Fold Change Calculation

### Current Method

```python
fold_change = cancer_mean / normal_mean
```

**Validity**: ⚠️ **NEEDS IMPROVEMENT**

**Issues**:
1. Division by zero when normal_mean = 0
2. Asymmetric (FC of 2 vs 0.5 treated differently)
3. Sensitive to small denominators

**Current Fix**: Returns `np.inf` if normal_mean = 0
**Problem**: Infinity is not statistically meaningful

**Better Approach**:
```python
# Log2 fold change with pseudocount
log2_fc = np.log2((cancer_mean + 1) / (normal_mean + 1))

# Advantages:
# - Symmetric: log2(2) = 1, log2(0.5) = -1
# - No division by zero
# - Standard in genomics/proteomics
```

---

## 6. Statistical Testing

### Volcano Plot: Mann-Whitney U Test + FDR

```python
statistic, p_value = stats.mannwhitneyu(cancer_nonzero, normal_nonzero)
_, fdr_values, _, _ = multipletests(p_values, method='fdr_bh')
```

**Validity**: ✅ **VALID**

**Scientific Basis**:
- **Mann-Whitney U**: Non-parametric test (doesn't assume normal distribution)
  - Appropriate for omics data (often non-normal)
  - Tests if two groups have different distributions
- **FDR Correction** (Benjamini-Hochberg):
  - Gold standard for multiple testing correction
  - Controls False Discovery Rate (expected proportion of false positives)
  - Better than Bonferroni (less conservative, more power)

**Sample Size Requirement**: n ≥ 3 per group
- **Current implementation**: Correctly checks this
- **Recommendation**: Warn if n < 5 (under-powered)

---

## 7. Visualization-Specific Validity

### PCA Plot

**Variance Explained Labels**: ✅ CORRECT
- Directly from PCA model (no calculations)
- Represents proportion of total variance captured by each PC

**95% Confidence Ellipses**: ✅ VALID
- Based on covariance matrix of each group
- Visualization aid, not statistical test

### Heatmap Clustering

**Current**: linkage='average', metric='euclidean'

**Validity**: ✅ VALID (but alternatives exist)
- **Average linkage**: Reasonable compromise between single and complete
- **Euclidean distance**: Assumes features on similar scales (true after normalization)
- **Alternative**: 'correlation' distance might be better for expression patterns

### Volcano Plot Thresholds

**Current**: FDR < 0.05, |FC| > 1.5

**Validity**: ✅ REASONABLE
- **FDR < 0.05**: Standard significance level
- **FC > 1.5**: Arbitrary but common (50% change)
- **Note**: Thresholds should be justified by biological context

---

## Priority Recommendations

### HIGH PRIORITY (Affects Scientific Validity)

#### 1. Fix Missing Data Handling
**Impact**: Systematic bias in all mean calculations, fold changes, and group comparisons

**Implementation**:
```python
# In utils.py, modify replace_empty_with_zero to:
def calculate_group_statistics(df, sample_cols):
    """Calculate statistics using only non-missing values"""
    values = df[sample_cols].replace('', np.nan).astype(float)

    return {
        'mean': values.mean(axis=1, skipna=True),
        'std': values.std(axis=1, skipna=True),
        'count': values.count(axis=1),  # Number of non-missing
        'missing_rate': values.isna().sum(axis=1) / len(sample_cols)
    }
```

#### 2. Improve Fold Change Calculation
**Impact**: More robust and interpretable

```python
def calculate_log2_fold_change(cancer_mean, normal_mean, pseudocount=1):
    """Log2 fold change with pseudocount"""
    return np.log2((cancer_mean + pseudocount) / (normal_mean + pseudocount))
```

### MEDIUM PRIORITY (Improves Robustness)

#### 3. Add VIP Cross-Validation
**Impact**: Identifies more stable biomarkers

```python
def validate_vip_scores(X, y, n_iterations=100):
    """Bootstrap validation of VIP scores"""
    vip_matrix = []
    for i in range(n_iterations):
        # Resample with replacement
        indices = np.random.choice(len(X), len(X), replace=True)
        X_boot = X[indices]
        y_boot = y[indices]

        # Fit PLS-DA
        pls = PLSRegression(n_components=2)
        pls.fit(X_boot, y_boot)

        # Calculate VIP
        vip = calculate_vip_scores(X_boot, y_boot, pls)
        vip_matrix.append(vip)

    # Return mean and CI
    vip_mean = np.mean(vip_matrix, axis=0)
    vip_ci_lower = np.percentile(vip_matrix, 2.5, axis=0)
    vip_ci_upper = np.percentile(vip_matrix, 97.5, axis=0)

    return vip_mean, vip_ci_lower, vip_ci_upper
```

#### 4. Add Sample Size Warnings
```python
def validate_statistical_power(cancer_n, normal_n, min_n=5):
    """Warn if sample size is insufficient"""
    if cancer_n < min_n:
        logger.warning(f"Cancer group has only {cancer_n} samples (< {min_n}). "
                      "Statistical tests may be under-powered.")
    if normal_n < min_n:
        logger.warning(f"Normal group has only {normal_n} samples (< {min_n}). "
                      "Statistical tests may be under-powered.")
```

### LOW PRIORITY (Nice to Have)

#### 5. Document All Assumptions
Add to README or methods document:
- TIC normalization assumes similar total glycoprotein abundance
- Log transform assumes multiplicative model
- Missing data handling strategy
- Statistical test assumptions

#### 6. Alternative Normalization Options
```python
# Quantile normalization
from sklearn.preprocessing import quantile_transform
normalized = quantile_transform(data, axis=0)

# Variance Stabilizing Normalization (VSN)
# More complex, requires external package
```

---

## Overall Assessment

### Strengths ✅

1. **Excellent Scaling Choice**: RobustScaler (better than standard)
2. **Proper Multiple Testing**: FDR correction (Benjamini-Hochberg)
3. **Appropriate Tests**: Mann-Whitney U (non-parametric)
4. **Standard Methods**: PCA, PLS-DA widely accepted
5. **Good Documentation**: Code is well-commented

### Weaknesses ⚠️

1. **Missing Data as Zero**: Systematic bias (HIGH PRIORITY)
2. **No VIP Validation**: May select unstable biomarkers
3. **Fold Change Method**: Could be more robust
4. **Limited Warnings**: No alerts for under-powered tests

### Recommended Actions

**Immediate** (Before publication):
- Fix missing data handling
- Switch to log2 fold change

**Short-term** (Strengthen analysis):
- Add VIP cross-validation
- Add sample size warnings

**Long-term** (Future versions):
- Alternative normalization options
- Comprehensive methods documentation

---

## Conclusion

The pGlyco Auto Combine pipeline uses **scientifically valid core methods** and makes some **excellent technical choices** (RobustScaler). However, the **treatment of missing data as zeros is a significant scientific concern** that could bias results.

**Recommendation**: Implement HIGH PRIORITY fixes before using for publication-quality analysis. The MEDIUM and LOW priority items would strengthen the analysis but are not critical for validity.

**Overall Scientific Validity**: **B+** (Good, with important caveats)

---

**Report Generated**: 2025-10-05
**Next Steps**: Implement recommended fixes based on priority
