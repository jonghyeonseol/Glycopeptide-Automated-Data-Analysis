# Materials and Methods

## Glycan Structure Annotation

Glycan compositions were parsed from pGlyco output using regular expression pattern matching to extract monosaccharide counts: H (hexose), N (N-acetylhexosamine/HexNAc), A (N-acetylneuraminic acid/NeuAc/sialic acid), and F (fucose).

**Classification Criteria**:
- **High Mannose (HM)**: H ≥ 5, N = 2, no sialic acid (A), fucose (F), or N-acetylgalactosamine (G)
- **Fucosylated (F)**: Presence of fucose without sialic acid
- **Sialylated (S)**: Presence of sialic acid without fucose
- **Sialofucosylated (SF)**: Presence of both sialic acid and fucose
- **Complex/Hybrid (C/H)**: All other glycan structures not meeting HM criteria

These classifications follow standard glycobiology nomenclature and enable functional interpretation of glycosylation changes.

## Data Processing and Quality Control

Glycoproteomics data from 24 cancer and 23 normal tissue samples were integrated using a custom Python pipeline (pGlyco Auto Combine). Raw pGlyco output files were merged based on peptide-glycan composition combinations, creating a wide-format data matrix with 6434 unique glycopeptides across 47 samples.

**Quality Control Filtering**: Glycopeptides were required to be detected in ≥30% of samples in at least one group (cancer or normal) to ensure reliable quantification. This filter removed 4120 (64.0%) low-frequency glycopeptides, retaining 2314 glycopeptides for downstream analysis.

**Data Normalization**: Intensity values were normalized using Total Ion Current (TIC) normalization to account for sample-to-sample variation in total signal intensity. TIC-normalized values were log2-transformed (log2(x+1)) to stabilize variance and approximate normal distribution.

**Missing Data Handling**: Missing values (non-detected glycopeptides) were handled using the 'skipna' method, where statistical calculations excluded missing values rather than imputing zeros. This approach is appropriate for missing-not-at-random (MNAR) data common in mass spectrometry-based proteomics.

## Statistical Analysis and Biomarker Validation

**Multivariate Analysis**: Principal Component Analysis (PCA) was performed on log2-transformed, TIC-normalized intensity data (2314 glycopeptides) using scikit-learn (v1.3.0). Data were scaled using RobustScaler (median and IQR-based scaling) prior to PCA to reduce the influence of outliers. The first two principal components explained 11.2% and 4.5% of variance, respectively. Statistical significance of cancer vs. normal separation was assessed using permutation testing (1,000 permutations), yielding p < 0.0001.

**Supervised Classification**: Partial Least Squares Discriminant Analysis (PLS-DA) was performed using scikit-learn with 2 components. Variable Importance in Projection (VIP) scores were calculated to identify glycopeptides contributing most to group discrimination, with VIP > 1.0 indicating importance.

**Biomarker Stability Assessment**: Bootstrap resampling validation (1,000 iterations with replacement) was performed to assess biomarker stability. For each iteration, PLS-DA was re-fit and VIP scores recalculated. Glycopeptides with VIP > 1.0 in ≥80% of bootstrap iterations were classified as "stable biomarkers", identifying 368 robust candidates.

**Cross-Validation**: 10-fold stratified cross-validation assessed model generalizability, achieving 98.0% ± 2.0% accuracy and 1.000 ± 0.000 ROC-AUC, indicating minimal overfitting.

**Differential Expression Analysis**: Mann-Whitney U tests (non-parametric) compared cancer vs. normal groups for each glycopeptide. P-values were adjusted for multiple testing using the Benjamini-Hochberg false discovery rate (FDR) correction. Glycopeptides with FDR < 0.05 and |fold change| > 1.5 were considered significantly differentially expressed.

**Effect Size Calculation**: Cohen's d effect sizes were calculated for all glycopeptides to quantify biological significance independent of sample size. Effect sizes were classified as small (0.2 ≤ |d| < 0.5), medium (0.5 ≤ |d| < 0.8), or large (|d| ≥ 0.8) following standard conventions.

## Data Visualization

All visualizations were generated using Python (v3.9+) with matplotlib (v3.7.0), seaborn (v0.12.0), and custom plotting modules. Figures were rendered at 300 DPI resolution suitable for publication.

**Color Scheme**: Glycan type colors were chosen to reflect biological significance: green (high mannose - simple structures), blue (complex/hybrid - mature structures), pink (sialylated - charged modifications), orange (sialofucosylated - dual modifications), and red (fucosylated - core modifications). Cancer vs. normal comparisons used red (#E41A1C) and blue (#377EB8), respectively.

**Sample Size Annotation**: All comparative plots include sample size annotations (n= values) to meet journal requirements for transparent reporting.

**Statistical Annotations**: Box plots display statistical significance using Mann-Whitney U tests with symbols: * (p < 0.05), ** (p < 0.01), *** (p < 0.001). Error bars represent standard deviation unless otherwise noted.

**Software and Reproducibility**: All analyses were performed using Python 3.9+ with pandas (v2.0.0), NumPy (v1.24.0), SciPy (v1.10.0), and scikit-learn (v1.3.0). Complete analysis code and parameters are available in the project repository to ensure reproducibility.

## Data Availability

Raw mass spectrometry data and processed glycopeptide intensity matrices are available upon request. Analysis code and pipeline documentation are available at [repository URL].

---

*Methods text auto-generated by pGlyco Auto Combine on 2025-10-06*