"""
Constants Module for pGlyco Auto Combine
Centralized configuration values and constants
"""

from typing import List

# ==============================================================================
# Column Names
# ==============================================================================

# Required columns for input CSV files
REQUIRED_INPUT_COLUMNS: List[str] = [
    'Peptide',
    'GlycanComposition',
    'IsotopeArea',
    'Proteins'
]

# Metadata columns (non-sample intensity columns)
METADATA_COLUMNS: List[str] = [
    'Peptide',
    'GlycanComposition',
    'Proteins',
    'Sialylation',
    'Fucosylation',
    'IsSialylated',
    'IsFucosylated',
    'SialylationCount',
    'FucosylationCount',
    'GlycanType',
    'HighMannose',
    'ComplexHybrid',
    'IsHighMannose',
    'IsComplexHybrid',
    'N_count',
    'PrimaryClassification',
    'SecondaryClassification',
    'GlycanTypeCategory'
]

# ==============================================================================
# Sample Groups
# ==============================================================================

CANCER_PREFIX: str = 'C'
NORMAL_PREFIX: str = 'N'

GROUP_CANCER: str = 'Cancer'
GROUP_NORMAL: str = 'Normal'

# ==============================================================================
# Monosaccharide Symbols
# ==============================================================================

MONOSACCHARIDE_H: str = 'H'  # Hexose
MONOSACCHARIDE_N: str = 'N'  # HexNAc
MONOSACCHARIDE_A: str = 'A'  # NeuAc (Sialic acid)
MONOSACCHARIDE_F: str = 'F'  # Fucose
MONOSACCHARIDE_G: str = 'G'  # NeuGc

# Default markers
DEFAULT_SIALYLATION_MARKER: str = 'A'
DEFAULT_FUCOSYLATION_MARKER: str = 'F'

# Monosaccharide sorting order for glycan composition
MONOSACCHARIDE_ORDER: dict = {
    'H': 0,
    'N': 1,
    'A': 2,
    'F': 3,
    'G': 4
}

# ==============================================================================
# Glycan Type Categories
# ==============================================================================

GLYCAN_TYPE_HM: str = 'HM'          # High-mannose
GLYCAN_TYPE_F: str = 'F'            # Fucosylated
GLYCAN_TYPE_S: str = 'S'            # Sialylated
GLYCAN_TYPE_SF: str = 'SF'          # Sialofucosylated
GLYCAN_TYPE_CH: str = 'C/H'         # Complex/Hybrid
GLYCAN_TYPE_UNKNOWN: str = 'Unknown'

# Glycan type category order (for consistent visualization)
GLYCAN_TYPE_ORDER: List[str] = [
    GLYCAN_TYPE_HM,
    GLYCAN_TYPE_F,
    GLYCAN_TYPE_S,
    GLYCAN_TYPE_SF,
    GLYCAN_TYPE_CH
]

# ==============================================================================
# Classification Labels
# ==============================================================================

# Primary classifications
PRIMARY_TRUNCATED: str = 'Truncated'
PRIMARY_HIGH_MANNOSE: str = 'High Mannose'
PRIMARY_COMPLEX_HYBRID: str = 'ComplexHybrid'
PRIMARY_OUTLIER: str = 'Outlier'
PRIMARY_UNKNOWN: str = 'Unknown'

# Secondary classifications
SECONDARY_TRUNCATED: str = 'Truncated'
SECONDARY_HIGH_MANNOSE: str = 'High Mannose'
SECONDARY_COMPLEX_HYBRID: str = 'Complex/Hybrid'
SECONDARY_FUCOSYLATED: str = 'Fucosylated'
SECONDARY_SIALYLATED: str = 'Sialylated'
SECONDARY_SIALOFUCOSYLATED: str = 'Sialofucosylated'
SECONDARY_OUTLIER: str = 'Outlier'
SECONDARY_UNKNOWN: str = 'Unknown'

# Legacy glycan types
LEGACY_NON: str = 'Non'
LEGACY_SIALYLATED: str = 'Sialylated'
LEGACY_FUCOSYLATED: str = 'Fucosylated'
LEGACY_BOTH: str = 'Both'

# ==============================================================================
# High-Mannose Criteria
# ==============================================================================

HIGH_MANNOSE_MIN_H: int = 5  # Minimum H count
HIGH_MANNOSE_EXACT_N: int = 2  # Exact N count required

# ==============================================================================
# Complex/Hybrid Criteria
# ==============================================================================

COMPLEX_HYBRID_MIN_N: int = 3  # Minimum N count

# ==============================================================================
# Analysis Parameters
# ==============================================================================

DEFAULT_PCA_COMPONENTS: int = 2
DEFAULT_PLSDA_COMPONENTS: int = 2
DEFAULT_LOG_TRANSFORM: bool = True
DEFAULT_SIGNIFICANCE_ALPHA: float = 0.05

# Log transformation pseudocount
LOG_TRANSFORM_PSEUDOCOUNT: float = 1.0

# ==============================================================================
# Visualization Parameters
# ==============================================================================

DEFAULT_DPI: int = 300

# Default figure sizes (width, height in inches)
FIGSIZE_PCA: tuple = (10, 8)
FIGSIZE_BOXPLOT: tuple = (12, 6)
FIGSIZE_HEATMAP: tuple = (14, 10)
FIGSIZE_COMPARISON_HEATMAP: tuple = (24, 16)

# Heatmap parameters
HEATMAP_TOP_N_DEFAULT: int = 50
COMPARISON_MAX_PEPTIDES_DEFAULT: int = 50
COMPARISON_MAX_GLYCANS_PER_TYPE_DEFAULT: int = 15

# VIP score parameters
VIP_TOP_N_DEFAULT: int = 10

# Volcano plot parameters
VOLCANO_PVALUE_THRESHOLD: float = 0.05
VOLCANO_FC_THRESHOLD: float = 2.0  # Linear fold change
VOLCANO_LOG2FC_THRESHOLD: float = 1.0  # log2(2) = 1

# ==============================================================================
# File Patterns
# ==============================================================================

CSV_FILE_PATTERN: str = '*.csv'
SAMPLE_ID_PATTERN: str = r'([CN])_(\d+)'

# ==============================================================================
# Output Filenames
# ==============================================================================

OUTPUT_INTEGRATED: str = 'integrated.csv'
OUTPUT_STATISTICS: str = 'glycan_type_statistics.csv'
OUTPUT_VIP_SCORES: str = 'vip_scores_all.csv'
OUTPUT_SUMMARY: str = 'analysis_summary.txt'

# Trace directory
TRACE_DIR: str = 'Trace'

# ==============================================================================
# Regex Patterns
# ==============================================================================

# Pattern to extract monosaccharide and count: e.g., H(5), N(4), A(2)
GLYCAN_COMPOSITION_PATTERN: str = r'([A-Z]+)\((\d+)\)'

# ==============================================================================
# Color Schemes
# ==============================================================================

# Default glycan type colors (for comparison heatmap)
DEFAULT_GLYCAN_TYPE_COLORS: dict = {
    GLYCAN_TYPE_HM: '#00CC00',   # Green
    GLYCAN_TYPE_F: '#FF0000',    # Red
    GLYCAN_TYPE_S: '#FF69B4',    # Pink
    GLYCAN_TYPE_SF: '#FFA500',   # Orange
    GLYCAN_TYPE_CH: '#0000FF'    # Blue
}

# Legacy colors (for backward compatibility)
DEFAULT_LEGACY_COLORS: dict = {
    LEGACY_NON: '#CCCCCC',
    LEGACY_SIALYLATED: '#E74C3C',
    LEGACY_FUCOSYLATED: '#3498DB',
    LEGACY_BOTH: '#9B59B6'
}

# Symbols for comparison heatmap
SYMBOL_CANCER: str = 'x'
SYMBOL_NORMAL: str = '+'

# ==============================================================================
# Statistical Parameters
# ==============================================================================

# Normalization
NORMALIZATION_METHOD_TIC: str = 'tic'  # Total Ion Current

# Scaling
SCALER_ROBUST: str = 'robust'
SCALER_STANDARD: str = 'standard'

# S0 parameter for t-test stability
# Prevents division by very small standard deviations which can lead to
# unreliable t-statistics. Inspired by SAM (Significance Analysis of Microarrays)
S0_DEFAULT: float = 0.05  # Default S0 parameter (5% of median standard deviation)
S0_MIN_SAMPLES: int = 3   # Minimum samples required per group for t-test

# ==============================================================================
# Logging
# ==============================================================================

LOG_FORMAT: str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
LOG_LEVEL_INFO: str = 'INFO'
LOG_LEVEL_DEBUG: str = 'DEBUG'
LOG_LEVEL_WARNING: str = 'WARNING'
LOG_LEVEL_ERROR: str = 'ERROR'

# ==============================================================================
# UI/Formatting Constants
# ==============================================================================

# Separator line for logging and reports
SEPARATOR_LINE: str = "=" * 80
SEPARATOR_LINE_SHORT: str = "-" * 80

# ==============================================================================
# Statistical Validation Constants
# ==============================================================================

# Bootstrap validation parameters
BOOTSTRAP_ITERATIONS: int = 1000
BOOTSTRAP_CONFIDENCE_LEVEL: float = 0.95

# Stability and effect size thresholds
STABILITY_THRESHOLD: float = 0.8  # 80% stability for biomarkers
COHENS_D_LARGE_EFFECT: float = 0.8  # Cohen's d threshold for large effects
VIP_THRESHOLD: float = 1.0  # VIP score threshold for important features

# Cross-validation parameters
CV_N_FOLDS: int = 10  # 10-fold cross-validation

# Permutation test parameters
PERMUTATION_N_ITERATIONS: int = 1000

# ==============================================================================
# File Name Constants
# ==============================================================================

# Output file names
FILENAME_INTEGRATED_RAW: str = 'integrated.csv'
FILENAME_INTEGRATED_FILTERED: str = 'integrated_filtered.csv'
FILENAME_FILTERING_REPORT: str = 'filtering_report.txt'
FILENAME_GLYCAN_STATS: str = 'glycan_type_statistics.csv'
FILENAME_VIP_SCORES: str = 'vip_scores_all.csv'
FILENAME_BOOTSTRAP_VALIDATION: str = 'vip_bootstrap_validation.csv'
FILENAME_STABLE_BIOMARKERS: str = 'stable_biomarkers.csv'
FILENAME_CV_RESULTS: str = 'plsda_cross_validation.txt'
FILENAME_COHENS_D: str = 'cohens_d_effect_sizes.csv'
FILENAME_PCA_PERMUTATION: str = 'pca_permutation_test.txt'
FILENAME_ANALYSIS_SUMMARY: str = 'analysis_summary.txt'
FILENAME_INPUT_MANIFEST: str = 'input_data_manifest.json'
FILENAME_OUTPUT_MANIFEST: str = 'output_data_manifest.json'
FILENAME_EXECUTION_METADATA: str = 'execution_metadata.json'


# Export all constants
__all__ = [
    # Column Names
    'REQUIRED_INPUT_COLUMNS',
    'METADATA_COLUMNS',
    # Sample Groups
    'CANCER_PREFIX',
    'NORMAL_PREFIX',
    'GROUP_CANCER',
    'GROUP_NORMAL',
    # Monosaccharide Symbols
    'MONOSACCHARIDE_H',
    'MONOSACCHARIDE_N',
    'MONOSACCHARIDE_A',
    'MONOSACCHARIDE_F',
    'MONOSACCHARIDE_G',
    'DEFAULT_SIALYLATION_MARKER',
    'DEFAULT_FUCOSYLATION_MARKER',
    'MONOSACCHARIDE_ORDER',
    # Glycan Type Categories
    'GLYCAN_TYPE_HM',
    'GLYCAN_TYPE_F',
    'GLYCAN_TYPE_S',
    'GLYCAN_TYPE_SF',
    'GLYCAN_TYPE_CH',
    'GLYCAN_TYPE_UNKNOWN',
    'GLYCAN_TYPE_ORDER',
    # Classification Labels
    'PRIMARY_TRUNCATED',
    'PRIMARY_HIGH_MANNOSE',
    'PRIMARY_COMPLEX_HYBRID',
    'PRIMARY_OUTLIER',
    'PRIMARY_UNKNOWN',
    'SECONDARY_TRUNCATED',
    'SECONDARY_HIGH_MANNOSE',
    'SECONDARY_COMPLEX_HYBRID',
    'SECONDARY_FUCOSYLATED',
    'SECONDARY_SIALYLATED',
    'SECONDARY_SIALOFUCOSYLATED',
    'SECONDARY_OUTLIER',
    'SECONDARY_UNKNOWN',
    'LEGACY_NON',
    'LEGACY_SIALYLATED',
    'LEGACY_FUCOSYLATED',
    'LEGACY_BOTH',
    # High-Mannose Criteria
    'HIGH_MANNOSE_MIN_H',
    'HIGH_MANNOSE_EXACT_N',
    # Complex/Hybrid Criteria
    'COMPLEX_HYBRID_MIN_N',
    # Analysis Parameters
    'DEFAULT_PCA_COMPONENTS',
    'DEFAULT_PLSDA_COMPONENTS',
    'DEFAULT_LOG_TRANSFORM',
    'DEFAULT_SIGNIFICANCE_ALPHA',
    'LOG_TRANSFORM_PSEUDOCOUNT',
    # Visualization Parameters
    'DEFAULT_DPI',
    'FIGSIZE_PCA',
    'FIGSIZE_BOXPLOT',
    'FIGSIZE_HEATMAP',
    'FIGSIZE_COMPARISON_HEATMAP',
    'HEATMAP_TOP_N_DEFAULT',
    'COMPARISON_MAX_PEPTIDES_DEFAULT',
    'COMPARISON_MAX_GLYCANS_PER_TYPE_DEFAULT',
    'VIP_TOP_N_DEFAULT',
    'VOLCANO_PVALUE_THRESHOLD',
    'VOLCANO_FC_THRESHOLD',
    'VOLCANO_LOG2FC_THRESHOLD',
    # File Patterns
    'CSV_FILE_PATTERN',
    'SAMPLE_ID_PATTERN',
    # Output Filenames
    'OUTPUT_INTEGRATED',
    'OUTPUT_STATISTICS',
    'OUTPUT_VIP_SCORES',
    'OUTPUT_SUMMARY',
    'TRACE_DIR',
    # Regex Patterns
    'GLYCAN_COMPOSITION_PATTERN',
    # Color Schemes
    'DEFAULT_GLYCAN_TYPE_COLORS',
    'DEFAULT_LEGACY_COLORS',
    'SYMBOL_CANCER',
    'SYMBOL_NORMAL',
    # Statistical Parameters
    'NORMALIZATION_METHOD_TIC',
    'SCALER_ROBUST',
    'SCALER_STANDARD',
    'S0_DEFAULT',
    'S0_MIN_SAMPLES',
    # Logging
    'LOG_FORMAT',
    'LOG_LEVEL_INFO',
    'LOG_LEVEL_DEBUG',
    'LOG_LEVEL_WARNING',
    'LOG_LEVEL_ERROR',
    # UI/Formatting
    'SEPARATOR_LINE',
    'SEPARATOR_LINE_SHORT',
    # Statistical Validation
    'BOOTSTRAP_ITERATIONS',
    'BOOTSTRAP_CONFIDENCE_LEVEL',
    'STABILITY_THRESHOLD',
    'COHENS_D_LARGE_EFFECT',
    'VIP_THRESHOLD',
    'CV_N_FOLDS',
    'PERMUTATION_N_ITERATIONS',
    # File Names
    'FILENAME_INTEGRATED_RAW',
    'FILENAME_INTEGRATED_FILTERED',
    'FILENAME_FILTERING_REPORT',
    'FILENAME_GLYCAN_STATS',
    'FILENAME_VIP_SCORES',
    'FILENAME_BOOTSTRAP_VALIDATION',
    'FILENAME_STABLE_BIOMARKERS',
    'FILENAME_CV_RESULTS',
    'FILENAME_COHENS_D',
    'FILENAME_PCA_PERMUTATION',
    'FILENAME_ANALYSIS_SUMMARY',
    'FILENAME_INPUT_MANIFEST',
    'FILENAME_OUTPUT_MANIFEST',
    'FILENAME_EXECUTION_METADATA',
]
