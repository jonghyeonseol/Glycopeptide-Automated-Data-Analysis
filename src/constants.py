"""
Constants Module for pGlyco Auto Combine
Centralized configuration values and constants
"""

from typing import List, Set

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

# ==============================================================================
# Logging
# ==============================================================================

LOG_FORMAT: str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
LOG_LEVEL_INFO: str = 'INFO'
LOG_LEVEL_DEBUG: str = 'DEBUG'
LOG_LEVEL_WARNING: str = 'WARNING'
LOG_LEVEL_ERROR: str = 'ERROR'
