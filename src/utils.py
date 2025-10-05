"""
Utility Functions for pGlyco Auto Combine
Handles common data operations with proper type handling
"""

import pandas as pd
import numpy as np
import logging
import re
from pathlib import Path
from typing import List, Tuple, Union, Optional, Dict
from functools import lru_cache

from .constants import (
    METADATA_COLUMNS,
    CANCER_PREFIX,
    NORMAL_PREFIX,
    GROUP_CANCER,
    GROUP_NORMAL,
    TRACE_DIR,
    SAMPLE_ID_PATTERN,
    LOG_TRANSFORM_PSEUDOCOUNT
)
from .exceptions import (
    FileOperationError,
    TraceDataSaveError,
    OutputDirectoryError,
    ValidationError
)

logger = logging.getLogger(__name__)


# ==============================================================================
# Data Type Conversion
# ==============================================================================

def replace_empty_with_zero(data: Union[pd.DataFrame, pd.Series]) -> Union[pd.DataFrame, pd.Series]:
    """
    Replace empty strings with 0 without FutureWarning

    Uses mask() instead of replace() to avoid pandas downcasting warnings.

    Args:
        data: pandas Series or DataFrame

    Returns:
        Data with empty strings replaced by 0
    """
    if isinstance(data, pd.DataFrame):
        # Use mask() to avoid replace() downcasting warnings
        result = data.copy()
        result = result.mask(result == '', np.nan)
        result = result.apply(pd.to_numeric, errors='coerce')
        result = result.fillna(0)
        return result
    elif isinstance(data, pd.Series):
        # Use mask() to avoid replace() downcasting warnings
        result = data.copy()
        result = result.mask(result == '', np.nan)
        result = pd.to_numeric(result, errors='coerce')
        result = result.fillna(0)
        return result
    else:
        return data


def to_numeric_safe(data: Union[pd.DataFrame, pd.Series]) -> Union[pd.DataFrame, pd.Series]:
    """
    Convert to numeric safely without downcasting warnings

    Args:
        data: pandas Series or DataFrame

    Returns:
        Numeric data
    """
    if isinstance(data, pd.DataFrame):
        return data.apply(pd.to_numeric, errors='coerce').fillna(0)
    elif isinstance(data, pd.Series):
        return pd.to_numeric(data, errors='coerce').fillna(0)
    else:
        return data


def convert_to_numeric_with_na(data: Union[pd.DataFrame, pd.Series]) -> Union[pd.DataFrame, pd.Series]:
    """
    Convert empty strings to NaN (NOT zero) for proper missing data handling

    SCIENTIFIC VALIDITY: Missing values should be treated as NaN, not zero.
    Zero implies "measured and found to be absent", while NaN means "not measured".

    Args:
        data: pandas Series or DataFrame with potential empty strings

    Returns:
        Data with empty strings converted to NaN (numeric dtype)
    """
    if isinstance(data, pd.DataFrame):
        result = data.copy()
        result = result.mask(result == '', np.nan)
        result = result.apply(pd.to_numeric, errors='coerce')
        return result  # Keep NaN, don't fill with zero
    elif isinstance(data, pd.Series):
        result = data.copy()
        result = result.mask(result == '', np.nan)
        result = pd.to_numeric(result, errors='coerce')
        return result  # Keep NaN, don't fill with zero
    else:
        return data


def calculate_group_statistics(df: pd.DataFrame, sample_cols: List[str]) -> Dict[str, pd.Series]:
    """
    Calculate statistics using only non-missing values (skipna=True)

    SCIENTIFIC VALIDITY: Using skipna=True ensures that missing values don't bias
    the mean downward. Only measured values are used in calculations.

    Args:
        df: DataFrame with intensity data
        sample_cols: List of sample column names

    Returns:
        Dictionary with statistics:
        - mean: Mean of non-missing values
        - std: Standard deviation of non-missing values
        - count: Number of non-missing values
        - min: Minimum of non-missing values
        - max: Maximum of non-missing values
        - missing_rate: Proportion of missing values
    """
    # Convert to numeric with NaN for missing
    values = convert_to_numeric_with_na(df[sample_cols])

    return {
        'mean': values.mean(axis=1, skipna=True),
        'std': values.std(axis=1, skipna=True),
        'count': values.count(axis=1),  # Count of non-NaN values
        'min': values.min(axis=1, skipna=True),
        'max': values.max(axis=1, skipna=True),
        'missing_rate': values.isna().sum(axis=1) / len(sample_cols)
    }


# ==============================================================================
# File Operations
# ==============================================================================

def ensure_directory(directory: Union[str, Path]) -> Path:
    """
    Ensure a directory exists, create if it doesn't

    Args:
        directory: Path to directory

    Returns:
        Path object to the directory

    Raises:
        OutputDirectoryError: If directory cannot be created
    """
    try:
        dir_path = Path(directory)
        dir_path.mkdir(parents=True, exist_ok=True)
        return dir_path
    except Exception as e:
        raise OutputDirectoryError(str(directory), str(e))


def ensure_trace_dir(output_dir: Union[str, Path]) -> Path:
    """
    Ensure Trace directory exists for visualization data exports

    Args:
        output_dir: Base output directory

    Returns:
        Path to Trace directory
    """
    trace_dir = Path(output_dir) / TRACE_DIR
    return ensure_directory(trace_dir)


def save_trace_data(data: pd.DataFrame, output_dir: Union[str, Path], filename: str) -> Path:
    """
    Save visualization source data to Trace folder

    Args:
        data: DataFrame to save
        output_dir: Base output directory
        filename: Name of the trace file (e.g., 'heatmap_data.csv')

    Returns:
        Path to saved file

    Raises:
        TraceDataSaveError: If trace data cannot be saved
    """
    try:
        trace_dir = ensure_trace_dir(output_dir)
        output_path = trace_dir / filename
        data.to_csv(output_path, index=False)
        logger.debug(f"Saved trace data to {output_path}")
        return output_path
    except Exception as e:
        raise TraceDataSaveError(filename, str(e))


# ==============================================================================
# Sample Column Extraction
# ==============================================================================

def get_sample_columns(df: pd.DataFrame,
                       exclude_metadata: bool = True) -> Tuple[List[str], List[str]]:
    """
    Extract cancer and normal sample columns from DataFrame

    Identifies columns that represent cancer samples (start with 'C' followed by digits)
    and normal samples (start with 'N' followed by digits).

    Args:
        df: DataFrame with sample columns
        exclude_metadata: If True, only return sample columns (no metadata)

    Returns:
        Tuple of (cancer_samples, normal_samples) where each is a list of column names
    """
    if exclude_metadata:
        # Get all columns except metadata
        sample_cols = [col for col in df.columns if col not in METADATA_COLUMNS]
    else:
        sample_cols = df.columns.tolist()

    cancer_samples = [col for col in sample_cols
                     if col.startswith(CANCER_PREFIX) and col[1:].isdigit()]
    normal_samples = [col for col in sample_cols
                     if col.startswith(NORMAL_PREFIX) and col[1:].isdigit()]

    return cancer_samples, normal_samples


def get_all_sample_columns(df: pd.DataFrame) -> List[str]:
    """
    Get all sample columns (both cancer and normal)

    Args:
        df: DataFrame with sample columns

    Returns:
        List of all sample column names sorted (C samples first, then N samples)
    """
    cancer_samples, normal_samples = get_sample_columns(df)
    return cancer_samples + normal_samples


def get_metadata_columns(df: pd.DataFrame) -> List[str]:
    """
    Get metadata columns from DataFrame

    Args:
        df: DataFrame

    Returns:
        List of metadata column names that exist in the DataFrame
    """
    return [col for col in df.columns if col in METADATA_COLUMNS]


# ==============================================================================
# Sample ID Utilities
# ==============================================================================

def extract_sample_id(filename: str) -> str:
    """
    Extract sample ID from filename (e.g., C_01.csv -> C1)

    Args:
        filename: Name of the CSV file

    Returns:
        Sample ID (e.g., 'C1', 'N10')
    """
    match = re.search(SAMPLE_ID_PATTERN, filename)
    if match:
        prefix = match.group(1)
        number = int(match.group(2))
        return f"{prefix}{number}"
    return filename


def get_sample_group(sample_name: str) -> str:
    """
    Determine if sample is Cancer or Normal based on name

    Args:
        sample_name: Sample column name (e.g., 'C1', 'N10')

    Returns:
        'Cancer' or 'Normal'
    """
    if sample_name.startswith(CANCER_PREFIX):
        return GROUP_CANCER
    elif sample_name.startswith(NORMAL_PREFIX):
        return GROUP_NORMAL
    else:
        return 'Unknown'


def is_cancer_sample(sample_name: str) -> bool:
    """Check if sample is a cancer sample"""
    return sample_name.startswith(CANCER_PREFIX)


def is_normal_sample(sample_name: str) -> bool:
    """Check if sample is a normal sample"""
    return sample_name.startswith(NORMAL_PREFIX)


# ==============================================================================
# Data Transformation
# ==============================================================================

def log_transform(data: Union[pd.DataFrame, pd.Series, np.ndarray],
                 pseudocount: float = LOG_TRANSFORM_PSEUDOCOUNT) -> Union[pd.DataFrame, pd.Series, np.ndarray]:
    """
    Apply log2 transformation with pseudocount

    Args:
        data: Data to transform
        pseudocount: Value to add before log transform (default: 1.0)

    Returns:
        Log2-transformed data
    """
    if isinstance(data, (pd.DataFrame, pd.Series)):
        return np.log2(data + pseudocount)
    else:
        return np.log2(data + pseudocount)


# ==============================================================================
# Statistical Utilities
# ==============================================================================

def validate_statistical_power(cancer_samples: List[str], normal_samples: List[str],
                               min_n: int = 5) -> None:
    """
    Validate that sample sizes are sufficient for statistical analysis

    SCIENTIFIC VALIDITY: Statistical tests require minimum sample sizes for adequate power.
    - n < 3: Invalid (cannot compute variance)
    - n < 5: Under-powered (high Type II error rate)
    - n >= 5: Acceptable for exploratory analysis
    - n >= 10: Good power for hypothesis testing

    Args:
        cancer_samples: List of cancer sample IDs
        normal_samples: List of normal sample IDs
        min_n: Minimum recommended sample size (default: 5)

    Raises:
        InsufficientDataError: If either group has < 3 samples (invalid)

    Warnings:
        Logs warning if either group has < min_n samples
    """
    import logging
    logger = logging.getLogger(__name__)

    cancer_n = len(cancer_samples)
    normal_n = len(normal_samples)

    # Critical: Cannot perform tests with n < 3
    if cancer_n < 3:
        from src.exceptions import InsufficientDataError
        raise InsufficientDataError(
            f"Cancer group has only {cancer_n} samples. Minimum 3 required for statistical tests."
        )

    if normal_n < 3:
        from src.exceptions import InsufficientDataError
        raise InsufficientDataError(
            f"Normal group has only {normal_n} samples. Minimum 3 required for statistical tests."
        )

    # Warning: Under-powered if n < min_n
    if cancer_n < min_n:
        logger.warning(
            f"⚠️  Cancer group has only {cancer_n} samples (< {min_n} recommended). "
            f"Statistical tests may be under-powered. Consider collecting more samples."
        )

    if normal_n < min_n:
        logger.warning(
            f"⚠️  Normal group has only {normal_n} samples (< {min_n} recommended). "
            f"Statistical tests may be under-powered. Consider collecting more samples."
        )

    # Info: Report sample sizes
    logger.info(f"Sample size validation: Cancer n={cancer_n}, Normal n={normal_n}")


def calculate_fold_change(cancer_mean: float, normal_mean: float,
                         log_scale: bool = False, pseudocount: float = 1.0) -> float:
    """
    Calculate fold change between cancer and normal with pseudocount

    SCIENTIFIC VALIDITY: Uses pseudocount to handle zeros and make log2 fold change symmetric.
    - Simple division: 2-fold up = 2.0, 2-fold down = 0.5 (asymmetric)
    - Log2 fold change: 2-fold up = +1.0, 2-fold down = -1.0 (symmetric)

    Args:
        cancer_mean: Mean intensity in cancer samples
        normal_mean: Mean intensity in normal samples
        log_scale: If True, return log2 fold change
        pseudocount: Small constant added to avoid division by zero (default: 1.0)

    Returns:
        Fold change value (linear or log2 scale)

    Examples:
        >>> calculate_fold_change(200, 100, log_scale=True)  # 2-fold increase
        1.0
        >>> calculate_fold_change(100, 200, log_scale=True)  # 2-fold decrease
        -1.0
    """
    if log_scale:
        # Log2 fold change with pseudocount (scientifically preferred)
        # Handles zeros gracefully and provides symmetric scale
        log2_fc = np.log2((cancer_mean + pseudocount) / (normal_mean + pseudocount))
        return log2_fc
    else:
        # Linear fold change (for backward compatibility)
        if normal_mean == 0:
            # Return inf if cancer > 0, else 0 (edge case)
            return np.inf if cancer_mean > 0 else 0

        fc = cancer_mean / normal_mean
        return fc


def calculate_statistics(values: Union[pd.Series, np.ndarray]) -> dict:
    """
    Calculate basic statistics for a set of values

    Args:
        values: Array of values

    Returns:
        Dictionary with mean, median, std, min, max
    """
    if isinstance(values, pd.Series):
        values = values.values

    # Remove NaN and zero values for statistics
    values_clean = values[~np.isnan(values) & (values > 0)]

    if len(values_clean) == 0:
        return {
            'mean': 0,
            'median': 0,
            'std': 0,
            'min': 0,
            'max': 0,
            'count': 0
        }

    return {
        'mean': float(np.mean(values_clean)),
        'median': float(np.median(values_clean)),
        'std': float(np.std(values_clean)),
        'min': float(np.min(values_clean)),
        'max': float(np.max(values_clean)),
        'count': len(values_clean)
    }


# ==============================================================================
# Validation Utilities
# ==============================================================================

def validate_sample_counts(df: pd.DataFrame,
                          min_cancer: int = 1,
                          min_normal: int = 1) -> Tuple[int, int]:
    """
    Validate sample counts in DataFrame

    Args:
        df: DataFrame with sample columns
        min_cancer: Minimum number of cancer samples required
        min_normal: Minimum number of normal samples required

    Returns:
        Tuple of (n_cancer, n_normal)

    Raises:
        ValidationError: If sample counts are insufficient
    """
    cancer_samples, normal_samples = get_sample_columns(df)
    n_cancer = len(cancer_samples)
    n_normal = len(normal_samples)

    if n_cancer < min_cancer:
        raise ValidationError(
            f"Insufficient cancer samples: {n_cancer} < {min_cancer}"
        )

    if n_normal < min_normal:
        raise ValidationError(
            f"Insufficient normal samples: {n_normal} < {min_normal}"
        )

    return n_cancer, n_normal


def validate_dataframe_not_empty(df: pd.DataFrame, context: str = "") -> None:
    """
    Validate that DataFrame is not empty

    Args:
        df: DataFrame to validate
        context: Context message for error

    Raises:
        ValidationError: If DataFrame is empty
    """
    if df is None or len(df) == 0:
        message = "DataFrame is empty"
        if context:
            message += f": {context}"
        raise ValidationError(message)


# ==============================================================================
# Formatting Utilities
# ==============================================================================

def format_percentage(value: float, decimals: int = 1) -> str:
    """Format value as percentage string"""
    return f"{value * 100:.{decimals}f}%"


def format_scientific(value: float, decimals: int = 2) -> str:
    """Format value in scientific notation"""
    return f"{value:.{decimals}e}"


def format_pvalue(pvalue: float) -> str:
    """
    Format p-value with appropriate precision

    Args:
        pvalue: P-value to format

    Returns:
        Formatted string (e.g., "0.001" or "< 0.001")
    """
    if pvalue < 0.001:
        return "< 0.001"
    elif pvalue < 0.01:
        return f"{pvalue:.4f}"
    else:
        return f"{pvalue:.3f}"


# ==============================================================================
# Caching Utilities
# ==============================================================================

@lru_cache(maxsize=128)
def get_cached_metadata_columns() -> Tuple[str, ...]:
    """
    Get metadata columns (cached for performance)

    Returns:
        Tuple of metadata column names
    """
    return tuple(METADATA_COLUMNS)
