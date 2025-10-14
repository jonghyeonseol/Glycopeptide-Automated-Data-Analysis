"""
Utility Functions for pGlyco Auto Combine
Handles common data operations with proper type handling
"""

import pandas as pd
import numpy as np
import logging
import re
from pathlib import Path
from typing import List, Tuple, Union, Dict, Optional
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

def clean_inf_nan(data: Union[pd.DataFrame, pd.Series, np.ndarray],
                  context: str = "transformation",
                  verbose: bool = True) -> Union[pd.DataFrame, pd.Series, np.ndarray]:
    """
    Clean infinity and NaN values from data after transformations

    This is CRITICAL for numerical stability after operations that can produce
    inf/nan (division, log, sqrt, etc.). Inspired by alphapeptstats' approach.

    SCIENTIFIC RATIONALE:
    - inf values corrupt all downstream calculations
    - They arise from: division by zero, log(0), overflow, etc.
    - Must be cleaned IMMEDIATELY after each transformation
    - Replacing with NaN allows proper skipna handling in statistics

    Args:
        data: Data to clean
        context: Description of where this is called (for logging)
        verbose: If True, log warnings when inf/nan are found

    Returns:
        Data with inf replaced by NaN
    """
    if isinstance(data, pd.DataFrame):
        inf_count_before = np.isinf(data.values).sum()
        nan_count_before = np.isnan(data.values).sum()

        # Replace inf with NaN
        data_clean = data.replace([np.inf, -np.inf], np.nan)

        if verbose and (inf_count_before > 0):
            logger.warning(
                f"⚠️  Cleaned {inf_count_before} inf values after {context} "
                f"(now NaN, will be handled by skipna)"
            )

        return data_clean

    elif isinstance(data, pd.Series):
        inf_count_before = np.isinf(data.values).sum()

        # Replace inf with NaN
        data_clean = data.replace([np.inf, -np.inf], np.nan)

        if verbose and (inf_count_before > 0):
            logger.warning(
                f"⚠️  Cleaned {inf_count_before} inf values after {context}"
            )

        return data_clean

    elif isinstance(data, np.ndarray):
        inf_count_before = np.isinf(data).sum()

        # Replace inf with NaN
        data_clean = data.copy()
        data_clean[np.isinf(data_clean)] = np.nan

        if verbose and (inf_count_before > 0):
            logger.warning(
                f"⚠️  Cleaned {inf_count_before} inf values after {context}"
            )

        return data_clean

    else:
        return data


def safe_log_transform(data: Union[pd.DataFrame, pd.Series, np.ndarray],
                       pseudocount: float = LOG_TRANSFORM_PSEUDOCOUNT,
                       base: int = 2) -> Union[pd.DataFrame, pd.Series, np.ndarray]:
    """
    Apply log transformation with automatic inf/nan cleanup

    CRITICAL: This function ensures numerical stability by:
    1. Adding pseudocount to avoid log(0) = -inf
    2. Cleaning any resulting inf/nan values
    3. Logging warnings if issues are found

    Args:
        data: Data to transform
        pseudocount: Value to add before log transform (default: 1.0)
        base: Logarithm base (2 for log2, 10 for log10, etc.)

    Returns:
        Log-transformed data with inf/nan cleaned
    """
    if base == 2:
        transformed = np.log2(data + pseudocount)
    elif base == 10:
        transformed = np.log10(data + pseudocount)
    elif base == np.e:
        transformed = np.log(data + pseudocount)
    else:
        transformed = np.log(data + pseudocount) / np.log(base)

    # CRITICAL: Clean inf/nan immediately after transformation
    return clean_inf_nan(transformed, context="log transformation", verbose=True)


def log_transform(data: Union[pd.DataFrame, pd.Series, np.ndarray],
                  pseudocount: float = LOG_TRANSFORM_PSEUDOCOUNT) -> Union[pd.DataFrame, pd.Series, np.ndarray]:
    """
    Apply log2 transformation with pseudocount

    DEPRECATED: Use safe_log_transform() for better numerical stability.
    This function is kept for backward compatibility.

    Args:
        data: Data to transform
        pseudocount: Value to add before log transform (default: 1.0)

    Returns:
        Log2-transformed data
    """
    return safe_log_transform(data, pseudocount=pseudocount, base=2)


def safe_division(numerator: Union[float, np.ndarray, pd.Series],
                  denominator: Union[float, np.ndarray, pd.Series],
                  fill_value: float = np.nan,
                  min_denominator: float = 1e-10) -> Union[float, np.ndarray, pd.Series]:
    """
    Perform division with protection against division by zero

    CRITICAL: Prevents inf values from corrupting downstream calculations.

    Args:
        numerator: Numerator values
        denominator: Denominator values
        fill_value: Value to use when denominator is too small (default: NaN)
        min_denominator: Minimum denominator value (default: 1e-10)

    Returns:
        Result of division with inf/nan handled

    Examples:
        >>> safe_division(10, 0)  # Returns NaN instead of inf
        nan
        >>> safe_division(10, 1e-12)  # Returns NaN (denominator too small)
        nan
        >>> safe_division(10, 5)  # Normal division
        2.0
    """
    if isinstance(denominator, (pd.Series, np.ndarray)):
        # Vectorized division
        result = np.where(
            np.abs(denominator) > min_denominator,
            numerator / denominator,
            fill_value
        )
        if isinstance(numerator, pd.Series):
            return pd.Series(result, index=numerator.index)
        return result
    else:
        # Scalar division
        if abs(denominator) > min_denominator:
            return numerator / denominator
        else:
            return fill_value


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
            "Statistical tests may be under-powered. Consider collecting more samples."
        )

    if normal_n < min_n:
        logger.warning(
            f"⚠️  Normal group has only {normal_n} samples (< {min_n} recommended). "
            "Statistical tests may be under-powered. Consider collecting more samples."
        )

    # Info: Report sample sizes
    logger.info(f"Sample size validation: Cancer n={cancer_n}, Normal n={normal_n}")


def calculate_fold_change(cancer_mean: Union[float, np.ndarray, pd.Series],
                          normal_mean: Union[float, np.ndarray, pd.Series],
                          log_scale: bool = False,
                          pseudocount: float = 1.0,
                          adaptive_pseudocount: bool = False) -> Union[float, np.ndarray, pd.Series]:
    """
    Calculate fold change between cancer and normal with robust zero handling

    ENHANCED VERSION with multiple improvements:
    1. Uses safe_division to avoid inf values
    2. Supports vectorized operations for pandas Series/arrays
    3. Optional adaptive pseudocount based on data scale
    4. Returns NaN for problematic cases instead of inf

    SCIENTIFIC VALIDITY: Uses pseudocount to handle zeros and make log2 fold change symmetric.
    - Simple division: 2-fold up = 2.0, 2-fold down = 0.5 (asymmetric)
    - Log2 fold change: 2-fold up = +1.0, 2-fold down = -1.0 (symmetric)

    Args:
        cancer_mean: Mean intensity in cancer samples (scalar, array, or Series)
        normal_mean: Mean intensity in normal samples (scalar, array, or Series)
        log_scale: If True, return log2 fold change
        pseudocount: Small constant added to avoid division by zero (default: 1.0)
        adaptive_pseudocount: If True, use 1% of median intensity as pseudocount (default: False)

    Returns:
        Fold change value (linear or log2 scale), same type as input

    Examples:
        >>> calculate_fold_change(200, 100, log_scale=True)  # 2-fold increase
        1.0
        >>> calculate_fold_change(100, 200, log_scale=True)  # 2-fold decrease
        -1.0
        >>> calculate_fold_change(0, 0, log_scale=False)  # Both zero
        nan
    """
    # Adaptive pseudocount: use 1% of median intensity (inspired by alphapeptstats)
    if adaptive_pseudocount:
        if isinstance(cancer_mean, (pd.Series, np.ndarray)) and isinstance(normal_mean, (pd.Series, np.ndarray)):
            # Vectorized: calculate median of all non-zero values
            all_values = np.concatenate([
                cancer_mean[cancer_mean > 0] if hasattr(cancer_mean, '__len__') else [cancer_mean],
                normal_mean[normal_mean > 0] if hasattr(normal_mean, '__len__') else [normal_mean]
            ])
            if len(all_values) > 0:
                pseudocount = np.median(all_values) * 0.01
            else:
                pseudocount = 1.0  # Fallback
        else:
            # Scalar: use 1% of average of the two values
            avg_intensity = (abs(cancer_mean) + abs(normal_mean)) / 2
            pseudocount = max(avg_intensity * 0.01, 1.0)

    if log_scale:
        # Log2 fold change with pseudocount (scientifically preferred)
        # Handles zeros gracefully and provides symmetric scale

        # Add pseudocount to both numerator and denominator
        cancer_adj = cancer_mean + pseudocount
        normal_adj = normal_mean + pseudocount

        # Use safe_division to avoid inf
        ratio = safe_division(cancer_adj, normal_adj, fill_value=np.nan, min_denominator=1e-10)

        # Apply log2 safely
        log2_fc = safe_log_transform(ratio, pseudocount=0, base=2)

        return log2_fc
    else:
        # Linear fold change with robust zero handling
        # Add pseudocount to both to handle zeros symmetrically
        cancer_adj = cancer_mean + pseudocount
        normal_adj = normal_mean + pseudocount

        # Use safe_division to avoid inf
        fc = safe_division(cancer_adj, normal_adj, fill_value=np.nan, min_denominator=1e-10)

        return fc


def detect_log_transform(data: Union[pd.DataFrame, pd.Series, np.ndarray],
                         sample_cols: Optional[List[str]] = None,
                         threshold_confidence: float = 0.7) -> Dict[str, any]:
    """
    Automatically detect if data is already log-transformed

    Uses multiple heuristics to determine if intensity data is on log scale or raw scale.
    Inspired by alphapeptstats' preprocessing tracking but adds automatic detection.

    HEURISTICS:
    1. Negative values: Log-transformed data can have negatives, raw cannot
    2. Value range: Log2 typically -10 to +30, raw is 0 to 1e9
    3. Median: Log2 median typically 5-20, raw median typically 1e4-1e7
    4. Max/Min ratio: Log2 has smaller ratio (~40), raw has huge ratio (>1e6)

    SCIENTIFIC RATIONALE:
    - Prevents double log-transformation (causes severe compression)
    - Ensures statistical tests use correct scale
    - Validates preprocessing assumptions

    Args:
        data: DataFrame, Series, or array with intensity values
        sample_cols: List of sample columns (for DataFrame). If None, uses all numeric columns
        threshold_confidence: Minimum confidence (0-1) to declare transformation state (default: 0.7)

    Returns:
        Dictionary with:
        - 'is_log_transformed': Boolean (True if likely log-transformed)
        - 'confidence': Float 0-1 (confidence in detection)
        - 'evidence': Dict with individual heuristic results
        - 'recommendation': String with human-readable recommendation

    Examples:
        >>> # Raw intensity data
        >>> detect_log_transform(raw_data)
        {'is_log_transformed': False, 'confidence': 0.95, ...}

        >>> # Log2-transformed data
        >>> detect_log_transform(log_data)
        {'is_log_transformed': True, 'confidence': 0.90, ...}
    """
    # Extract numeric values
    if isinstance(data, pd.DataFrame):
        if sample_cols is not None:
            values = data[sample_cols].values.flatten()
        else:
            # Use all numeric columns
            numeric_cols = data.select_dtypes(include=[np.number]).columns
            values = data[numeric_cols].values.flatten()
    elif isinstance(data, pd.Series):
        values = data.values
    else:
        values = np.asarray(data).flatten()

    # Remove NaN and zero values for analysis
    values_clean = values[~np.isnan(values) & (values != 0)]

    if len(values_clean) == 0:
        return {
            'is_log_transformed': None,
            'confidence': 0.0,
            'evidence': {},
            'recommendation': 'Cannot detect: No valid values found'
        }

    # Calculate statistics
    has_negatives = np.any(values_clean < 0)
    min_val = np.min(values_clean)
    max_val = np.max(values_clean)
    median_val = np.median(values_clean)
    mean_val = np.mean(values_clean)
    value_range = max_val - min_val
    dynamic_range = max_val / abs(min_val) if min_val != 0 else np.inf

    # Evidence dictionary
    evidence = {}
    log_score = 0.0  # Score 0-1, higher = more likely log-transformed
    weight_total = 0.0

    # Heuristic 1: Negative values (STRONG EVIDENCE WHEN PRESENT)
    # Note: Absence of negatives is NEUTRAL, not evidence against log-transform
    if has_negatives:
        weight = 0.3
        evidence['negative_values'] = {
            'result': 'HAS_NEGATIVES',
            'score': 1.0,
            'weight': weight,
            'interpretation': 'Negative values PROVE log-transformed (raw intensities cannot be negative)'
        }
        log_score += weight * 1.0
        weight_total += weight
    else:
        # Don't count absence of negatives - it's uninformative
        # Many log2-transformed datasets have no negatives (if all raw values > 1)
        weight = 0.0
        evidence['negative_values'] = {
            'result': 'NO_NEGATIVES',
            'score': None,
            'weight': 0.0,
            'interpretation': 'Absence of negatives is neutral (consistent with both raw and log-transformed)'
        }
        # Don't add to scores since it's uninformative

    # Heuristic 2: Maximum value (STRONG EVIDENCE)
    weight = 0.25
    if max_val > 1000:
        # Very high max suggests raw data
        score = 0.0
        evidence['max_value'] = {
            'value': float(max_val),
            'result': 'VERY_HIGH',
            'score': score,
            'weight': weight,
            'interpretation': f'Max value {max_val:.2e} suggests raw intensity'
        }
    elif max_val > 100:
        # Moderate max, could be either
        score = 0.3
        evidence['max_value'] = {
            'value': float(max_val),
            'result': 'MODERATE',
            'score': score,
            'weight': weight,
            'interpretation': f'Max value {max_val:.2f} is ambiguous'
        }
    elif max_val > 30:
        # Typical log2 max is ~20-30
        score = 0.7
        evidence['max_value'] = {
            'value': float(max_val),
            'result': 'TYPICAL_LOG',
            'score': score,
            'weight': weight,
            'interpretation': f'Max value {max_val:.2f} suggests log-transformed'
        }
    else:
        # Low max strongly suggests log2
        score = 1.0
        evidence['max_value'] = {
            'value': float(max_val),
            'result': 'LOW',
            'score': score,
            'weight': weight,
            'interpretation': f'Max value {max_val:.2f} strongly suggests log-transformed'
        }
    log_score += weight * score
    weight_total += weight

    # Heuristic 3: Median value (MODERATE EVIDENCE)
    weight = 0.2
    if median_val > 1000:
        score = 0.0
        evidence['median_value'] = {
            'value': float(median_val),
            'result': 'VERY_HIGH',
            'score': score,
            'weight': weight,
            'interpretation': f'Median {median_val:.2e} suggests raw intensity'
        }
    elif median_val > 100:
        score = 0.2
        evidence['median_value'] = {
            'value': float(median_val),
            'result': 'HIGH',
            'score': score,
            'weight': weight,
            'interpretation': f'Median {median_val:.2f} likely raw intensity'
        }
    elif median_val > 30:
        score = 0.5
        evidence['median_value'] = {
            'value': float(median_val),
            'result': 'MODERATE',
            'score': score,
            'weight': weight,
            'interpretation': f'Median {median_val:.2f} is ambiguous'
        }
    else:
        score = 0.9
        evidence['median_value'] = {
            'value': float(median_val),
            'result': 'LOW',
            'score': score,
            'weight': weight,
            'interpretation': f'Median {median_val:.2f} suggests log-transformed'
        }
    log_score += weight * score
    weight_total += weight

    # Heuristic 4: Value range (MODERATE EVIDENCE)
    weight = 0.15
    if value_range > 10000:
        score = 0.0
        evidence['value_range'] = {
            'range': float(value_range),
            'result': 'VERY_WIDE',
            'score': score,
            'weight': weight,
            'interpretation': f'Range {value_range:.2e} suggests raw intensity'
        }
    elif value_range > 1000:
        score = 0.2
        evidence['value_range'] = {
            'range': float(value_range),
            'result': 'WIDE',
            'score': score,
            'weight': weight,
            'interpretation': f'Range {value_range:.2f} likely raw intensity'
        }
    elif value_range > 100:
        score = 0.5
        evidence['value_range'] = {
            'range': float(value_range),
            'result': 'MODERATE',
            'score': score,
            'weight': weight,
            'interpretation': f'Range {value_range:.2f} is ambiguous'
        }
    else:
        score = 0.9
        evidence['value_range'] = {
            'range': float(value_range),
            'result': 'NARROW',
            'score': score,
            'weight': weight,
            'interpretation': f'Range {value_range:.2f} suggests log-transformed'
        }
    log_score += weight * score
    weight_total += weight

    # Heuristic 5: Minimum value (WEAK EVIDENCE)
    weight = 0.1
    if min_val < 0:
        score = 1.0
        evidence['min_value'] = {
            'value': float(min_val),
            'result': 'NEGATIVE',
            'score': score,
            'weight': weight,
            'interpretation': f'Min value {min_val:.2f} proves log-transformed'
        }
    elif min_val < 1:
        score = 0.8
        evidence['min_value'] = {
            'value': float(min_val),
            'result': 'SMALL_POSITIVE',
            'score': score,
            'weight': weight,
            'interpretation': f'Min value {min_val:.4f} suggests log-transformed'
        }
    elif min_val < 100:
        # Log2(100) ≈ 6.6, so min < 100 is consistent with log2
        score = 0.6
        evidence['min_value'] = {
            'value': float(min_val),
            'result': 'MODERATE',
            'score': score,
            'weight': weight,
            'interpretation': f'Min value {min_val:.2f} consistent with log-transformed'
        }
    else:
        score = 0.1
        evidence['min_value'] = {
            'value': float(min_val),
            'result': 'LARGE_POSITIVE',
            'score': score,
            'weight': weight,
            'interpretation': f'Min value {min_val:.2f} weakly suggests raw'
        }
    log_score += weight * score
    weight_total += weight

    # Calculate final confidence
    confidence = log_score / weight_total if weight_total > 0 else 0.5

    # Make decision
    if confidence >= threshold_confidence:
        is_log = True
        recommendation = f"Data appears LOG-TRANSFORMED (confidence: {confidence:.2%}). Do NOT apply log transformation."
    elif confidence <= (1 - threshold_confidence):
        is_log = False
        recommendation = f"Data appears RAW (confidence: {1-confidence:.2%}). Consider log2 transformation for normality."
    else:
        is_log = None
        recommendation = f"UNCERTAIN (confidence: {max(confidence, 1-confidence):.2%}). Manual inspection recommended."

    return {
        'is_log_transformed': is_log,
        'confidence': float(confidence),
        'evidence': evidence,
        'recommendation': recommendation,
        'statistics': {
            'n_values': len(values_clean),
            'min': float(min_val),
            'max': float(max_val),
            'median': float(median_val),
            'mean': float(mean_val),
            'range': float(value_range),
            'has_negatives': bool(has_negatives)
        }
    }


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
