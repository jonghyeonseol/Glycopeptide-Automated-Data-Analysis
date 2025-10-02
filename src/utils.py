"""
Utility Functions for pGlyco Auto Combine
Handles common data operations with proper type handling
"""

import pandas as pd
import numpy as np
from pathlib import Path


def replace_empty_with_zero(data):
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


def to_numeric_safe(data):
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


def ensure_trace_dir(output_dir: Path) -> Path:
    """
    Ensure Trace directory exists for visualization data exports

    Args:
        output_dir: Base output directory

    Returns:
        Path to Trace directory
    """
    trace_dir = Path(output_dir) / 'Trace'
    trace_dir.mkdir(parents=True, exist_ok=True)
    return trace_dir


def save_trace_data(data: pd.DataFrame, output_dir: Path, filename: str):
    """
    Save visualization source data to Trace folder

    Args:
        data: DataFrame to save
        output_dir: Base output directory
        filename: Name of the trace file (e.g., 'heatmap_data.csv')
    """
    trace_dir = ensure_trace_dir(output_dir)
    output_path = trace_dir / filename
    data.to_csv(output_path, index=False)
    return output_path


def get_sample_columns(df: pd.DataFrame) -> tuple:
    """
    Extract cancer and normal sample columns from DataFrame

    Identifies columns that represent cancer samples (start with 'C' followed by digits)
    and normal samples (start with 'N' followed by digits).

    Args:
        df: DataFrame with sample columns

    Returns:
        Tuple of (cancer_samples, normal_samples) where each is a list of column names
    """
    cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
    normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]
    return cancer_samples, normal_samples
