"""
Centralized Data Preparation Module for pGlyco Auto Combine
Ensures consistent data filtering and processing across all visualizations

CRITICAL: All visualization modules MUST use these standardized functions
to ensure data consistency and scientific reproducibility.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
from .logger_config import get_logger
from .utils import get_sample_columns

logger = get_logger(__name__)


class DataPreparationConfig:
    """Configuration for standardized data preparation"""

    def __init__(self, min_detection_pct: float = 0.30,
                 min_samples: int = 5,
                 missing_data_method: str = 'skipna'):
        """
        Initialize data preparation configuration

        Args:
            min_detection_pct: Minimum detection % in at least one group (default: 0.30 = 30%)
            min_samples: Minimum number of detected samples for statistical tests (default: 5)
            missing_data_method: Method for handling missing data ('skipna' or 'replace_zero')
        """
        self.min_detection_pct = min_detection_pct
        self.min_samples = min_samples
        self.missing_data_method = missing_data_method

        # Validate configuration
        if not 0 < min_detection_pct <= 1:
            raise ValueError(f"min_detection_pct must be between 0 and 1, got {min_detection_pct}")
        if min_samples < 1:
            raise ValueError(f"min_samples must be >= 1, got {min_samples}")
        if missing_data_method not in ['skipna', 'replace_zero']:
            raise ValueError(f"missing_data_method must be 'skipna' or 'replace_zero', got {missing_data_method}")

        logger.info(f"Data preparation config: detection_pct={min_detection_pct*100:.0f}%, "
                   f"min_samples={min_samples}, missing_data={missing_data_method}")


def calculate_detection_statistics(df: pd.DataFrame,
                                   sample_cols: List[str],
                                   group_name: str = "Group") -> Dict[str, pd.Series]:
    """
    Calculate detection statistics for a group of samples

    STANDARDIZED METHOD: Used by all visualizations to ensure consistency

    Args:
        df: DataFrame with sample columns
        sample_cols: List of sample column names
        group_name: Name of the group (for logging)

    Returns:
        Dictionary with:
        - 'count': Number of detected samples per glycopeptide
        - 'detection_pct': Detection percentage per glycopeptide
        - 'has_detection': Boolean mask for glycopeptides with any detection
    """
    # Count non-empty, non-zero values
    # Empty strings, NaN, and 0 are all considered "not detected"
    detection_count = pd.Series(0, index=df.index)

    for col in sample_cols:
        # Check if value exists, is not empty string, and is not zero
        detected = (df[col].notna()) & (df[col] != '') & (pd.to_numeric(df[col], errors='coerce') > 0)
        detection_count += detected.astype(int)

    detection_pct = detection_count / len(sample_cols)
    has_detection = detection_count > 0

    logger.debug(f"{group_name}: {has_detection.sum()}/{len(df)} glycopeptides with detection, "
                f"mean detection rate: {detection_pct.mean()*100:.1f}%")

    return {
        'count': detection_count,
        'detection_pct': detection_pct,
        'has_detection': has_detection
    }


def calculate_group_statistics_standardized(df: pd.DataFrame,
                                            sample_cols: List[str],
                                            method: str = 'skipna',
                                            group_name: str = "Group") -> Dict[str, pd.Series]:
    """
    Calculate standardized statistics for a group of samples

    CRITICAL: This is the SINGLE SOURCE OF TRUTH for mean calculations
    All visualizations MUST use this function to ensure consistency

    Args:
        df: DataFrame with sample columns
        sample_cols: List of sample column names
        method: 'skipna' (exclude missing) or 'replace_zero' (include as zero)
        group_name: Name of the group (for logging)

    Returns:
        Dictionary with:
        - 'mean': Mean intensity (NaN if no detection)
        - 'std': Standard deviation
        - 'min': Minimum value
        - 'max': Maximum value
        - 'count': Number of detected samples
        - 'sum': Sum of intensities
        - 'detection_pct': Detection percentage
    """
    result = {
        'mean': pd.Series(dtype=float, index=df.index),
        'std': pd.Series(dtype=float, index=df.index),
        'min': pd.Series(dtype=float, index=df.index),
        'max': pd.Series(dtype=float, index=df.index),
        'count': pd.Series(dtype=int, index=df.index),
        'sum': pd.Series(dtype=float, index=df.index),
        'detection_pct': pd.Series(dtype=float, index=df.index)
    }

    for idx in df.index:
        row = df.loc[idx, sample_cols]

        # Convert to numeric, handling empty strings and non-numeric values
        numeric_values = pd.to_numeric(row, errors='coerce')

        if method == 'skipna':
            # RECOMMENDED: Exclude missing values (NaN, empty, 0) from calculations
            # This is scientifically correct for MNAR (Missing Not At Random) data
            valid_values = numeric_values[numeric_values > 0]

            if len(valid_values) > 0:
                result['mean'][idx] = valid_values.mean()
                result['std'][idx] = valid_values.std()
                result['min'][idx] = valid_values.min()
                result['max'][idx] = valid_values.max()
                result['count'][idx] = len(valid_values)
                result['sum'][idx] = valid_values.sum()
            else:
                result['mean'][idx] = np.nan
                result['std'][idx] = np.nan
                result['min'][idx] = np.nan
                result['max'][idx] = np.nan
                result['count'][idx] = 0
                result['sum'][idx] = 0.0

        elif method == 'replace_zero':
            # LEGACY: Replace missing with 0 and include in calculations
            # WARNING: This underestimates intensity for low-detection glycopeptides
            filled_values = numeric_values.fillna(0)

            result['mean'][idx] = filled_values.mean()
            result['std'][idx] = filled_values.std()
            result['min'][idx] = filled_values.min()
            result['max'][idx] = filled_values.max()
            result['count'][idx] = (filled_values > 0).sum()
            result['sum'][idx] = filled_values.sum()

        result['detection_pct'][idx] = result['count'][idx] / len(sample_cols)

    logger.debug(f"{group_name} stats: mean_intensity={result['mean'].mean():.2e}, "
                f"mean_detection={result['detection_pct'].mean()*100:.1f}%")

    return result


def filter_by_detection_frequency(df: pd.DataFrame,
                                  config: DataPreparationConfig,
                                  cancer_samples: List[str] = None,
                                  normal_samples: List[str] = None,
                                  log_prefix: str = "") -> pd.DataFrame:
    """
    Apply standardized detection frequency filter

    CRITICAL: This is the SINGLE SOURCE OF TRUTH for detection filtering
    All visualizations MUST use this function to ensure consistency

    Args:
        df: DataFrame to filter
        config: DataPreparationConfig with filtering parameters
        cancer_samples: List of cancer sample column names (auto-detected if None)
        normal_samples: List of normal sample column names (auto-detected if None)
        log_prefix: Prefix for log messages

    Returns:
        Filtered DataFrame
    """
    # Auto-detect sample columns if not provided
    if cancer_samples is None or normal_samples is None:
        cancer_samples, normal_samples = get_sample_columns(df)

    total_before = len(df)

    # Calculate detection statistics for both groups
    cancer_stats = calculate_detection_statistics(df, cancer_samples, "Cancer")
    normal_stats = calculate_detection_statistics(df, normal_samples, "Normal")

    # Maximum detection across both groups
    max_detection_pct = pd.concat([cancer_stats['detection_pct'],
                                   normal_stats['detection_pct']], axis=1).max(axis=1)

    # Maximum detection count across both groups
    max_detection_count = pd.concat([cancer_stats['count'],
                                    normal_stats['count']], axis=1).max(axis=1)

    # Apply filter: require min_detection_pct OR min_samples in at least one group
    # Using OR (not AND) to be more inclusive while still maintaining quality
    filter_mask = (max_detection_pct >= config.min_detection_pct) | \
                  (max_detection_count >= config.min_samples)

    df_filtered = df[filter_mask].copy()
    total_after = len(df_filtered)
    removed = total_before - total_after

    logger.info(f"{log_prefix}Detection filter (≥{config.min_detection_pct*100:.0f}% OR "
               f"≥{config.min_samples} samples in at least one group):")
    logger.info(f"  Before: {total_before} glycopeptides")
    logger.info(f"  After: {total_after} glycopeptides")
    logger.info(f"  Removed: {removed} ({removed/total_before*100:.1f}%)")

    if total_after == 0:
        logger.warning(f"{log_prefix}WARNING: No glycopeptides pass detection filter!")

    return df_filtered


def prepare_visualization_data(df: pd.DataFrame,
                               config: DataPreparationConfig,
                               vip_scores: Optional[pd.DataFrame] = None,
                               merge_method: str = 'left',
                               apply_detection_filter: bool = True,
                               log_prefix: str = "") -> pd.DataFrame:
    """
    Centralized data preparation for visualizations

    CRITICAL: This is the SINGLE SOURCE OF TRUTH for visualization data preparation

    Pipeline:
    1. Merge with VIP scores (if provided)
    2. Calculate standardized statistics for Cancer and Normal groups
    3. Apply detection frequency filter (if enabled)
    4. Add derived metrics (fold change, log2 fold change)

    Args:
        df: Annotated DataFrame with sample columns
        config: DataPreparationConfig
        vip_scores: Optional VIP scores DataFrame to merge
        merge_method: 'left' (keep all) or 'inner' (keep only VIP-scored)
        apply_detection_filter: Whether to apply detection filter
        log_prefix: Prefix for log messages

    Returns:
        Prepared DataFrame with standardized statistics
    """
    logger.info(f"{log_prefix}Preparing visualization data...")
    logger.info(f"  Input: {len(df)} glycopeptides")

    df_prep = df.copy()

    # Step 1: Merge with VIP scores if provided
    if vip_scores is not None:
        logger.info(f"  Merging with VIP scores ({merge_method} join)...")
        merge_cols = ['Peptide', 'GlycanComposition']

        if merge_method == 'left':
            # Keep all glycopeptides from df, add VIP scores where available
            df_prep = df_prep.merge(vip_scores, on=merge_cols, how='left')
            # Fill missing VIP scores with 0
            if 'VIP_Score' in df_prep.columns:
                df_prep['VIP_Score'] = df_prep['VIP_Score'].fillna(0)
            logger.info(f"    After merge: {len(df_prep)} glycopeptides (all kept)")
        elif merge_method == 'inner':
            # Keep only glycopeptides with VIP scores
            before_merge = len(df_prep)
            df_prep = df_prep.merge(vip_scores, on=merge_cols, how='inner')
            logger.info(f"    After merge: {len(df_prep)} glycopeptides "
                       f"({before_merge - len(df_prep)} removed without VIP scores)")
        else:
            raise ValueError(f"merge_method must be 'left' or 'inner', got {merge_method}")

    # Step 2: Get sample columns
    cancer_samples, normal_samples = get_sample_columns(df_prep)
    logger.info(f"  Cancer samples: {len(cancer_samples)}, Normal samples: {len(normal_samples)}")

    # Step 3: Calculate standardized statistics
    logger.info(f"  Calculating statistics (method: {config.missing_data_method})...")

    cancer_stats = calculate_group_statistics_standardized(
        df_prep, cancer_samples, method=config.missing_data_method, group_name="Cancer"
    )
    normal_stats = calculate_group_statistics_standardized(
        df_prep, normal_samples, method=config.missing_data_method, group_name="Normal"
    )

    # Add statistics to dataframe
    df_prep['Cancer_Mean'] = cancer_stats['mean']
    df_prep['Cancer_StdDev'] = cancer_stats['std']
    df_prep['Cancer_Min'] = cancer_stats['min']
    df_prep['Cancer_Max'] = cancer_stats['max']
    df_prep['Cancer_SampleCount'] = cancer_stats['count']
    df_prep['Cancer_Sum'] = cancer_stats['sum']
    df_prep['Cancer_Detection_Pct'] = cancer_stats['detection_pct']

    df_prep['Normal_Mean'] = normal_stats['mean']
    df_prep['Normal_StdDev'] = normal_stats['std']
    df_prep['Normal_Min'] = normal_stats['min']
    df_prep['Normal_Max'] = normal_stats['max']
    df_prep['Normal_SampleCount'] = normal_stats['count']
    df_prep['Normal_Sum'] = normal_stats['sum']
    df_prep['Normal_Detection_Pct'] = normal_stats['detection_pct']

    # Maximum detection across both groups
    df_prep['Max_Detection_Pct'] = df_prep[['Cancer_Detection_Pct', 'Normal_Detection_Pct']].max(axis=1)
    df_prep['Max_SampleCount'] = df_prep[['Cancer_SampleCount', 'Normal_SampleCount']].max(axis=1)

    # Step 4: Apply detection filter if enabled
    if apply_detection_filter:
        df_prep = filter_by_detection_frequency(
            df_prep, config, cancer_samples, normal_samples, log_prefix=log_prefix
        )

    # Step 5: Calculate derived metrics (fold change)
    logger.info(f"  Calculating fold change...")

    # Linear fold change (Cancer / Normal)
    # Handle division by zero: use pseudocount of 1
    df_prep['Fold_Change'] = (df_prep['Cancer_Mean'] + 1) / (df_prep['Normal_Mean'] + 1)

    # Log2 fold change (symmetric, better for visualization)
    df_prep['Log2_Fold_Change'] = np.log2(df_prep['Fold_Change'])

    # Regulation direction
    df_prep['Regulation_Direction'] = 'Unchanged'
    df_prep.loc[df_prep['Log2_Fold_Change'] > 0, 'Regulation_Direction'] = 'Up in Cancer'
    df_prep.loc[df_prep['Log2_Fold_Change'] < 0, 'Regulation_Direction'] = 'Down in Cancer'

    logger.info(f"{log_prefix}Data preparation complete:")
    logger.info(f"  Final: {len(df_prep)} glycopeptides")
    logger.info(f"  Mean Cancer intensity: {df_prep['Cancer_Mean'].mean():.2e}")
    logger.info(f"  Mean Normal intensity: {df_prep['Normal_Mean'].mean():.2e}")
    logger.info(f"  Mean detection rate: {df_prep['Max_Detection_Pct'].mean()*100:.1f}%")

    return df_prep


def get_standard_config_from_dict(config_dict: dict) -> DataPreparationConfig:
    """
    Create DataPreparationConfig from configuration dictionary

    Args:
        config_dict: Configuration dictionary (from config.yaml)

    Returns:
        DataPreparationConfig instance
    """
    # Extract parameters with defaults
    detection_filter = config_dict.get('analysis', {}).get('detection_filter', {})
    missing_data = config_dict.get('analysis', {}).get('missing_data_handling', {})

    min_detection_pct = detection_filter.get('min_detection_pct', 0.30)
    min_samples = detection_filter.get('min_samples', 5)
    method = missing_data.get('method', 'skipna')

    return DataPreparationConfig(
        min_detection_pct=min_detection_pct,
        min_samples=min_samples,
        missing_data_method=method
    )


def calculate_statistical_significance(df_prep: pd.DataFrame,
                                       cancer_samples: List[str],
                                       normal_samples: List[str],
                                       method: str = 'mannwhitneyu',
                                       fdr_correction: bool = True) -> pd.DataFrame:
    """
    Calculate statistical significance for differential expression

    Args:
        df_prep: Prepared DataFrame from prepare_visualization_data()
        cancer_samples: List of cancer sample columns
        normal_samples: List of normal sample columns
        method: Statistical test ('mannwhitneyu' or 'ttest')
        fdr_correction: Apply Benjamini-Hochberg FDR correction

    Returns:
        DataFrame with p-values and FDR values
    """
    from scipy import stats as scipy_stats

    logger.info(f"Calculating statistical significance ({method})...")

    p_values = []

    for idx, row in df_prep.iterrows():
        # Get values for both groups (excluding zeros/missing)
        cancer_vals = pd.to_numeric(row[cancer_samples], errors='coerce')
        normal_vals = pd.to_numeric(row[normal_samples], errors='coerce')

        cancer_nonzero = cancer_vals[cancer_vals > 0].dropna()
        normal_nonzero = normal_vals[normal_vals > 0].dropna()

        # Require at least 3 samples in each group for valid test
        if len(cancer_nonzero) >= 3 and len(normal_nonzero) >= 3:
            try:
                if method == 'mannwhitneyu':
                    stat, p_val = scipy_stats.mannwhitneyu(
                        cancer_nonzero, normal_nonzero, alternative='two-sided'
                    )
                elif method == 'ttest':
                    stat, p_val = scipy_stats.ttest_ind(
                        cancer_nonzero, normal_nonzero, equal_var=False
                    )
                else:
                    raise ValueError(f"Unknown method: {method}")

                p_values.append(p_val)
            except Exception as e:
                logger.warning(f"Statistical test failed for {row.get('Peptide', 'unknown')}: {e}")
                p_values.append(np.nan)
        else:
            p_values.append(np.nan)

    df_prep['P_Value'] = p_values
    df_prep['-Log10_P_Value'] = -np.log10(df_prep['P_Value'].replace(0, 1e-300))

    # FDR correction
    if fdr_correction:
        from statsmodels.stats.multitest import multipletests

        valid_p = df_prep['P_Value'].dropna()
        if len(valid_p) > 0:
            _, fdr_values, _, _ = multipletests(valid_p.values, method='fdr_bh')

            # Map back to original dataframe
            df_prep['FDR'] = np.nan
            df_prep.loc[valid_p.index, 'FDR'] = fdr_values
            df_prep['-Log10_FDR'] = -np.log10(df_prep['FDR'].replace(0, 1e-300))

            logger.info(f"  Significant (FDR < 0.05): {(df_prep['FDR'] < 0.05).sum()} glycopeptides")
        else:
            df_prep['FDR'] = np.nan
            df_prep['-Log10_FDR'] = np.nan

    logger.info(f"  Valid p-values: {(~df_prep['P_Value'].isna()).sum()}/{len(df_prep)}")

    return df_prep
