"""
Permutation-Based FDR Calculation for pGlyco Auto Combine

Implements permutation testing to estimate False Discovery Rate (FDR) empirically.
More accurate than parametric methods for small samples or non-normal data.

Inspired by alphapeptstats' permutation approach and SAM (Significance Analysis of Microarrays).

Author: pGlyco Auto Combine Pipeline
Created: 2025-10-14
"""

import numpy as np
import pandas as pd
from typing import Tuple, List, Optional, Callable
from scipy import stats
import logging

# Try to import tqdm for progress bars, but make it optional
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    tqdm = lambda x, **kwargs: x  # Fallback: no progress bar

from .constants import S0_DEFAULT, S0_MIN_SAMPLES
from .data_preparation import calculate_ttest_with_s0

logger = logging.getLogger(__name__)


def permutation_fdr(
    cancer_data: pd.DataFrame,
    normal_data: pd.DataFrame,
    n_permutations: int = 100,
    test_statistic: str = 'ttest_s0',
    s0: float = S0_DEFAULT,
    fdr_threshold: float = 0.05,
    min_samples: int = S0_MIN_SAMPLES,
    random_seed: Optional[int] = None,
    show_progress: bool = True
) -> pd.DataFrame:
    """
    Calculate FDR using permutation testing

    SCIENTIFIC RATIONALE:
    - Permutation testing makes NO parametric assumptions
    - Empirically estimates null distribution from data
    - More accurate for small samples than BH correction
    - Controls FDR by comparing real vs. permuted statistics

    ALGORITHM:
    1. Calculate real test statistics for all features
    2. For each permutation:
       a. Randomly shuffle group labels
       b. Recalculate test statistics
       c. Count how many exceed each real statistic
    3. FDR(t) = E[#permuted > t] / #real > t

    Args:
        cancer_data: DataFrame with cancer samples (features × samples)
        normal_data: DataFrame with normal samples (features × samples)
        n_permutations: Number of random permutations (default: 100)
                       More permutations = more accurate but slower
        test_statistic: 'ttest', 'ttest_s0', or 'fold_change' (default: 'ttest_s0')
        s0: S0 parameter for t-test stability (default: 0.05)
        fdr_threshold: FDR cutoff for significance calling (default: 0.05)
        min_samples: Minimum samples per group (default: 3)
        random_seed: Random seed for reproducibility (default: None)
        show_progress: Show progress bar (default: True)

    Returns:
        DataFrame with columns:
        - 'feature_index': Feature index
        - 'test_statistic': Real test statistic
        - 'p_value': Parametric p-value
        - 'p_value_permutation': Permutation-based p-value
        - 'fdr_permutation': Permutation-based FDR
        - 'significant': Boolean (FDR < fdr_threshold)

    Examples:
        >>> # Perform permutation FDR
        >>> results = permutation_fdr(cancer_df, normal_df, n_permutations=100)
        >>> sig_features = results[results['significant']]
        >>> print(f"Found {len(sig_features)} significant features")
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    # Validate inputs
    if cancer_data.shape[0] != normal_data.shape[0]:
        raise ValueError("Cancer and normal data must have same number of features")

    n_features = cancer_data.shape[0]
    n_cancer = cancer_data.shape[1]
    n_normal = normal_data.shape[1]

    logger.info(f"Permutation FDR: {n_features} features, {n_cancer} cancer, {n_normal} normal samples")
    logger.info(f"Running {n_permutations} permutations with test: {test_statistic}")

    # Combine data for permutation
    combined_data = pd.concat([cancer_data, normal_data], axis=1)

    # Calculate real test statistics
    logger.info("Calculating real test statistics...")
    real_stats = []
    real_pvalues = []

    for idx in range(n_features):
        cancer_vals = cancer_data.iloc[idx, :].values
        normal_vals = normal_data.iloc[idx, :].values

        # Remove NaN
        cancer_clean = cancer_vals[~np.isnan(cancer_vals)]
        normal_clean = normal_vals[~np.isnan(normal_vals)]

        if len(cancer_clean) >= min_samples and len(normal_clean) >= min_samples:
            if test_statistic == 'ttest_s0':
                _, _, t_stat_s0, p_val_s0 = calculate_ttest_with_s0(
                    cancer_clean, normal_clean, s0=s0
                )
                real_stats.append(abs(t_stat_s0))  # Use absolute value
                real_pvalues.append(p_val_s0)
            elif test_statistic == 'ttest':
                t_stat, p_val = stats.ttest_ind(cancer_clean, normal_clean, equal_var=False)
                real_stats.append(abs(t_stat))
                real_pvalues.append(p_val)
            elif test_statistic == 'fold_change':
                fc = np.mean(cancer_clean) / (np.mean(normal_clean) + 1e-10)
                log2fc = np.log2(fc + 1e-10)
                real_stats.append(abs(log2fc))
                real_pvalues.append(np.nan)  # No p-value for FC
            else:
                raise ValueError(f"Unknown test_statistic: {test_statistic}")
        else:
            real_stats.append(0.0)
            real_pvalues.append(1.0)

    real_stats = np.array(real_stats)
    real_pvalues = np.array(real_pvalues)

    # Permutation testing
    logger.info("Running permutations...")
    perm_stats_matrix = np.zeros((n_permutations, n_features))

    iterator = tqdm(range(n_permutations), desc="Permutations") if (show_progress and HAS_TQDM) else range(n_permutations)

    for perm_idx in iterator:
        # Shuffle column labels
        shuffled_indices = np.random.permutation(n_cancer + n_normal)

        # Split into permuted groups
        perm_cancer = combined_data.iloc[:, shuffled_indices[:n_cancer]]
        perm_normal = combined_data.iloc[:, shuffled_indices[n_cancer:]]

        # Calculate test statistics for permuted data
        for feature_idx in range(n_features):
            cancer_vals = perm_cancer.iloc[feature_idx, :].values
            normal_vals = perm_normal.iloc[feature_idx, :].values

            # Remove NaN
            cancer_clean = cancer_vals[~np.isnan(cancer_vals)]
            normal_clean = normal_vals[~np.isnan(normal_vals)]

            if len(cancer_clean) >= min_samples and len(normal_clean) >= min_samples:
                if test_statistic == 'ttest_s0':
                    _, _, t_stat_s0, _ = calculate_ttest_with_s0(
                        cancer_clean, normal_clean, s0=s0
                    )
                    perm_stats_matrix[perm_idx, feature_idx] = abs(t_stat_s0)
                elif test_statistic == 'ttest':
                    t_stat, _ = stats.ttest_ind(cancer_clean, normal_clean, equal_var=False)
                    perm_stats_matrix[perm_idx, feature_idx] = abs(t_stat)
                elif test_statistic == 'fold_change':
                    fc = np.mean(cancer_clean) / (np.mean(normal_clean) + 1e-10)
                    log2fc = np.log2(fc + 1e-10)
                    perm_stats_matrix[perm_idx, feature_idx] = abs(log2fc)

    # Calculate permutation-based p-values
    logger.info("Calculating permutation p-values...")
    perm_pvalues = np.zeros(n_features)

    for feature_idx in range(n_features):
        # Count how many permuted stats >= real stat
        count_exceeding = np.sum(perm_stats_matrix[:, feature_idx] >= real_stats[feature_idx])
        perm_pvalues[feature_idx] = (count_exceeding + 1) / (n_permutations + 1)

    # Calculate permutation-based FDR
    logger.info("Calculating permutation FDR...")
    perm_fdr = calculate_permutation_fdr(real_stats, perm_stats_matrix)

    # Create results DataFrame
    results = pd.DataFrame({
        'feature_index': range(n_features),
        'test_statistic': real_stats,
        'p_value': real_pvalues,
        'p_value_permutation': perm_pvalues,
        'fdr_permutation': perm_fdr,
        'significant': perm_fdr < fdr_threshold
    })

    n_significant = (perm_fdr < fdr_threshold).sum()
    logger.info(f"Permutation FDR complete: {n_significant}/{n_features} features significant at FDR < {fdr_threshold}")

    return results


def calculate_permutation_fdr(real_stats: np.ndarray, perm_stats_matrix: np.ndarray) -> np.ndarray:
    """
    Calculate FDR from permutation statistics

    FDR(t) = E[V(t)] / R(t)
    where:
    - V(t) = number of false positives (estimated from permutations)
    - R(t) = number of features called significant at threshold t

    Args:
        real_stats: Array of real test statistics (n_features,)
        perm_stats_matrix: Matrix of permuted statistics (n_permutations, n_features)

    Returns:
        Array of FDR values for each feature (n_features,)
    """
    n_permutations, n_features = perm_stats_matrix.shape

    # Sort real statistics in descending order
    sorted_indices = np.argsort(-real_stats)  # Descending
    sorted_stats = real_stats[sorted_indices]

    # Calculate FDR for each threshold
    fdr_values = np.zeros(n_features)

    for i, threshold in enumerate(sorted_stats):
        # Count real features above threshold
        n_real_above = np.sum(real_stats >= threshold)

        if n_real_above == 0:
            fdr_values[sorted_indices[i]] = 1.0
            continue

        # Count permuted features above threshold (across all permutations)
        n_perm_above = np.sum(perm_stats_matrix >= threshold)

        # Expected number of false positives
        expected_fp = n_perm_above / n_permutations

        # FDR = E[FP] / #positives
        fdr = expected_fp / n_real_above

        # FDR cannot exceed 1
        fdr = min(fdr, 1.0)

        fdr_values[sorted_indices[i]] = fdr

    # Enforce monotonicity (FDR cannot decrease as threshold relaxes)
    # This is the Benjamini-Hochberg step-up procedure adapted for permutations
    for i in range(n_features - 2, -1, -1):
        idx_current = sorted_indices[i]
        idx_next = sorted_indices[i + 1]
        fdr_values[idx_current] = min(fdr_values[idx_current], fdr_values[idx_next])

    return fdr_values


def compare_fdr_methods(
    cancer_data: pd.DataFrame,
    normal_data: pd.DataFrame,
    n_permutations: int = 100,
    s0: float = S0_DEFAULT,
    fdr_threshold: float = 0.05,
    random_seed: Optional[int] = None
) -> pd.DataFrame:
    """
    Compare Benjamini-Hochberg FDR vs. Permutation FDR

    Useful for understanding differences between parametric and non-parametric FDR.

    Args:
        cancer_data: DataFrame with cancer samples
        normal_data: DataFrame with normal samples
        n_permutations: Number of permutations
        s0: S0 parameter for t-test
        fdr_threshold: FDR cutoff
        random_seed: Random seed

    Returns:
        DataFrame comparing both methods
    """
    from statsmodels.stats.multitest import multipletests

    # Get permutation FDR
    perm_results = permutation_fdr(
        cancer_data, normal_data,
        n_permutations=n_permutations,
        test_statistic='ttest_s0',
        s0=s0,
        fdr_threshold=fdr_threshold,
        random_seed=random_seed,
        show_progress=True
    )

    # Calculate BH FDR
    valid_pvals = perm_results['p_value'].dropna()
    if len(valid_pvals) > 0:
        _, fdr_bh, _, _ = multipletests(valid_pvals.values, method='fdr_bh')
        perm_results.loc[valid_pvals.index, 'fdr_bh'] = fdr_bh
        perm_results['significant_bh'] = perm_results['fdr_bh'] < fdr_threshold
    else:
        perm_results['fdr_bh'] = np.nan
        perm_results['significant_bh'] = False

    # Comparison statistics
    n_perm_sig = perm_results['significant'].sum()
    n_bh_sig = perm_results['significant_bh'].sum()
    n_both_sig = ((perm_results['significant']) & (perm_results['significant_bh'])).sum()
    n_only_perm = ((perm_results['significant']) & (~perm_results['significant_bh'])).sum()
    n_only_bh = ((~perm_results['significant']) & (perm_results['significant_bh'])).sum()

    logger.info(f"\nFDR Method Comparison (threshold: {fdr_threshold}):")
    logger.info(f"  Permutation FDR: {n_perm_sig} significant")
    logger.info(f"  Benjamini-Hochberg: {n_bh_sig} significant")
    logger.info(f"  Agreement: {n_both_sig} features")
    logger.info(f"  Only permutation: {n_only_perm} features")
    logger.info(f"  Only BH: {n_only_bh} features")

    return perm_results
