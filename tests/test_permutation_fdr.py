"""
Test Permutation-Based FDR Calculation

Verifies that permutation testing correctly estimates FDR and
identifies truly differential features.

Author: pGlyco Auto Combine Pipeline
Created: 2025-10-14
"""

import numpy as np
import pandas as pd
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.permutation_fdr import permutation_fdr, calculate_permutation_fdr, compare_fdr_methods


def create_test_data(n_features=100, n_cancer=10, n_normal=10, n_differential=20, effect_size=2.0, noise_std=1.0):
    """
    Create test data with known differential features

    Args:
        n_features: Total number of features
        n_cancer: Number of cancer samples
        n_normal: Number of normal samples
        n_differential: Number of truly differential features
        effect_size: Fold change for differential features
        noise_std: Standard deviation of noise

    Returns:
        Tuple of (cancer_df, normal_df, true_differential_indices)
    """
    np.random.seed(42)

    # Create baseline data (all features)
    cancer_data = np.random.normal(loc=10, scale=noise_std, size=(n_features, n_cancer))
    normal_data = np.random.normal(loc=10, scale=noise_std, size=(n_features, n_normal))

    # Add differential expression to first n_differential features
    cancer_data[:n_differential, :] += effect_size

    # Create DataFrames
    cancer_df = pd.DataFrame(cancer_data, columns=[f'C{i+1}' for i in range(n_cancer)])
    normal_df = pd.DataFrame(normal_data, columns=[f'N{i+1}' for i in range(n_normal)])

    true_differential = set(range(n_differential))

    return cancer_df, normal_df, true_differential


def test_basic_permutation_fdr():
    """Test basic permutation FDR calculation"""
    print("\n" + "="*80)
    print("TEST 1: Basic Permutation FDR")
    print("="*80)

    # Create test data with 20 differential features out of 100
    cancer_df, normal_df, true_diff = create_test_data(
        n_features=100, n_cancer=10, n_normal=10,
        n_differential=20, effect_size=3.0
    )

    print(f"\nTest data:")
    print(f"  Features: {len(cancer_df)}")
    print(f"  Cancer samples: {len(cancer_df.columns)}")
    print(f"  Normal samples: {len(normal_df.columns)}")
    print(f"  True differential: {len(true_diff)}")

    # Run permutation FDR
    results = permutation_fdr(
        cancer_df, normal_df,
        n_permutations=50,  # Use fewer for speed in testing
        test_statistic='ttest_s0',
        fdr_threshold=0.05,
        random_seed=42,
        show_progress=False
    )

    print(f"\nResults:")
    print(f"  Significant features (FDR < 0.05): {results['significant'].sum()}")
    print(f"  Mean test statistic (significant): {results[results['significant']]['test_statistic'].mean():.2f}")
    print(f"  Mean test statistic (non-significant): {results[~results['significant']]['test_statistic'].mean():.2f}")

    # Check recovery of true differential features
    detected_indices = set(results[results['significant']]['feature_index'].values)
    true_positives = len(detected_indices & true_diff)
    false_positives = len(detected_indices - true_diff)
    false_negatives = len(true_diff - detected_indices)

    sensitivity = true_positives / len(true_diff) if len(true_diff) > 0 else 0
    precision = true_positives / len(detected_indices) if len(detected_indices) > 0 else 0

    print(f"\nPerformance:")
    print(f"  True positives: {true_positives}/{len(true_diff)}")
    print(f"  False positives: {false_positives}")
    print(f"  False negatives: {false_negatives}")
    print(f"  Sensitivity (recall): {sensitivity:.2%}")
    print(f"  Precision: {precision:.2%}")

    # Assertions
    assert results['significant'].sum() > 0, "Should find some significant features"
    assert sensitivity > 0.5, "Should recover at least 50% of true differential features"
    assert precision > 0.5, "Precision should be > 50%"

    print("\n✅ TEST 1 PASSED: Basic permutation FDR works correctly")


def test_no_differential_features():
    """Test that permutation FDR doesn't find false positives when there are none"""
    print("\n" + "="*80)
    print("TEST 2: No Differential Features (Null Case)")
    print("="*80)

    # Create data with NO differential features
    np.random.seed(42)
    n_features = 50
    n_samples = 10

    cancer_df = pd.DataFrame(
        np.random.normal(loc=10, scale=1, size=(n_features, n_samples)),
        columns=[f'C{i+1}' for i in range(n_samples)]
    )
    normal_df = pd.DataFrame(
        np.random.normal(loc=10, scale=1, size=(n_features, n_samples)),
        columns=[f'N{i+1}' for i in range(n_samples)]
    )

    print(f"\nNull data (no true differences):")
    print(f"  Features: {n_features}")
    print(f"  Samples per group: {n_samples}")

    # Run permutation FDR
    results = permutation_fdr(
        cancer_df, normal_df,
        n_permutations=50,
        fdr_threshold=0.05,
        random_seed=42,
        show_progress=False
    )

    n_significant = results['significant'].sum()
    false_positive_rate = n_significant / n_features

    print(f"\nResults:")
    print(f"  Significant features: {n_significant}/{n_features}")
    print(f"  False positive rate: {false_positive_rate:.2%}")

    # With proper FDR control, we expect ~5% false positives at FDR 0.05
    assert false_positive_rate <= 0.10, "False positive rate should be ≤ 10% (allowing some variance)"

    print("\n✅ TEST 2 PASSED: Correctly controls false positives in null case")


def test_all_differential_features():
    """Test with all features being differential"""
    print("\n" + "="*80)
    print("TEST 3: All Features Differential")
    print("="*80)

    # Create data where ALL features are differential
    cancer_df, normal_df, true_diff = create_test_data(
        n_features=50, n_cancer=10, n_normal=10,
        n_differential=50, effect_size=3.0  # ALL features differential
    )

    print(f"\nTest data (all differential):")
    print(f"  Total features: {len(cancer_df)}")
    print(f"  Differential features: {len(true_diff)}")

    # Run permutation FDR
    results = permutation_fdr(
        cancer_df, normal_df,
        n_permutations=50,
        fdr_threshold=0.05,
        random_seed=42,
        show_progress=False
    )

    n_significant = results['significant'].sum()
    sensitivity = n_significant / len(true_diff)

    print(f"\nResults:")
    print(f"  Significant features: {n_significant}/{len(true_diff)}")
    print(f"  Sensitivity: {sensitivity:.2%}")

    # Should recover most or all truly differential features
    assert sensitivity > 0.8, "Should recover > 80% when all features are differential"

    print("\n✅ TEST 3 PASSED: High sensitivity when all features differential")


def test_different_effect_sizes():
    """Test sensitivity to different effect sizes"""
    print("\n" + "="*80)
    print("TEST 4: Different Effect Sizes")
    print("="*80)

    effect_sizes = [1.0, 2.0, 3.0, 4.0]
    sensitivities = []

    for effect_size in effect_sizes:
        cancer_df, normal_df, true_diff = create_test_data(
            n_features=100, n_cancer=10, n_normal=10,
            n_differential=20, effect_size=effect_size
        )

        results = permutation_fdr(
            cancer_df, normal_df,
            n_permutations=50,
            fdr_threshold=0.05,
            random_seed=42,
            show_progress=False
        )

        detected = set(results[results['significant']]['feature_index'].values)
        tp = len(detected & true_diff)
        sensitivity = tp / len(true_diff)
        sensitivities.append(sensitivity)

        print(f"  Effect size {effect_size:.1f}: Sensitivity = {sensitivity:.2%}")

    # Sensitivity should increase with effect size
    for i in range(len(sensitivities) - 1):
        assert sensitivities[i+1] >= sensitivities[i] - 0.1, \
            "Sensitivity should generally increase with effect size"

    print("\n✅ TEST 4 PASSED: Sensitivity increases with effect size")


def test_permutation_p_values():
    """Test that permutation p-values are reasonable"""
    print("\n" + "="*80)
    print("TEST 5: Permutation P-Values")
    print("="*80)

    cancer_df, normal_df, true_diff = create_test_data(
        n_features=100, n_cancer=10, n_normal=10,
        n_differential=20, effect_size=3.0
    )

    results = permutation_fdr(
        cancer_df, normal_df,
        n_permutations=50,
        fdr_threshold=0.05,
        random_seed=42,
        show_progress=False
    )

    # Check p-value properties
    perm_pvals = results['p_value_permutation'].values
    param_pvals = results['p_value'].values

    print(f"\nPermutation p-values:")
    print(f"  Min: {np.min(perm_pvals):.4f}")
    print(f"  Max: {np.max(perm_pvals):.4f}")
    print(f"  Median: {np.median(perm_pvals):.4f}")

    print(f"\nParametric p-values:")
    print(f"  Min: {np.nanmin(param_pvals):.4f}")
    print(f"  Max: {np.nanmax(param_pvals):.4f}")
    print(f"  Median: {np.nanmedian(param_pvals):.4f}")

    # P-values should be in [0, 1]
    assert np.all((perm_pvals >= 0) & (perm_pvals <= 1)), "P-values must be in [0, 1]"

    # Correlation between permutation and parametric p-values should be positive
    valid_mask = ~np.isnan(param_pvals)
    if np.sum(valid_mask) > 10:
        correlation = np.corrcoef(perm_pvals[valid_mask], param_pvals[valid_mask])[0, 1]
        print(f"\nCorrelation (perm vs param p-values): {correlation:.3f}")
        assert correlation > 0.5, "Permutation and parametric p-values should be correlated"

    print("\n✅ TEST 5 PASSED: Permutation p-values are reasonable")


def test_fdr_comparison():
    """Test comparison between BH FDR and Permutation FDR"""
    print("\n" + "="*80)
    print("TEST 6: FDR Method Comparison")
    print("="*80)

    cancer_df, normal_df, true_diff = create_test_data(
        n_features=100, n_cancer=10, n_normal=10,
        n_differential=20, effect_size=2.5
    )

    print(f"\nComparing FDR methods...")

    # Compare methods
    comparison = compare_fdr_methods(
        cancer_df, normal_df,
        n_permutations=50,
        fdr_threshold=0.05,
        random_seed=42
    )

    n_perm_sig = comparison['significant'].sum()
    n_bh_sig = comparison['significant_bh'].sum()

    print(f"\nComparison:")
    print(f"  Permutation FDR: {n_perm_sig} significant")
    print(f"  Benjamini-Hochberg: {n_bh_sig} significant")

    # Both methods should find some significant features
    assert n_perm_sig > 0, "Permutation FDR should find some significant features"
    assert n_bh_sig > 0, "BH FDR should find some significant features"

    # Results shouldn't be wildly different
    ratio = max(n_perm_sig, n_bh_sig) / max(min(n_perm_sig, n_bh_sig), 1)
    print(f"  Ratio: {ratio:.2f}")
    assert ratio < 3.0, "Methods shouldn't differ by more than 3x"

    print("\n✅ TEST 6 PASSED: FDR methods produce comparable results")


def test_reproducibility():
    """Test that results are reproducible with same random seed"""
    print("\n" + "="*80)
    print("TEST 7: Reproducibility")
    print("="*80)

    cancer_df, normal_df, _ = create_test_data(
        n_features=50, n_cancer=10, n_normal=10,
        n_differential=10, effect_size=2.0
    )

    # Run twice with same seed
    results1 = permutation_fdr(
        cancer_df, normal_df,
        n_permutations=30,
        random_seed=123,
        show_progress=False
    )

    results2 = permutation_fdr(
        cancer_df, normal_df,
        n_permutations=30,
        random_seed=123,
        show_progress=False
    )

    # Results should be identical
    assert np.allclose(results1['test_statistic'].values, results2['test_statistic'].values), \
        "Test statistics should be identical with same seed"
    assert np.allclose(results1['fdr_permutation'].values, results2['fdr_permutation'].values), \
        "FDR values should be identical with same seed"
    assert (results1['significant'] == results2['significant']).all(), \
        "Significance calls should be identical with same seed"

    print("\n✅ TEST 7 PASSED: Results are reproducible with same random seed")


def run_all_tests():
    """Run all permutation FDR tests"""
    print("\n" + "="*80)
    print("PERMUTATION-BASED FDR CALCULATION TESTS")
    print("Testing: permutation_fdr() with empirical null distribution")
    print("="*80)

    try:
        test_basic_permutation_fdr()
        test_no_differential_features()
        test_all_differential_features()
        test_different_effect_sizes()
        test_permutation_p_values()
        test_fdr_comparison()
        test_reproducibility()

        print("\n" + "="*80)
        print("✅ ALL TESTS PASSED")
        print("="*80)
        print("\nPermutation-based FDR calculation is working correctly!")
        print("The function successfully:")
        print("  1. Identifies truly differential features with good sensitivity")
        print("  2. Controls false positive rate in null case")
        print("  3. Scales sensitivity with effect size")
        print("  4. Generates reasonable permutation p-values")
        print("  5. Produces results comparable to Benjamini-Hochberg")
        print("  6. Provides reproducible results with random seed")
        print("  7. Handles edge cases (no diff, all diff)")
        print("\nPermutation FDR is more accurate than parametric methods for:")
        print("  - Small sample sizes")
        print("  - Non-normal distributions")
        print("  - Unknown dependence structures")
        print("\n" + "="*80)
        return True

    except AssertionError as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
