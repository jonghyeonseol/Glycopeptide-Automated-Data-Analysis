"""
Test S0 Parameter Implementation

Verifies that the S0 parameter correctly stabilizes t-statistics
for features with small variances.

Author: pGlyco Auto Combine Pipeline
Created: 2025-10-14
"""

import numpy as np
import pandas as pd
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.data_preparation import calculate_ttest_with_s0
from src.constants import S0_DEFAULT


def test_s0_basic_functionality():
    """Test that S0 parameter produces valid results"""
    print("\n" + "="*80)
    print("TEST 1: Basic S0 Functionality")
    print("="*80)

    # Create two groups with different means
    group1 = np.array([10.0, 12.0, 11.0, 10.5, 11.5])  # Cancer
    group2 = np.array([5.0, 6.0, 5.5, 5.8, 6.2])       # Normal

    t_stat, p_val, t_stat_s0, p_val_s0 = calculate_ttest_with_s0(
        group1, group2, s0=S0_DEFAULT
    )

    print(f"\nGroup 1 (Cancer): mean={np.mean(group1):.2f}, std={np.std(group1, ddof=1):.2f}")
    print(f"Group 2 (Normal): mean={np.mean(group2):.2f}, std={np.std(group2, ddof=1):.2f}")
    print(f"\nRegular t-test:")
    print(f"  t-statistic: {t_stat:.4f}")
    print(f"  p-value: {p_val:.6f}")
    print(f"\nS0-adjusted t-test:")
    print(f"  t-statistic (S0): {t_stat_s0:.4f}")
    print(f"  p-value (S0): {p_val_s0:.6f}")

    # Assertions
    assert t_stat != 0, "Regular t-statistic should not be zero"
    assert t_stat_s0 != 0, "S0-adjusted t-statistic should not be zero"
    assert p_val < 1.0, "P-value should be < 1.0"
    assert p_val_s0 < 1.0, "S0-adjusted p-value should be < 1.0"

    print("\n✅ TEST 1 PASSED: S0 produces valid results")


def test_s0_stabilizes_small_variance():
    """
    Test that S0 parameter stabilizes t-statistics for small variance cases

    CRITICAL TEST: This demonstrates the main benefit of S0 parameter.
    When variance is very small, regular t-test can produce artificially large
    t-statistics. S0 parameter prevents this inflation.
    """
    print("\n" + "="*80)
    print("TEST 2: S0 Stabilization for Small Variance")
    print("="*80)

    # Create two groups with SMALL variance but moderate difference
    group1 = np.array([10.01, 10.02, 10.01, 10.02, 10.01])  # Very small variance
    group2 = np.array([9.50, 9.51, 9.50, 9.51, 9.50])      # Very small variance

    t_stat, p_val, t_stat_s0, p_val_s0 = calculate_ttest_with_s0(
        group1, group2, s0=S0_DEFAULT
    )

    print(f"\nGroup 1 (Cancer): mean={np.mean(group1):.4f}, std={np.std(group1, ddof=1):.6f} (VERY SMALL)")
    print(f"Group 2 (Normal): mean={np.mean(group2):.4f}, std={np.std(group2, ddof=1):.6f} (VERY SMALL)")
    print(f"Mean difference: {np.mean(group1) - np.mean(group2):.4f}")

    print(f"\nRegular t-test:")
    print(f"  t-statistic: {t_stat:.4f} (INFLATED due to small variance)")
    print(f"  p-value: {p_val:.6f}")
    print(f"\nS0-adjusted t-test:")
    print(f"  t-statistic (S0): {t_stat_s0:.4f} (STABILIZED)")
    print(f"  p-value (S0): {p_val_s0:.6f}")

    # The key assertion: S0-adjusted should be smaller (more conservative)
    print(f"\nAssertion: S0-adjusted t-statistic should be smaller (more conservative)")
    print(f"  |t_stat| = {abs(t_stat):.4f}")
    print(f"  |t_stat_s0| = {abs(t_stat_s0):.4f}")
    print(f"  Ratio: {abs(t_stat_s0)/abs(t_stat):.4f}")

    assert abs(t_stat_s0) < abs(t_stat), "S0-adjusted should be more conservative for small variance"

    print("\n✅ TEST 2 PASSED: S0 successfully stabilizes small variance cases")


def test_s0_maintains_large_effects():
    """
    Test that S0 parameter does NOT penalize features with genuinely large effects

    When there's a large effect size with reasonable variance, S0 should have
    minimal impact on the t-statistic.
    """
    print("\n" + "="*80)
    print("TEST 3: S0 Preserves Large Effect Sizes")
    print("="*80)

    # Create two groups with LARGE difference and normal variance
    group1 = np.array([20.0, 22.0, 21.0, 19.5, 20.5])  # Cancer (high intensity)
    group2 = np.array([5.0, 6.0, 5.5, 4.8, 5.2])       # Normal (low intensity)

    t_stat, p_val, t_stat_s0, p_val_s0 = calculate_ttest_with_s0(
        group1, group2, s0=S0_DEFAULT
    )

    print(f"\nGroup 1 (Cancer): mean={np.mean(group1):.2f}, std={np.std(group1, ddof=1):.2f}")
    print(f"Group 2 (Normal): mean={np.mean(group2):.2f}, std={np.std(group2, ddof=1):.2f}")
    print(f"Mean difference: {np.mean(group1) - np.mean(group2):.2f} (LARGE EFFECT)")

    print(f"\nRegular t-test:")
    print(f"  t-statistic: {t_stat:.4f}")
    print(f"  p-value: {p_val:.6f}")
    print(f"\nS0-adjusted t-test:")
    print(f"  t-statistic (S0): {t_stat_s0:.4f}")
    print(f"  p-value (S0): {p_val_s0:.6f}")

    # For large effects, S0 should have minimal impact
    impact_ratio = abs(t_stat_s0) / abs(t_stat)
    print(f"\nImpact of S0 on large effect:")
    print(f"  Ratio: {impact_ratio:.4f} (should be close to 1.0)")

    assert impact_ratio > 0.9, "S0 should have minimal impact on large effects"
    assert p_val < 0.05, "Regular p-value should be significant"
    assert p_val_s0 < 0.05, "S0-adjusted p-value should also be significant"

    print("\n✅ TEST 3 PASSED: S0 preserves large effect sizes")


def test_s0_minimum_samples():
    """Test that function handles insufficient samples correctly"""
    print("\n" + "="*80)
    print("TEST 4: Minimum Sample Requirement")
    print("="*80)

    # Create groups with insufficient samples
    group1 = np.array([10.0, 12.0])  # Only 2 samples (< 3 required)
    group2 = np.array([5.0, 6.0])    # Only 2 samples (< 3 required)

    t_stat, p_val, t_stat_s0, p_val_s0 = calculate_ttest_with_s0(
        group1, group2, s0=S0_DEFAULT
    )

    print(f"\nGroup 1: {len(group1)} samples (< 3 required)")
    print(f"Group 2: {len(group2)} samples (< 3 required)")
    print(f"\nResults:")
    print(f"  t-statistic: {t_stat:.4f}")
    print(f"  p-value: {p_val:.4f}")
    print(f"  t-statistic (S0): {t_stat_s0:.4f}")
    print(f"  p-value (S0): {p_val_s0:.4f}")

    # Should return default values indicating insufficient data
    assert t_stat == 0.0, "t-statistic should be 0.0 for insufficient samples"
    assert p_val == 1.0, "p-value should be 1.0 for insufficient samples"

    print("\n✅ TEST 4 PASSED: Handles insufficient samples correctly")


def run_all_tests():
    """Run all S0 parameter tests"""
    print("\n" + "="*80)
    print("S0 PARAMETER IMPLEMENTATION TESTS")
    print("Testing: calculate_ttest_with_s0() function")
    print("="*80)

    try:
        test_s0_basic_functionality()
        test_s0_stabilizes_small_variance()
        test_s0_maintains_large_effects()
        test_s0_minimum_samples()

        print("\n" + "="*80)
        print("✅ ALL TESTS PASSED")
        print("="*80)
        print("\nS0 parameter implementation is working correctly!")
        print("The function successfully:")
        print("  1. Produces valid statistical results")
        print("  2. Stabilizes t-statistics for small variance cases")
        print("  3. Preserves large effect sizes")
        print("  4. Handles edge cases (insufficient samples)")
        print("\n" + "="*80)
        return True

    except AssertionError as e:
        print(f"\n❌ TEST FAILED: {e}")
        return False


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
