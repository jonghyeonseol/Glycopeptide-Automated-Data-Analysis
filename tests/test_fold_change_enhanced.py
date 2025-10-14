"""
Test Enhanced Fold Change Calculation

Verifies robust zero handling, vectorized operations, and numerical stability
in fold change calculations. Enhanced version prevents inf values and
handles edge cases gracefully.

Author: pGlyco Auto Combine Pipeline
Created: 2025-10-14
"""

import numpy as np
import pandas as pd
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.utils import calculate_fold_change


def test_basic_fold_change():
    """Test basic fold change calculations"""
    print("\n" + "="*80)
    print("TEST 1: Basic Fold Change Calculations")
    print("="*80)

    # Test 1: Simple 2-fold increase
    fc = calculate_fold_change(200, 100, log_scale=False, pseudocount=0)
    log2fc = calculate_fold_change(200, 100, log_scale=True, pseudocount=0)

    print(f"\n200 / 100:")
    print(f"  Linear FC: {fc:.3f} (expected: 2.0)")
    print(f"  Log2 FC: {log2fc:.3f} (expected: 1.0)")

    assert np.isclose(fc, 2.0), "Linear FC should be 2.0"
    assert np.isclose(log2fc, 1.0), "Log2 FC should be 1.0"

    # Test 2: Simple 2-fold decrease
    fc = calculate_fold_change(100, 200, log_scale=False, pseudocount=0)
    log2fc = calculate_fold_change(100, 200, log_scale=True, pseudocount=0)

    print(f"\n100 / 200:")
    print(f"  Linear FC: {fc:.3f} (expected: 0.5)")
    print(f"  Log2 FC: {log2fc:.3f} (expected: -1.0)")

    assert np.isclose(fc, 0.5), "Linear FC should be 0.5"
    assert np.isclose(log2fc, -1.0), "Log2 FC should be -1.0"

    print("\n✅ TEST 1 PASSED: Basic fold change calculations work correctly")


def test_zero_handling():
    """Test robust zero handling"""
    print("\n" + "="*80)
    print("TEST 2: Zero Handling")
    print("="*80)

    # Test 1: Division by zero (both zero)
    fc = calculate_fold_change(0, 0, log_scale=False, pseudocount=1.0)
    log2fc = calculate_fold_change(0, 0, log_scale=True, pseudocount=1.0)

    print(f"\n0 / 0 (with pseudocount=1):")
    print(f"  Linear FC: {fc:.3f} (expected: 1.0, no change)")
    print(f"  Log2 FC: {log2fc:.3f} (expected: 0.0, no change)")
    print(f"  Contains inf: {np.isinf(fc) or np.isinf(log2fc)}")

    assert not np.isinf(fc), "Linear FC should not be inf"
    assert not np.isinf(log2fc), "Log2 FC should not be inf"
    assert np.isclose(fc, 1.0), "0/0 with pseudocount should be 1.0"
    assert np.isclose(log2fc, 0.0), "Log2(0/0) with pseudocount should be 0.0"

    # Test 2: Numerator zero, denominator non-zero
    fc = calculate_fold_change(0, 100, log_scale=False, pseudocount=1.0)
    log2fc = calculate_fold_change(0, 100, log_scale=True, pseudocount=1.0)

    print(f"\n0 / 100 (with pseudocount=1):")
    print(f"  Linear FC: {fc:.3f} (expected: ~0.01)")
    print(f"  Log2 FC: {log2fc:.3f} (expected: negative)")
    print(f"  Contains inf: {np.isinf(fc) or np.isinf(log2fc)}")

    assert not np.isinf(fc), "Linear FC should not be inf"
    assert not np.isinf(log2fc), "Log2 FC should not be inf"
    assert fc < 1.0, "FC should be < 1 (decrease)"
    assert log2fc < 0, "Log2 FC should be negative (decrease)"

    # Test 3: Denominator zero, numerator non-zero
    fc = calculate_fold_change(100, 0, log_scale=False, pseudocount=1.0)
    log2fc = calculate_fold_change(100, 0, log_scale=True, pseudocount=1.0)

    print(f"\n100 / 0 (with pseudocount=1):")
    print(f"  Linear FC: {fc:.3f} (expected: large, ~101)")
    print(f"  Log2 FC: {log2fc:.3f} (expected: positive)")
    print(f"  Contains inf: {np.isinf(fc) or np.isinf(log2fc)}")

    assert not np.isinf(fc), "Linear FC should not be inf"
    assert not np.isinf(log2fc), "Log2 FC should not be inf"
    assert fc > 1.0, "FC should be > 1 (increase)"
    assert log2fc > 0, "Log2 FC should be positive (increase)"

    print("\n✅ TEST 2 PASSED: Zero handling prevents inf values")


def test_vectorized_operations():
    """Test vectorized fold change calculations"""
    print("\n" + "="*80)
    print("TEST 3: Vectorized Operations")
    print("="*80)

    # Create test data
    cancer_means = pd.Series([200, 100, 0, 150, 50])
    normal_means = pd.Series([100, 200, 0, 100, 100])

    print("\nInput data:")
    print(f"  Cancer means: {cancer_means.values}")
    print(f"  Normal means: {normal_means.values}")

    # Calculate vectorized fold changes
    fc_linear = calculate_fold_change(cancer_means, normal_means, log_scale=False, pseudocount=1.0)
    fc_log2 = calculate_fold_change(cancer_means, normal_means, log_scale=True, pseudocount=1.0)

    print("\nLinear fold changes:")
    print(f"  {fc_linear.values}")
    print(f"\nLog2 fold changes:")
    print(f"  {fc_log2.values}")

    # Verify no inf values
    assert not np.any(np.isinf(fc_linear)), "Linear FC should not contain inf"
    assert not np.any(np.isinf(fc_log2)), "Log2 FC should not contain inf"

    # Verify expected patterns
    assert fc_linear.iloc[0] > 1.0, "First FC should be > 1 (increase)"
    assert fc_linear.iloc[1] < 1.0, "Second FC should be < 1 (decrease)"
    assert np.isclose(fc_linear.iloc[2], 1.0), "Third FC should be 1 (both zero)"

    print("\n✅ TEST 3 PASSED: Vectorized operations work correctly")


def test_adaptive_pseudocount():
    """Test adaptive pseudocount based on data scale"""
    print("\n" + "="*80)
    print("TEST 4: Adaptive Pseudocount")
    print("="*80)

    # Test with high-intensity data
    cancer_high = pd.Series([10000, 20000, 30000])
    normal_high = pd.Series([8000, 15000, 25000])

    # Fixed pseudocount
    fc_fixed = calculate_fold_change(
        cancer_high, normal_high,
        log_scale=True,
        pseudocount=1.0,
        adaptive_pseudocount=False
    )

    # Adaptive pseudocount (1% of median)
    fc_adaptive = calculate_fold_change(
        cancer_high, normal_high,
        log_scale=True,
        pseudocount=1.0,  # Will be overridden
        adaptive_pseudocount=True
    )

    print("\nHigh-intensity data:")
    print(f"  Cancer: {cancer_high.values}")
    print(f"  Normal: {normal_high.values}")
    print(f"  Fixed pseudocount FC: {fc_fixed.values}")
    print(f"  Adaptive pseudocount FC: {fc_adaptive.values}")

    # Verify no inf values
    assert not np.any(np.isinf(fc_fixed)), "Fixed pseudocount should not produce inf"
    assert not np.any(np.isinf(fc_adaptive)), "Adaptive pseudocount should not produce inf"

    # Adaptive pseudocount should scale with data
    median_intensity = np.median(np.concatenate([cancer_high.values, normal_high.values]))
    adaptive_pseudo = median_intensity * 0.01
    print(f"\nAdaptive pseudocount: {adaptive_pseudo:.2f} (1% of median {median_intensity:.2f})")

    print("\n✅ TEST 4 PASSED: Adaptive pseudocount works correctly")


def test_edge_cases():
    """Test various edge cases"""
    print("\n" + "="*80)
    print("TEST 5: Edge Cases")
    print("="*80)

    # Test 1: Very small values
    fc = calculate_fold_change(1e-10, 1e-10, log_scale=False, pseudocount=1.0)
    print(f"\n1e-10 / 1e-10:")
    print(f"  Linear FC: {fc:.6f}")
    print(f"  Is inf: {np.isinf(fc)}")
    assert not np.isinf(fc), "Very small values should not produce inf"

    # Test 2: Very large values
    fc = calculate_fold_change(1e10, 1e10, log_scale=False, pseudocount=1.0)
    print(f"\n1e10 / 1e10:")
    print(f"  Linear FC: {fc:.6f}")
    print(f"  Is inf: {np.isinf(fc)}")
    assert not np.isinf(fc), "Very large values should not produce inf"

    # Test 3: Negative values (should still work with pseudocount)
    fc = calculate_fold_change(-100, -50, log_scale=False, pseudocount=1.0)
    print(f"\n-100 / -50 (with pseudocount):")
    print(f"  Linear FC: {fc:.6f}")
    print(f"  Is inf: {np.isinf(fc)}")
    assert not np.isinf(fc), "Negative values with pseudocount should not produce inf"

    # Test 4: NaN in Series
    cancer_with_nan = pd.Series([100, np.nan, 200])
    normal_with_nan = pd.Series([50, 100, np.nan])
    fc = calculate_fold_change(cancer_with_nan, normal_with_nan, log_scale=False, pseudocount=1.0)
    print(f"\nWith NaN values:")
    print(f"  Input cancer: {cancer_with_nan.values}")
    print(f"  Input normal: {normal_with_nan.values}")
    print(f"  Output FC: {fc.values}")
    # NaN propagates, which is expected behavior
    assert np.isnan(fc.iloc[1]) or np.isnan(fc.iloc[2]), "NaN should propagate"

    print("\n✅ TEST 5 PASSED: Edge cases handled correctly")


def test_symmetry():
    """Test that log2 fold change is symmetric"""
    print("\n" + "="*80)
    print("TEST 6: Log2 Fold Change Symmetry")
    print("="*80)

    # 2-fold increase
    fc_up = calculate_fold_change(200, 100, log_scale=True, pseudocount=0)

    # 2-fold decrease (reciprocal)
    fc_down = calculate_fold_change(100, 200, log_scale=True, pseudocount=0)

    print(f"\n2-fold increase (200/100): Log2 FC = {fc_up:.3f}")
    print(f"2-fold decrease (100/200): Log2 FC = {fc_down:.3f}")
    print(f"Symmetric: {fc_up:.3f} = -{fc_down:.3f}?")

    # Should be symmetric (equal magnitude, opposite sign)
    assert np.isclose(fc_up, -fc_down), "Log2 FC should be symmetric"
    assert np.isclose(fc_up, 1.0), "Log2(2) should be 1.0"
    assert np.isclose(fc_down, -1.0), "Log2(0.5) should be -1.0"

    print("\n✅ TEST 6 PASSED: Log2 fold change is symmetric")


def test_no_inf_in_real_data():
    """Test that realistic glycoproteomics data produces no inf values"""
    print("\n" + "="*80)
    print("TEST 7: Realistic Glycoproteomics Data")
    print("="*80)

    # Simulate realistic data with various patterns
    np.random.seed(42)
    n = 100

    cancer_means = np.random.lognormal(mean=5, sigma=2, size=n)
    normal_means = np.random.lognormal(mean=4.8, sigma=2, size=n)

    # Add some zeros (common in MS data)
    cancer_means[::10] = 0
    normal_means[::7] = 0

    cancer_series = pd.Series(cancer_means)
    normal_series = pd.Series(normal_means)

    print(f"\nSimulated data:")
    print(f"  N = {n} glycopeptides")
    print(f"  Cancer zeros: {(cancer_series == 0).sum()}")
    print(f"  Normal zeros: {(normal_series == 0).sum()}")
    print(f"  Both zero: {((cancer_series == 0) & (normal_series == 0)).sum()}")

    # Calculate fold changes
    fc_linear = calculate_fold_change(cancer_series, normal_series, log_scale=False, pseudocount=1.0)
    fc_log2 = calculate_fold_change(cancer_series, normal_series, log_scale=True, pseudocount=1.0)

    print(f"\nResults:")
    print(f"  Linear FC range: {fc_linear.min():.3f} - {fc_linear.max():.3f}")
    print(f"  Log2 FC range: {fc_log2.min():.3f} - {fc_log2.max():.3f}")
    print(f"  Linear FC inf count: {np.isinf(fc_linear).sum()}")
    print(f"  Log2 FC inf count: {np.isinf(fc_log2).sum()}")

    # Critical assertions
    assert not np.any(np.isinf(fc_linear)), "Linear FC should not contain inf"
    assert not np.any(np.isinf(fc_log2)), "Log2 FC should not contain inf"

    # Most values should be valid (non-NaN)
    valid_pct = (~np.isnan(fc_log2)).sum() / len(fc_log2) * 100
    print(f"  Valid FC values: {valid_pct:.1f}%")
    assert valid_pct > 80, "At least 80% of values should be valid"

    print("\n✅ TEST 7 PASSED: Realistic data produces no inf values")


def run_all_tests():
    """Run all fold change enhancement tests"""
    print("\n" + "="*80)
    print("ENHANCED FOLD CHANGE CALCULATION TESTS")
    print("Testing: calculate_fold_change() with robust zero handling")
    print("="*80)

    try:
        test_basic_fold_change()
        test_zero_handling()
        test_vectorized_operations()
        test_adaptive_pseudocount()
        test_edge_cases()
        test_symmetry()
        test_no_inf_in_real_data()

        print("\n" + "="*80)
        print("✅ ALL TESTS PASSED")
        print("="*80)
        print("\nEnhanced fold change calculation is working correctly!")
        print("The function successfully:")
        print("  1. Handles zero values without producing inf")
        print("  2. Supports vectorized operations for pandas Series")
        print("  3. Provides adaptive pseudocount based on data scale")
        print("  4. Handles edge cases gracefully")
        print("  5. Maintains log2 fold change symmetry")
        print("  6. Works correctly with realistic MS data")
        print("  7. Prevents numerical instability in downstream analyses")
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
