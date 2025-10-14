"""
Test Systematic Infinity/NaN Handling

Verifies that all transformation steps properly clean inf/nan values
to prevent them from corrupting downstream analyses.

Author: pGlyco Auto Combine Pipeline
Created: 2025-10-14
"""

import numpy as np
import pandas as pd
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.utils import clean_inf_nan, safe_log_transform, safe_division


def test_clean_inf_nan_dataframe():
    """Test that clean_inf_nan removes inf values from DataFrame"""
    print("\n" + "="*80)
    print("TEST 1: clean_inf_nan() on DataFrame")
    print("="*80)

    # Create DataFrame with inf values
    df = pd.DataFrame({
        'A': [1.0, np.inf, 3.0, -np.inf],
        'B': [10.0, 20.0, np.inf, 40.0],
        'C': [100.0, 200.0, 300.0, 400.0]
    })

    print("\nOriginal DataFrame:")
    print(df)
    print(f"Inf count: {np.isinf(df.values).sum()}")

    # Clean inf/nan
    df_clean = clean_inf_nan(df, context="test", verbose=False)

    print("\nCleaned DataFrame:")
    print(df_clean)
    print(f"Inf count after cleaning: {np.isinf(df_clean.values).sum()}")
    print(f"NaN count after cleaning: {np.isnan(df_clean.values).sum()}")

    # Assertions
    assert np.isinf(df_clean.values).sum() == 0, "All inf values should be removed"
    assert np.isnan(df_clean.values).sum() == 3, "Inf values should become NaN"

    print("\n✅ TEST 1 PASSED: clean_inf_nan successfully removes inf values")


def test_safe_log_transform():
    """Test that safe_log_transform handles edge cases"""
    print("\n" + "="*80)
    print("TEST 2: safe_log_transform() Edge Cases")
    print("="*80)

    # Create data with potential problems
    data = pd.Series([0.0, 1.0, 10.0, 100.0, 1000.0])

    print("\nOriginal data:")
    print(data.values)

    # Apply safe log transform
    log_data = safe_log_transform(data, pseudocount=1.0, base=2)

    print("\nLog2-transformed data:")
    print(log_data.values)
    print(f"Inf count: {np.isinf(log_data.values).sum()}")
    print(f"NaN count: {np.isnan(log_data.values).sum()}")

    # Assertions
    assert np.isinf(log_data.values).sum() == 0, "No inf values should remain"
    assert not np.isnan(log_data.values[0]), "log2(0+1) = 0 should not be NaN"

    # Check specific values
    assert np.isclose(log_data.iloc[0], 0.0), "log2(0+1) = log2(1) = 0"
    assert np.isclose(log_data.iloc[1], 1.0), "log2(1+1) = log2(2) = 1"

    print("\n✅ TEST 2 PASSED: safe_log_transform handles edge cases correctly")


def test_safe_division():
    """Test that safe_division prevents division by zero"""
    print("\n" + "="*80)
    print("TEST 3: safe_division() Division by Zero Protection")
    print("="*80)

    # Test scalar division
    result1 = safe_division(10, 0, fill_value=np.nan)
    result2 = safe_division(10, 1e-15, fill_value=np.nan)  # Too small
    result3 = safe_division(10, 2, fill_value=np.nan)      # Normal

    print(f"\n10 / 0 = {result1} (expected: nan)")
    print(f"10 / 1e-15 = {result2} (expected: nan, denominator too small)")
    print(f"10 / 2 = {result3} (expected: 5.0)")

    assert np.isnan(result1), "Division by zero should return NaN"
    assert np.isnan(result2), "Division by very small number should return NaN"
    assert result3 == 5.0, "Normal division should work"

    # Test vectorized division
    numerators = np.array([10.0, 20.0, 30.0, 40.0])
    denominators = np.array([2.0, 0.0, 1e-15, 5.0])

    results = safe_division(numerators, denominators, fill_value=np.nan)

    print(f"\nVectorized division results:")
    print(f"  10 / 2 = {results[0]} (expected: 5.0)")
    print(f"  20 / 0 = {results[1]} (expected: nan)")
    print(f"  30 / 1e-15 = {results[2]} (expected: nan)")
    print(f"  40 / 5 = {results[3]} (expected: 8.0)")

    assert results[0] == 5.0, "First division should be 5.0"
    assert np.isnan(results[1]), "Second division should be NaN"
    assert np.isnan(results[2]), "Third division should be NaN"
    assert results[3] == 8.0, "Fourth division should be 8.0"

    print("\n✅ TEST 3 PASSED: safe_division prevents division by zero")


def test_transformation_pipeline():
    """Test that a full transformation pipeline cleans inf/nan at each step"""
    print("\n" + "="*80)
    print("TEST 4: Full Transformation Pipeline")
    print("="*80)

    # Simulate a problematic dataset
    df = pd.DataFrame({
        'Sample1': [100.0, 0.0, 500.0, 1e-20],
        'Sample2': [200.0, 300.0, 0.0, 600.0],
        'Sample3': [150.0, 250.0, 400.0, 0.0]
    })

    print("\nOriginal data (rows=glycopeptides, cols=samples):")
    print(df)

    # Step 1: TIC normalization (division)
    print("\nStep 1: TIC Normalization...")
    sample_sums = df.sum(axis=1)
    median_sum = sample_sums.median()
    sample_sums_safe = sample_sums.replace(0, 1)
    df_norm = df.div(sample_sums_safe, axis=0) * median_sum
    df_norm = clean_inf_nan(df_norm, context="TIC normalization", verbose=False)

    print(f"  After normalization - Inf count: {np.isinf(df_norm.values).sum()}")
    assert np.isinf(df_norm.values).sum() == 0, "No inf after normalization"

    # Step 2: Log2 transform
    print("\nStep 2: Log2 Transformation...")
    df_log = safe_log_transform(df_norm, pseudocount=1.0, base=2)

    print(f"  After log transform - Inf count: {np.isinf(df_log.values).sum()}")
    assert np.isinf(df_log.values).sum() == 0, "No inf after log transform"

    # Step 3: Fold change calculation
    print("\nStep 3: Fold Change Calculation...")
    mean1 = df_log['Sample1'].mean()
    mean2 = df_log['Sample2'].mean()

    # Division that might produce inf
    fc = safe_division(mean1, mean2, fill_value=np.nan)
    log2fc = safe_log_transform(pd.Series([fc]), pseudocount=0, base=2).iloc[0]

    print(f"  Fold change: {fc}")
    print(f"  Log2 fold change: {log2fc}")
    print(f"  Is inf: {np.isinf(log2fc)}")

    assert not np.isinf(log2fc), "Fold change should not be inf"

    print("\n✅ TEST 4 PASSED: Full pipeline maintains numerical stability")


def run_all_tests():
    """Run all inf/nan handling tests"""
    print("\n" + "="*80)
    print("SYSTEMATIC INF/NAN HANDLING TESTS")
    print("Testing: clean_inf_nan(), safe_log_transform(), safe_division()")
    print("="*80)

    try:
        test_clean_inf_nan_dataframe()
        test_safe_log_transform()
        test_safe_division()
        test_transformation_pipeline()

        print("\n" + "="*80)
        print("✅ ALL TESTS PASSED")
        print("="*80)
        print("\nSystematic inf/nan handling is working correctly!")
        print("The functions successfully:")
        print("  1. Remove inf values from DataFrames")
        print("  2. Handle log transformation edge cases")
        print("  3. Prevent division by zero")
        print("  4. Maintain numerical stability throughout pipeline")
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
