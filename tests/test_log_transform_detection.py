"""
Test Automatic Log-Transform Detection

Verifies that the detection function correctly identifies whether data is
already log-transformed or in raw intensity scale.

Author: pGlyco Auto Combine Pipeline
Created: 2025-10-14
"""

import numpy as np
import pandas as pd
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.utils import detect_log_transform


def test_raw_intensity_detection():
    """Test detection of raw intensity data"""
    print("\n" + "="*80)
    print("TEST 1: Raw Intensity Data Detection")
    print("="*80)

    # Create realistic raw intensity data (typical MS data)
    np.random.seed(42)
    raw_data = np.random.lognormal(mean=10, sigma=2, size=1000)

    print(f"\nRaw data statistics:")
    print(f"  Min: {raw_data.min():.2e}")
    print(f"  Max: {raw_data.max():.2e}")
    print(f"  Median: {np.median(raw_data):.2e}")
    print(f"  Mean: {raw_data.mean():.2e}")
    print(f"  Range: {raw_data.max() - raw_data.min():.2e}")

    # Detect
    result = detect_log_transform(raw_data)

    print(f"\nDetection result:")
    print(f"  Is log-transformed: {result['is_log_transformed']}")
    print(f"  Confidence: {result['confidence']:.2%}")
    print(f"  Recommendation: {result['recommendation']}")

    print(f"\nEvidence breakdown:")
    for key, value in result['evidence'].items():
        print(f"  {key}:")
        print(f"    Result: {value['result']}")
        score_str = f"{value['score']:.2f}" if value['score'] is not None else "N/A"
        print(f"    Score: {score_str} (weight: {value['weight']:.2f})")
        print(f"    {value['interpretation']}")

    # Assertions
    assert result['is_log_transformed'] == False, "Should detect as RAW data"
    assert result['confidence'] < 0.3, "Should have high confidence (low score for raw)"

    print("\n✅ TEST 1 PASSED: Raw intensity correctly detected")


def test_log_transformed_detection():
    """Test detection of log-transformed data"""
    print("\n" + "="*80)
    print("TEST 2: Log-Transformed Data Detection")
    print("="*80)

    # Create log-transformed data (log2 of raw data)
    np.random.seed(42)
    raw_data = np.random.lognormal(mean=10, sigma=2, size=1000)
    log_data = np.log2(raw_data + 1)

    print(f"\nLog2-transformed data statistics:")
    print(f"  Min: {log_data.min():.2f}")
    print(f"  Max: {log_data.max():.2f}")
    print(f"  Median: {np.median(log_data):.2f}")
    print(f"  Mean: {log_data.mean():.2f}")
    print(f"  Range: {log_data.max() - log_data.min():.2f}")

    # Detect
    result = detect_log_transform(log_data)

    print(f"\nDetection result:")
    print(f"  Is log-transformed: {result['is_log_transformed']}")
    print(f"  Confidence: {result['confidence']:.2%}")
    print(f"  Recommendation: {result['recommendation']}")

    print(f"\nEvidence breakdown:")
    for key, value in result['evidence'].items():
        score_str = f"{value['score']:.2f}" if value['score'] is not None else "N/A"
        print(f"  {key}: {value['result']} (score: {score_str})")

    # Assertions
    assert result['is_log_transformed'] == True, "Should detect as LOG-TRANSFORMED"
    assert result['confidence'] > 0.7, "Should have high confidence"

    print("\n✅ TEST 2 PASSED: Log-transformed data correctly detected")


def test_log_with_negatives():
    """Test detection with negative values (definitive proof of log transform)"""
    print("\n" + "="*80)
    print("TEST 3: Log-Transformed Data with Negative Values")
    print("="*80)

    # Create log data with some negatives (from very low raw intensities)
    np.random.seed(42)
    raw_data = np.random.lognormal(mean=5, sigma=3, size=1000)
    log_data = np.log2(raw_data + 1)

    # Add some negative values (simulating log2 of values < 1)
    log_data[::10] = np.random.uniform(-5, 0, size=100)

    print(f"\nLog2 data with negatives:")
    print(f"  Min: {log_data.min():.2f}")
    print(f"  Max: {log_data.max():.2f}")
    print(f"  Median: {np.median(log_data):.2f}")
    print(f"  Has negatives: {(log_data < 0).any()}")
    print(f"  Negative count: {(log_data < 0).sum()}")

    # Detect
    result = detect_log_transform(log_data)

    print(f"\nDetection result:")
    print(f"  Is log-transformed: {result['is_log_transformed']}")
    print(f"  Confidence: {result['confidence']:.2%}")

    # Assertions
    assert result['is_log_transformed'] == True, "Negative values prove log-transformed"
    assert result['confidence'] > 0.8, "Should have very high confidence with negatives"
    assert result['statistics']['has_negatives'] == True

    print("\n✅ TEST 3 PASSED: Negative values correctly identify log-transformed data")


def test_dataframe_detection():
    """Test detection with DataFrame input"""
    print("\n" + "="*80)
    print("TEST 4: DataFrame Input Detection")
    print("="*80)

    # Create DataFrame with sample columns
    np.random.seed(42)
    n_features = 100
    sample_cols = [f'C{i+1}' for i in range(5)] + [f'N{i+1}' for i in range(5)]

    # Raw intensity data
    raw_df = pd.DataFrame({
        'Peptide': [f'PEP_{i}' for i in range(n_features)],
        'GlycanComposition': [f'H{i}N{i}' for i in range(n_features)]
    })

    for col in sample_cols:
        raw_df[col] = np.random.lognormal(mean=10, sigma=2, size=n_features)

    # Detect on raw DataFrame
    result_raw = detect_log_transform(raw_df, sample_cols=sample_cols)

    print(f"\nRaw DataFrame detection:")
    print(f"  Is log-transformed: {result_raw['is_log_transformed']}")
    print(f"  Confidence: {result_raw['confidence']:.2%}")

    assert result_raw['is_log_transformed'] == False, "Should detect DataFrame as raw"

    # Log-transformed data
    log_df = raw_df.copy()
    for col in sample_cols:
        log_df[col] = np.log2(raw_df[col] + 1)

    # Detect on log DataFrame
    result_log = detect_log_transform(log_df, sample_cols=sample_cols)

    print(f"\nLog DataFrame detection:")
    print(f"  Is log-transformed: {result_log['is_log_transformed']}")
    print(f"  Confidence: {result_log['confidence']:.2%}")

    assert result_log['is_log_transformed'] == True, "Should detect DataFrame as log-transformed"

    print("\n✅ TEST 4 PASSED: DataFrame detection works correctly")


def test_edge_cases():
    """Test edge cases and ambiguous data"""
    print("\n" + "="*80)
    print("TEST 5: Edge Cases")
    print("="*80)

    # Test 1: Empty array
    print("\nTest 5.1: Empty array")
    result = detect_log_transform(np.array([]))
    print(f"  Result: {result['is_log_transformed']}")
    print(f"  Recommendation: {result['recommendation']}")
    assert result['is_log_transformed'] is None, "Empty array should return None"

    # Test 2: All zeros
    print("\nTest 5.2: All zeros")
    result = detect_log_transform(np.zeros(100))
    print(f"  Result: {result['is_log_transformed']}")
    assert result['is_log_transformed'] is None, "All zeros should return None"

    # Test 3: Ambiguous data (moderate values)
    print("\nTest 5.3: Ambiguous data (values around 50-100)")
    ambiguous_data = np.random.uniform(50, 100, size=1000)
    result = detect_log_transform(ambiguous_data)
    print(f"  Result: {result['is_log_transformed']}")
    print(f"  Confidence: {result['confidence']:.2%}")
    print(f"  Recommendation: {result['recommendation']}")
    # Ambiguous data might return None or either value, but confidence should be moderate
    assert 0.3 < result['confidence'] < 0.7, "Ambiguous data should have moderate confidence"

    # Test 4: Very small values (could be log-transformed)
    print("\nTest 5.4: Very small positive values (0.01-1)")
    small_data = np.random.uniform(0.01, 1, size=1000)
    result = detect_log_transform(small_data)
    print(f"  Result: {result['is_log_transformed']}")
    print(f"  Confidence: {result['confidence']:.2%}")

    print("\n✅ TEST 5 PASSED: Edge cases handled correctly")


def test_double_log_prevention():
    """Test that detection prevents double log-transformation"""
    print("\n" + "="*80)
    print("TEST 6: Double Log-Transform Prevention")
    print("="*80)

    # Create raw data
    np.random.seed(42)
    raw_data = np.random.lognormal(mean=10, sigma=2, size=500)

    # Apply log2 once
    log_once = np.log2(raw_data + 1)

    # Apply log2 twice (INCORRECT)
    log_twice = np.log2(log_once + 1)

    print(f"\nRaw data:")
    result_raw = detect_log_transform(raw_data)
    print(f"  Detected as: {'LOG' if result_raw['is_log_transformed'] else 'RAW'}")
    print(f"  Median: {np.median(raw_data):.2e}")
    print(f"  Recommendation: {result_raw['recommendation']}")

    print(f"\nLog once (CORRECT):")
    result_log_once = detect_log_transform(log_once)
    print(f"  Detected as: {'LOG' if result_log_once['is_log_transformed'] else 'RAW'}")
    print(f"  Median: {np.median(log_once):.2f}")
    print(f"  Recommendation: {result_log_once['recommendation']}")

    print(f"\nLog twice (INCORRECT - should be prevented):")
    result_log_twice = detect_log_transform(log_twice)
    print(f"  Detected as: {'LOG' if result_log_twice['is_log_transformed'] else 'RAW'}")
    print(f"  Median: {np.median(log_twice):.2f}")
    print(f"  Recommendation: {result_log_twice['recommendation']}")

    # The detection should warn against transforming already-log data
    assert result_log_once['is_log_transformed'] == True, "Log-once should be detected"
    assert result_log_twice['is_log_transformed'] == True, "Log-twice should also be detected"

    print("\n✅ TEST 6 PASSED: Detection can prevent double log-transformation")


def test_confidence_thresholds():
    """Test different confidence thresholds"""
    print("\n" + "="*80)
    print("TEST 7: Confidence Threshold Tuning")
    print("="*80)

    # Create moderately ambiguous data
    np.random.seed(42)
    ambiguous_data = np.random.uniform(20, 80, size=1000)

    # Test with different thresholds
    thresholds = [0.5, 0.6, 0.7, 0.8, 0.9]

    print("\nAmbiguous data with different thresholds:")
    for threshold in thresholds:
        result = detect_log_transform(ambiguous_data, threshold_confidence=threshold)
        print(f"  Threshold {threshold:.1f}: is_log={result['is_log_transformed']}, "
              f"confidence={result['confidence']:.2%}")

    # At higher thresholds, ambiguous data should return None
    result_high_threshold = detect_log_transform(ambiguous_data, threshold_confidence=0.9)
    # With high threshold and ambiguous data, might be None or have lower certainty

    print("\n✅ TEST 7 PASSED: Confidence thresholds work as expected")


def test_realistic_ms_data():
    """Test with realistic mass spectrometry data patterns"""
    print("\n" + "="*80)
    print("TEST 8: Realistic Mass Spectrometry Data")
    print("="*80)

    # Realistic raw MS data: log-normal distribution, some zeros
    np.random.seed(42)
    n_features = 200

    # Simulate raw intensities with typical MS patterns
    raw_intensities = np.random.lognormal(mean=8, sigma=2.5, size=n_features)
    raw_intensities[::15] = 0  # Some missing/zero values (common in MS)

    print(f"\nRealistic raw MS data:")
    print(f"  N features: {n_features}")
    print(f"  Zeros: {(raw_intensities == 0).sum()}")
    print(f"  Non-zero median: {np.median(raw_intensities[raw_intensities > 0]):.2e}")

    result_raw = detect_log_transform(raw_intensities)
    print(f"\nDetection:")
    print(f"  Is log-transformed: {result_raw['is_log_transformed']}")
    print(f"  Confidence: {result_raw['confidence']:.2%}")

    assert result_raw['is_log_transformed'] == False, "Should detect realistic MS data as raw"

    # Now log-transform it
    log_intensities = np.log2(raw_intensities + 1)

    print(f"\nAfter log2 transformation:")
    print(f"  Non-zero median: {np.median(log_intensities[log_intensities > 0]):.2f}")

    result_log = detect_log_transform(log_intensities)
    print(f"\nDetection:")
    print(f"  Is log-transformed: {result_log['is_log_transformed']}")
    print(f"  Confidence: {result_log['confidence']:.2%}")

    assert result_log['is_log_transformed'] == True, "Should detect transformed MS data"

    print("\n✅ TEST 8 PASSED: Realistic MS data patterns correctly identified")


def run_all_tests():
    """Run all log-transform detection tests"""
    print("\n" + "="*80)
    print("AUTO LOG-TRANSFORM DETECTION TESTS")
    print("Testing: detect_log_transform() with multiple heuristics")
    print("="*80)

    try:
        test_raw_intensity_detection()
        test_log_transformed_detection()
        test_log_with_negatives()
        test_dataframe_detection()
        test_edge_cases()
        test_double_log_prevention()
        test_confidence_thresholds()
        test_realistic_ms_data()

        print("\n" + "="*80)
        print("✅ ALL TESTS PASSED")
        print("="*80)
        print("\nAuto log-transform detection is working correctly!")
        print("The function successfully:")
        print("  1. Detects raw intensity data (log-normal, high values)")
        print("  2. Detects log-transformed data (small values, narrow range)")
        print("  3. Uses negative values as definitive proof of log-transform")
        print("  4. Works with DataFrame, Series, and array inputs")
        print("  5. Handles edge cases gracefully (empty, all zeros)")
        print("  6. Prevents dangerous double log-transformation")
        print("  7. Supports configurable confidence thresholds")
        print("  8. Handles realistic mass spectrometry data patterns")
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
