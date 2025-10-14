"""
Test PreprocessingTracker Integration

Verifies that the PreprocessingTracker correctly records all transformation
steps during analysis and can save/load state for reproducibility.

Author: pGlyco Auto Combine Pipeline
Created: 2025-10-14
"""

import numpy as np
import pandas as pd
import sys
import tempfile
import json
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.analyzer import GlycanAnalyzer
from src.preprocessing_tracker import PreprocessingState


def create_test_dataset():
    """Create a simple test dataset for preprocessing tracking"""
    # Create 20 samples (10 cancer, 10 normal)
    sample_cols = [f'C{i+1}' for i in range(10)] + [f'N{i+1}' for i in range(10)]

    # Create 50 glycopeptides with realistic intensity values
    np.random.seed(42)
    n_features = 50

    data = {
        'Peptide': [f'PEPTIDE_{i}' for i in range(n_features)],
        'GlycanComposition': [f'H{i}N{i}' for i in range(n_features)]
    }

    # Add intensity columns with some variation
    for sample in sample_cols:
        # Cancer samples have higher mean intensity
        if sample.startswith('C'):
            data[sample] = np.random.lognormal(mean=5, sigma=1, size=n_features)
        else:
            data[sample] = np.random.lognormal(mean=4.5, sigma=1, size=n_features)

    return pd.DataFrame(data)


def test_preprocessing_tracker_basic():
    """Test that preprocessing tracker records transformations"""
    print("\n" + "="*80)
    print("TEST 1: Basic Preprocessing Tracking")
    print("="*80)

    # Create test dataset
    df = create_test_dataset()
    print(f"\nCreated test dataset: {df.shape}")

    # Initialize analyzer with tracking enabled
    analyzer = GlycanAnalyzer(n_components=2, log_transform=True, track_preprocessing=True)

    # Run PCA (this will trigger all preprocessing steps)
    print("\nRunning PCA with preprocessing tracking...")
    pca_results = analyzer.perform_pca(df)

    # Verify tracker exists and has recorded steps
    assert analyzer.preprocessing_tracker is not None, "Tracker should exist"

    state = analyzer.preprocessing_tracker.state
    assert state.tic_normalized, "TIC normalization should be marked"
    assert state.log2_transformed, "Log2 transformation should be marked"
    assert state.scaled, "Scaling should be marked"

    # Check transformation history
    history = state.transformation_history
    print(f"\nRecorded {len(history)} transformation steps:")
    for step in history:
        print(f"  {step.step_number}. {step.name}")

    assert len(history) >= 3, "Should have at least 3 transformation steps"

    # Check that transformation names are recorded
    step_names = [step.name for step in history]
    assert "TIC Normalization" in step_names, "TIC Normalization should be recorded"
    assert "Log2 Transformation" in step_names, "Log2 Transformation should be recorded"
    assert "RobustScaler" in step_names, "RobustScaler should be recorded"

    # Check sample counts
    assert state.n_samples_cancer == 10, "Should have 10 cancer samples"
    assert state.n_samples_normal == 10, "Should have 10 normal samples"
    assert state.n_samples_total == 20, "Should have 20 total samples"

    # Check glycopeptide counts
    assert state.n_glycopeptides_original == 50, "Should have 50 glycopeptides"

    print("\n✅ TEST 1 PASSED: Preprocessing tracking records all transformations")


def test_preprocessing_tracker_save_load():
    """Test that preprocessing state can be saved and loaded"""
    print("\n" + "="*80)
    print("TEST 2: Save and Load Preprocessing State")
    print("="*80)

    # Create test dataset
    df = create_test_dataset()

    # Initialize analyzer and run PCA
    analyzer = GlycanAnalyzer(n_components=2, log_transform=True, track_preprocessing=True)
    pca_results = analyzer.perform_pca(df)

    # Save preprocessing state to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as tmp:
        tmp_path = tmp.name

    print(f"\nSaving preprocessing state to: {tmp_path}")
    analyzer.save_preprocessing_state(tmp_path)

    # Load and verify JSON
    with open(tmp_path, 'r') as f:
        saved_data = json.load(f)

    print("\nLoaded JSON keys:", list(saved_data.keys()))

    # Verify structure
    assert 'transformations' in saved_data, "Should have transformations section"
    assert 'parameters' in saved_data, "Should have parameters section"
    assert 'statistics' in saved_data, "Should have statistics section"
    assert 'history' in saved_data, "Should have history section"
    assert 'metadata' in saved_data, "Should have metadata section"

    # Verify transformations
    trans = saved_data['transformations']
    assert trans['tic_normalized'] is True
    assert trans['log2_transformed'] is True
    assert trans['scaled'] is True

    # Verify parameters
    params = saved_data['parameters']
    assert params['log_transform_pseudocount'] == 1.0
    assert params['log_transform_base'] == 2
    assert params['scaler_type'] == 'RobustScaler'
    assert params['normalization_method'] == 'tic'

    # Verify statistics
    stats = saved_data['statistics']
    assert stats['n_samples_cancer'] == 10
    assert stats['n_samples_normal'] == 10
    assert stats['n_samples_total'] == 20
    assert stats['n_glycopeptides_original'] == 50

    # Verify history
    history = saved_data['history']
    assert len(history) >= 3, "Should have at least 3 transformation steps"

    # Load state from file
    print("\nLoading preprocessing state from file...")
    loaded_state = PreprocessingState.load_from_file(tmp_path)

    assert loaded_state.tic_normalized is True
    assert loaded_state.log2_transformed is True
    assert loaded_state.scaled is True
    assert loaded_state.n_samples_total == 20
    assert len(loaded_state.transformation_history) >= 3

    # Cleanup
    Path(tmp_path).unlink()

    print("\n✅ TEST 2 PASSED: Preprocessing state can be saved and loaded correctly")


def test_alphapeptstats_compatibility():
    """Test that preprocessing state can be exported in alphapeptstats format"""
    print("\n" + "="*80)
    print("TEST 3: Alphapeptstats Compatibility")
    print("="*80)

    # Create test dataset
    df = create_test_dataset()

    # Initialize analyzer and run PCA
    analyzer = GlycanAnalyzer(n_components=2, log_transform=True, track_preprocessing=True)
    pca_results = analyzer.perform_pca(df)

    # Get alphapeptstats-style dictionary
    alphapeptstats_dict = analyzer.get_alphapeptstats_preprocessing_dict()

    print("\nAlphapeptstats-style dictionary:")
    for key, value in alphapeptstats_dict.items():
        print(f"  {key}: {value}")

    # Verify expected keys (alphapeptstats compatible)
    expected_keys = [
        'log2_transformed',
        'normalized',
        'scaled',
        'filtered',
        'imputed',
        'normalization_method',
        'scaler',
        'log_pseudocount',
        'n_features_original',
        'n_features_filtered',
        'n_samples',
        'missing_data_method'
    ]

    for key in expected_keys:
        assert key in alphapeptstats_dict, f"Key '{key}' should be in alphapeptstats dict"

    # Verify values
    assert alphapeptstats_dict['log2_transformed'] is True
    assert alphapeptstats_dict['normalized'] is True
    assert alphapeptstats_dict['scaled'] is True
    assert alphapeptstats_dict['normalization_method'] == 'tic'
    assert alphapeptstats_dict['scaler'] == 'RobustScaler'
    assert alphapeptstats_dict['log_pseudocount'] == 1.0
    assert alphapeptstats_dict['n_samples'] == 20

    print("\n✅ TEST 3 PASSED: Alphapeptstats-compatible export works correctly")


def test_preprocessing_summary():
    """Test that preprocessing summary can be generated"""
    print("\n" + "="*80)
    print("TEST 4: Preprocessing Summary Generation")
    print("="*80)

    # Create test dataset
    df = create_test_dataset()

    # Initialize analyzer and run PCA
    analyzer = GlycanAnalyzer(n_components=2, log_transform=True, track_preprocessing=True)
    pca_results = analyzer.perform_pca(df)

    # Print summary
    print("\n--- PREPROCESSING SUMMARY ---")
    analyzer.print_preprocessing_summary()
    print("--- END SUMMARY ---")

    # Get summary as string
    summary_str = analyzer.preprocessing_tracker.state.summary()

    # Verify summary contains key information
    assert "TIC Normalization" in summary_str
    assert "Log2 Transformation" in summary_str
    assert "RobustScaler" in summary_str
    assert "Cancer samples: 10" in summary_str
    assert "Normal samples: 10" in summary_str
    assert "Total samples: 20" in summary_str

    print("\n✅ TEST 4 PASSED: Preprocessing summary generation works correctly")


def test_tracking_disabled():
    """Test that analyzer works correctly with tracking disabled"""
    print("\n" + "="*80)
    print("TEST 5: Tracking Disabled")
    print("="*80)

    # Create test dataset
    df = create_test_dataset()

    # Initialize analyzer with tracking DISABLED
    analyzer = GlycanAnalyzer(n_components=2, log_transform=True, track_preprocessing=False)

    # Verify tracker is None
    assert analyzer.preprocessing_tracker is None, "Tracker should be None when disabled"

    # Run PCA (should work without tracking)
    print("\nRunning PCA with tracking disabled...")
    pca_results = analyzer.perform_pca(df)

    # Verify PCA still works
    assert 'pca_df' in pca_results
    assert len(pca_results['pca_df']) == 20, "Should still process 20 samples"

    # Get preprocessing dict should return empty dict
    preprocessing_dict = analyzer.get_preprocessing_dict()
    assert preprocessing_dict == {}, "Should return empty dict when tracking disabled"

    print("\n✅ TEST 5 PASSED: Analyzer works correctly with tracking disabled")


def run_all_tests():
    """Run all preprocessing tracker tests"""
    print("\n" + "="*80)
    print("PREPROCESSING TRACKER INTEGRATION TESTS")
    print("Testing: PreprocessingTracker + GlycanAnalyzer integration")
    print("="*80)

    try:
        test_preprocessing_tracker_basic()
        test_preprocessing_tracker_save_load()
        test_alphapeptstats_compatibility()
        test_preprocessing_summary()
        test_tracking_disabled()

        print("\n" + "="*80)
        print("✅ ALL TESTS PASSED")
        print("="*80)
        print("\nPreprocessing tracking is working correctly!")
        print("The system successfully:")
        print("  1. Records all transformation steps during analysis")
        print("  2. Saves and loads preprocessing state to/from JSON")
        print("  3. Exports in alphapeptstats-compatible format")
        print("  4. Generates human-readable summaries")
        print("  5. Works correctly with tracking disabled")
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
