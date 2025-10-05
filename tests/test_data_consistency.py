"""
Unit Tests for Data Consistency
Validates that centralized data preparation ensures consistency across visualizations
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from src.data_preparation import (
    DataPreparationConfig,
    calculate_group_statistics_standardized,
    prepare_visualization_data,
    filter_by_detection_frequency
)
from src.data_validator import DataConsistencyValidator, quick_consistency_check


class TestDataPreparationConfig:
    """Test DataPreparationConfig initialization and validation"""

    def test_default_config(self):
        """Test default configuration"""
        config = DataPreparationConfig()
        assert config.min_detection_pct == 0.30
        assert config.min_samples == 5
        assert config.missing_data_method == 'skipna'

    def test_custom_config(self):
        """Test custom configuration"""
        config = DataPreparationConfig(
            min_detection_pct=0.50,
            min_samples=10,
            missing_data_method='replace_zero'
        )
        assert config.min_detection_pct == 0.50
        assert config.min_samples == 10
        assert config.missing_data_method == 'replace_zero'

    def test_invalid_detection_pct(self):
        """Test that invalid detection percentage raises error"""
        with pytest.raises(ValueError):
            DataPreparationConfig(min_detection_pct=1.5)

    def test_invalid_method(self):
        """Test that invalid missing data method raises error"""
        with pytest.raises(ValueError):
            DataPreparationConfig(missing_data_method='invalid')


class TestGroupStatistics:
    """Test group statistics calculation"""

    @pytest.fixture
    def sample_data(self):
        """Create sample data for testing"""
        data = {
            'Peptide': ['PEP1', 'PEP2', 'PEP3'],
            'GlycanComposition': ['H(5)N(4)', 'H(6)N(5)A(1)', 'H(5)N(4)F(1)'],
            'C1': [100, 200, ''],  # Empty string (missing)
            'C2': [150, '', 300],
            'C3': [120, 250, 280],
            'N1': [80, 180, 250],
            'N2': ['', 190, 260],
            'N3': [90, 200, '']
        }
        return pd.DataFrame(data)

    def test_skipna_method(self, sample_data):
        """Test skipna method excludes missing values"""
        cancer_samples = ['C1', 'C2', 'C3']
        config = DataPreparationConfig(missing_data_method='skipna')

        stats = calculate_group_statistics_standardized(
            sample_data, cancer_samples, method=config.missing_data_method
        )

        # PEP1: mean of [100, 150, 120] = 123.33
        assert abs(stats['mean'][0] - 123.33) < 0.1

        # PEP2: mean of [200, 250] = 225 (C2 is missing)
        assert abs(stats['mean'][1] - 225.0) < 0.1

        # Detection counts should be correct
        assert stats['count'][0] == 3  # All 3 detected
        assert stats['count'][1] == 2  # 2 detected (C2 missing)

    def test_replace_zero_method(self, sample_data):
        """Test replace_zero method includes missing as 0"""
        cancer_samples = ['C1', 'C2', 'C3']
        config = DataPreparationConfig(missing_data_method='replace_zero')

        stats = calculate_group_statistics_standardized(
            sample_data, cancer_samples, method=config.missing_data_method
        )

        # PEP1: mean of [100, 150, 120] = 123.33
        assert abs(stats['mean'][0] - 123.33) < 0.1

        # PEP2: mean of [200, 0, 250] = 150 (C2 counted as 0)
        assert abs(stats['mean'][1] - 150.0) < 0.1

        # But detection counts still exclude zeros
        assert stats['count'][1] == 2  # Still 2 (zeros not counted as detection)


class TestDetectionFiltering:
    """Test detection frequency filtering"""

    @pytest.fixture
    def test_data(self):
        """Create test data with varying detection rates"""
        data = {
            'Peptide': ['PEP_HIGH', 'PEP_MED', 'PEP_LOW'],
            'GlycanComposition': ['H(5)N(4)', 'H(6)N(5)', 'H(7)N(6)'],
            # PEP_HIGH: 80% detection in cancer
            'C1': [100, 50, ''],
            'C2': [110, '', ''],
            'C3': [105, 60, ''],
            'C4': [120, '', ''],
            'C5': [115, 55, 10],
            # PEP_MED: 40% detection in cancer
            'N1': [80, 45, ''],
            'N2': ['', '', ''],
            'N3': [85, 50, ''],
            'N4': ['', '', ''],
            'N5': [90, 48, 8]
        }
        return pd.DataFrame(data)

    def test_30_percent_filter(self, test_data):
        """Test 30% detection filter"""
        config = DataPreparationConfig(min_detection_pct=0.30)
        filtered = filter_by_detection_frequency(test_data, config)

        # PEP_HIGH: 80% in C, 60% in N → PASS
        # PEP_MED: 40% in C, 60% in N → PASS
        # PEP_LOW: 20% in C, 20% in N → FAIL
        assert len(filtered) == 2
        assert 'PEP_HIGH' in filtered['Peptide'].values
        assert 'PEP_MED' in filtered['Peptide'].values

    def test_50_percent_filter(self, test_data):
        """Test 50% detection filter"""
        config = DataPreparationConfig(min_detection_pct=0.50)
        filtered = filter_by_detection_frequency(test_data, config)

        # PEP_HIGH: 80% in C, 60% in N → PASS
        # PEP_MED: 40% in C, 60% in N → PASS (60% in N)
        # PEP_LOW: 20% in C, 20% in N → FAIL
        assert len(filtered) == 2


class TestConsistencyValidation:
    """Test data consistency validation"""

    @pytest.fixture
    def consistent_datasets(self):
        """Create two consistent datasets"""
        data1 = {
            'Peptide': ['PEP1', 'PEP2'],
            'GlycanComposition': ['H(5)N(4)', 'H(6)N(5)'],
            'Cancer_Mean': [150.0, 200.0],
            'Normal_Mean': [100.0, 180.0],
            'Cancer_Detection_Pct': [0.8, 0.6],
            'Normal_Detection_Pct': [0.7, 0.5]
        }
        data2 = {
            'Peptide': ['PEP1', 'PEP2'],
            'GlycanComposition': ['H(5)N(4)', 'H(6)N(5)'],
            'Cancer_Mean': [150.0, 200.0],  # Same values
            'Normal_Mean': [100.0, 180.0],
            'Cancer_Detection_Pct': [0.8, 0.6],
            'Normal_Detection_Pct': [0.7, 0.5]
        }
        return pd.DataFrame(data1), pd.DataFrame(data2)

    @pytest.fixture
    def inconsistent_datasets(self):
        """Create two inconsistent datasets"""
        data1 = {
            'Peptide': ['PEP1', 'PEP2'],
            'GlycanComposition': ['H(5)N(4)', 'H(6)N(5)'],
            'Cancer_Mean': [150.0, 200.0],
            'Normal_Mean': [100.0, 180.0]
        }
        data2 = {
            'Peptide': ['PEP1', 'PEP2'],
            'GlycanComposition': ['H(5)N(4)', 'H(6)N(5)'],
            'Cancer_Mean': [120.0, 200.0],  # Different value for PEP1
            'Normal_Mean': [100.0, 180.0]
        }
        return pd.DataFrame(data1), pd.DataFrame(data2)

    def test_consistent_data(self, consistent_datasets):
        """Test that consistent data passes validation"""
        df1, df2 = consistent_datasets
        result = quick_consistency_check(df1, df2, "Dataset1", "Dataset2")
        assert result is True

    def test_inconsistent_data(self, inconsistent_datasets):
        """Test that inconsistent data fails validation"""
        df1, df2 = inconsistent_datasets

        # This should fail because means differ
        with pytest.raises(Exception):  # Will raise ValidationError
            validator = DataConsistencyValidator()
            validator.validate_intensity_consistency(
                df1, df2, 'Cancer_Mean', 'Dataset1', 'Dataset2', max_diff_pct=0.01
            )


class TestPieChartData:
    """Test pie chart data preparation"""

    @pytest.fixture
    def pie_chart_data(self):
        """Create sample data for pie chart testing"""
        data = {
            'Peptide': ['PEP1', 'PEP2', 'PEP3', 'PEP4'],
            'GlycanComposition': ['H(5)N(4)', 'H(6)N(5)A(1)', 'H(5)N(4)F(1)', 'H(6)N(5)A(1)F(1)'],
            'GlycanType': ['Non', 'Sialylated', 'Fucosylated', 'Both'],
            'PrimaryClassification': ['High Mannose', 'ComplexHybrid', 'ComplexHybrid', 'ComplexHybrid'],
            'SecondaryClassification': ['High Mannose', 'Sialylated', 'Fucosylated', 'Sialofucosylated'],
            'C1': [100, 200, 300, 400],
            'C2': [110, 210, 310, 410],
            'N1': [90, 180, 270, 360],
            'N2': [95, 190, 280, 370]
        }
        return pd.DataFrame(data)

    def test_glycan_type_totals(self, pie_chart_data):
        """Test glycan type total calculations"""
        cancer_samples = ['C1', 'C2']

        # Simulate pie chart calculation
        glycan_types = ['Non', 'Sialylated', 'Fucosylated', 'Both']
        totals = {}

        for glycan_type in glycan_types:
            mask = pie_chart_data['GlycanType'] == glycan_type
            total = pie_chart_data[mask][cancer_samples].sum().sum()
            totals[glycan_type] = total

        # Verify totals
        assert totals['Non'] == 210  # PEP1: 100+110
        assert totals['Sialylated'] == 410  # PEP2: 200+210
        assert totals['Fucosylated'] == 610  # PEP3: 300+310
        assert totals['Both'] == 810  # PEP4: 400+410


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])
