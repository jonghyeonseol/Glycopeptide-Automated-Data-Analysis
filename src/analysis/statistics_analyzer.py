"""
Statistics Analyzer Module
Handles statistical calculations and boxplot data preparation
"""

import pandas as pd
import numpy as np

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.utils import (
    replace_empty_with_zero,
    to_numeric_safe,
    get_sample_columns,
    get_all_sample_columns,
    get_sample_group,
    calculate_fold_change
)
from src.logger_config import get_logger
from .base_analyzer import BaseAnalyzer

logger = get_logger(__name__)


class StatisticsAnalyzer(BaseAnalyzer):
    """
    Statistical Analysis and Data Preparation

    Handles:
    - Statistics calculation by glycan type
    - Boxplot data preparation
    - Fold change calculations
    """

    def __init__(self, log_transform: bool = True):
        """
        Initialize StatisticsAnalyzer

        Args:
            log_transform: Whether to apply log transformation for boxplots
        """
        super().__init__(log_transform)

    def calculate_statistics_by_glycan_type(
        self,
        df: pd.DataFrame,
        group_col: str = 'GlycanType'
    ) -> pd.DataFrame:
        """
        Calculate statistics for each glycan type

        Args:
            df: Annotated DataFrame
            group_col: Column to group by (default: 'GlycanType')

        Returns:
            DataFrame with statistics per glycan type
        """
        # Get sample columns
        c_samples, n_samples = get_sample_columns(df)

        # Calculate statistics for each type
        stats_list = []

        for glycan_type in df[group_col].unique():
            subset = df[df[group_col] == glycan_type]

            # Convert to numeric and replace empty with 0
            c_values = replace_empty_with_zero(subset[c_samples])
            n_values = replace_empty_with_zero(subset[n_samples])

            # Calculate mean intensity for each sample
            c_mean = c_values.sum(axis=0).mean()
            n_mean = n_values.sum(axis=0).mean()

            # Calculate total intensity
            c_total = c_values.sum().sum()
            n_total = n_values.sum().sum()

            # Count glycopeptides
            count = len(subset)

            stats_list.append({
                'GlycanType': glycan_type,
                'Count': count,
                'Cancer_Mean': c_mean,
                'Normal_Mean': n_mean,
                'Cancer_Total': c_total,
                'Normal_Total': n_total,
                'Fold_Change': calculate_fold_change(c_mean, n_mean, log_scale=False),
                'Log2_Fold_Change': calculate_fold_change(c_mean, n_mean, log_scale=True)
            })

        stats_df = pd.DataFrame(stats_list)

        logger.info(f"\nStatistics by {group_col}:")
        logger.info(stats_df.to_string(index=False))

        return stats_df

    def prepare_boxplot_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Prepare data for boxplot visualization

        Args:
            df: Annotated DataFrame

        Returns:
            Long-format DataFrame for plotting
        """
        # Get sample columns
        sample_cols = get_all_sample_columns(df)

        # Melt to long format
        df_long = df.melt(
            id_vars=['GlycanType'],
            value_vars=sample_cols,
            var_name='Sample',
            value_name='Intensity'
        )

        # Convert intensity to numeric
        df_long['Intensity'] = to_numeric_safe(
            df_long['Intensity'].replace('', np.nan)
        )

        # Add group information
        df_long['Group'] = df_long['Sample'].apply(get_sample_group)

        # Log transform if specified
        if self.log_transform:
            df_long['Intensity'] = np.log2(df_long['Intensity'] + 1)

        # Remove zero intensities for better visualization
        df_long = df_long[df_long['Intensity'] > 0]

        logger.info(f"Boxplot data prepared: {len(df_long)} observations")

        return df_long

    def prepare_boxplot_data_extended(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Prepare data for extended boxplot with 5 categories

        Args:
            df: Annotated DataFrame

        Returns:
            Long-format DataFrame with extended categories
        """
        # Get sample columns
        sample_cols = get_all_sample_columns(df)

        # Determine extended category
        def determine_extended_category(row):
            if row['IsHighMannose']:
                return 'HM'
            elif row['IsComplexHybrid']:
                return 'C/H'
            elif row['IsSialylated'] and row['IsFucosylated']:
                return 'Sialofucosylated'
            elif row['IsSialylated']:
                return 'Sialylated'
            elif row['IsFucosylated']:
                return 'Fucosylated'
            else:
                return 'Other'

        df_with_extended = df.copy()
        df_with_extended['ExtendedCategory'] = df_with_extended.apply(
            determine_extended_category,
            axis=1
        )

        # Melt to long format
        df_long = df_with_extended.melt(
            id_vars=['ExtendedCategory'],
            value_vars=sample_cols,
            var_name='Sample',
            value_name='Intensity'
        )

        # Convert to numeric
        df_long['Intensity'] = to_numeric_safe(
            df_long['Intensity'].replace('', np.nan)
        )

        # Add group
        df_long['Group'] = df_long['Sample'].apply(get_sample_group)

        # Log transform if specified
        if self.log_transform:
            df_long['Intensity'] = np.log2(df_long['Intensity'] + 1)

        # Remove zeros
        df_long = df_long[df_long['Intensity'] > 0]

        logger.info(f"Extended boxplot data prepared: {len(df_long)} observations")

        return df_long
