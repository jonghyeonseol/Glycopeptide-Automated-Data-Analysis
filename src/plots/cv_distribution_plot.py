"""
CV Distribution Plot Module for pGlyco Auto Combine
Visualizes coefficient of variation for quality control

UPDATED: Now uses centralized data preparation for consistency
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
from ..utils import save_trace_data
from ..data_preparation import (
    DataPreparationConfig,
    calculate_group_statistics_standardized
)

logger = logging.getLogger(__name__)


class CVDistributionPlotMixin:
    """Mixin class for CV distribution visualization"""

    def plot_cv_distribution(self, df: pd.DataFrame, figsize: tuple = (14, 6)):
        """
        Create CV (Coefficient of Variation) distribution plots for Cancer and Normal samples

        Args:
            df: Annotated DataFrame with intensity data
            figsize: Figure size (width, height)
        """
        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]

        # STANDARDIZED: Calculate CV using centralized statistics
        config = DataPreparationConfig(missing_data_method='skipna')
        cv_data = []

        for idx, row in df.iterrows():
            peptide = row['Peptide']
            glycan_comp = row['GlycanComposition']
            glycopeptide = f"{peptide}_{glycan_comp}"

            # Get single row as DataFrame for standardized function
            glycopeptide_row = df.loc[[idx]]

            # Use standardized statistics calculation
            cancer_stats = calculate_group_statistics_standardized(
                glycopeptide_row, cancer_samples, method=config.missing_data_method
            )
            normal_stats = calculate_group_statistics_standardized(
                glycopeptide_row, normal_samples, method=config.missing_data_method
            )

            # Calculate CV (coefficient of variation = std / mean * 100)
            # Only calculate if we have at least 3 samples
            if cancer_stats['count'].iloc[0] >= 3:
                cancer_mean = cancer_stats['mean'].iloc[0]
                cancer_std = cancer_stats['std'].iloc[0]
                cancer_cv = (cancer_std / cancer_mean * 100) if cancer_mean > 0 else np.nan
            else:
                cancer_cv = np.nan

            if normal_stats['count'].iloc[0] >= 3:
                normal_mean = normal_stats['mean'].iloc[0]
                normal_std = normal_stats['std'].iloc[0]
                normal_cv = (normal_std / normal_mean * 100) if normal_mean > 0 else np.nan
            else:
                normal_cv = np.nan

            cv_data.append({
                'Glycopeptide': glycopeptide,
                'Peptide': peptide,
                'GlycanComposition': glycan_comp,
                'Cancer_CV': cancer_cv,
                'Normal_CV': normal_cv
            })

        # Create DataFrame
        cv_df = pd.DataFrame(cv_data)

        # Remove NaN values
        cv_df_clean = cv_df.dropna(subset=['Cancer_CV', 'Normal_CV'])

        # Create figure with two subplots
        fig, axes = plt.subplots(1, 2, figsize=figsize)

        # Define colors
        cancer_color = '#E74C3C'
        normal_color = '#3498DB'

        # Plot Cancer CV distribution
        ax1 = axes[0]
        ax1.hist(cv_df_clean['Cancer_CV'], bins=50, color=cancer_color,
                 alpha=0.7, edgecolor='black', linewidth=0.5)

        # Add median line
        cancer_median = cv_df_clean['Cancer_CV'].median()
        ax1.axvline(cancer_median, color='darkred', linestyle='--',
                    linewidth=2, label=f'Median: {cancer_median:.1f}%')

        ax1.set_xlabel('CV (%)', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Count', fontsize=12, fontweight='bold')
        ax1.set_title('CV Distribution: Cancer Samples', fontsize=13, fontweight='bold')
        ax1.legend(loc='upper right', frameon=True, fontsize=10)
        ax1.grid(True, alpha=0.3)

        # Add statistics text
        stats_text = f"n = {len(cv_df_clean)}\n"
        stats_text += f"Mean: {cv_df_clean['Cancer_CV'].mean():.1f}%\n"
        stats_text += f"Median: {cancer_median:.1f}%\n"
        q1 = cv_df_clean['Cancer_CV'].quantile(0.25)
        q3 = cv_df_clean['Cancer_CV'].quantile(0.75)
        stats_text += f"Q1-Q3: {q1:.1f}%-{q3:.1f}%"

        ax1.text(0.98, 0.98, stats_text, transform=ax1.transAxes,
                 fontsize=9, verticalalignment='top', horizontalalignment='right',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # Plot Normal CV distribution
        ax2 = axes[1]
        ax2.hist(cv_df_clean['Normal_CV'], bins=50, color=normal_color,
                 alpha=0.7, edgecolor='black', linewidth=0.5)

        # Add median line
        normal_median = cv_df_clean['Normal_CV'].median()
        ax2.axvline(normal_median, color='darkblue', linestyle='--',
                    linewidth=2, label=f'Median: {normal_median:.1f}%')

        ax2.set_xlabel('CV (%)', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Count', fontsize=12, fontweight='bold')
        ax2.set_title('CV Distribution: Normal Samples', fontsize=13, fontweight='bold')
        ax2.legend(loc='upper right', frameon=True, fontsize=10)
        ax2.grid(True, alpha=0.3)

        # Add statistics text
        stats_text = f"n = {len(cv_df_clean)}\n"
        stats_text += f"Mean: {cv_df_clean['Normal_CV'].mean():.1f}%\n"
        stats_text += f"Median: {normal_median:.1f}%\n"
        q1 = cv_df_clean['Normal_CV'].quantile(0.25)
        q3 = cv_df_clean['Normal_CV'].quantile(0.75)
        stats_text += f"Q1-Q3: {q1:.1f}%-{q3:.1f}%"

        ax2.text(0.98, 0.98, stats_text, transform=ax2.transAxes,
                 fontsize=9, verticalalignment='top', horizontalalignment='right',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # Main title
        fig.suptitle('Coefficient of Variation (CV) Distribution Analysis',
                     fontsize=15, fontweight='bold', y=1.02)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'cv_distribution.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved CV distribution plot to {output_file}")

        # Save trace data
        save_trace_data(cv_df_clean, self.output_dir, 'cv_distribution_data.csv')

        plt.close()

        return cv_df_clean
