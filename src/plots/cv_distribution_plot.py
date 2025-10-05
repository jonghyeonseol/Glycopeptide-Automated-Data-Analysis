"""
CV Distribution Plot Module for pGlyco Auto Combine
Visualizes coefficient of variation for quality control
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from ..utils import replace_empty_with_zero, save_trace_data

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

        # Calculate CV for each glycopeptide
        cv_data = []

        for idx, row in df.iterrows():
            peptide = row['Peptide']
            glycan_comp = row['GlycanComposition']
            glycopeptide = f"{peptide}_{glycan_comp}"

            # Cancer CV
            cancer_values = replace_empty_with_zero(row[cancer_samples]).values.astype(float)
            cancer_nonzero = cancer_values[cancer_values > 0]

            if len(cancer_nonzero) >= 3:
                cancer_mean = np.mean(cancer_nonzero)
                cancer_std = np.std(cancer_nonzero, ddof=1)
                cancer_cv = (cancer_std / cancer_mean * 100) if cancer_mean > 0 else np.nan
            else:
                cancer_cv = np.nan

            # Normal CV
            normal_values = replace_empty_with_zero(row[normal_samples]).values.astype(float)
            normal_nonzero = normal_values[normal_values > 0]

            if len(normal_nonzero) >= 3:
                normal_mean = np.mean(normal_nonzero)
                normal_std = np.std(normal_nonzero, ddof=1)
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
        stats_text += f"Q1-Q3: {cv_df_clean['Cancer_CV'].quantile(0.25):.1f}%-{cv_df_clean['Cancer_CV'].quantile(0.75):.1f}%"

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
        stats_text += f"Q1-Q3: {cv_df_clean['Normal_CV'].quantile(0.25):.1f}%-{cv_df_clean['Normal_CV'].quantile(0.75):.1f}%"

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
