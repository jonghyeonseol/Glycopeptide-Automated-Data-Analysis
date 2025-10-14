"""
CV Distribution Plot Module for pGlyco Auto Combine
Visualizes coefficient of variation for quality control

Dependencies:
    External:
        - pandas: Data manipulation
        - numpy: Numerical computations
        - matplotlib: Plotting backend

    Internal:
        - src.utils: save_trace_data, get_sample_columns
        - src.data_preparation: DataPreparationConfig, calculate_group_statistics_standardized
        - src.plots.plot_config: save_publication_figure

UPDATED: Now uses centralized data preparation for consistency
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
from ..utils import save_trace_data, get_sample_columns
from ..data_preparation import (
    DataPreparationConfig,
    calculate_group_statistics_standardized
)
from .plot_config import (
    save_publication_figure, COLOR_CANCER, COLOR_NORMAL, DPI_MAIN,
    TITLE_SIZE, AXIS_LABEL_SIZE, LEGEND_SIZE, ANNOTATION_SIZE,  # Font constants
    EDGE_LINEWIDTH_THIN, PLOT_LINE_LINEWIDTH,  # Linewidth constants
    ALPHA_HIGH, ALPHA_VERY_HIGH, ALPHA_MEDIUM_LIGHT,  # Alpha constants,
    EDGE_COLOR_BLACK,  # Edge color standardization
    # Linestyle constants (Phase 10.3.8)
    LINESTYLE_DASHED
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
        logger.info("Creating CV distribution plot...")

        # Get sample columns using centralized function
        cancer_samples, normal_samples = get_sample_columns(df)
        logger.debug(f"  Analyzing {len(cancer_samples)} Cancer and {len(normal_samples)} Normal samples")

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
        cancer_color = COLOR_CANCER
        normal_color = COLOR_NORMAL

        # Plot Cancer CV distribution
        ax1 = axes[0]
        ax1.hist(cv_df_clean['Cancer_CV'], bins=50, color=cancer_color,
                 alpha=ALPHA_HIGH, edgecolor=EDGE_COLOR_BLACK, linewidth=EDGE_LINEWIDTH_THIN)

        # Add median line
        cancer_median = cv_df_clean['Cancer_CV'].median()
        ax1.axvline(cancer_median, color='darkred', linestyle=LINESTYLE_DASHED,
                    linewidth=PLOT_LINE_LINEWIDTH, label=f'Median: {cancer_median:.1f}%')

        ax1.set_xlabel('CV (%)', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax1.set_ylabel('Count', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax1.set_title('CV Distribution: Cancer Samples', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax1.legend(loc='upper right', frameon=True, fontsize=LEGEND_SIZE)
        ax1.grid(True, alpha=ALPHA_MEDIUM_LIGHT)

        # Add statistics text
        stats_text = f"n = {len(cv_df_clean)}\n"
        stats_text += f"Mean: {cv_df_clean['Cancer_CV'].mean():.1f}%\n"
        stats_text += f"Median: {cancer_median:.1f}%\n"
        q1 = cv_df_clean['Cancer_CV'].quantile(0.25)
        q3 = cv_df_clean['Cancer_CV'].quantile(0.75)
        stats_text += f"Q1-Q3: {q1:.1f}%-{q3:.1f}%"

        ax1.text(0.98, 0.98, stats_text, transform=ax1.transAxes,
                 fontsize=ANNOTATION_SIZE, verticalalignment='top', horizontalalignment='right',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=ALPHA_VERY_HIGH))

        # Plot Normal CV distribution
        ax2 = axes[1]
        ax2.hist(cv_df_clean['Normal_CV'], bins=50, color=normal_color,
                 alpha=ALPHA_HIGH, edgecolor=EDGE_COLOR_BLACK, linewidth=EDGE_LINEWIDTH_THIN)

        # Add median line
        normal_median = cv_df_clean['Normal_CV'].median()
        ax2.axvline(normal_median, color='darkblue', linestyle=LINESTYLE_DASHED,
                    linewidth=PLOT_LINE_LINEWIDTH, label=f'Median: {normal_median:.1f}%')

        ax2.set_xlabel('CV (%)', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax2.set_ylabel('Count', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax2.set_title('CV Distribution: Normal Samples', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax2.legend(loc='upper right', frameon=True, fontsize=LEGEND_SIZE)
        ax2.grid(True, alpha=ALPHA_MEDIUM_LIGHT)

        # Add statistics text
        stats_text = f"n = {len(cv_df_clean)}\n"
        stats_text += f"Mean: {cv_df_clean['Normal_CV'].mean():.1f}%\n"
        stats_text += f"Median: {normal_median:.1f}%\n"
        q1 = cv_df_clean['Normal_CV'].quantile(0.25)
        q3 = cv_df_clean['Normal_CV'].quantile(0.75)
        stats_text += f"Q1-Q3: {q1:.1f}%-{q3:.1f}%"

        ax2.text(0.98, 0.98, stats_text, transform=ax2.transAxes,
                 fontsize=ANNOTATION_SIZE, verticalalignment='top', horizontalalignment='right',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=ALPHA_VERY_HIGH))

        # Main title
        fig.suptitle('Coefficient of Variation (CV) Distribution Analysis',
                     fontsize=TITLE_SIZE, fontweight='bold', y=1.02)

        plt.tight_layout()

        # Save plot using standardized function
        output_file = self.output_dir / 'cv_distribution.png'
        save_publication_figure(fig, output_file, dpi=DPI_MAIN)
        logger.info(f"Saved CV distribution plot to {output_file} (optimized, {DPI_MAIN} DPI)")

        # Save trace data
        save_trace_data(cv_df_clean, self.output_dir, 'cv_distribution_data.csv')

        plt.close()

        logger.info(f"✓ CV distribution analysis complete ({len(cv_df_clean)} glycopeptides)")
        return cv_df_clean
