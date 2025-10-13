"""
Radar Chart Plot Module for pGlyco Auto Combine
Visualizes glycan type profiles as radar charts
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
from ..utils import replace_empty_with_zero, save_trace_data, get_sample_columns
from .plot_config import (
    save_publication_figure, DPI_MAIN, COLOR_CANCER, COLOR_NORMAL,
    TITLE_SIZE, LEGEND_SIZE, ANNOTATION_SIZE,
    PLOT_LINE_LINEWIDTH,  # Linewidth constant
    ALPHA_LIGHT, ALPHA_MEDIUM_LIGHT  # Alpha constants
)

logger = logging.getLogger(__name__)


class RadarChartPlotMixin:
    """Mixin class for radar/spider chart visualization"""

    def plot_radar_chart(self, df: pd.DataFrame, figsize: tuple = (12, 10)):
        """
        Create radar chart comparing glycan profiles between Cancer and Normal

        Uses mutually exclusive GlycanTypeCategory for accurate percentage distribution.

        Args:
            df: Annotated DataFrame with intensity data and GlycanTypeCategory
            figsize: Figure size (width, height)
        """
        # Get sample columns using centralized function
        cancer_samples, normal_samples = get_sample_columns(df)

        # Check if GlycanTypeCategory exists
        if 'GlycanTypeCategory' not in df.columns:
            logger.error("GlycanTypeCategory column not found in DataFrame")
            return None

        # Define mutually exclusive categories (5-category system)
        category_order = ['S', 'SF', 'F', 'HM', 'C/H']
        category_labels = {
            'S': 'Sialylated\n(only)',
            'SF': 'Sialofucosylated\n(both)',
            'F': 'Fucosylated\n(only)',
            'HM': 'High Mannose',
            'C/H': 'Complex/Hybrid\n(no modifications)'
        }

        # Calculate total intensities per category
        categories = []
        cancer_values = []
        normal_values = []

        # Calculate total intensity for each group
        cancer_total = replace_empty_with_zero(df[cancer_samples]).sum().sum()
        normal_total = replace_empty_with_zero(df[normal_samples]).sum().sum()

        # Calculate percentage for each mutually exclusive category
        for cat_code in category_order:
            cat_mask = df['GlycanTypeCategory'] == cat_code

            # Cancer group
            cancer_cat_total = replace_empty_with_zero(df[cat_mask][cancer_samples]).sum().sum()
            cancer_cat_pct = (cancer_cat_total / cancer_total * 100) if cancer_total > 0 else 0

            # Normal group
            normal_cat_total = replace_empty_with_zero(df[cat_mask][normal_samples]).sum().sum()
            normal_cat_pct = (normal_cat_total / normal_total * 100) if normal_total > 0 else 0

            categories.append(category_labels[cat_code])
            cancer_values.append(cancer_cat_pct)
            normal_values.append(normal_cat_pct)

            logger.info(f"{cat_code}: Cancer={cancer_cat_pct:.2f}%, Normal={normal_cat_pct:.2f}%")

        # Number of variables
        num_vars = len(categories)

        # Compute angle for each axis
        angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()

        # Complete the circle
        cancer_values += cancer_values[:1]
        normal_values += normal_values[:1]
        angles += angles[:1]

        # Create plot
        fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(projection='polar'))

        # Plot data
        ax.plot(angles, cancer_values, 'o-', linewidth=PLOT_LINE_LINEWIDTH, label='Cancer', color=COLOR_CANCER)
        ax.fill(angles, cancer_values, alpha=ALPHA_LIGHT, color=COLOR_CANCER)

        ax.plot(angles, normal_values, 'o-', linewidth=PLOT_LINE_LINEWIDTH, label='Normal', color=COLOR_NORMAL)
        ax.fill(angles, normal_values, alpha=ALPHA_LIGHT, color=COLOR_NORMAL)

        # Fix axis to go in the right order and start at 12 o'clock
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)

        # Set category labels
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories, fontsize=ANNOTATION_SIZE, fontweight='bold')

        # Set y-axis limits (0-100%)
        ax.set_ylim(0, 100)
        ax.set_yticks([20, 40, 60, 80, 100])
        ax.set_yticklabels(['20%', '40%', '60%', '80%', '100%'], fontsize=ANNOTATION_SIZE)

        # Add grid
        ax.grid(True, alpha=ALPHA_MEDIUM_LIGHT)

        # Add legend
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=LEGEND_SIZE)

        # Title
        plt.title('Glycan Type Distribution: Cancer vs Normal\n(Mutually Exclusive Categories)',
                  fontsize=TITLE_SIZE, fontweight='bold', pad=20)

        # Log totals for validation
        cancer_sum = sum(cancer_values[:-1])  # Exclude duplicated first element
        normal_sum = sum(normal_values[:-1])
        logger.info(f"Radar chart percentage totals - Cancer: {cancer_sum:.2f}%, Normal: {normal_sum:.2f}%")

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'radar_chart_glycan_profile.png'
        save_publication_figure(fig, output_file, dpi=DPI_MAIN)
        logger.info(f"Saved radar chart to {output_file} (optimized, {DPI_MAIN} DPI)")

        # Prepare trace data (before completing the circle)
        radar_data = pd.DataFrame({
            'Category': categories,
            'Cancer_Percentage': cancer_values[:-1],  # Remove duplicate first element
            'Normal_Percentage': normal_values[:-1]   # Remove duplicate first element
        })

        save_trace_data(radar_data, self.output_dir, 'radar_chart_glycan_profile_data.csv')

        plt.close()

        return radar_data
