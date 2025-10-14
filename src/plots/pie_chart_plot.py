"""
Pie Chart Plot Module for pGlyco Auto Combine
Visualizes glycan type and classification distributions as pie charts

Dependencies:
    External:
        - pandas: Data manipulation
        - matplotlib: Plotting backend

    Internal:
        - src.utils: replace_empty_with_zero, save_trace_data, get_sample_columns
        - src.plots.plot_config: LEGACY_GLYCAN_COLORS, EXTENDED_CATEGORY_COLORS, save_publication_figure

NEW MODULE: Created for comprehensive glycan distribution visualization
"""

import pandas as pd
import matplotlib.pyplot as plt
import logging
from ..utils import replace_empty_with_zero, save_trace_data, get_sample_columns
from .plot_config import (
    LEGACY_GLYCAN_COLORS,
    EXTENDED_CATEGORY_COLORS,
    DEFAULT_FALLBACK_COLOR,
    save_publication_figure,
    DPI_MAIN,
    TITLE_SIZE, AXIS_LABEL_SIZE  # Font constants
)

logger = logging.getLogger(__name__)


class PieChartPlotMixin:
    """Mixin class for pie chart visualizations"""

    def _plot_dual_pie_chart(self, df: pd.DataFrame, classification_column: str,
                             categories: list, color_palette: dict,
                             output_filename: str, chart_title: str,
                             trace_column_name: str, figsize: tuple = (16, 8),
                             label_fontsize: int = 11):
        """
        Unified base method for creating side-by-side pie charts (Cancer vs Normal)

        This method consolidates the logic shared across all three pie chart methods,
        eliminating ~180 lines of code duplication.

        Args:
            df: Annotated DataFrame with intensity data
            classification_column: Column name for classification (e.g., 'GlycanType')
            categories: List of category names to plot
            color_palette: Dictionary mapping category names to colors
            output_filename: Output PNG filename (e.g., 'pie_chart_glycan_types.png')
            chart_title: Main title for the figure
            trace_column_name: Column name for trace data (e.g., 'GlycanType')
            figsize: Figure size (width, height)
            label_fontsize: Font size for pie chart labels

        Returns:
            None (saves figure and trace data)
        """
        # Get sample columns using centralized function
        cancer_samples, normal_samples = get_sample_columns(df)

        # Aggregate data by classification
        cancer_data = {}
        normal_data = {}

        for category in categories:
            mask = df[classification_column] == category

            # Sum intensities across all samples in group
            cancer_total = replace_empty_with_zero(df[mask][cancer_samples]).sum().sum()
            normal_total = replace_empty_with_zero(df[mask][normal_samples]).sum().sum()

            cancer_data[category] = cancer_total
            normal_data[category] = normal_total

        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

        # Get colors from palette
        colors = [color_palette.get(cat, DEFAULT_FALLBACK_COLOR) for cat in categories]

        # Cancer pie chart
        cancer_values = [cancer_data[cat] for cat in categories]
        cancer_total_sum = sum(cancer_values)

        # Safe autopct lambda with zero-denominator guard
        wedges1, texts1, autotexts1 = ax1.pie(
            cancer_values,
            labels=categories,
            colors=colors,
            autopct=lambda pct: f'{pct:.1f}%\n({pct * cancer_total_sum / 100:.2e})' if pct > 0 and cancer_total_sum > 0 else '',
            startangle=90,
            textprops={'fontsize': label_fontsize, 'weight': 'bold'}
        )

        # Extract classification type from column name for title
        classification_type = classification_column.replace('Classification', '').replace('Glycan', 'Glycan ')
        ax1.set_title(f'Cancer Group\n{classification_type}', fontsize=AXIS_LABEL_SIZE, weight='bold', pad=20)

        # Normal pie chart
        normal_values = [normal_data[cat] for cat in categories]
        normal_total_sum = sum(normal_values)

        # Safe autopct lambda with zero-denominator guard
        wedges2, texts2, autotexts2 = ax2.pie(
            normal_values,
            labels=categories,
            colors=colors,
            autopct=lambda pct: f'{pct:.1f}%\n({pct * normal_total_sum / 100:.2e})' if pct > 0 and normal_total_sum > 0 else '',
            startangle=90,
            textprops={'fontsize': label_fontsize, 'weight': 'bold'}
        )

        ax2.set_title(f'Normal Group\n{classification_type}', fontsize=AXIS_LABEL_SIZE, weight='bold', pad=20)

        # Overall title
        fig.suptitle(chart_title, fontsize=TITLE_SIZE, weight='bold', y=0.98)

        plt.tight_layout()

        # Save plot using standardized function
        output_file = self.output_dir / output_filename
        save_publication_figure(fig, output_file, dpi=DPI_MAIN)
        logger.info(f"Saved {classification_type.lower()} pie charts to {output_file} (optimized, {DPI_MAIN} DPI)")

        # Save trace data with zero-denominator guards
        trace_data = pd.DataFrame({
            trace_column_name: categories,
            'Cancer_Intensity': cancer_values,
            'Cancer_Percentage': [v / cancer_total_sum * 100 if cancer_total_sum > 0 else 0 for v in cancer_values],
            'Normal_Intensity': normal_values,
            'Normal_Percentage': [v / normal_total_sum * 100 if normal_total_sum > 0 else 0 for v in normal_values]
        })
        trace_filename = output_filename.replace('.png', '_data.csv')
        save_trace_data(trace_data, self.output_dir, trace_filename)

        plt.close()

    def plot_pie_chart_glycan_types(self, df: pd.DataFrame, figsize: tuple = (16, 8)):
        """
        Create side-by-side pie charts showing glycan type distribution (Cancer vs Normal)

        Args:
            df: Annotated DataFrame with intensity data
            figsize: Figure size (width, height)
        """
        logger.info("Creating glycan type distribution pie charts...")

        # Use unified helper method
        self._plot_dual_pie_chart(
            df=df,
            classification_column='GlycanType',
            categories=['Non', 'Sialylated', 'Fucosylated', 'Both'],
            color_palette=LEGACY_GLYCAN_COLORS,
            output_filename='pie_chart_glycan_types.png',
            chart_title='Glycan Type Distribution by Group (Based on Total Intensity)',
            trace_column_name='GlycanType',
            figsize=figsize,
            label_fontsize=11
        )

    def plot_pie_chart_primary_classification(self, df: pd.DataFrame, figsize: tuple = (16, 8)):
        """
        Create side-by-side pie charts showing primary classification distribution (Cancer vs Normal)

        Args:
            df: Annotated DataFrame with intensity data
            figsize: Figure size (width, height)
        """
        logger.info("Creating primary classification distribution pie charts...")

        # Use unified helper method
        self._plot_dual_pie_chart(
            df=df,
            classification_column='PrimaryClassification',
            categories=['Truncated', 'High Mannose', 'ComplexHybrid', 'Outlier'],
            color_palette=EXTENDED_CATEGORY_COLORS,
            output_filename='pie_chart_primary_classification.png',
            chart_title='Primary Classification Distribution by Group (Based on Total Intensity)',
            trace_column_name='PrimaryClassification',
            figsize=figsize,
            label_fontsize=11
        )

    def plot_pie_chart_secondary_classification(self, df: pd.DataFrame, figsize: tuple = (18, 8)):
        """
        Create side-by-side pie charts showing secondary classification distribution (Cancer vs Normal)

        Args:
            df: Annotated DataFrame with intensity data
            figsize: Figure size (width, height)
        """
        logger.info("Creating secondary classification distribution pie charts...")

        # Use unified helper method
        self._plot_dual_pie_chart(
            df=df,
            classification_column='SecondaryClassification',
            categories=['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated'],
            color_palette=EXTENDED_CATEGORY_COLORS,
            output_filename='pie_chart_secondary_classification.png',
            chart_title='Secondary Classification Distribution by Group (Based on Total Intensity)',
            trace_column_name='SecondaryClassification',
            figsize=figsize,
            label_fontsize=10  # Slightly smaller font for 5 categories
        )
