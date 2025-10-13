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

    def plot_pie_chart_glycan_types(self, df: pd.DataFrame, figsize: tuple = (16, 8)):
        """
        Create side-by-side pie charts showing glycan type distribution (Cancer vs Normal)

        Args:
            df: Annotated DataFrame with intensity data
            figsize: Figure size (width, height)
        """
        logger.info("Creating glycan type distribution pie charts...")

        # Get sample columns using centralized function
        cancer_samples, normal_samples = get_sample_columns(df)

        # Calculate total intensities per glycan type for each group
        # Using total intensity (sum across all samples) as the metric
        glycan_types = ['Non', 'Sialylated', 'Fucosylated', 'Both']

        cancer_data = {}
        normal_data = {}

        for glycan_type in glycan_types:
            mask = df['GlycanType'] == glycan_type

            # Sum intensities across all samples in group
            cancer_total = replace_empty_with_zero(df[mask][cancer_samples]).sum().sum()
            normal_total = replace_empty_with_zero(df[mask][normal_samples]).sum().sum()

            cancer_data[glycan_type] = cancer_total
            normal_data[glycan_type] = normal_total

        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

        # Define colors using standard palette
        colors = [LEGACY_GLYCAN_COLORS.get(gt, DEFAULT_FALLBACK_COLOR) for gt in glycan_types]

        # Cancer pie chart
        cancer_values = [cancer_data[gt] for gt in glycan_types]
        cancer_total_sum = sum(cancer_values)

        # Safe autopct lambda with zero-denominator guard
        wedges1, texts1, autotexts1 = ax1.pie(
            cancer_values,
            labels=glycan_types,
            colors=colors,
            autopct =lambda pct: f'{pct:.1f}%\n({pct * cancer_total_sum / 100:.2e})' if pct > 0 and cancer_total_sum > 0 else '',
            startangle=90,
            textprops={'fontsize': 11, 'weight': 'bold'}
        )

        ax1.set_title('Cancer Group\nGlycan Type Distribution', fontsize=AXIS_LABEL_SIZE, weight='bold', pad=20)

        # Normal pie chart
        normal_values = [normal_data[gt] for gt in glycan_types]
        normal_total_sum = sum(normal_values)

        # Safe autopct lambda with zero-denominator guard
        wedges2, texts2, autotexts2 = ax2.pie(
            normal_values,
            labels=glycan_types,
            colors=colors,
            autopct =lambda pct: f'{pct:.1f}%\n({pct * normal_total_sum / 100:.2e})' if pct > 0 and normal_total_sum > 0 else '',
            startangle=90,
            textprops={'fontsize': 11, 'weight': 'bold'}
        )

        ax2.set_title('Normal Group\nGlycan Type Distribution', fontsize=AXIS_LABEL_SIZE, weight='bold', pad=20)

        # Overall title
        fig.suptitle('Glycan Type Distribution by Group (Based on Total Intensity)',
                     fontsize=TITLE_SIZE, weight='bold', y=0.98)

        plt.tight_layout()

        # Save plot using standardized function
        output_file = self.output_dir / 'pie_chart_glycan_types.png'
        save_publication_figure(fig, output_file, dpi=DPI_MAIN)
        logger.info(f"Saved glycan type pie charts to {output_file} (optimized, {DPI_MAIN} DPI)")

        # Save trace data with zero-denominator guards
        trace_data = pd.DataFrame({
            'GlycanType': glycan_types,
            'Cancer_Intensity': cancer_values,
            'Cancer_Percentage': [v / cancer_total_sum * 100 if cancer_total_sum > 0 else 0 for v in cancer_values],
            'Normal_Intensity': normal_values,
            'Normal_Percentage': [v / normal_total_sum * 100 if normal_total_sum > 0 else 0 for v in normal_values]
        })
        save_trace_data(trace_data, self.output_dir, 'pie_chart_glycan_types_data.csv')

        plt.close()

    def plot_pie_chart_primary_classification(self, df: pd.DataFrame, figsize: tuple = (16, 8)):
        """
        Create side-by-side pie charts showing primary classification distribution (Cancer vs Normal)

        Args:
            df: Annotated DataFrame with intensity data
            figsize: Figure size (width, height)
        """
        logger.info("Creating primary classification distribution pie charts...")

        # Get sample columns using centralized function
        cancer_samples, normal_samples = get_sample_columns(df)

        # Primary classification categories
        primary_categories = ['Truncated', 'High Mannose', 'ComplexHybrid', 'Outlier']

        cancer_data = {}
        normal_data = {}

        for category in primary_categories:
            mask = df['PrimaryClassification'] == category

            # Sum intensities across all samples in group
            cancer_total = replace_empty_with_zero(df[mask][cancer_samples]).sum().sum()
            normal_total = replace_empty_with_zero(df[mask][normal_samples]).sum().sum()

            cancer_data[category] = cancer_total
            normal_data[category] = normal_total

        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

        # Use centralized colors from EXTENDED_CATEGORY_COLORS
        colors = [EXTENDED_CATEGORY_COLORS.get(cat, DEFAULT_FALLBACK_COLOR) for cat in primary_categories]

        # Cancer pie chart
        cancer_values = [cancer_data[cat] for cat in primary_categories]
        cancer_total_sum = sum(cancer_values)

        # Safe autopct lambda with zero-denominator guard
        wedges1, texts1, autotexts1 = ax1.pie(
            cancer_values,
            labels=primary_categories,
            colors=colors,
            autopct =lambda pct: f'{pct:.1f}%\n({pct * cancer_total_sum / 100:.2e})' if pct > 0 and cancer_total_sum > 0 else '',
            startangle=90,
            textprops={'fontsize': 11, 'weight': 'bold'}
        )

        ax1.set_title('Cancer Group\nPrimary Classification', fontsize=AXIS_LABEL_SIZE, weight='bold', pad=20)

        # Normal pie chart
        normal_values = [normal_data[cat] for cat in primary_categories]
        normal_total_sum = sum(normal_values)

        # Safe autopct lambda with zero-denominator guard
        wedges2, texts2, autotexts2 = ax2.pie(
            normal_values,
            labels=primary_categories,
            colors=colors,
            autopct =lambda pct: f'{pct:.1f}%\n({pct * normal_total_sum / 100:.2e})' if pct > 0 and normal_total_sum > 0 else '',
            startangle=90,
            textprops={'fontsize': 11, 'weight': 'bold'}
        )

        ax2.set_title('Normal Group\nPrimary Classification', fontsize=AXIS_LABEL_SIZE, weight='bold', pad=20)

        # Overall title
        fig.suptitle('Primary Classification Distribution by Group (Based on Total Intensity)',
                     fontsize=TITLE_SIZE, weight='bold', y=0.98)

        plt.tight_layout()

        # Save plot using standardized function
        output_file = self.output_dir / 'pie_chart_primary_classification.png'
        save_publication_figure(fig, output_file, dpi=DPI_MAIN)
        logger.info(f"Saved primary classification pie charts to {output_file} (optimized, {DPI_MAIN} DPI)")

        # Save trace data
        trace_data = pd.DataFrame({
            'PrimaryClassification': primary_categories,
            'Cancer_Intensity': cancer_values,
            'Cancer_Percentage': [v / cancer_total_sum * 100 if cancer_total_sum > 0 else 0 for v in cancer_values],
            'Normal_Intensity': normal_values,
            'Normal_Percentage': [v / normal_total_sum * 100 if normal_total_sum > 0 else 0 for v in normal_values]
        })
        save_trace_data(trace_data, self.output_dir, 'pie_chart_primary_classification_data.csv')

        plt.close()

    def plot_pie_chart_secondary_classification(self, df: pd.DataFrame, figsize: tuple = (18, 8)):
        """
        Create side-by-side pie charts showing secondary classification distribution (Cancer vs Normal)

        Args:
            df: Annotated DataFrame with intensity data
            figsize: Figure size (width, height)
        """
        logger.info("Creating secondary classification distribution pie charts...")

        # Get sample columns using centralized function
        cancer_samples, normal_samples = get_sample_columns(df)

        # Secondary classification categories (5 categories)
        secondary_categories = ['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        cancer_data = {}
        normal_data = {}

        for category in secondary_categories:
            mask = df['SecondaryClassification'] == category

            # Sum intensities across all samples in group
            cancer_total = replace_empty_with_zero(df[mask][cancer_samples]).sum().sum()
            normal_total = replace_empty_with_zero(df[mask][normal_samples]).sum().sum()

            cancer_data[category] = cancer_total
            normal_data[category] = normal_total

        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

        # Use centralized colors from EXTENDED_CATEGORY_COLORS
        colors = [EXTENDED_CATEGORY_COLORS.get(cat, DEFAULT_FALLBACK_COLOR) for cat in secondary_categories]

        # Cancer pie chart
        cancer_values = [cancer_data[cat] for cat in secondary_categories]
        cancer_total_sum = sum(cancer_values)

        # Safe autopct lambda with zero-denominator guard
        wedges1, texts1, autotexts1 = ax1.pie(
            cancer_values,
            labels=secondary_categories,
            colors=colors,
            autopct =lambda pct: f'{pct:.1f}%\n({pct * cancer_total_sum / 100:.2e})' if pct > 0 and cancer_total_sum > 0 else '',
            startangle=90,
            textprops={'fontsize': 10, 'weight': 'bold'}
        )

        ax1.set_title('Cancer Group\nSecondary Classification', fontsize=AXIS_LABEL_SIZE, weight='bold', pad=20)

        # Normal pie chart
        normal_values = [normal_data[cat] for cat in secondary_categories]
        normal_total_sum = sum(normal_values)

        # Safe autopct lambda with zero-denominator guard
        wedges2, texts2, autotexts2 = ax2.pie(
            normal_values,
            labels=secondary_categories,
            colors=colors,
            autopct =lambda pct: f'{pct:.1f}%\n({pct * normal_total_sum / 100:.2e})' if pct > 0 and normal_total_sum > 0 else '',
            startangle=90,
            textprops={'fontsize': 10, 'weight': 'bold'}
        )

        ax2.set_title('Normal Group\nSecondary Classification', fontsize=AXIS_LABEL_SIZE, weight='bold', pad=20)

        # Overall title
        fig.suptitle('Secondary Classification Distribution by Group (Based on Total Intensity)',
                     fontsize=TITLE_SIZE, weight='bold', y=0.98)

        plt.tight_layout()

        # Save plot using standardized function
        output_file = self.output_dir / 'pie_chart_secondary_classification.png'
        save_publication_figure(fig, output_file, dpi=DPI_MAIN)
        logger.info(f"Saved secondary classification pie charts to {output_file} (optimized, {DPI_MAIN} DPI)")

        # Save trace data
        trace_data = pd.DataFrame({
            'SecondaryClassification': secondary_categories,
            'Cancer_Intensity': cancer_values,
            'Cancer_Percentage': [v / cancer_total_sum * 100 if cancer_total_sum > 0 else 0 for v in cancer_values],
            'Normal_Intensity': normal_values,
            'Normal_Percentage': [v / normal_total_sum * 100 if normal_total_sum > 0 else 0 for v in normal_values]
        })
        save_trace_data(trace_data, self.output_dir, 'pie_chart_secondary_classification_data.csv')

        plt.close()
