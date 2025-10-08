"""
Pie Chart Plot Module for pGlyco Auto Combine
Visualizes glycan type and classification distributions as pie charts

NEW MODULE: Created for comprehensive glycan distribution visualization
"""

import pandas as pd
import matplotlib.pyplot as plt
import logging
from ..utils import replace_empty_with_zero, save_trace_data
from .plot_config import (
    LEGACY_GLYCAN_COLORS,
    EXTENDED_CATEGORY_COLORS
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

        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]

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
        colors = [LEGACY_GLYCAN_COLORS.get(gt, '#CCCCCC') for gt in glycan_types]

        # Cancer pie chart
        cancer_values = [cancer_data[gt] for gt in glycan_types]
        cancer_total_sum = sum(cancer_values)

        wedges1, texts1, autotexts1 = ax1.pie(
            cancer_values,
            labels=glycan_types,
            colors=colors,
            autopct =lambda pct: f'{pct:.1f}%\n({pct * cancer_total_sum / 100:.2e})' if pct > 0 else '',
            startangle=90,
            textprops={'fontsize': 11, 'weight': 'bold'}
        )

        ax1.set_title('Cancer Group\nGlycan Type Distribution', fontsize=14, weight='bold', pad=20)

        # Normal pie chart
        normal_values = [normal_data[gt] for gt in glycan_types]
        normal_total_sum = sum(normal_values)

        wedges2, texts2, autotexts2 = ax2.pie(
            normal_values,
            labels=glycan_types,
            colors=colors,
            autopct =lambda pct: f'{pct:.1f}%\n({pct * normal_total_sum / 100:.2e})' if pct > 0 else '',
            startangle=90,
            textprops={'fontsize': 11, 'weight': 'bold'}
        )

        ax2.set_title('Normal Group\nGlycan Type Distribution', fontsize=14, weight='bold', pad=20)

        # Overall title
        fig.suptitle('Glycan Type Distribution by Group (Based on Total Intensity)',
                     fontsize=16, weight='bold', y=0.98)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'pie_chart_glycan_types.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved glycan type pie charts to {output_file}")

        # Save trace data
        trace_data = pd.DataFrame({
            'GlycanType': glycan_types,
            'Cancer_Intensity': cancer_values,
            'Cancer_Percentage': [v / cancer_total_sum * 100 for v in cancer_values],
            'Normal_Intensity': normal_values,
            'Normal_Percentage': [v / normal_total_sum * 100 for v in normal_values]
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

        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]

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

        # Define colors
        colors_primary = {
            'Truncated': '#CCCCCC',       # Gray
            'High Mannose': '#2ECC71',    # Green
            'ComplexHybrid': '#3498DB',   # Blue
            'Outlier': '#95A5A6'          # Light gray
        }
        colors = [colors_primary.get(cat, '#CCCCCC') for cat in primary_categories]

        # Cancer pie chart
        cancer_values = [cancer_data[cat] for cat in primary_categories]
        cancer_total_sum = sum(cancer_values)

        wedges1, texts1, autotexts1 = ax1.pie(
            cancer_values,
            labels=primary_categories,
            colors=colors,
            autopct =lambda pct: f'{pct:.1f}%\n({pct * cancer_total_sum / 100:.2e})' if pct > 0 else '',
            startangle=90,
            textprops={'fontsize': 11, 'weight': 'bold'}
        )

        ax1.set_title('Cancer Group\nPrimary Classification', fontsize=14, weight='bold', pad=20)

        # Normal pie chart
        normal_values = [normal_data[cat] for cat in primary_categories]
        normal_total_sum = sum(normal_values)

        wedges2, texts2, autotexts2 = ax2.pie(
            normal_values,
            labels=primary_categories,
            colors=colors,
            autopct =lambda pct: f'{pct:.1f}%\n({pct * normal_total_sum / 100:.2e})' if pct > 0 else '',
            startangle=90,
            textprops={'fontsize': 11, 'weight': 'bold'}
        )

        ax2.set_title('Normal Group\nPrimary Classification', fontsize=14, weight='bold', pad=20)

        # Overall title
        fig.suptitle('Primary Classification Distribution by Group (Based on Total Intensity)',
                     fontsize=16, weight='bold', y=0.98)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'pie_chart_primary_classification.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved primary classification pie charts to {output_file}")

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

        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]

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

        # Define colors (using standardized palette)
        colors = [EXTENDED_CATEGORY_COLORS.get(cat, '#CCCCCC') for cat in secondary_categories]

        # Cancer pie chart
        cancer_values = [cancer_data[cat] for cat in secondary_categories]
        cancer_total_sum = sum(cancer_values)

        wedges1, texts1, autotexts1 = ax1.pie(
            cancer_values,
            labels=secondary_categories,
            colors=colors,
            autopct =lambda pct: f'{pct:.1f}%\n({pct * cancer_total_sum / 100:.2e})' if pct > 0 else '',
            startangle=90,
            textprops={'fontsize': 10, 'weight': 'bold'}
        )

        ax1.set_title('Cancer Group\nSecondary Classification', fontsize=14, weight='bold', pad=20)

        # Normal pie chart
        normal_values = [normal_data[cat] for cat in secondary_categories]
        normal_total_sum = sum(normal_values)

        wedges2, texts2, autotexts2 = ax2.pie(
            normal_values,
            labels=secondary_categories,
            colors=colors,
            autopct =lambda pct: f'{pct:.1f}%\n({pct * normal_total_sum / 100:.2e})' if pct > 0 else '',
            startangle=90,
            textprops={'fontsize': 10, 'weight': 'bold'}
        )

        ax2.set_title('Normal Group\nSecondary Classification', fontsize=14, weight='bold', pad=20)

        # Overall title
        fig.suptitle('Secondary Classification Distribution by Group (Based on Total Intensity)',
                     fontsize=16, weight='bold', y=0.98)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'pie_chart_secondary_classification.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved secondary classification pie charts to {output_file}")

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
