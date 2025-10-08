"""
Histogram Plot Module for pGlyco Auto Combine
Handles histogram visualizations
"""

import pandas as pd
import matplotlib.pyplot as plt
import logging
from ..utils import replace_empty_with_zero, save_trace_data, get_sample_columns
from .plot_config import (
    HISTOGRAM_X_ROTATION, HISTOGRAM_X_HA, EXTENDED_CATEGORY_COLORS,
    apply_standard_axis_style,
    apply_standard_legend
)

logger = logging.getLogger(__name__)


class HistogramMixin:
    """Mixin class for histogram-related plots"""

    def plot_histogram_normalized(self, df: pd.DataFrame, figsize: tuple = (20, 12)):
        """
        Create histogram showing glycan type intensities per sample (TIC normalized data)

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # Identify sample columns
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Get intensity matrix
        intensity_matrix = replace_empty_with_zero(df[sample_cols])
        # Apply TIC normalization (same as analyzer.py)
        sample_sums = intensity_matrix.sum(axis=0)
        median_sum = sample_sums.median()
        sample_sums_safe = sample_sums.replace(0, 1)
        intensity_normalized = intensity_matrix / sample_sums_safe * median_sum

        # Calculate total intensity per sample per category (normalized)
        data_for_plot = []

        for i, sample in enumerate(sample_cols):
            sample_data = {'Sample': sample}

            # Get normalized intensity column
            intensity_col = intensity_normalized.iloc[:, i]

            # High mannose (priority 1)
            high_mannose_mask = df['IsHighMannose']
            sample_data['High mannose'] = intensity_col[high_mannose_mask].sum()

            # C/H (priority 2)
            ch_mask = df['IsComplexHybrid']
            sample_data['C/H'] = intensity_col[ch_mask].sum()

            # Both (Sialylated AND Fucosylated) (priority 3)
            both_mask = df['IsSialylated'] & df['IsFucosylated']
            sample_data['Both'] = intensity_col[both_mask].sum()

            # Sialylated only (priority 4)
            sia_only_mask = df['IsSialylated'] & ~df['IsFucosylated'] & ~df['IsHighMannose'] & ~df['IsComplexHybrid']
            sample_data['Sialylated'] = intensity_col[sia_only_mask].sum()

            # Fucosylated only (priority 5)
            fuc_only_mask = df['IsFucosylated'] & ~df['IsSialylated'] & ~df['IsHighMannose'] & ~df['IsComplexHybrid']
            sample_data['Fucosylated'] = intensity_col[fuc_only_mask].sum()

            data_for_plot.append(sample_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot)
        plot_df = plot_df.set_index('Sample')

        # Reorder columns
        column_order = ['High mannose', 'C/H', 'Fucosylated', 'Sialylated', 'Both']
        plot_df = plot_df[column_order]

        # Create stacked bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Define colors for each category (using standardized palette)
        colors = {
            'High mannose': EXTENDED_CATEGORY_COLORS['HM'],
            'C/H': EXTENDED_CATEGORY_COLORS['C/H'],
            'Fucosylated': EXTENDED_CATEGORY_COLORS['Fucosylated'],
            'Sialylated': EXTENDED_CATEGORY_COLORS['Sialylated'],
            'Both': EXTENDED_CATEGORY_COLORS['Sialofucosylated']
        }

        plot_df.plot(
            kind='bar',
            stacked=True,
            ax=ax,
            color=[colors[col] for col in column_order],
            width=0.8,
            edgecolor='black',
            linewidth=0.5
        )

        # Apply standardized styling
        apply_standard_axis_style(
            ax,
            xlabel='Sample',
            ylabel='Total Signal Intensity (TIC Normalized)',
            title='Glycan Type Distribution by Sample (TIC Normalized Data)',
            grid=False  # Grid not needed for histograms
        )

        # Apply standardized legend (positioned outside plot area)
        apply_standard_legend(ax, title='Glycan Type')

        # Rotate x-axis labels using standardized settings
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=HISTOGRAM_X_ROTATION, ha=HISTOGRAM_X_HA)

        # Use scientific notation for y-axis
        ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'histogram_glycan_types_by_sample_normalized.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved normalized histogram to {output_file}")

        # Save trace data
        save_trace_data(plot_df.reset_index(), self.output_dir, 'histogram_glycan_types_by_sample_normalized_data.csv')

        plt.close()

    def plot_histogram_primary_classification(self, df: pd.DataFrame, normalization: str = 'raw',
                                              figsize: tuple = (20, 10)):
        """
        Create histogram for primary classification (Truncated, High Mannose, ComplexHybrid)
        across all samples

        Args:
            df: Annotated DataFrame
            normalization: 'raw' (normalize raw data then sum) or 'aggregated' (sum then normalize)
            figsize: Figure size
        """
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Get intensity matrix
        replace_empty_with_zero(df[sample_cols])
        # Primary classification categories
        primary_categories = ['Outlier', 'High Mannose', 'ComplexHybrid']

        data_for_plot = []

        for sample in sample_cols:
            sample_data = {'Sample': sample}

            if normalization == 'raw':
                # Normalize raw data using min-max scaling
                intensity_col = replace_empty_with_zero(df[sample])
                # Min-max normalization
                min_val = intensity_col.min()
                max_val = intensity_col.max()
                if max_val > min_val:
                    intensity_col = (intensity_col - min_val) / (max_val - min_val)
                else:
                    intensity_col = intensity_col * 0  # All zeros if no variation
            else:
                # Use raw intensities for aggregation
                intensity_col = replace_empty_with_zero(df[sample])
            # Sum by primary classification
            for category in primary_categories:
                mask = df['PrimaryClassification'] == category
                sample_data[category] = intensity_col[mask].sum()

            data_for_plot.append(sample_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot)
        plot_df = plot_df.set_index('Sample')

        # Apply normalization after aggregation if needed
        if normalization == 'aggregated':
            # Within-sample proportional normalization: show proportion of each type within each sample
            # This allows visual comparison of HM vs CH ratios within samples
            row_sums = plot_df[primary_categories].sum(axis=1)
            for col in primary_categories:
                plot_df[col] = plot_df[col] / row_sums
            # Replace any NaN with 0 (if row sum was 0)
            plot_df = plot_df.fillna(0)

        # Create grouped bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Define colors
        colors_primary = {
            'Outlier': '#95A5A6',  # Gray - Outlier glycans
            'High Mannose': '#2ECC71',  # Green
            'ComplexHybrid': '#3498DB'  # Blue
        }

        plot_df[primary_categories].plot(
            kind='bar',
            ax=ax,
            color=[colors_primary[cat] for cat in primary_categories],
            width=0.8,
            edgecolor='black',
            linewidth=0.5
        )

        if normalization == 'aggregated':
            norm_text = 'Proportional (Within-Sample)'
            ylabel_text = 'Proportion of Total Intensity'
        else:
            norm_text = 'Raw Normalized then Summed'
            ylabel_text = 'Normalized Intensity Sum'

        # Apply standardized styling
        apply_standard_axis_style(
            ax,
            xlabel='Sample',
            ylabel=ylabel_text,
            title=f'Primary Classification Distribution by Sample ({norm_text})',
            grid=False
        )

        # Apply standardized legend (positioned outside plot area)
        apply_standard_legend(ax, title='Primary Classification')

        plt.setp(ax.xaxis.get_majorticklabels(), rotation=HISTOGRAM_X_ROTATION, ha=HISTOGRAM_X_HA)
        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / f'histogram_primary_{normalization}_normalized.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved primary classification histogram to {output_file}")

        # Save trace data
        save_trace_data(plot_df.reset_index(), self.output_dir,
                        f'histogram_primary_{normalization}_normalized_data.csv')

        plt.close()

    def plot_histogram_secondary_classification(self, df: pd.DataFrame, normalization: str = 'raw',
                                                figsize: tuple = (20, 10)):
        """
        Create histogram for secondary classification across all samples
        (excludes Truncated and Outlier)

        Args:
            df: Annotated DataFrame
            normalization: 'raw' (normalize raw data then sum) or 'aggregated' (sum then normalize)
            figsize: Figure size
        """
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Secondary classification categories (5 categories, exclude Truncated and Outlier)
        secondary_categories = ['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        data_for_plot = []

        for sample in sample_cols:
            sample_data = {'Sample': sample}

            if normalization == 'raw':
                # Normalize raw data using min-max scaling
                intensity_col = replace_empty_with_zero(df[sample])
                # Min-max normalization
                min_val = intensity_col.min()
                max_val = intensity_col.max()
                if max_val > min_val:
                    intensity_col = (intensity_col - min_val) / (max_val - min_val)
                else:
                    intensity_col = intensity_col * 0
            else:
                # Use raw intensities for aggregation
                intensity_col = replace_empty_with_zero(df[sample])
            # Sum by secondary classification
            for category in secondary_categories:
                mask = df['SecondaryClassification'] == category
                sample_data[category] = intensity_col[mask].sum()

            data_for_plot.append(sample_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot)
        plot_df = plot_df.set_index('Sample')

        # Apply normalization after aggregation if needed
        if normalization == 'aggregated':
            # Within-sample proportional normalization: show proportion of each type within each sample
            # This allows visual comparison of ratios within samples
            row_sums = plot_df[secondary_categories].sum(axis=1)
            for col in secondary_categories:
                plot_df[col] = plot_df[col] / row_sums
            # Replace any NaN with 0 (if row sum was 0)
            plot_df = plot_df.fillna(0)

        # Create grouped bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Define colors (using standardized palette to avoid conflicts)
        colors_secondary = {
            'High Mannose': EXTENDED_CATEGORY_COLORS['HM'],
            'Complex/Hybrid': EXTENDED_CATEGORY_COLORS['C/H'],
            'Fucosylated': EXTENDED_CATEGORY_COLORS['Fucosylated'],
            'Sialylated': EXTENDED_CATEGORY_COLORS['Sialylated'],
            'Sialofucosylated': EXTENDED_CATEGORY_COLORS['Sialofucosylated']
        }

        plot_df[secondary_categories].plot(
            kind='bar',
            ax=ax,
            color=[colors_secondary[cat] for cat in secondary_categories],
            width=0.8,
            edgecolor='black',
            linewidth=0.5
        )

        if normalization == 'aggregated':
            norm_text = 'Proportional (Within-Sample)'
            ylabel_text = 'Proportion of Total Intensity'
        else:
            norm_text = 'Raw Normalized then Summed'
            ylabel_text = 'Normalized Intensity Sum'

        # Apply standardized styling
        apply_standard_axis_style(
            ax,
            xlabel='Sample',
            ylabel=ylabel_text,
            title=f'Secondary Classification Distribution by Sample ({norm_text})',
            grid=False
        )

        # Apply standardized legend (positioned outside plot area)
        apply_standard_legend(ax, title='Secondary Classification')

        plt.setp(ax.xaxis.get_majorticklabels(), rotation=HISTOGRAM_X_ROTATION, ha=HISTOGRAM_X_HA)
        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / f'histogram_secondary_{normalization}_normalized.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved secondary classification histogram to {output_file}")

        # Save trace data
        save_trace_data(plot_df.reset_index(), self.output_dir,
                        f'histogram_secondary_{normalization}_normalized_data.csv')

        plt.close()

    def plot_histogram_cancer_vs_normal_primary(self, df: pd.DataFrame, figsize: tuple = (12, 8)):
        """
        Create histogram comparing Cancer vs Normal for primary classification
        (Aggregated with Log2 scaling)

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Separate Cancer and Normal samples
        cancer_samples = [col for col in sample_cols if col.startswith('C')]
        normal_samples = [col for col in sample_cols if col.startswith('N')]

        # Primary classification categories
        primary_categories = ['Outlier', 'High Mannose', 'ComplexHybrid']

        # Calculate sums for each group
        data_for_plot = []

        for group_name, group_samples in [('Cancer', cancer_samples), ('Normal', normal_samples)]:
            group_data = {'Group': group_name}

            # Get intensity matrix for this group
            intensity_matrix = replace_empty_with_zero(df[group_samples])
            # Sum across all samples in the group, then by category
            for category in primary_categories:
                mask = df['PrimaryClassification'] == category
                group_data[category] = intensity_matrix[mask].sum().sum()

            data_for_plot.append(group_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot)
        plot_df = plot_df.set_index('Group')

        # Calculate proportions within each group (Cancer and Normal separately)
        # This shows what % of each group's total is each glycan type
        row_sums = plot_df[primary_categories].sum(axis=1)
        for col in primary_categories:
            plot_df[col] = plot_df[col] / row_sums
        # Replace any NaN with 0
        plot_df = plot_df.fillna(0)

        # Create grouped bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Transpose for proper grouping (categories as x-axis)
        plot_df_t = plot_df.T

        # Define colors
        group_colors = {'Cancer': '#E74C3C', 'Normal': '#3498DB'}

        plot_df_t.plot(
            kind='bar',
            ax=ax,
            color=[group_colors[g] for g in plot_df_t.columns],
            width=0.7,
            edgecolor='black',
            linewidth=1
        )

        # Apply standardized styling
        apply_standard_axis_style(
            ax,
            xlabel='Primary Classification',
            ylabel='Proportion of Total Intensity',
            title='Primary Classification: Cancer vs Normal (Proportional)',
            grid=True  # Grid helpful for comparison plots
        )

        # Apply standardized legend (positioned outside plot area)
        apply_standard_legend(ax, title='Group')

        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'histogram_primary_cancer_vs_normal.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved primary Cancer vs Normal histogram to {output_file}")

        # Save trace data
        save_trace_data(plot_df.reset_index(), self.output_dir, 'histogram_primary_cancer_vs_normal_data.csv')

        plt.close()

    def plot_histogram_cancer_vs_normal_secondary(self, df: pd.DataFrame, figsize: tuple = (12, 8)):
        """
        Create histogram comparing Cancer vs Normal for secondary classification
        (Aggregated with Log2 scaling)

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Separate Cancer and Normal samples
        cancer_samples = [col for col in sample_cols if col.startswith('C')]
        normal_samples = [col for col in sample_cols if col.startswith('N')]

        # Secondary classification categories (5 categories)
        secondary_categories = ['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        # Calculate sums for each group
        data_for_plot = []

        for group_name, group_samples in [('Cancer', cancer_samples), ('Normal', normal_samples)]:
            group_data = {'Group': group_name}

            # Get intensity matrix for this group
            intensity_matrix = replace_empty_with_zero(df[group_samples])
            # Sum across all samples in the group, then by category
            for category in secondary_categories:
                mask = df['SecondaryClassification'] == category
                group_data[category] = intensity_matrix[mask].sum().sum()

            data_for_plot.append(group_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot)
        plot_df = plot_df.set_index('Group')

        # Calculate proportions within each group (Cancer and Normal separately)
        # This shows what % of each group's total is each glycan type
        row_sums = plot_df[secondary_categories].sum(axis=1)
        for col in secondary_categories:
            plot_df[col] = plot_df[col] / row_sums
        # Replace any NaN with 0
        plot_df = plot_df.fillna(0)

        # Create grouped bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Transpose for proper grouping (categories as x-axis)
        plot_df_t = plot_df.T

        # Define colors
        group_colors = {'Cancer': '#E74C3C', 'Normal': '#3498DB'}

        plot_df_t.plot(
            kind='bar',
            ax=ax,
            color=[group_colors[g] for g in plot_df_t.columns],
            width=0.7,
            edgecolor='black',
            linewidth=1
        )

        # Apply standardized styling
        apply_standard_axis_style(
            ax,
            xlabel='Secondary Classification',
            ylabel='Proportion of Total Intensity',
            title='Secondary Classification: Cancer vs Normal (Proportional)',
            grid=True  # Grid helpful for comparison plots
        )

        # Apply standardized legend (positioned outside plot area)
        apply_standard_legend(ax, title='Group')

        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'histogram_secondary_cancer_vs_normal.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved secondary Cancer vs Normal histogram to {output_file}")

        # Save trace data
        save_trace_data(plot_df.reset_index(), self.output_dir, 'histogram_secondary_cancer_vs_normal_data.csv')

        plt.close()
