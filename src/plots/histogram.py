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
    GROUP_PALETTE,  # ✨ Centralized group colors
    apply_standard_axis_style,
    apply_standard_legend,
    save_publication_figure, HISTOGRAM_DPI,
    apply_publication_theme,  # ✨ Enhanced styling
    DESIGN_SYSTEM_AVAILABLE,
    EDGE_LINEWIDTH_THIN, EDGE_LINEWIDTH_NORMAL,  # Centralized linewidth constants
    EDGE_COLOR_BLACK  # Edge color standardization
)

# Import premium design system if available
if DESIGN_SYSTEM_AVAILABLE:
    from .design_system import VisualEffects

logger = logging.getLogger(__name__)


class HistogramMixin:
    """Mixin class for histogram-related plots"""

    @staticmethod
    def _normalize_intensity_column(intensity_col: pd.Series, normalization: str) -> pd.Series:
        """
        Normalize a single intensity column using min-max scaling or keep raw.

        Args:
            intensity_col: Intensity column to normalize
            normalization: 'raw' for min-max normalization, else return as-is

        Returns:
            Normalized intensity column
        """
        if normalization == 'raw':
            min_val = intensity_col.min()
            max_val = intensity_col.max()
            if max_val > min_val:
                return (intensity_col - min_val) / (max_val - min_val)
            else:
                return intensity_col * 0  # All zeros if no variation
        else:
            return intensity_col

    @staticmethod
    def _apply_proportional_normalization(plot_df: pd.DataFrame, categories: list) -> pd.DataFrame:
        """
        Apply within-row proportional normalization (convert to proportions).

        Args:
            plot_df: DataFrame to normalize
            categories: List of category columns

        Returns:
            Normalized DataFrame
        """
        row_sums = plot_df[categories].sum(axis=1)
        for col in categories:
            plot_df[col] = plot_df[col] / row_sums
        return plot_df.fillna(0)

    def _aggregate_by_classification_samplewise(
        self,
        df: pd.DataFrame,
        sample_cols: list,
        categories: list,
        classification_col: str,
        normalization: str
    ) -> pd.DataFrame:
        """
        Aggregate intensities by classification across individual samples.

        Args:
            df: Annotated DataFrame
            sample_cols: List of sample column names
            categories: List of classification categories
            classification_col: Name of classification column
            normalization: 'raw' or 'aggregated'

        Returns:
            DataFrame with aggregated data (samples as rows, categories as columns)
        """
        data_for_plot = []

        for sample in sample_cols:
            sample_data = {'Sample': sample}

            # Get intensity column and apply normalization
            intensity_col = replace_empty_with_zero(df[sample])
            intensity_col = self._normalize_intensity_column(intensity_col, normalization)

            # Sum by classification
            for category in categories:
                mask = df[classification_col] == category
                sample_data[category] = intensity_col[mask].sum()

            data_for_plot.append(sample_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot).set_index('Sample')

        # Apply proportional normalization if requested
        if normalization == 'aggregated':
            plot_df = self._apply_proportional_normalization(plot_df, categories)

        return plot_df

    def _aggregate_by_classification_groupwise(
        self,
        df: pd.DataFrame,
        cancer_samples: list,
        normal_samples: list,
        categories: list,
        classification_col: str
    ) -> pd.DataFrame:
        """
        Aggregate intensities by classification for cancer vs normal groups.

        Args:
            df: Annotated DataFrame
            cancer_samples: List of cancer sample columns
            normal_samples: List of normal sample columns
            categories: List of classification categories
            classification_col: Name of classification column

        Returns:
            DataFrame with aggregated data (groups as rows, categories as columns)
        """
        data_for_plot = []

        for group_name, group_samples in [('Cancer', cancer_samples), ('Normal', normal_samples)]:
            group_data = {'Group': group_name}

            # Get intensity matrix for this group
            intensity_matrix = replace_empty_with_zero(df[group_samples])

            # Sum across all samples in the group, then by category
            for category in categories:
                mask = df[classification_col] == category
                group_data[category] = intensity_matrix[mask].sum().sum()

            data_for_plot.append(group_data)

        # Create DataFrame and apply proportional normalization
        plot_df = pd.DataFrame(data_for_plot).set_index('Group')
        plot_df = self._apply_proportional_normalization(plot_df, categories)

        return plot_df

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
            edgecolor=EDGE_COLOR_BLACK,
            linewidth=EDGE_LINEWIDTH_THIN
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

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'histogram_glycan_types_by_sample_normalized.png'
        save_publication_figure(fig, output_file, dpi=HISTOGRAM_DPI)
        logger.info(f"Saved normalized histogram to {output_file} (optimized, {HISTOGRAM_DPI} DPI)")

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

        # Primary classification categories
        primary_categories = ['Outlier', 'High Mannose', 'ComplexHybrid']

        # Aggregate data using unified helper
        plot_df = self._aggregate_by_classification_samplewise(
            df=df,
            sample_cols=sample_cols,
            categories=primary_categories,
            classification_col='PrimaryClassification',
            normalization=normalization
        )

        # Create grouped bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Use centralized colors from EXTENDED_CATEGORY_COLORS
        colors_primary = {
            'Outlier': EXTENDED_CATEGORY_COLORS['Outlier'],
            'High Mannose': EXTENDED_CATEGORY_COLORS['High Mannose'],
            'ComplexHybrid': EXTENDED_CATEGORY_COLORS['ComplexHybrid']
        }

        plot_df[primary_categories].plot(
            kind='bar',
            ax=ax,
            color=[colors_primary[cat] for cat in primary_categories],
            width=0.8,
            edgecolor=EDGE_COLOR_BLACK,
            linewidth=EDGE_LINEWIDTH_THIN
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
        save_publication_figure(fig, output_file, dpi=HISTOGRAM_DPI)
        logger.info(f"Saved primary classification histogram to {output_file} (optimized, {HISTOGRAM_DPI} DPI)")

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

        # Aggregate data using unified helper
        plot_df = self._aggregate_by_classification_samplewise(
            df=df,
            sample_cols=sample_cols,
            categories=secondary_categories,
            classification_col='SecondaryClassification',
            normalization=normalization
        )

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
            edgecolor=EDGE_COLOR_BLACK,
            linewidth=EDGE_LINEWIDTH_THIN
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
        save_publication_figure(fig, output_file, dpi=HISTOGRAM_DPI)
        logger.info(f"Saved secondary classification histogram to {output_file} (optimized, {HISTOGRAM_DPI} DPI)")

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

        # Primary classification categories
        primary_categories = ['Outlier', 'High Mannose', 'ComplexHybrid']

        # Aggregate data using unified helper
        plot_df = self._aggregate_by_classification_groupwise(
            df=df,
            cancer_samples=cancer_samples,
            normal_samples=normal_samples,
            categories=primary_categories,
            classification_col='PrimaryClassification'
        )

        # Create grouped bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Transpose for proper grouping (categories as x-axis)
        plot_df_t = plot_df.T

        # Use centralized group colors from GROUP_PALETTE
        plot_df_t.plot(
            kind='bar',
            ax=ax,
            color=[GROUP_PALETTE[g] for g in plot_df_t.columns],
            width=0.7,
            edgecolor=EDGE_COLOR_BLACK,
            linewidth=EDGE_LINEWIDTH_NORMAL
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
        save_publication_figure(fig, output_file, dpi=HISTOGRAM_DPI)
        logger.info(f"Saved primary Cancer vs Normal histogram to {output_file} (optimized, {HISTOGRAM_DPI} DPI)")

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

        # Secondary classification categories (5 categories)
        secondary_categories = ['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        # Aggregate data using unified helper
        plot_df = self._aggregate_by_classification_groupwise(
            df=df,
            cancer_samples=cancer_samples,
            normal_samples=normal_samples,
            categories=secondary_categories,
            classification_col='SecondaryClassification'
        )

        # Create grouped bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Transpose for proper grouping (categories as x-axis)
        plot_df_t = plot_df.T

        # Use centralized group colors from GROUP_PALETTE
        plot_df_t.plot(
            kind='bar',
            ax=ax,
            color=[GROUP_PALETTE[g] for g in plot_df_t.columns],
            width=0.7,
            edgecolor=EDGE_COLOR_BLACK,
            linewidth=EDGE_LINEWIDTH_NORMAL
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
        save_publication_figure(fig, output_file, dpi=HISTOGRAM_DPI)
        logger.info(f"Saved secondary Cancer vs Normal histogram to {output_file} (optimized, {HISTOGRAM_DPI} DPI)")

        # Save trace data
        save_trace_data(plot_df.reset_index(), self.output_dir, 'histogram_secondary_cancer_vs_normal_data.csv')

        plt.close()
