"""
Boxplot Module for pGlyco Auto Combine
Handles all boxplot visualizations
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from scipy import stats

logger = logging.getLogger(__name__)


class BoxplotMixin:
    """Mixin class for boxplot-related visualizations"""

    def plot_boxplot(self, boxplot_data: pd.DataFrame, figsize: tuple = (12, 6)):
        """
        Create boxplot comparing glycan types between groups with statistical significance

        Args:
            boxplot_data: Long-format DataFrame from analyzer
            figsize: Figure size
        """
        fig, ax = plt.subplots(figsize=figsize)

        # Define fixed order for glycan types
        glycan_order = ['Non', 'Sialylated', 'Fucosylated', 'Both']

        # Filter to only include existing glycan types
        existing_types = [gt for gt in glycan_order if gt in boxplot_data['GlycanType'].unique()]

        # Create boxplot with ordered categories
        sns.boxplot(
            data=boxplot_data,
            x='GlycanType',
            y='Intensity',
            hue='Group',
            order=existing_types,
            hue_order=['Normal', 'Cancer'],
            palette={'Cancer': '#E74C3C', 'Normal': '#3498DB'},
            width=0.6,
            ax=ax
        )

        # Perform statistical tests and add significance markers
        y_max = boxplot_data['Intensity'].max()
        y_range = boxplot_data['Intensity'].max() - boxplot_data['Intensity'].min()

        for i, glycan_type in enumerate(existing_types):
            # Get data for each group
            cancer_data = boxplot_data[
                (boxplot_data['GlycanType'] == glycan_type) &
                (boxplot_data['Group'] == 'Cancer')
            ]['Intensity'].values

            normal_data = boxplot_data[
                (boxplot_data['GlycanType'] == glycan_type) &
                (boxplot_data['Group'] == 'Normal')
            ]['Intensity'].values

            # Skip if either group has insufficient data
            if len(cancer_data) < 3 or len(normal_data) < 3:
                continue

            # Perform Mann-Whitney U test (non-parametric)
            try:
                statistic, p_value = stats.mannwhitneyu(cancer_data, normal_data, alternative='two-sided')

                # Determine significance level
                if p_value < 0.001:
                    sig_marker = '***'
                elif p_value < 0.01:
                    sig_marker = '**'
                elif p_value < 0.05:
                    sig_marker = '*'
                else:
                    sig_marker = 'ns'

                # Add significance marker if significant
                if sig_marker != 'ns':
                    # Position marker inside the plot area, near the top
                    y_position = y_max - y_range * 0.02 - (y_range * 0.03 * (i % 2))
                    ax.text(
                        i,
                        y_position,
                        sig_marker,
                        ha='center',
                        va='top',
                        fontsize=14,
                        fontweight='bold'
                    )

                    logger.info(f"{glycan_type}: p={p_value:.4f} ({sig_marker})")

            except Exception as e:
                logger.warning(f"Statistical test failed for {glycan_type}: {str(e)}")

        ax.set_xlabel('Glycan Type')
        ax.set_ylabel('Log2(Intensity + 1)')
        ax.set_title('Glycan Intensity Distribution by Type and Group')
        ax.legend(title='Group', loc='best')

        plt.tight_layout()

        output_file = self.output_dir / 'boxplot_glycan_types.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved boxplot to {output_file}")

        plt.close()

    def plot_boxplot_extended(self, boxplot_data: pd.DataFrame, figsize: tuple = (14, 6)):
        """
        Create extended boxplot comparing 5 glycan categories between groups

        Args:
            boxplot_data: Long-format DataFrame from analyzer (extended)
            figsize: Figure size
        """
        fig, ax = plt.subplots(figsize=figsize)

        # Define fixed order for extended categories
        category_order = ['HM', 'C/H', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        # Filter to only include existing categories
        existing_categories = [cat for cat in category_order if cat in boxplot_data['ExtendedCategory'].unique()]

        # Create boxplot with ordered categories
        sns.boxplot(
            data=boxplot_data,
            x='ExtendedCategory',
            y='Intensity',
            hue='Group',
            order=existing_categories,
            hue_order=['Normal', 'Cancer'],
            palette={'Cancer': '#E74C3C', 'Normal': '#3498DB'},
            width=0.6,
            ax=ax
        )

        # Perform statistical tests and add significance markers
        y_max = boxplot_data['Intensity'].max()
        y_range = boxplot_data['Intensity'].max() - boxplot_data['Intensity'].min()

        for i, category in enumerate(existing_categories):
            # Get data for each group
            cancer_data = boxplot_data[
                (boxplot_data['ExtendedCategory'] == category) &
                (boxplot_data['Group'] == 'Cancer')
            ]['Intensity'].values

            normal_data = boxplot_data[
                (boxplot_data['ExtendedCategory'] == category) &
                (boxplot_data['Group'] == 'Normal')
            ]['Intensity'].values

            # Skip if either group has insufficient data
            if len(cancer_data) < 3 or len(normal_data) < 3:
                continue

            # Perform Mann-Whitney U test (non-parametric)
            try:
                statistic, p_value = stats.mannwhitneyu(cancer_data, normal_data, alternative='two-sided')

                # Determine significance level
                if p_value < 0.001:
                    sig_marker = '***'
                elif p_value < 0.01:
                    sig_marker = '**'
                elif p_value < 0.05:
                    sig_marker = '*'
                else:
                    sig_marker = 'ns'

                # Add significance marker if significant
                if sig_marker != 'ns':
                    # Position marker inside the plot area, near the top
                    y_position = y_max - y_range * 0.02 - (y_range * 0.03 * (i % 2))
                    ax.text(
                        i,
                        y_position,
                        sig_marker,
                        ha='center',
                        va='top',
                        fontsize=14,
                        fontweight='bold'
                    )

                    logger.info(f"{category}: p={p_value:.4f} ({sig_marker})")

            except Exception as e:
                logger.warning(f"Statistical test failed for {category}: {str(e)}")

        ax.set_xlabel('Glycan Category', fontsize=12)
        ax.set_ylabel('Log2(Intensity + 1)', fontsize=12)
        ax.set_title('Extended Glycan Category Distribution (HM, C/H, Fucosylated, Sialylated, Sialofucosylated)', fontsize=14)
        ax.legend(title='Group', loc='upper left', bbox_to_anchor=(1, 1), frameon=True)

        plt.tight_layout()

        output_file = self.output_dir / 'boxplot_extended_categories.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved extended boxplot to {output_file}")

        plt.close()

    def plot_boxplot_primary_classification(self, df: pd.DataFrame, normalization: str = 'raw', figsize: tuple = (12, 8)):
        """
        Create boxplot for primary classification

        Args:
            df: Annotated DataFrame
            normalization: 'raw' or 'aggregated'
            figsize: Figure size
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']
        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Primary classification categories
        primary_categories = ['Truncated', 'High Mannose', 'ComplexHybrid']

        # Prepare long format data
        data_for_plot = []

        for sample in sample_cols:
            if normalization == 'raw':
                # Min-max normalize per sample
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)
                min_val, max_val = intensity_col.min(), intensity_col.max()
                if max_val > min_val:
                    intensity_col = (intensity_col - min_val) / (max_val - min_val)
            else:
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

            group = 'Cancer' if sample.startswith('C') else 'Normal'

            for idx, row in df.iterrows():
                if row['PrimaryClassification'] in primary_categories:
                    data_for_plot.append({
                        'Sample': sample,
                        'Group': group,
                        'Classification': row['PrimaryClassification'],
                        'Intensity': intensity_col.iloc[idx] if normalization == 'raw' else intensity_col.iloc[idx]
                    })

        plot_df = pd.DataFrame(data_for_plot)

        # Apply aggregated normalization if needed
        if normalization == 'aggregated':
            plot_df['Intensity'] = np.log2(plot_df['Intensity'] + 1)

        fig, ax = plt.subplots(figsize=figsize)

        sns.boxplot(data=plot_df, x='Classification', y='Intensity', hue='Group',
                   palette={'Cancer': '#E74C3C', 'Normal': '#3498DB'}, ax=ax)

        norm_text = 'Raw Normalized' if normalization == 'raw' else 'Log2 Scaled'
        ax.set_xlabel('Primary Classification', fontsize=12)
        ax.set_ylabel(f'Intensity ({norm_text})', fontsize=12)
        ax.set_title(f'Primary Classification Distribution ({norm_text})', fontsize=14)
        ax.legend(title='Group', loc='upper right')

        plt.tight_layout()

        output_file = self.output_dir / f'boxplot_primary_{normalization}_normalized.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved primary boxplot to {output_file}")

        plt.close()

    def plot_boxplot_secondary_classification(self, df: pd.DataFrame, normalization: str = 'raw', figsize: tuple = (14, 8)):
        """
        Create boxplot for secondary classification

        Args:
            df: Annotated DataFrame
            normalization: 'raw' or 'aggregated'
            figsize: Figure size
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']
        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Secondary classification categories
        secondary_categories = ['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        # Prepare long format data
        data_for_plot = []

        for sample in sample_cols:
            if normalization == 'raw':
                # Min-max normalize per sample
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)
                min_val, max_val = intensity_col.min(), intensity_col.max()
                if max_val > min_val:
                    intensity_col = (intensity_col - min_val) / (max_val - min_val)
            else:
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

            group = 'Cancer' if sample.startswith('C') else 'Normal'

            for idx, row in df.iterrows():
                if row['SecondaryClassification'] in secondary_categories:
                    data_for_plot.append({
                        'Sample': sample,
                        'Group': group,
                        'Classification': row['SecondaryClassification'],
                        'Intensity': intensity_col.iloc[idx] if normalization == 'raw' else intensity_col.iloc[idx]
                    })

        plot_df = pd.DataFrame(data_for_plot)

        # Apply aggregated normalization if needed
        if normalization == 'aggregated':
            plot_df['Intensity'] = np.log2(plot_df['Intensity'] + 1)

        fig, ax = plt.subplots(figsize=figsize)

        sns.boxplot(data=plot_df, x='Classification', y='Intensity', hue='Group',
                   order=secondary_categories,
                   palette={'Cancer': '#E74C3C', 'Normal': '#3498DB'}, ax=ax)

        norm_text = 'Raw Normalized' if normalization == 'raw' else 'Log2 Scaled'
        ax.set_xlabel('Secondary Classification', fontsize=12)
        ax.set_ylabel(f'Intensity ({norm_text})', fontsize=12)
        ax.set_title(f'Secondary Classification Distribution ({norm_text})', fontsize=14)
        ax.legend(title='Group', loc='upper right')
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

        plt.tight_layout()

        output_file = self.output_dir / f'boxplot_secondary_{normalization}_normalized.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved secondary boxplot to {output_file}")

        plt.close()

    def plot_boxplot_cancer_vs_normal_primary(self, df: pd.DataFrame, figsize: tuple = (10, 6)):
        """
        Create boxplot for Cancer vs Normal comparison (Primary classification)

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']
        sample_cols = [col for col in df.columns if col not in metadata_cols]

        primary_categories = ['Truncated', 'High Mannose', 'ComplexHybrid']

        # Prepare data
        data_for_plot = []
        for sample in sample_cols:
            intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)
            group = 'Cancer' if sample.startswith('C') else 'Normal'

            for idx, row in df.iterrows():
                if row['PrimaryClassification'] in primary_categories:
                    data_for_plot.append({
                        'Group': group,
                        'Classification': row['PrimaryClassification'],
                        'Intensity': np.log2(intensity_col.iloc[idx] + 1)
                    })

        plot_df = pd.DataFrame(data_for_plot)

        fig, ax = plt.subplots(figsize=figsize)

        sns.boxplot(data=plot_df, x='Classification', y='Intensity', hue='Group',
                   palette={'Cancer': '#E74C3C', 'Normal': '#3498DB'}, ax=ax)

        ax.set_xlabel('Primary Classification', fontsize=12)
        ax.set_ylabel('Log2(Intensity + 1)', fontsize=12)
        ax.set_title('Primary Classification: Cancer vs Normal', fontsize=14)
        ax.legend(title='Group')

        plt.tight_layout()

        output_file = self.output_dir / 'boxplot_primary_cancer_vs_normal.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved primary Cancer vs Normal boxplot to {output_file}")

        plt.close()

    def plot_boxplot_cancer_vs_normal_secondary(self, df: pd.DataFrame, figsize: tuple = (12, 6)):
        """
        Create boxplot for Cancer vs Normal comparison (Secondary classification)

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']
        sample_cols = [col for col in df.columns if col not in metadata_cols]

        secondary_categories = ['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        # Prepare data
        data_for_plot = []
        for sample in sample_cols:
            intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)
            group = 'Cancer' if sample.startswith('C') else 'Normal'

            for idx, row in df.iterrows():
                if row['SecondaryClassification'] in secondary_categories:
                    data_for_plot.append({
                        'Group': group,
                        'Classification': row['SecondaryClassification'],
                        'Intensity': np.log2(intensity_col.iloc[idx] + 1)
                    })

        plot_df = pd.DataFrame(data_for_plot)

        fig, ax = plt.subplots(figsize=figsize)

        sns.boxplot(data=plot_df, x='Classification', y='Intensity', hue='Group',
                   order=secondary_categories,
                   palette={'Cancer': '#E74C3C', 'Normal': '#3498DB'}, ax=ax)

        ax.set_xlabel('Secondary Classification', fontsize=12)
        ax.set_ylabel('Log2(Intensity + 1)', fontsize=12)
        ax.set_title('Secondary Classification: Cancer vs Normal', fontsize=14)
        ax.legend(title='Group')
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

        plt.tight_layout()

        output_file = self.output_dir / 'boxplot_secondary_cancer_vs_normal.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved secondary Cancer vs Normal boxplot to {output_file}")

        plt.close()
