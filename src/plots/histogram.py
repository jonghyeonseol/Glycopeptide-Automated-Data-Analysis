"""
Histogram Plot Module for pGlyco Auto Combine
Handles histogram visualizations
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

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
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid']

        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Get intensity matrix
        intensity_matrix = df[sample_cols].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

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

        # Define colors for each category
        colors = {
            'High mannose': '#2ECC71',  # Green
            'C/H': "#1500FF",           # Blue
            'Fucosylated': "#FF0000",   # Red
            'Sialylated': "#FF00FB",    # Pink
            'Both': "#FF9D00"           # Orange
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

        ax.set_xlabel('Sample', fontsize=12)
        ax.set_ylabel('Total Signal Intensity (TIC Normalized)', fontsize=12)
        ax.set_title('Glycan Type Distribution by Sample (TIC Normalized Data)', fontsize=14, fontweight='bold')
        ax.legend(title='Glycan Type', loc='upper left', bbox_to_anchor=(1, 1), frameon=True)

        # Rotate x-axis labels
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=90, ha='right')

        # Use scientific notation for y-axis
        ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

        plt.tight_layout()

        output_file = self.output_dir / 'histogram_glycan_types_by_sample_normalized.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved normalized histogram to {output_file}")

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
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']

        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Get intensity matrix
        intensity_matrix = df[sample_cols].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

        # Primary classification categories (exclude Outlier for visualization)
        primary_categories = ['Truncated', 'High Mannose', 'ComplexHybrid']

        data_for_plot = []

        for sample in sample_cols:
            sample_data = {'Sample': sample}

            if normalization == 'raw':
                # Normalize raw data using min-max scaling
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

                # Min-max normalization
                min_val = intensity_col.min()
                max_val = intensity_col.max()
                if max_val > min_val:
                    intensity_col = (intensity_col - min_val) / (max_val - min_val)
                else:
                    intensity_col = intensity_col * 0  # All zeros if no variation
            else:
                # Use raw intensities for aggregation
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

            # Sum by primary classification
            for category in primary_categories:
                mask = df['PrimaryClassification'] == category
                sample_data[category] = intensity_col[mask].sum()

            data_for_plot.append(sample_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot)
        plot_df = plot_df.set_index('Sample')

        # Apply min-max normalization after aggregation if needed
        if normalization == 'aggregated':
            for col in primary_categories:
                min_val = plot_df[col].min()
                max_val = plot_df[col].max()
                if max_val > min_val:
                    plot_df[col] = (plot_df[col] - min_val) / (max_val - min_val)

        # Create grouped bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Define colors
        colors_primary = {
            'Truncated': '#CCCCCC',
            'High Mannose': '#2ECC71',
            'ComplexHybrid': '#3498DB'
        }

        plot_df[primary_categories].plot(
            kind='bar',
            ax=ax,
            color=[colors_primary[cat] for cat in primary_categories],
            width=0.8,
            edgecolor='black',
            linewidth=0.5
        )

        norm_text = 'Raw Normalized then Summed' if normalization == 'raw' else 'Summed then Normalized'
        ax.set_xlabel('Sample', fontsize=12)
        ax.set_ylabel('Normalized Intensity Sum', fontsize=12)
        ax.set_title(f'Primary Classification Distribution by Sample ({norm_text})', fontsize=14)
        ax.legend(title='Primary Classification', loc='upper left', bbox_to_anchor=(1, 1), frameon=True)

        plt.setp(ax.xaxis.get_majorticklabels(), rotation=90, ha='right')
        plt.tight_layout()

        output_file = self.output_dir / f'histogram_primary_{normalization}_normalized.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved primary classification histogram to {output_file}")

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
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']

        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Secondary classification categories (5 categories, exclude Truncated and Outlier)
        secondary_categories = ['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        data_for_plot = []

        for sample in sample_cols:
            sample_data = {'Sample': sample}

            if normalization == 'raw':
                # Normalize raw data using min-max scaling
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

                # Min-max normalization
                min_val = intensity_col.min()
                max_val = intensity_col.max()
                if max_val > min_val:
                    intensity_col = (intensity_col - min_val) / (max_val - min_val)
                else:
                    intensity_col = intensity_col * 0
            else:
                # Use raw intensities for aggregation
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

            # Sum by secondary classification
            for category in secondary_categories:
                mask = df['SecondaryClassification'] == category
                sample_data[category] = intensity_col[mask].sum()

            data_for_plot.append(sample_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot)
        plot_df = plot_df.set_index('Sample')

        # Apply min-max normalization after aggregation if needed
        if normalization == 'aggregated':
            for col in secondary_categories:
                min_val = plot_df[col].min()
                max_val = plot_df[col].max()
                if max_val > min_val:
                    plot_df[col] = (plot_df[col] - min_val) / (max_val - min_val)

        # Create grouped bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Define colors
        colors_secondary = {
            'High Mannose': '#2ECC71',
            'Complex/Hybrid': '#3498DB',
            'Fucosylated': '#E74C3C',
            'Sialylated': '#9B59B6',
            'Sialofucosylated': '#F39C12'
        }

        plot_df[secondary_categories].plot(
            kind='bar',
            ax=ax,
            color=[colors_secondary[cat] for cat in secondary_categories],
            width=0.8,
            edgecolor='black',
            linewidth=0.5
        )

        norm_text = 'Raw Normalized then Summed' if normalization == 'raw' else 'Summed then Normalized'
        ax.set_xlabel('Sample', fontsize=12)
        ax.set_ylabel('Normalized Intensity Sum', fontsize=12)
        ax.set_title(f'Secondary Classification Distribution by Sample ({norm_text})', fontsize=14)
        ax.legend(title='Secondary Classification', loc='upper left', bbox_to_anchor=(1, 1), frameon=True)

        plt.setp(ax.xaxis.get_majorticklabels(), rotation=90, ha='right')
        plt.tight_layout()

        output_file = self.output_dir / f'histogram_secondary_{normalization}_normalized.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved secondary classification histogram to {output_file}")

        plt.close()

    def plot_histogram_cancer_vs_normal_primary(self, df: pd.DataFrame, figsize: tuple = (12, 8)):
        """
        Create histogram comparing Cancer vs Normal for primary classification
        (Aggregated with Log2 scaling)

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

        # Separate Cancer and Normal samples
        cancer_samples = [col for col in sample_cols if col.startswith('C')]
        normal_samples = [col for col in sample_cols if col.startswith('N')]

        # Primary classification categories (exclude Outlier)
        primary_categories = ['Truncated', 'High Mannose', 'ComplexHybrid']

        # Calculate sums for each group
        data_for_plot = []

        for group_name, group_samples in [('Cancer', cancer_samples), ('Normal', normal_samples)]:
            group_data = {'Group': group_name}

            # Get intensity matrix for this group
            intensity_matrix = df[group_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

            # Sum across all samples in the group, then by category
            for category in primary_categories:
                mask = df['PrimaryClassification'] == category
                group_data[category] = intensity_matrix[mask].sum().sum()

            data_for_plot.append(group_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot)
        plot_df = plot_df.set_index('Group')

        # Apply Log2 transformation (no min-max normalization)
        for col in primary_categories:
            plot_df[col] = np.log2(plot_df[col] + 1)

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

        ax.set_xlabel('Primary Classification', fontsize=12)
        ax.set_ylabel('Log2(Total Intensity + 1)', fontsize=12)
        ax.set_title('Primary Classification: Cancer vs Normal (Log2 Scaled)', fontsize=14)
        ax.legend(title='Group', loc='upper right', frameon=True)

        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
        plt.tight_layout()

        output_file = self.output_dir / 'histogram_primary_cancer_vs_normal.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved primary Cancer vs Normal histogram to {output_file}")

        plt.close()

    def plot_histogram_cancer_vs_normal_secondary(self, df: pd.DataFrame, figsize: tuple = (12, 8)):
        """
        Create histogram comparing Cancer vs Normal for secondary classification
        (Aggregated with Log2 scaling)

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
            intensity_matrix = df[group_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

            # Sum across all samples in the group, then by category
            for category in secondary_categories:
                mask = df['SecondaryClassification'] == category
                group_data[category] = intensity_matrix[mask].sum().sum()

            data_for_plot.append(group_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot)
        plot_df = plot_df.set_index('Group')

        # Apply Log2 transformation (no min-max normalization)
        for col in secondary_categories:
            plot_df[col] = np.log2(plot_df[col] + 1)

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

        ax.set_xlabel('Secondary Classification', fontsize=12)
        ax.set_ylabel('Log2(Total Intensity + 1)', fontsize=12)
        ax.set_title('Secondary Classification: Cancer vs Normal (Log2 Scaled)', fontsize=14)
        ax.legend(title='Group', loc='upper right', frameon=True)

        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
        plt.tight_layout()

        output_file = self.output_dir / 'histogram_secondary_cancer_vs_normal.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved secondary Cancer vs Normal histogram to {output_file}")

        plt.close()
