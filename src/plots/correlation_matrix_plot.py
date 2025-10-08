"""
Correlation Matrix Plot Module for pGlyco Auto Combine
Visualizes sample-to-sample correlations
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from ..utils import replace_empty_with_zero, save_trace_data

logger = logging.getLogger(__name__)


class CorrelationMatrixPlotMixin:
    """Mixin class for correlation matrix visualization"""

    def plot_correlation_matrix(self, df: pd.DataFrame, figsize: tuple = (20, 18)):
        """
        Create correlation matrix heatmap for samples

        Args:
            df: Annotated DataFrame with intensity data
            figsize: Figure size (width, height)
        """
        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]

        # Create two separate correlation matrices
        self._plot_single_correlation_matrix(df, cancer_samples, 'Cancer', '#E74C3C')
        self._plot_single_correlation_matrix(df, normal_samples, 'Normal', '#3498DB')

    def _plot_single_correlation_matrix(self, df: pd.DataFrame, samples: list,
                                        group_name: str, color: str):
        """
        Create a single correlation matrix for a group of samples

        Pipeline: TIC Normalization → Log2 Transform → Pearson Correlation

        Args:
            df: Annotated DataFrame
            samples: List of sample columns
            group_name: Name of the group (Cancer/Normal)
            color: Color for the group
        """
        # Prepare intensity matrix
        intensity_data = replace_empty_with_zero(df[samples])

        # Step 1: TIC (Total Ion Current) Normalization
        sample_sums = intensity_data.sum(axis=0)
        median_sum = sample_sums.median()
        sample_sums_safe = sample_sums.replace(0, 1)
        intensity_normalized = intensity_data / sample_sums_safe * median_sum

        # Step 2: Log2 transform
        intensity_log = np.log2(intensity_normalized + 1)

        # Step 3: Calculate correlation matrix
        corr_matrix = intensity_log.corr(method='pearson')

        # Create figure
        fig, ax = plt.subplots(figsize=(14, 12))

        # Create mask for upper triangle (we'll show values there)
        mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)

        # Plot heatmap
        sns.heatmap(corr_matrix, ax=ax,
                    mask=~mask,  # Show upper triangle only
                    cmap='RdYlBu_r', center=0.8, vmin=0.5, vmax=1.0,
                    square=True, linewidths=0.5, linecolor='white',
                    cbar_kws={'label': 'Pearson Correlation', 'shrink': 0.8},
                    annot=True, fmt='.2f', annot_kws={'size': 8})

        # Plot lower triangle as colored squares without numbers
        sns.heatmap(corr_matrix, ax=ax,
                    mask=mask,  # Show lower triangle
                    cmap='RdYlBu_r', center=0.8, vmin=0.5, vmax=1.0,
                    square=True, linewidths=0.5, linecolor='white',
                    cbar=False, annot=False)

        ax.set_title(f'Sample Correlation Matrix: {group_name} Samples\n'
                     '(Upper triangle: values, Lower triangle: heatmap)',
                     fontsize=14, fontweight='bold', pad=20)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / f'correlation_matrix_{group_name.lower()}.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved {group_name} correlation matrix to {output_file}")

        # Save trace data
        corr_trace = corr_matrix.copy()
        corr_trace['Sample'] = corr_trace.index
        save_trace_data(corr_trace, self.output_dir, f'correlation_matrix_{group_name.lower()}_data.csv')

        plt.close()

    def plot_correlation_clustermap(self, df: pd.DataFrame):
        """
        Create hierarchical clustered correlation heatmaps for Cancer and Normal

        Args:
            df: Annotated DataFrame with intensity data
        """
        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]

        # Create clustermaps
        self._plot_single_clustermap(df, cancer_samples, 'Cancer')
        self._plot_single_clustermap(df, normal_samples, 'Normal')

    def _plot_single_clustermap(self, df: pd.DataFrame, samples: list, group_name: str):
        """
        Create a single correlation clustermap

        Pipeline: TIC Normalization → Log2 Transform → Pearson Correlation → Hierarchical Clustering

        Args:
            df: Annotated DataFrame
            samples: List of sample columns
            group_name: Name of the group (Cancer/Normal)
        """
        # Prepare intensity matrix
        intensity_data = replace_empty_with_zero(df[samples])

        # Step 1: TIC (Total Ion Current) Normalization
        sample_sums = intensity_data.sum(axis=0)
        median_sum = sample_sums.median()
        sample_sums_safe = sample_sums.replace(0, 1)
        intensity_normalized = intensity_data / sample_sums_safe * median_sum

        # Step 2: Log2 transform
        intensity_log = np.log2(intensity_normalized + 1)

        # Step 3: Calculate correlation matrix
        corr_matrix = intensity_log.corr(method='pearson')

        # Create clustermap
        g = sns.clustermap(corr_matrix,
                           cmap='RdYlBu_r', center=0.8, vmin=0.5, vmax=1.0,
                           linewidths=0.5, linecolor='white',
                           cbar_kws={'label': 'Pearson Correlation'},
                           figsize=(12, 10),
                           dendrogram_ratio=0.15)

        g.fig.suptitle(f'Hierarchical Clustering: {group_name} Sample Correlation',
                       fontsize=14, fontweight='bold', y=0.98)

        # Save plot
        output_file = self.output_dir / f'correlation_clustermap_{group_name.lower()}.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved {group_name} correlation clustermap to {output_file}")

        # Get dendrogram linkage for trace
        linkage_matrix = g.dendrogram_row.linkage
        linkage_df = pd.DataFrame(
            linkage_matrix,
            columns=['cluster1', 'cluster2', 'distance', 'n_samples']
        )
        save_trace_data(linkage_df, self.output_dir, f'correlation_clustermap_{group_name.lower()}_linkage_data.csv')

        plt.close()
