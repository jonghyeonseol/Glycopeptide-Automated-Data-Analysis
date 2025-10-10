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
from .plot_config import (
    save_publication_figure, DPI_COMPLEX,
    HEATMAP_CMAP_CORRELATION, enhance_heatmap_colorbar, apply_publication_theme  # ✨ Enhanced styling
)

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

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / f'correlation_matrix_{group_name.lower()}.png'
        save_publication_figure(fig, output_file, dpi=DPI_COMPLEX)
        logger.info(f"Saved {group_name} correlation matrix to {output_file} (optimized, {DPI_COMPLEX} DPI)")

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
        save_publication_figure(g.fig, output_file, dpi=DPI_COMPLEX)
        logger.info(f"Saved {group_name} correlation clustermap to {output_file} (optimized, {DPI_COMPLEX} DPI)")

        # Get dendrogram linkage for trace
        linkage_matrix = g.dendrogram_row.linkage
        linkage_df = pd.DataFrame(
            linkage_matrix,
            columns=['cluster1', 'cluster2', 'distance', 'n_samples']
        )
        save_trace_data(linkage_df, self.output_dir, f'correlation_clustermap_{group_name.lower()}_linkage_data.csv')

        plt.close()

    def plot_correlation_matrix_combined(self, df: pd.DataFrame):
        """
        Create combined correlation matrix for all samples (Cancer + Normal)

        Shows 4 quadrants:
        - Top-left: Cancer-Cancer correlations
        - Bottom-right: Normal-Normal correlations
        - Top-right & Bottom-left: Cancer-Normal cross-correlations (most informative)

        Args:
            df: Annotated DataFrame with intensity data
        """
        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]
        all_samples = cancer_samples + normal_samples

        # Prepare intensity matrix
        intensity_data = replace_empty_with_zero(df[all_samples])

        # Step 1: TIC (Total Ion Current) Normalization
        sample_sums = intensity_data.sum(axis=0)
        median_sum = sample_sums.median()
        sample_sums_safe = sample_sums.replace(0, 1)
        intensity_normalized = intensity_data / sample_sums_safe * median_sum

        # Step 2: Log2 transform
        intensity_log = np.log2(intensity_normalized + 1)

        # Step 3: Calculate correlation matrix
        corr_matrix = intensity_log.corr(method='pearson')

        # Create figure with larger size to accommodate all samples
        fig, ax = plt.subplots(figsize=(18, 16))

        # Create custom colormap annotation for group boundaries
        # Add boundary lines to separate Cancer/Normal quadrants
        n_cancer = len(cancer_samples)

        # Plot heatmap
        sns.heatmap(corr_matrix, ax=ax,
                    cmap='RdYlBu_r', center=0.8, vmin=0.5, vmax=1.0,
                    square=True, linewidths=0.5, linecolor='white',
                    cbar_kws={'label': 'Pearson Correlation', 'shrink': 0.6},
                    annot=False)  # Too many samples for annotations

        # Add quadrant boundary lines
        ax.axhline(y=n_cancer, color='black', linewidth=2, linestyle='--', alpha=0.5)
        ax.axvline(x=n_cancer, color='black', linewidth=2, linestyle='--', alpha=0.5)

        # Add group labels
        cancer_center = n_cancer / 2
        normal_center = n_cancer + len(normal_samples) / 2

        ax.text(cancer_center, -1, 'Cancer', ha='center', va='top', fontsize=14, fontweight='bold', color='#E41A1C')
        ax.text(normal_center, -1, 'Normal', ha='center', va='top', fontsize=14, fontweight='bold', color='#377EB8')
        ax.text(-1, cancer_center, 'Cancer', ha='right', va='center', fontsize=14, fontweight='bold', color='#E41A1C', rotation=90)
        ax.text(-1, normal_center, 'Normal', ha='right', va='center', fontsize=14, fontweight='bold', color='#377EB8', rotation=90)

        # Add quadrant annotations
        ax.text(cancer_center, cancer_center, 'Cancer-Cancer\n(Within-Group)',
                ha='center', va='center', fontsize=10, alpha=0.3, fontweight='bold')
        ax.text(normal_center, normal_center, 'Normal-Normal\n(Within-Group)',
                ha='center', va='center', fontsize=10, alpha=0.3, fontweight='bold')
        ax.text(normal_center, cancer_center, 'Cancer-Normal\n(Cross-Group)',
                ha='center', va='center', fontsize=10, alpha=0.3, fontweight='bold')

        ax.set_title('Sample Correlation Matrix: Combined Analysis (Cancer + Normal)\n'
                     'Cross-group correlations reveal biological vs technical variation',
                     fontsize=14, fontweight='bold', pad=20)

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'correlation_matrix_combined.png'
        save_publication_figure(fig, output_file, dpi=DPI_COMPLEX)
        logger.info(f"Saved combined correlation matrix to {output_file} (optimized, {DPI_COMPLEX} DPI)")

        # Save trace data
        corr_trace = corr_matrix.copy()
        corr_trace['Sample'] = corr_trace.index
        save_trace_data(corr_trace, self.output_dir, 'correlation_matrix_combined_data.csv')

        plt.close()

    def plot_correlation_cross_group(self, df: pd.DataFrame):
        """
        Create Cancer vs Normal cross-correlation heatmap

        Rectangular matrix: Cancer samples (rows) × Normal samples (columns)
        Shows only cross-group correlations (most informative for group comparison)

        Args:
            df: Annotated DataFrame with intensity data
        """
        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]
        all_samples = cancer_samples + normal_samples

        # Prepare intensity matrix
        intensity_data = replace_empty_with_zero(df[all_samples])

        # Step 1: TIC normalization
        sample_sums = intensity_data.sum(axis=0)
        median_sum = sample_sums.median()
        sample_sums_safe = sample_sums.replace(0, 1)
        intensity_normalized = intensity_data / sample_sums_safe * median_sum

        # Step 2: Log2 transform
        intensity_log = np.log2(intensity_normalized + 1)

        # Step 3: Calculate full correlation matrix
        corr_matrix = intensity_log.corr(method='pearson')

        # Extract cross-group correlations only
        cross_corr = corr_matrix.loc[cancer_samples, normal_samples]

        # Create figure
        fig, ax = plt.subplots(figsize=(14, 12))

        # Plot heatmap
        sns.heatmap(cross_corr, ax=ax,
                    cmap='RdYlBu_r', center=0.7, vmin=0.4, vmax=1.0,
                    square=False, linewidths=0.5, linecolor='white',
                    cbar_kws={'label': 'Pearson Correlation', 'shrink': 0.8},
                    annot=False)

        ax.set_title('Cancer vs Normal Cross-Group Correlation\n'
                     'Higher values indicate similar glycosylation patterns',
                     fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('Normal Samples', fontsize=12, fontweight='bold')
        ax.set_ylabel('Cancer Samples', fontsize=12, fontweight='bold')

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'correlation_cross_group.png'
        save_publication_figure(fig, output_file, dpi=DPI_COMPLEX)
        logger.info(f"Saved cross-group correlation heatmap to {output_file} (optimized, {DPI_COMPLEX} DPI)")

        # Save trace data
        cross_corr_trace = cross_corr.copy()
        cross_corr_trace['Cancer_Sample'] = cross_corr_trace.index
        save_trace_data(cross_corr_trace, self.output_dir, 'correlation_cross_group_data.csv')

        # Calculate statistics
        mean_corr = cross_corr.values.mean()
        std_corr = cross_corr.values.std()
        min_corr = cross_corr.values.min()
        max_corr = cross_corr.values.max()

        logger.info(f"Cross-group correlation statistics:")
        logger.info(f"  Mean: {mean_corr:.3f}")
        logger.info(f"  Std: {std_corr:.3f}")
        logger.info(f"  Range: {min_corr:.3f} - {max_corr:.3f}")

        plt.close()

    def plot_correlation_clustermap_combined(self, df: pd.DataFrame):
        """
        Create hierarchical clustered correlation heatmap for all samples

        Shows whether samples cluster by disease state (Cancer vs Normal)
        or by other factors (batch effects, individual variation)

        Args:
            df: Annotated DataFrame with intensity data
        """
        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]
        all_samples = cancer_samples + normal_samples

        # Prepare intensity matrix
        intensity_data = replace_empty_with_zero(df[all_samples])

        # Step 1: TIC normalization
        sample_sums = intensity_data.sum(axis=0)
        median_sum = sample_sums.median()
        sample_sums_safe = sample_sums.replace(0, 1)
        intensity_normalized = intensity_data / sample_sums_safe * median_sum

        # Step 2: Log2 transform
        intensity_log = np.log2(intensity_normalized + 1)

        # Step 3: Calculate correlation matrix
        corr_matrix = intensity_log.corr(method='pearson')

        # Create color annotation for sample groups
        sample_colors = []
        for sample in all_samples:
            if sample in cancer_samples:
                sample_colors.append('#E41A1C')  # Cancer = Red
            else:
                sample_colors.append('#377EB8')  # Normal = Blue

        # Create clustermap
        g = sns.clustermap(corr_matrix,
                           cmap='RdYlBu_r', center=0.8, vmin=0.5, vmax=1.0,
                           linewidths=0.5, linecolor='white',
                           cbar_kws={'label': 'Pearson Correlation'},
                           figsize=(16, 14),
                           dendrogram_ratio=0.1,
                           row_colors=sample_colors,
                           col_colors=sample_colors)

        g.fig.suptitle('Hierarchical Clustering: All Samples (Cancer + Normal)\n'
                       'Color bar: Red=Cancer, Blue=Normal',
                       fontsize=14, fontweight='bold', y=0.98)

        # Save plot
        output_file = self.output_dir / 'correlation_clustermap_combined.png'
        save_publication_figure(g.fig, output_file, dpi=DPI_COMPLEX)
        logger.info(f"Saved combined correlation clustermap to {output_file} (optimized, {DPI_COMPLEX} DPI)")

        # Get dendrogram linkage for trace
        linkage_matrix = g.dendrogram_row.linkage
        linkage_df = pd.DataFrame(
            linkage_matrix,
            columns=['cluster1', 'cluster2', 'distance', 'n_samples']
        )
        save_trace_data(linkage_df, self.output_dir, 'correlation_clustermap_combined_linkage_data.csv')

        # Get reordered sample order for trace
        reordered_indices = g.dendrogram_row.reordered_ind
        reordered_samples = [all_samples[i] for i in reordered_indices]
        reorder_df = pd.DataFrame({
            'Original_Order': range(len(all_samples)),
            'Sample': all_samples,
            'Clustered_Order': reordered_indices,
            'Clustered_Sample': reordered_samples
        })
        save_trace_data(reorder_df, self.output_dir, 'correlation_clustermap_combined_sample_order.csv')

        plt.close()
