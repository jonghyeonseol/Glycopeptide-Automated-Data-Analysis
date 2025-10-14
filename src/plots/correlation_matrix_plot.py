"""
Correlation Matrix Plot Module for pGlyco Auto Combine
Visualizes sample-to-sample correlations with dynamic centering

Dependencies:
    External:
        - pandas: Data manipulation
        - numpy: Numerical computations
        - matplotlib: Plotting backend
        - seaborn: Statistical visualization (heatmap, clustermap)

    Internal:
        - src.utils: replace_empty_with_zero, save_trace_data, get_sample_columns
        - src.plots.plot_config: CORR_* constants, save_publication_figure
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging

from ..utils import replace_empty_with_zero, save_trace_data, get_sample_columns
from .plot_config import (
    CORR_FIGSIZE, CORR_VMIN, CORR_VMAX, CORR_CENTER_AUTO, CORR_CENTER_FIXED,
    CORR_CMAP, CORR_LINEWIDTH, CORR_LINECOLOR, CORR_ANNOT_FONTSIZE,
    CORR_SQUARE, CORR_DPI, HEATMAP_CMAP_CORRELATION,
    save_publication_figure, apply_publication_theme, enhance_heatmap_colorbar,
    COLOR_CANCER, COLOR_NORMAL,
    TITLE_SIZE, AXIS_LABEL_SIZE, ANNOTATION_SIZE,  # Font constants
    PLOT_LINE_LINEWIDTH,  # Phase 10.3.3: Linewidth unification
    ALPHA_MEDIUM_HIGH, ALPHA_MEDIUM_LIGHT,  # Phase 10.3.4: Alpha unification
    # Linestyle constants (Phase 10.3.8)
    THRESHOLD_LINESTYLE
)

logger = logging.getLogger(__name__)


class CorrelationMatrixPlotMixin:
    """Mixin class for correlation matrix visualization"""

    @staticmethod
    def _prepare_correlation_matrix(df: pd.DataFrame, samples: list, method: str = 'pearson') -> tuple:
        """
        Unified correlation matrix preparation pipeline.

        Pipeline: Extract Data → TIC Normalization → Log2 Transform → Correlation

        Args:
            df: Annotated DataFrame
            samples: List of sample columns
            method: Correlation method (default='pearson')

        Returns:
            Tuple of (corr_matrix, intensity_log)
        """
        # Step 1: Extract intensity data
        intensity_data = replace_empty_with_zero(df[samples])

        # Step 2: TIC (Total Ion Current) Normalization
        sample_sums = intensity_data.sum(axis=0)
        median_sum = sample_sums.median()
        sample_sums_safe = sample_sums.replace(0, 1)
        intensity_normalized = intensity_data / sample_sums_safe * median_sum

        # Step 3: Log2 transform
        intensity_log = np.log2(intensity_normalized + 1)

        # Step 4: Calculate correlation matrix
        corr_matrix = intensity_log.corr(method=method)

        return corr_matrix, intensity_log

    @staticmethod
    def _get_correlation_center(corr_matrix: pd.DataFrame, default_center: float = None) -> float:
        """
        Calculate correlation center (dynamic or fixed).

        Args:
            corr_matrix: Correlation matrix
            default_center: Optional default center for special cases

        Returns:
            Center value for heatmap color scale
        """
        if CORR_CENTER_AUTO:
            center = np.median(corr_matrix.values)
            logger.debug(f"  Using dynamic center: {center:.3f}")
        else:
            center = default_center if default_center is not None else CORR_CENTER_FIXED
            logger.debug(f"  Using fixed center: {center:.3f}")
        return center

    def plot_correlation_matrix(self, df: pd.DataFrame):
        """
        Create correlation matrix heatmaps for Cancer and Normal samples

        Args:
            df: Annotated DataFrame with intensity data
        """
        logger.info("Creating correlation matrices...")

        # Use centralized sample column extraction
        cancer_samples, normal_samples = get_sample_columns(df)

        # Create two separate correlation matrices
        self._plot_single_correlation_matrix(df, cancer_samples, 'Cancer', COLOR_CANCER)
        self._plot_single_correlation_matrix(df, normal_samples, 'Normal', COLOR_NORMAL)

        logger.info("✓ Correlation matrices created")

    def _plot_single_correlation_matrix(self, df: pd.DataFrame, samples: list,
                                        group_name: str, color: str):
        """
        Create a single correlation matrix for a group of samples

        Pipeline: TIC Normalization → Log2 Transform → Pearson Correlation

        Args:
            df: Annotated DataFrame
            samples: List of sample columns
            group_name: Name of the group (Cancer/Normal)
            color: Color for the group (unused, for backward compatibility)
        """
        logger.info(f"  Creating {group_name} correlation matrix...")

        # Prepare correlation matrix using unified helper
        corr_matrix, intensity_log = self._prepare_correlation_matrix(df, samples)

        logger.debug(f"  Correlation range: [{corr_matrix.values.min():.3f}, "
                    f"{corr_matrix.values.max():.3f}]")

        # Get correlation center using unified helper
        corr_center = self._get_correlation_center(corr_matrix)

        # Create figure
        fig, ax = plt.subplots(figsize=CORR_FIGSIZE)

        # Create mask for upper triangle (we'll show values there)
        mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)

        # ========================================
        # HEATMAP with centralized constants
        # ========================================
        # Plot upper triangle with annotations
        sns.heatmap(corr_matrix, ax=ax,
                    mask=~mask,  # Show upper triangle only
                    cmap=CORR_CMAP,
                    center=corr_center,
                    vmin=CORR_VMIN,
                    vmax=CORR_VMAX,
                    square=CORR_SQUARE,
                    linewidths=CORR_LINEWIDTH,
                    linecolor=CORR_LINECOLOR,
                    cbar_kws={'label': 'Pearson Correlation', 'shrink': 0.8},
                    annot=True,
                    fmt='.2f',
                    annot_kws={'size': CORR_ANNOT_FONTSIZE})

        # Plot lower triangle as colored squares without numbers
        sns.heatmap(corr_matrix, ax=ax,
                    mask=mask,  # Show lower triangle
                    cmap=CORR_CMAP,
                    center=corr_center,
                    vmin=CORR_VMIN,
                    vmax=CORR_VMAX,
                    square=CORR_SQUARE,
                    linewidths=CORR_LINEWIDTH,
                    linecolor=CORR_LINECOLOR,
                    cbar=False,
                    annot=False)

        center_text = "dynamic" if CORR_CENTER_AUTO else f"fixed={CORR_CENTER_FIXED}"
        ax.set_title(f'Sample Correlation Matrix: {group_name} Samples\n'
                     f'(Upper: values, Lower: heatmap, center={center_text})',
                     fontsize=TITLE_SIZE, fontweight='bold', pad=20)

        # Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        # Save with standardized function
        output_file = self.output_dir / f'correlation_matrix_{group_name.lower()}.png'
        save_publication_figure(fig, output_file, dpi=CORR_DPI)
        logger.info(f"  Saved {group_name} correlation matrix to {output_file} "
                   f"(optimized, {CORR_DPI} DPI)")

        # Save trace data
        corr_trace = corr_matrix.copy()
        corr_trace['Sample'] = corr_trace.index
        save_trace_data(corr_trace, self.output_dir,
                       f'correlation_matrix_{group_name.lower()}_data.csv')

        plt.close()

    def plot_correlation_clustermap(self, df: pd.DataFrame):
        """
        Create hierarchical clustered correlation heatmaps for Cancer and Normal

        Args:
            df: Annotated DataFrame with intensity data
        """
        logger.info("Creating correlation clustermaps...")

        # Use centralized sample column extraction
        cancer_samples, normal_samples = get_sample_columns(df)

        # Create clustermaps
        self._plot_single_clustermap(df, cancer_samples, 'Cancer')
        self._plot_single_clustermap(df, normal_samples, 'Normal')

        logger.info("✓ Correlation clustermaps created")

    def _plot_single_clustermap(self, df: pd.DataFrame, samples: list, group_name: str):
        """
        Create a single correlation clustermap

        Pipeline: TIC Normalization → Log2 Transform → Pearson Correlation → Hierarchical Clustering

        Args:
            df: Annotated DataFrame
            samples: List of sample columns
            group_name: Name of the group (Cancer/Normal)
        """
        logger.info(f"  Creating {group_name} correlation clustermap...")

        # Prepare correlation matrix using unified helper
        corr_matrix, intensity_log = self._prepare_correlation_matrix(df, samples)

        # Get correlation center using unified helper
        corr_center = self._get_correlation_center(corr_matrix)

        # Create clustermap with centralized constants
        g = sns.clustermap(corr_matrix,
                           cmap=CORR_CMAP,
                           center=corr_center,
                           vmin=CORR_VMIN,
                           vmax=CORR_VMAX,
                           linewidths=CORR_LINEWIDTH,
                           linecolor=CORR_LINECOLOR,
                           cbar_kws={'label': 'Pearson Correlation'},
                           figsize=(12, 10),
                           dendrogram_ratio=0.15)

        center_text = "dynamic" if CORR_CENTER_AUTO else f"fixed={CORR_CENTER_FIXED}"
        g.fig.suptitle(f'Hierarchical Clustering: {group_name} Sample Correlation\n'
                       f'(center={center_text})',
                       fontsize=TITLE_SIZE, fontweight='bold', y=0.98)

        # Save plot
        output_file = self.output_dir / f'correlation_clustermap_{group_name.lower()}.png'
        save_publication_figure(g.fig, output_file, dpi=CORR_DPI)
        logger.info(f"  Saved {group_name} correlation clustermap to {output_file} "
                   f"(optimized, {CORR_DPI} DPI)")

        # Get dendrogram linkage for trace
        linkage_matrix = g.dendrogram_row.linkage
        linkage_df = pd.DataFrame(
            linkage_matrix,
            columns=['cluster1', 'cluster2', 'distance', 'n_samples']
        )
        save_trace_data(linkage_df, self.output_dir,
                       f'correlation_clustermap_{group_name.lower()}_linkage_data.csv')

        plt.close()

    def plot_correlation_matrix_combined(self, df: pd.DataFrame):
        """
        Create combined correlation matrix for all samples (Cancer + Normal)

        Shows 4 quadrants:
        - Top-left: Cancer-Cancer correlations
        - Bottom-right: Normal-Normal correlations
        - Top-right & Bottom-left: Cancer-Normal cross-correlations

        Args:
            df: Annotated DataFrame with intensity data
        """
        logger.info("Creating combined correlation matrix...")

        # Use centralized sample column extraction
        cancer_samples, normal_samples = get_sample_columns(df)
        all_samples = cancer_samples + normal_samples

        logger.debug(f"  Total samples: {len(all_samples)} "
                    f"(Cancer: {len(cancer_samples)}, Normal: {len(normal_samples)})")

        # Prepare correlation matrix using unified helper
        corr_matrix, intensity_log = self._prepare_correlation_matrix(df, all_samples)

        # Get correlation center using unified helper
        corr_center = self._get_correlation_center(corr_matrix)

        # Create figure with larger size to accommodate all samples
        fig, ax = plt.subplots(figsize=(18, 16))

        # Create custom colormap annotation for group boundaries
        # Add boundary lines to separate Cancer/Normal quadrants
        n_cancer = len(cancer_samples)

        # Plot heatmap with centralized constants
        sns.heatmap(corr_matrix, ax=ax,
                    cmap=CORR_CMAP,
                    center=corr_center,
                    vmin=CORR_VMIN,
                    vmax=CORR_VMAX,
                    square=CORR_SQUARE,
                    linewidths=CORR_LINEWIDTH,
                    linecolor=CORR_LINECOLOR,
                    cbar_kws={'label': 'Pearson Correlation', 'shrink': 0.6},
                    annot=False)  # Too many samples for annotations

        # Add quadrant boundary lines
        ax.axhline(y=n_cancer, color='black', linewidth=PLOT_LINE_LINEWIDTH, linestyle=THRESHOLD_LINESTYLE, alpha=ALPHA_MEDIUM_HIGH)
        ax.axvline(x=n_cancer, color='black', linewidth=PLOT_LINE_LINEWIDTH, linestyle=THRESHOLD_LINESTYLE, alpha=ALPHA_MEDIUM_HIGH)

        # Add group labels
        cancer_center = n_cancer / 2
        normal_center = n_cancer + len(normal_samples) / 2

        ax.text(cancer_center, -1, 'Cancer', ha='center', va='top',
               fontsize=AXIS_LABEL_SIZE, fontweight='bold', color='#E41A1C')
        ax.text(normal_center, -1, 'Normal', ha='center', va='top',
               fontsize=AXIS_LABEL_SIZE, fontweight='bold', color='#377EB8')
        ax.text(-1, cancer_center, 'Cancer', ha='right', va='center',
               fontsize=AXIS_LABEL_SIZE, fontweight='bold', color='#E41A1C', rotation=90)
        ax.text(-1, normal_center, 'Normal', ha='right', va='center',
               fontsize=AXIS_LABEL_SIZE, fontweight='bold', color='#377EB8', rotation=90)

        # Add quadrant annotations
        ax.text(cancer_center, cancer_center, 'Cancer-Cancer\n(Within-Group)',
                ha='center', va='center', fontsize=ANNOTATION_SIZE, alpha=ALPHA_MEDIUM_LIGHT, fontweight='bold')
        ax.text(normal_center, normal_center, 'Normal-Normal\n(Within-Group)',
                ha='center', va='center', fontsize=ANNOTATION_SIZE, alpha=ALPHA_MEDIUM_LIGHT, fontweight='bold')
        ax.text(normal_center, cancer_center, 'Cancer-Normal\n(Cross-Group)',
                ha='center', va='center', fontsize=ANNOTATION_SIZE, alpha=ALPHA_MEDIUM_LIGHT, fontweight='bold')

        center_text = "dynamic" if CORR_CENTER_AUTO else f"fixed={CORR_CENTER_FIXED}"
        ax.set_title(f'Sample Correlation Matrix: Combined Analysis (Cancer + Normal)\n'
                     f'Cross-group correlations reveal biological vs technical variation\n'
                     f'(center={center_text})',
                     fontsize=TITLE_SIZE, fontweight='bold', pad=20)

        # Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'correlation_matrix_combined.png'
        save_publication_figure(fig, output_file, dpi=CORR_DPI)
        logger.info(f"Saved combined correlation matrix to {output_file} "
                   f"(optimized, {CORR_DPI} DPI)")

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
        logger.info("Creating cross-group correlation heatmap...")

        # Use centralized sample column extraction
        cancer_samples, normal_samples = get_sample_columns(df)
        all_samples = cancer_samples + normal_samples

        # Prepare correlation matrix using unified helper
        corr_matrix, intensity_log = self._prepare_correlation_matrix(df, all_samples)

        # Extract cross-group correlations only
        cross_corr = corr_matrix.loc[cancer_samples, normal_samples]

        # Get correlation center using unified helper (adjusted for cross-group)
        corr_center = self._get_correlation_center(cross_corr, default_center=0.7)

        # Create figure
        fig, ax = plt.subplots(figsize=CORR_FIGSIZE)

        # Plot heatmap with adjusted vmin for cross-group (typically lower correlation)
        sns.heatmap(cross_corr, ax=ax,
                    cmap=CORR_CMAP,
                    center=corr_center,
                    vmin=0.4,  # Lower vmin for cross-group
                    vmax=CORR_VMAX,
                    square=False,  # Rectangular for cross-group
                    linewidths=CORR_LINEWIDTH,
                    linecolor=CORR_LINECOLOR,
                    cbar_kws={'label': 'Pearson Correlation', 'shrink': 0.8},
                    annot=False)

        center_text = "dynamic" if CORR_CENTER_AUTO else f"adjusted={corr_center:.1f}"
        ax.set_title(f'Cancer vs Normal Cross-Group Correlation\n'
                     f'Higher values indicate similar glycosylation patterns\n'
                     f'(center={center_text})',
                     fontsize=TITLE_SIZE, fontweight='bold', pad=20)
        ax.set_xlabel('Normal Samples', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax.set_ylabel('Cancer Samples', fontsize=AXIS_LABEL_SIZE, fontweight='bold')

        # Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'correlation_cross_group.png'
        save_publication_figure(fig, output_file, dpi=CORR_DPI)
        logger.info(f"Saved cross-group correlation heatmap to {output_file} "
                   f"(optimized, {CORR_DPI} DPI)")

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
        logger.info(f"  Mean: {mean_corr:.3f}, Std: {std_corr:.3f}, "
                   f"Range: [{min_corr:.3f}, {max_corr:.3f}]")

        plt.close()

    def plot_correlation_clustermap_combined(self, df: pd.DataFrame):
        """
        Create hierarchical clustered correlation heatmap for all samples

        Shows whether samples cluster by disease state (Cancer vs Normal)
        or by other factors (batch effects, individual variation)

        Args:
            df: Annotated DataFrame with intensity data
        """
        logger.info("Creating combined correlation clustermap...")

        # Use centralized sample column extraction
        cancer_samples, normal_samples = get_sample_columns(df)
        all_samples = cancer_samples + normal_samples

        # Prepare correlation matrix using unified helper
        corr_matrix, intensity_log = self._prepare_correlation_matrix(df, all_samples)

        # Get correlation center using unified helper
        corr_center = self._get_correlation_center(corr_matrix)

        # Create color annotation for sample groups
        sample_colors = []
        for sample in all_samples:
            if sample in cancer_samples:
                sample_colors.append('#E41A1C')  # Cancer = Red
            else:
                sample_colors.append('#377EB8')  # Normal = Blue

        # Create clustermap with centralized constants
        g = sns.clustermap(corr_matrix,
                           cmap=CORR_CMAP,
                           center=corr_center,
                           vmin=CORR_VMIN,
                           vmax=CORR_VMAX,
                           linewidths=CORR_LINEWIDTH,
                           linecolor=CORR_LINECOLOR,
                           cbar_kws={'label': 'Pearson Correlation'},
                           figsize=(16, 14),
                           dendrogram_ratio=0.1,
                           row_colors=sample_colors,
                           col_colors=sample_colors)

        center_text = "dynamic" if CORR_CENTER_AUTO else f"fixed={CORR_CENTER_FIXED}"
        g.fig.suptitle(f'Hierarchical Clustering: All Samples (Cancer + Normal)\n'
                       f'Color bar: Red=Cancer, Blue=Normal (center={center_text})',
                       fontsize=TITLE_SIZE, fontweight='bold', y=0.98)

        # Save plot
        output_file = self.output_dir / 'correlation_clustermap_combined.png'
        save_publication_figure(g.fig, output_file, dpi=CORR_DPI)
        logger.info(f"Saved combined correlation clustermap to {output_file} "
                   f"(optimized, {CORR_DPI} DPI)")

        # Get dendrogram linkage for trace
        linkage_matrix = g.dendrogram_row.linkage
        linkage_df = pd.DataFrame(
            linkage_matrix,
            columns=['cluster1', 'cluster2', 'distance', 'n_samples']
        )
        save_trace_data(linkage_df, self.output_dir,
                       'correlation_clustermap_combined_linkage_data.csv')

        # Get reordered sample order for trace
        reordered_indices = g.dendrogram_row.reordered_ind
        reordered_samples = [all_samples[i] for i in reordered_indices]
        reorder_df = pd.DataFrame({
            'Original_Order': range(len(all_samples)),
            'Sample': all_samples,
            'Clustered_Order': reordered_indices,
            'Clustered_Sample': reordered_samples
        })
        save_trace_data(reorder_df, self.output_dir,
                       'correlation_clustermap_combined_sample_order.csv')

        logger.info("✓ Combined correlation clustermap complete")
        plt.close()
