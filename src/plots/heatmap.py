"""
Heatmap Plot Module for pGlyco Auto Combine
Handles heatmap visualizations
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from ..utils import replace_empty_with_zero, save_trace_data, get_sample_columns
from .plot_config import (
    HEATMAP_DPI, HEATMAP_FIGSIZE, save_publication_figure,
    HEATMAP_CMAP_INTENSITY, enhance_heatmap_colorbar, apply_publication_theme,
    COLOR_CANCER, COLOR_NORMAL,  # Premium colors
    DESIGN_SYSTEM_AVAILABLE,
    TITLE_SIZE, AXIS_LABEL_SIZE, LEGEND_SIZE, ANNOTATION_SIZE,  # Font constants
    LINE_THIN, EDGE_LINEWIDTH_THIN,  # Linewidth constants
    EDGE_COLOR_BLACK  # Edge color standardization
)

# Import premium design system if available
if DESIGN_SYSTEM_AVAILABLE:
    from .design_system import ColorSystem

logger = logging.getLogger(__name__)


class HeatmapMixin:
    """Mixin class for heatmap-related plots"""

    @staticmethod
    def _normalize_and_transform(intensity_matrix: pd.DataFrame) -> pd.DataFrame:
        """
        Unified TIC normalization and log2 transformation pipeline

        Eliminates ~14 lines of duplicated normalization logic.

        Pipeline:
            1. TIC (Total Ion Current) Normalization - median-based
            2. Log2 Transform - log2(x + 1)

        Args:
            intensity_matrix: Raw intensity DataFrame (glycopeptides × samples)

        Returns:
            Normalized and log2-transformed DataFrame

        Pattern Used:
            Helper Extraction - consolidates repeated normalization pipeline
        """
        # Step 1: TIC (Total Ion Current) Normalization
        sample_sums = intensity_matrix.sum(axis=0)
        median_sum = sample_sums.median()
        sample_sums_safe = sample_sums.replace(0, 1)
        intensity_normalized = intensity_matrix / sample_sums_safe * median_sum

        # Step 2: Log2 transform
        intensity_log2 = np.log2(intensity_normalized + 1)

        return intensity_log2

    @staticmethod
    def _create_sample_color_annotation(columns, cancer_samples, normal_samples):
        """
        Create color list for sample annotation in heatmap

        Eliminates ~6 lines of duplicated color assignment logic.

        Args:
            columns: Sample column names
            cancer_samples: List of cancer sample names
            normal_samples: List of normal sample names

        Returns:
            List of colors (red for cancer, blue for normal)

        Pattern Used:
            Helper Extraction - consolidates repeated color mapping
        """
        sample_colors = []
        for col in columns:
            if col in cancer_samples:
                sample_colors.append(COLOR_CANCER)  # Cancer - premium red
            else:
                sample_colors.append(COLOR_NORMAL)  # Normal - premium blue
        return sample_colors

    def _plot_heatmap_base(
        self, df: pd.DataFrame, heatmap_data: pd.DataFrame,
        row_labels: list, cancer_samples: list, normal_samples: list,
        figsize: tuple, title: str, output_file_path,
        trace_filename: str, dendrogram_ratio: float = 0.15,
        cbar_pos: tuple = (0.02, 0.83, 0.03, 0.15),
        linewidths: float = 0.3, yticklabels: bool = True,
        label_fontsize: int = 13, tick_fontsize: int = 11
    ):
        """
        Unified base method for heatmap visualization with hierarchical clustering

        Eliminates ~150 lines of code duplication across 2 heatmap methods.

        Args:
            df: Original annotated DataFrame (for reference)
            heatmap_data: Pre-processed intensity matrix (log2-transformed)
            row_labels: Y-axis labels for glycopeptides
            cancer_samples: List of cancer sample column names
            normal_samples: List of normal sample column names
            figsize: Figure size tuple
            title: Plot title
            output_file_path: Full path for output PNG file
            trace_filename: Filename for trace CSV
            dendrogram_ratio: Ratio of dendrogram to heatmap (default: 0.15)
            cbar_pos: Colorbar position tuple (left, bottom, width, height)
            linewidths: Cell border width (default: 0.3)
            yticklabels: Show y-axis tick labels (default: True)
            label_fontsize: Font size for axis labels
            tick_fontsize: Font size for tick labels

        Returns:
            None (saves figure to output directory)

        Pattern Used:
            Template Method Pattern - consolidates common visualization logic
        """
        # Set row labels
        heatmap_data_copy = heatmap_data.copy()
        heatmap_data_copy.index = row_labels

        # Create color annotation
        sample_colors = self._create_sample_color_annotation(
            heatmap_data_copy.columns, cancer_samples, normal_samples
        )

        # Create clustermap with publication-quality styling
        g = sns.clustermap(
            heatmap_data_copy,
            cmap=HEATMAP_CMAP_INTENSITY,
            figsize=figsize,
            dendrogram_ratio=dendrogram_ratio,
            cbar_pos=cbar_pos,
            cbar_kws={'label': 'Log2(Intensity + 1)'},
            col_cluster=True,
            row_cluster=True,
            xticklabels=True,
            yticklabels=yticklabels,
            linewidths=linewidths,
            linecolor='white',
            col_colors=sample_colors,
            method='average',
            metric='euclidean'
        )

        # Apply publication theme
        apply_publication_theme(g.fig)

        # Enhance colorbar styling
        enhance_heatmap_colorbar(g.cax, label='Log2(Intensity + 1)', fontsize=ANNOTATION_SIZE)

        # Adjust labels
        g.ax_heatmap.set_xlabel('Sample', fontsize=label_fontsize)
        g.ax_heatmap.set_ylabel('Glycopeptide', fontsize=label_fontsize)

        # Position title at the top
        g.fig.suptitle(title, fontsize=TITLE_SIZE, y=1.00, fontweight='bold')

        # Rotate sample labels
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90, fontsize=tick_fontsize)
        if yticklabels:
            plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize=tick_fontsize)

        # Enhanced legend with premium colors
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=COLOR_CANCER, label='Cancer', edgecolor=EDGE_COLOR_BLACK, linewidth=EDGE_LINEWIDTH_THIN),
            Patch(facecolor=COLOR_NORMAL, label='Normal', edgecolor=EDGE_COLOR_BLACK, linewidth=EDGE_LINEWIDTH_THIN)
        ]
        g.ax_heatmap.legend(
            handles=legend_elements,
            loc='lower left',
            bbox_to_anchor=(1.05, 0.5),
            frameon=True,
            fontsize=LEGEND_SIZE
        )

        # Save plot
        save_publication_figure(plt.gcf(), output_file_path, dpi=HEATMAP_DPI)
        logger.info(f"Saved clustered heatmap to {output_file_path} (optimized, 150 DPI)")

        # Save trace data
        trace_data = heatmap_data_copy.copy()
        trace_data.insert(0, 'Glycopeptide', row_labels)
        save_trace_data(trace_data, self.output_dir, trace_filename)

        plt.close()

    def plot_heatmap(self, df: pd.DataFrame, figsize: tuple = (16, 12), top_n: int = 20,
                     output_suffix: str = 'main'):
        """
        Create clustered heatmap of top glycopeptides with hierarchical clustering

        Pipeline: TIC Normalization → Log2 Transform → Hierarchical Clustering

        Phase 1.2 Enhancement: Default changed to top 20 for main figures
        - Top 20: Excellent readability (0.45"/label at 9" height)
        - Top 50: Use for supplementary materials

        Args:
            df: Annotated DataFrame
            figsize: Figure size (adjust based on top_n for optimal readability)
            top_n: Number of top glycopeptides to show (default: 20)
            output_suffix: Filename suffix ('main' or 'supplementary')
        """
        # Get sample columns
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Get intensity matrix and normalize
        intensity_matrix = replace_empty_with_zero(df[sample_cols])
        intensity_log2 = self._normalize_and_transform(intensity_matrix)

        # Calculate mean intensity for ranking (use intermediate normalized data before log2)
        # Re-compute just for ranking purposes
        sample_sums = intensity_matrix.sum(axis=0)
        median_sum = sample_sums.median()
        sample_sums_safe = sample_sums.replace(0, 1)
        intensity_normalized = intensity_matrix / sample_sums_safe * median_sum

        df_copy = df.copy()
        df_copy['MeanIntensity'] = intensity_normalized.mean(axis=1)

        # Select top N glycopeptides
        top_glycopeptides = df_copy.nlargest(top_n, 'MeanIntensity')

        # Extract top N from log2-transformed data
        heatmap_data = intensity_log2.loc[top_glycopeptides.index]

        # Create row labels
        row_labels = [
            f"{row['Peptide']}_{row['GlycanComposition']}_{row['GlycanType']}"
            for _, row in top_glycopeptides.iterrows()
        ]

        # Build title and output paths
        title = f'Top {top_n} Glycopeptides Heatmap with Hierarchical Clustering'
        if output_suffix == 'supplementary':
            title += ' (Supplementary)'
            output_file = self.output_dir / f'heatmap_top{top_n}_supplementary.png'
            trace_file = f'heatmap_top{top_n}_supplementary_data.csv'
        else:
            output_file = self.output_dir / f'heatmap_top{top_n}_main.png'
            trace_file = f'heatmap_top{top_n}_main_data.csv'

        # Determine font sizes based on top_n
        label_fontsize = 13 if top_n <= 20 else 12
        tick_fontsize = 11 if top_n <= 20 else 9

        # Call unified base method
        self._plot_heatmap_base(
            df=df,
            heatmap_data=heatmap_data,
            row_labels=row_labels,
            cancer_samples=cancer_samples,
            normal_samples=normal_samples,
            figsize=figsize,
            title=title,
            output_file_path=output_file,
            trace_filename=trace_file,
            dendrogram_ratio=0.15,
            cbar_pos=(0.02, 0.83, 0.03, 0.15),
            linewidths=0.3,
            yticklabels=True,
            label_fontsize=label_fontsize,
            tick_fontsize=tick_fontsize
        )

    def plot_heatmap_full_profile(self, df: pd.DataFrame, figsize: tuple = (18, 14)):
        """
        Create clustered heatmap of ALL glycopeptides (full glycan profile per sample)

        Pipeline: TIC Normalization → Log2 Transform → Hierarchical Clustering

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # Get sample columns
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Get intensity matrix and normalize
        intensity_matrix = replace_empty_with_zero(df[sample_cols])
        intensity_log2 = self._normalize_and_transform(intensity_matrix)

        # Filter out glycopeptides with all zeros (use normalized data before log2 for filtering)
        sample_sums = intensity_matrix.sum(axis=0)
        median_sum = sample_sums.median()
        sample_sums_safe = sample_sums.replace(0, 1)
        intensity_normalized = intensity_matrix / sample_sums_safe * median_sum
        row_sums = intensity_normalized.sum(axis=1)

        # Apply filter
        df_filtered = df[row_sums > 0].copy()
        heatmap_data = intensity_log2[row_sums > 0].copy()

        logger.info(f"Full profile heatmap: {len(heatmap_data)} glycopeptides across {len(sample_cols)} samples")

        # Create row labels
        row_labels = [
            f"{row['Peptide']}_{row['GlycanComposition']}_{row['GlycanType']}"
            for _, row in df_filtered.iterrows()
        ]

        # Build title and output paths
        title = f'Full Glycan Profile Heatmap ({len(heatmap_data)} glycopeptides)'
        output_file = self.output_dir / 'heatmap_full_glycan_profile.png'
        trace_file = 'heatmap_full_glycan_profile_data.csv'

        # Call unified base method with parameters for full profile
        self._plot_heatmap_base(
            df=df,
            heatmap_data=heatmap_data,
            row_labels=row_labels,
            cancer_samples=cancer_samples,
            normal_samples=normal_samples,
            figsize=figsize,
            title=title,
            output_file_path=output_file,
            trace_filename=trace_file,
            dendrogram_ratio=0.1,
            cbar_pos=(0.02, 0.85, 0.03, 0.12),
            linewidths=0,  # No borders for dense heatmap
            yticklabels=False,  # Too many glycopeptides to show labels
            label_fontsize=AXIS_LABEL_SIZE,
            tick_fontsize=11
        )
