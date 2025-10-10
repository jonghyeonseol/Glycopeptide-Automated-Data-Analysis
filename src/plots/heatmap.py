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
    HEATMAP_CMAP_INTENSITY, enhance_heatmap_colorbar, apply_publication_theme
)

logger = logging.getLogger(__name__)


class HeatmapMixin:
    """Mixin class for heatmap-related plots"""

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
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Get intensity matrix
        intensity_matrix = replace_empty_with_zero(df[sample_cols])

        # Step 1: TIC (Total Ion Current) Normalization
        sample_sums = intensity_matrix.sum(axis=0)
        median_sum = sample_sums.median()
        sample_sums_safe = sample_sums.replace(0, 1)
        intensity_normalized = intensity_matrix / sample_sums_safe * median_sum

        # Calculate mean intensity across all samples (TIC normalized)
        df_copy = df.copy()
        df_copy['MeanIntensity'] = intensity_normalized.mean(axis=1)

        # Select top N glycopeptides
        top_glycopeptides = df_copy.nlargest(top_n, 'MeanIntensity')

        # Prepare data for heatmap (extract top N from normalized matrix)
        heatmap_data = intensity_normalized.loc[top_glycopeptides.index]

        # Step 2: Log2 transform
        heatmap_data = np.log2(heatmap_data + 1)

        # Transpose for sample clustering (samples as rows)
        heatmap_data.T

        # Create row labels (Peptide_GlycanComposition_GlycanType)
        row_labels = [
            f"{row['Peptide']}_{row['GlycanComposition']}_{row['GlycanType']}"
            for _, row in top_glycopeptides.iterrows()
        ]

        heatmap_data.index = row_labels

        # Create color annotation for sample groups
        sample_colors = []
        for col in heatmap_data.columns:
            if col.startswith('C'):
                sample_colors.append('#E74C3C')  # Cancer - red
            else:
                sample_colors.append('#3498DB')  # Normal - blue

        # ENHANCED: Create clustermap with publication-quality styling
        g = sns.clustermap(
            heatmap_data,
            cmap=HEATMAP_CMAP_INTENSITY,  # ✨ Perceptually uniform viridis
            figsize=figsize,
            dendrogram_ratio=0.15,
            cbar_pos=(0.02, 0.83, 0.03, 0.15),
            cbar_kws={'label': 'Log2(Intensity + 1)'},
            col_cluster=True,  # Cluster samples
            row_cluster=True,  # Cluster glycopeptides
            xticklabels=True,
            yticklabels=True,
            linewidths=0.3,  # ✨ Subtle cell borders
            linecolor='white',  # ✨ White borders for cleaner look
            col_colors=sample_colors,
            method='average',  # Linkage method
            metric='euclidean'  # Distance metric
        )

        # ✨ ENHANCED: Apply publication theme to entire figure
        apply_publication_theme(g.fig)

        # ✨ ENHANCED: Improve colorbar styling
        enhance_heatmap_colorbar(g.cax, label='Log2(Intensity + 1)', fontsize=11)

        # Adjust labels (Phase 1.2: larger fonts for top 20)
        label_fontsize = 13 if top_n <= 20 else 12
        g.ax_heatmap.set_xlabel('Sample', fontsize=label_fontsize)
        g.ax_heatmap.set_ylabel('Glycopeptide', fontsize=label_fontsize)

        # Position title at the top with more space to avoid overlap
        title_text = f'Top {top_n} Glycopeptides Heatmap with Hierarchical Clustering'
        if output_suffix == 'supplementary':
            title_text += ' (Supplementary)'
        g.fig.suptitle(title_text, fontsize=14, y=1.00, fontweight='bold')

        # Rotate sample labels (Phase 1.2: larger ticks for top 20)
        tick_fontsize = 11 if top_n <= 20 else 9
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90, fontsize=tick_fontsize)
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize=tick_fontsize)

        # Add legend for sample colors - position to avoid clustering line
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#E74C3C', label='Cancer'),
            Patch(facecolor='#3498DB', label='Normal')
        ]
        g.ax_heatmap.legend(
            handles=legend_elements,
            loc='lower left',
            bbox_to_anchor=(1.05, 0.5),
            frameon=True,
            fontsize=10
        )

        # Save plot (Phase 1.2: filename based on suffix)
        if output_suffix == 'main':
            output_file = self.output_dir / f'heatmap_top{top_n}_main.png'
            trace_file = f'heatmap_top{top_n}_main_data.csv'
        else:
            output_file = self.output_dir / f'heatmap_top{top_n}_supplementary.png'
            trace_file = f'heatmap_top{top_n}_supplementary_data.csv'

        save_publication_figure(plt.gcf(), output_file, dpi=HEATMAP_DPI)
        logger.info(f"Saved clustered heatmap to {output_file} (optimized, 150 DPI)")

        # Save trace data
        trace_data = heatmap_data.copy()
        trace_data.insert(0, 'Glycopeptide', row_labels)
        save_trace_data(trace_data, self.output_dir, trace_file)

        plt.close()

    def plot_heatmap_full_profile(self, df: pd.DataFrame, figsize: tuple = (18, 14)):
        """
        Create clustered heatmap of ALL glycopeptides (full glycan profile per sample)

        Pipeline: TIC Normalization → Log2 Transform → Hierarchical Clustering

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Get intensity matrix for ALL glycopeptides
        intensity_matrix = replace_empty_with_zero(df[sample_cols])

        # Step 1: TIC (Total Ion Current) Normalization
        sample_sums = intensity_matrix.sum(axis=0)
        median_sum = sample_sums.median()
        sample_sums_safe = sample_sums.replace(0, 1)
        intensity_normalized = intensity_matrix / sample_sums_safe * median_sum

        # Filter out glycopeptides with all zeros (not detected in any sample)
        row_sums = intensity_normalized.sum(axis=1)
        df_filtered = df[row_sums > 0].copy()
        heatmap_data = intensity_normalized[row_sums > 0].copy()

        logger.info(f"Full profile heatmap: {len(heatmap_data)} glycopeptides across {len(sample_cols)} samples")

        # Step 2: Log2 transform
        heatmap_data = np.log2(heatmap_data + 1)

        # Create row labels (Peptide_GlycanComposition_GlycanType)
        row_labels = [
            f"{row['Peptide']}_{row['GlycanComposition']}_{row['GlycanType']}"
            for _, row in df_filtered.iterrows()
        ]

        heatmap_data.index = row_labels

        # Create color annotation for sample groups
        sample_colors = []
        for col in heatmap_data.columns:
            if col.startswith('C'):
                sample_colors.append('#E74C3C')  # Cancer - red
            else:
                sample_colors.append('#3498DB')  # Normal - blue

        # ENHANCED: Create clustermap with publication-quality styling
        g = sns.clustermap(
            heatmap_data,
            cmap=HEATMAP_CMAP_INTENSITY,  # ✨ Perceptually uniform viridis
            figsize=figsize,
            dendrogram_ratio=0.1,
            cbar_pos=(0.02, 0.85, 0.03, 0.12),
            cbar_kws={'label': 'Log2(Intensity + 1)'},
            col_cluster=True,  # Cluster samples
            row_cluster=True,  # Cluster glycopeptides
            xticklabels=True,
            yticklabels=False,  # Too many glycopeptides to show labels
            linewidths=0,  # No borders for dense heatmap
            col_colors=sample_colors,
            method='average',  # Linkage method
            metric='euclidean'  # Distance metric
        )

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(g.fig)

        # ✨ ENHANCED: Improve colorbar styling
        enhance_heatmap_colorbar(g.cax, label='Log2(Intensity + 1)', fontsize=11)

        # Adjust labels
        g.ax_heatmap.set_xlabel('Sample', fontsize=12)
        g.ax_heatmap.set_ylabel('Glycopeptide', fontsize=12)

        # Position title at the top with more space to avoid overlap
        g.fig.suptitle(f'Full Glycan Profile Heatmap ({len(heatmap_data)} glycopeptides)',
                       fontsize=14, y=1.00, fontweight='bold')

        # Rotate sample labels
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

        # Add legend for sample colors - position to avoid clustering line
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#E74C3C', label='Cancer'),
            Patch(facecolor='#3498DB', label='Normal')
        ]
        g.ax_heatmap.legend(
            handles=legend_elements,
            loc='lower left',
            bbox_to_anchor=(1.05, 0.5),
            frameon=True,
            fontsize=10
        )

        # Save plot with optimized settings (150 DPI for complex heatmap)
        output_file = self.output_dir / 'heatmap_full_glycan_profile.png'
        save_publication_figure(plt.gcf(), output_file, dpi=HEATMAP_DPI)
        logger.info(f"Saved full glycan profile heatmap to {output_file} (optimized, 150 DPI)")

        # Save trace data
        trace_data = heatmap_data.copy()
        trace_data.insert(0, 'Glycopeptide', row_labels)
        save_trace_data(trace_data, self.output_dir, 'heatmap_full_glycan_profile_data.csv')

        plt.close()
