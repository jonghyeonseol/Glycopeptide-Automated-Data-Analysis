"""
Heatmap Plot Module for pGlyco Auto Combine
Handles heatmap visualizations
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from utils import replace_empty_with_zero, save_trace_data, get_sample_columns

logger = logging.getLogger(__name__)


class HeatmapMixin:
    """Mixin class for heatmap-related plots"""

    def plot_heatmap(self, df: pd.DataFrame, figsize: tuple = (16, 12), top_n: int = 50):
        """
        Create clustered heatmap of top glycopeptides with hierarchical clustering

        Args:
            df: Annotated DataFrame
            figsize: Figure size
            top_n: Number of top glycopeptides to show
        """
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Get intensity matrix
        intensity_matrix = replace_empty_with_zero(df[sample_cols])
        # Calculate mean intensity across all samples
        df_copy = df.copy()
        df_copy['MeanIntensity'] = intensity_matrix.mean(axis=1)

        # Select top N glycopeptides
        top_glycopeptides = df_copy.nlargest(top_n, 'MeanIntensity')

        # Prepare data for heatmap
        heatmap_data = replace_empty_with_zero(top_glycopeptides[sample_cols])
        # Log transform
        heatmap_data = np.log2(heatmap_data + 1)

        # Transpose for sample clustering (samples as rows)
        heatmap_data_t = heatmap_data.T

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

        # Create clustermap with hierarchical clustering on both axes
        g = sns.clustermap(
            heatmap_data,
            cmap='YlOrRd',
            figsize=figsize,
            dendrogram_ratio=0.15,
            cbar_pos=(0.02, 0.83, 0.03, 0.15),
            cbar_kws={'label': 'Log2(Intensity + 1)'},
            col_cluster=True,  # Cluster samples
            row_cluster=True,  # Cluster glycopeptides
            xticklabels=True,
            yticklabels=True,
            linewidths=0.5,
            linecolor='lightgray',
            col_colors=sample_colors,
            method='average',  # Linkage method
            metric='euclidean'  # Distance metric
        )

        # Adjust labels
        g.ax_heatmap.set_xlabel('Sample', fontsize=12)
        g.ax_heatmap.set_ylabel('Glycopeptide', fontsize=12)

        # Position title at the top with more space to avoid overlap
        g.fig.suptitle(f'Top {top_n} Glycopeptides Heatmap with Hierarchical Clustering',
                       fontsize=14, y=1.00, fontweight='bold')

        # Rotate sample labels
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

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

        # Save plot
        output_file = self.output_dir / 'heatmap_top_glycopeptides.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved clustered heatmap to {output_file}")

        # Save trace data
        trace_data = heatmap_data.copy()
        trace_data.insert(0, 'Glycopeptide', row_labels)
        save_trace_data(trace_data, self.output_dir, 'heatmap_top_glycopeptides_data.csv')

        plt.close()

    def plot_heatmap_full_profile(self, df: pd.DataFrame, figsize: tuple = (18, 14)):
        """
        Create clustered heatmap of ALL glycopeptides (full glycan profile per sample)

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Get intensity matrix for ALL glycopeptides
        intensity_matrix = replace_empty_with_zero(df[sample_cols])
        # Filter out glycopeptides with all zeros (not detected in any sample)
        row_sums = intensity_matrix.sum(axis=1)
        df_filtered = df[row_sums > 0].copy()
        heatmap_data = intensity_matrix[row_sums > 0].copy()

        logger.info(f"Full profile heatmap: {len(heatmap_data)} glycopeptides across {len(sample_cols)} samples")

        # Log transform
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

        # Create clustermap with hierarchical clustering on both axes
        g = sns.clustermap(
            heatmap_data,
            cmap='YlOrRd',
            figsize=figsize,
            dendrogram_ratio=0.1,
            cbar_pos=(0.02, 0.85, 0.03, 0.12),
            cbar_kws={'label': 'Log2(Intensity + 1)'},
            col_cluster=True,  # Cluster samples
            row_cluster=True,  # Cluster glycopeptides
            xticklabels=True,
            yticklabels=False,  # Too many glycopeptides to show labels
            linewidths=0,
            col_colors=sample_colors,
            method='average',  # Linkage method
            metric='euclidean'  # Distance metric
        )

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

        # Save plot
        output_file = self.output_dir / 'heatmap_full_glycan_profile.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved full glycan profile heatmap to {output_file}")

        # Save trace data
        trace_data = heatmap_data.copy()
        trace_data.insert(0, 'Glycopeptide', row_labels)
        save_trace_data(trace_data, self.output_dir, 'heatmap_full_glycan_profile_data.csv')

        plt.close()
