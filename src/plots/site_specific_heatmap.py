"""
Site-Specific Glycosylation Heatmap Module for pGlyco Auto Combine
Visualizes glycan compositions per peptide site
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from ..utils import replace_empty_with_zero, save_trace_data

logger = logging.getLogger(__name__)


class SiteSpecificHeatmapMixin:
    """Mixin class for site-specific glycosylation heatmap"""

    def plot_site_specific_heatmap(self, df: pd.DataFrame, vip_df: pd.DataFrame,
                                   top_n_peptides: int = 20,
                                   figsize: tuple = (16, 12)):
        """
        Create heatmap showing glycan compositions for top peptides

        Args:
            df: Annotated DataFrame with intensity data
            vip_df: DataFrame with VIP scores
            top_n_peptides: Number of top peptides to show
            figsize: Figure size (width, height)
        """
        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]

        # Get top peptides by max VIP score
        top_peptides = vip_df.groupby('Peptide')['VIP_Score'].max().nlargest(top_n_peptides).index.tolist()

        # Prepare heatmap data
        heatmap_data = []
        row_labels = []
        glycan_types = []

        for peptide in top_peptides:
            # Get all glycan compositions for this peptide
            peptide_data = df[df['Peptide'] == peptide].copy()

            for idx, row in peptide_data.iterrows():
                glycan_comp = row['GlycanComposition']

                # Calculate mean intensities
                cancer_values = replace_empty_with_zero(row[cancer_samples]).values.astype(float)
                normal_values = replace_empty_with_zero(row[normal_samples]).values.astype(float)

                cancer_mean = np.mean(cancer_values[cancer_values > 0]) if np.any(cancer_values > 0) else 0
                normal_mean = np.mean(normal_values[normal_values > 0]) if np.any(normal_values > 0) else 0

                # Calculate log2 fold change
                if normal_mean > 0 and cancer_mean > 0:
                    log2_fc = np.log2(cancer_mean / normal_mean)
                else:
                    log2_fc = 0

                heatmap_data.append(log2_fc)
                row_labels.append(f"{peptide}_{glycan_comp}")

                # Get glycan type
                glycan_type = row.get('SecondaryClassification', 'Unknown')
                glycan_types.append(glycan_type)

        # Create DataFrame for heatmap
        heatmap_df = pd.DataFrame({
            'Glycopeptide': row_labels,
            'Log2FC': heatmap_data,
            'GlycanType': glycan_types
        })

        # Reshape for heatmap (each row is a glycopeptide, single column for log2FC)
        heatmap_matrix = heatmap_df[['Log2FC']].T

        # Create figure with subplots for heatmap and annotation
        fig = plt.figure(figsize=figsize)

        # Create grid spec: main heatmap (90% width) + annotation bar (10% width)
        gs = fig.add_gridspec(1, 2, width_ratios=[0.95, 0.05], wspace=0.02)

        ax_main = fig.add_subplot(gs[0])
        ax_anno = fig.add_subplot(gs[1])

        # Color mapping for glycan types
        glycan_type_colors = {
            'High Mannose': '#8E44AD',
            'Sialylated': '#E74C3C',
            'Fucosylated': '#3498DB',
            'Sialofucosylated': '#16A085',
            'Complex/Hybrid': '#F39C12',
            'Outlier': '#95A5A6',
            'Unknown': '#BDC3C7'
        }

        # Create annotation color map: map each type to a numeric index
        unique_types = list(set(glycan_types))
        type_to_idx = {gt: i for i, gt in enumerate(unique_types)}
        anno_indices = np.array([type_to_idx[gt] for gt in glycan_types]).reshape(1, -1)

        # Create colormap from unique colors
        type_colors = [glycan_type_colors.get(gt, '#BDC3C7') for gt in unique_types]
        from matplotlib.colors import ListedColormap
        cmap_anno = ListedColormap(type_colors)

        # Plot main heatmap
        sns.heatmap(heatmap_matrix, ax=ax_main,
                   cmap='RdBu_r', center=0,
                   cbar_kws={'label': 'Log2 Fold Change\n(Cancer / Normal)', 'shrink': 0.8},
                   yticklabels=['Log2FC'],
                   xticklabels=row_labels,
                   linewidths=0.5, linecolor='white')

        ax_main.set_xlabel('Glycopeptide (Peptide_GlycanComposition)', fontsize=11, fontweight='bold')
        ax_main.set_ylabel('')
        ax_main.set_title(f'Site-Specific Glycosylation: Top {top_n_peptides} Peptides by VIP Score',
                         fontsize=14, fontweight='bold', pad=20)

        # Rotate x-axis labels
        ax_main.set_xticklabels(ax_main.get_xticklabels(), rotation=90, ha='right', fontsize=7)

        # Plot annotation bar
        ax_anno.imshow(anno_indices, aspect='auto', cmap=cmap_anno, vmin=0, vmax=len(unique_types)-1)
        ax_anno.set_xticks([])
        ax_anno.set_yticks([0])
        ax_anno.set_yticklabels(['Glycan Type'], fontsize=10)
        ax_anno.set_title('Type', fontsize=10, fontweight='bold')

        # Create custom legend for glycan types
        from matplotlib.patches import Patch
        unique_types = sorted(set(glycan_types))
        legend_elements = [Patch(facecolor=glycan_type_colors.get(gt, '#BDC3C7'),
                                 label=gt) for gt in unique_types]

        ax_main.legend(handles=legend_elements, loc='upper left',
                      bbox_to_anchor=(0, -0.15), ncol=len(unique_types),
                      frameon=True, fontsize=9, title='Glycan Type')

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'site_specific_glycosylation_heatmap.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved site-specific heatmap to {output_file}")

        # Save trace data
        trace_data = heatmap_df.copy()
        save_trace_data(trace_data, self.output_dir, 'site_specific_glycosylation_heatmap_data.csv')

        plt.close()

        return heatmap_df
