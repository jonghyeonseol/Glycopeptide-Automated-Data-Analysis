"""
Site-Specific Glycosylation Heatmap Module for pGlyco Auto Combine
Visualizes glycan compositions per peptide site

UPDATED: Now uses centralized data preparation for consistency
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from ..utils import save_trace_data, get_sample_columns
from ..data_preparation import (
    DataPreparationConfig,
    calculate_group_statistics_standardized
)
from .plot_config import (
    EXTENDED_CATEGORY_COLORS,
    DEFAULT_FALLBACK_COLOR,
    save_publication_figure,
    DPI_COMPLEX,
    TITLE_SIZE, AXIS_LABEL_SIZE, TICK_LABEL_SIZE, LEGEND_SIZE, ANNOTATION_SIZE,  # Font constants
    ALPHA_VERY_HIGH,  # Alpha constants
    EDGE_COLOR_BLACK  # Edge color standardization
)

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
        # Get sample columns using centralized function
        cancer_samples, normal_samples = get_sample_columns(df)

        # Get top peptides by max VIP score
        top_peptides = vip_df.groupby('Peptide')['VIP_Score'].max().nlargest(top_n_peptides).index.tolist()

        # Prepare heatmap data
        heatmap_data = []
        row_labels = []
        glycan_types = []

        # STANDARDIZED: Use centralized statistics calculation
        config = DataPreparationConfig(missing_data_method='skipna')

        for peptide in top_peptides:
            # Get all glycan compositions for this peptide
            peptide_data = df[df['Peptide'] == peptide].copy()

            for idx, row in peptide_data.iterrows():
                glycan_comp = row['GlycanComposition']

                # Use standardized statistics calculation
                glycopeptide_row = peptide_data[peptide_data.index == idx]

                cancer_stats = calculate_group_statistics_standardized(
                    glycopeptide_row, cancer_samples, method=config.missing_data_method
                )
                normal_stats = calculate_group_statistics_standardized(
                    glycopeptide_row, normal_samples, method=config.missing_data_method
                )

                cancer_mean = cancer_stats['mean'].iloc[0] if not cancer_stats['mean'].isna().all() else 0
                normal_mean = normal_stats['mean'].iloc[0] if not normal_stats['mean'].isna().all() else 0

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

        # Use centralized color scheme from EXTENDED_CATEGORY_COLORS
        # All colors including fallbacks are now centralized
        glycan_type_colors = {
            'High Mannose': EXTENDED_CATEGORY_COLORS['High Mannose'],  # Green
            'Sialylated': EXTENDED_CATEGORY_COLORS['Sialylated'],  # Pink
            'Fucosylated': EXTENDED_CATEGORY_COLORS['Fucosylated'],  # Red
            'Sialofucosylated': EXTENDED_CATEGORY_COLORS['Sialofucosylated'],  # Orange
            'Complex/Hybrid': EXTENDED_CATEGORY_COLORS['Complex/Hybrid'],  # Blue
            'Outlier': EXTENDED_CATEGORY_COLORS['Outlier'],  # Gray (now centralized)
            'Unknown': EXTENDED_CATEGORY_COLORS['Unknown']  # Light gray (now centralized)
        }

        # Create annotation color map: map each type to a numeric index
        unique_types = list(set(glycan_types))
        type_to_idx = {gt: i for i, gt in enumerate(unique_types)}
        anno_indices = np.array([type_to_idx[gt] for gt in glycan_types]).reshape(1, -1)

        # Create colormap from unique colors (use centralized fallback)
        type_colors = [glycan_type_colors.get(gt, DEFAULT_FALLBACK_COLOR) for gt in unique_types]
        from matplotlib.colors import ListedColormap
        cmap_anno = ListedColormap(type_colors)

        # Plot main heatmap
        sns.heatmap(heatmap_matrix, ax=ax_main,
                    cmap='RdBu_r', center=0,
                    cbar_kws={'label': 'Log2 Fold Change\n(Cancer / Normal)', 'shrink': 0.8},
                    yticklabels=['Log2FC'],
                    xticklabels=row_labels,
                    linewidths=0.5, linecolor='white')

        ax_main.set_xlabel('Glycopeptide (Peptide_GlycanComposition)', fontsize=ANNOTATION_SIZE, fontweight='bold')
        ax_main.set_ylabel('')
        ax_main.set_title(f'Site-Specific Glycosylation: Top {top_n_peptides} Peptides by VIP Score',
                          fontsize=TITLE_SIZE, fontweight='bold', pad=20)

        # Rotate x-axis labels
        ax_main.set_xticklabels(ax_main.get_xticklabels(), rotation=90, ha='right', fontsize=TICK_LABEL_SIZE)

        # Plot annotation bar
        ax_anno.imshow(anno_indices, aspect='auto', cmap=cmap_anno, vmin=0, vmax=len(unique_types) - 1)
        ax_anno.set_xticks([])
        ax_anno.set_yticks([0])
        ax_anno.set_yticklabels(['Glycan Type'], fontsize=ANNOTATION_SIZE)
        ax_anno.set_title('Type', fontsize=ANNOTATION_SIZE, fontweight='bold')

        # Create custom legend for glycan types (use centralized fallback)
        from matplotlib.patches import Patch
        unique_types = sorted(set(glycan_types))
        legend_elements = [Patch(facecolor=glycan_type_colors.get(gt, DEFAULT_FALLBACK_COLOR),
                                 label=gt) for gt in unique_types]

        # Position legend on the RIGHT side (outside plot area, non-interruptive)
        ax_main.legend(handles=legend_elements, loc='center left',
                       bbox_to_anchor=(1.12, 0.5), ncol=1,
                       frameon=True, fontsize=LEGEND_SIZE, title='Glycan Type',
                       title_fontsize=ANNOTATION_SIZE, framealpha=ALPHA_VERY_HIGH, edgecolor=EDGE_COLOR_BLACK)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'site_specific_glycosylation_heatmap.png'
        save_publication_figure(fig, output_file, dpi=DPI_COMPLEX)
        logger.info(f"Saved site-specific heatmap to {output_file} (optimized, {DPI_COMPLEX} DPI)")

        # Save trace data
        trace_data = heatmap_df.copy()
        save_trace_data(trace_data, self.output_dir, 'site_specific_glycosylation_heatmap_data.csv')

        plt.close()

        return heatmap_df
