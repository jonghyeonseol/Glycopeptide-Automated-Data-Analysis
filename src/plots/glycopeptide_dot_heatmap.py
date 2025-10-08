"""
Glycopeptide Dot Heatmap Plot Module for pGlyco Auto Combine
Handles dot-based heatmap visualization with VIP-sorted peptides
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import logging
from ..utils import replace_empty_with_zero, save_trace_data
from .plot_config import GLYCAN_COLORS

logger = logging.getLogger(__name__)


class GlycopeptideDotHeatmapMixin:
    """Mixin class for glycopeptide dot heatmap visualization"""

    def plot_glycopeptide_dot_heatmap(
        self,
        df: pd.DataFrame,
        vip_scores: pd.DataFrame,
        sample_name: str,
        figsize: tuple = (20, 14),
        max_peptides: int = 50,
        max_glycans_per_type: int = 20
    ):
        """
        Create dot-based heatmap showing peptide-glycan combinations

        Features:
        - Y-axis: Peptides sorted by VIP score (descending)
        - X-axis: Glycan compositions grouped by type (HM, F, S, SF, C/H)
        - Dots: Color by glycan type, transparency by intensity
        - Top panel: Line plot of aggregated intensity per glycan

        Args:
            df: Annotated DataFrame with all samples
            vip_scores: VIP scores DataFrame from PLS-DA
            sample_name: Sample to visualize (e.g., 'C01')
            figsize: Figure size
            max_peptides: Maximum number of peptides to show
            max_glycans_per_type: Maximum glycans per type to show
        """
        logger.info(f"Creating glycopeptide dot heatmap for sample {sample_name}...")

        # Define glycan type order and use standardized colors from plot_config
        glycan_type_order = ['HM', 'F', 'S', 'SF', 'C/H']
        glycan_type_colors = GLYCAN_COLORS  # Standardized colors (no conflicts with Cancer/Normal)

        # Get top peptides by VIP score
        top_peptides = vip_scores.nlargest(max_peptides * 5, 'VIP_Score')['Peptide'].unique()[:max_peptides]

        # Filter data to top peptides
        df_filtered = df[df['Peptide'].isin(top_peptides)].copy()

        # Add VIP scores to filtered data
        peptide_vip_map = vip_scores.groupby('Peptide')['VIP_Score'].mean().to_dict()
        df_filtered['PeptideVIP'] = df_filtered['Peptide'].map(peptide_vip_map)

        # Select glycans by type (limit per type to avoid overcrowding)
        selected_glycans = []
        for glycan_type in glycan_type_order:
            type_glycans = df_filtered[df_filtered['GlycanTypeCategory'] == glycan_type]

            if len(type_glycans) > 0:
                # Sort by mean intensity for this sample
                type_glycans = type_glycans.copy()
                type_glycans['SampleIntensity'] = replace_empty_with_zero(
                    type_glycans[[sample_name]]
                ).values.flatten()

                # Get top N glycans of this type
                top_type_glycans = type_glycans.nlargest(
                    min(max_glycans_per_type, len(type_glycans)),
                    'SampleIntensity'
                )
                selected_glycans.append(top_type_glycans)

        if len(selected_glycans) == 0:
            logger.warning(f"No glycopeptides found for sample {sample_name}")
            return

        df_plot = pd.concat(selected_glycans, ignore_index=False)

        # Create peptide order (by VIP score, descending)
        peptide_order = df_plot.groupby('Peptide')['PeptideVIP'].mean().sort_values(ascending=False).index.tolist()

        # Create glycan order (grouped by type)
        glycan_order = []
        for glycan_type in glycan_type_order:
            type_glycans = df_plot[df_plot['GlycanTypeCategory'] == glycan_type]['GlycanComposition'].unique()
            glycan_order.extend(type_glycans)

        # Prepare data matrix
        peptide_to_idx = {p: i for i, p in enumerate(peptide_order)}
        glycan_to_idx = {g: i for i, g in enumerate(glycan_order)}

        # Create figure with two panels
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(2, 1, height_ratios=[1, 4], hspace=0.05)

        ax_top = fig.add_subplot(gs[0])
        ax_main = fig.add_subplot(gs[1])

        # === TOP PANEL: Line plot of aggregated intensity ===
        glycan_intensities = []
        for glycan in glycan_order:
            glycan_data = df_plot[df_plot['GlycanComposition'] == glycan]
            intensity = replace_empty_with_zero(glycan_data[[sample_name]]).sum()
            glycan_intensities.append(intensity)

        ax_top.plot(range(len(glycan_order)), glycan_intensities,
                    color='#333333', linewidth=2, marker='o', markersize=4)
        ax_top.set_xlim(-0.5, len(glycan_order) - 0.5)
        ax_top.set_ylabel('Aggregated\nIntensity', fontsize=10)
        ax_top.set_xticks([])
        ax_top.grid(axis='y', alpha=0.3)
        ax_top.spines['bottom'].set_visible(False)
        ax_top.spines['top'].set_visible(False)
        ax_top.spines['right'].set_visible(False)

        # === MAIN PANEL: Dot heatmap ===

        # Normalize intensities for transparency (0-1 scale per sample)
        sample_intensities = replace_empty_with_zero(df_plot[[sample_name]]).values.flatten()
        if sample_intensities.max() > 0:
            sample_intensities / sample_intensities.max()
        else:
            pass

        # Plot dots
        for idx, row in df_plot.iterrows():
            peptide = row['Peptide']
            glycan = row['GlycanComposition']
            glycan_type = row['GlycanTypeCategory']

            if peptide not in peptide_to_idx or glycan not in glycan_to_idx:
                continue

            y_pos = peptide_to_idx[peptide]
            x_pos = glycan_to_idx[glycan]

            # Get intensity
            intensity_raw = replace_empty_with_zero(pd.DataFrame({sample_name: [row[sample_name]]})).values[0, 0]

            if intensity_raw > 0:
                # Normalize for alpha
                alpha = min(0.3 + (intensity_raw / sample_intensities.max()) * 0.7, 1.0)

                # Plot dot
                color = glycan_type_colors.get(glycan_type, '#CCCCCC')
                ax_main.scatter(x_pos, y_pos, s=200, c=color, alpha=alpha,
                                edgecolors='black', linewidths=0.5, zorder=3)

        # Add background shading for glycan type groups
        x_pos = 0
        for glycan_type in glycan_type_order:
            type_count = sum(
                1 for g in glycan_order
                if len(df_plot[df_plot['GlycanComposition'] == g]) > 0
                if df_plot[df_plot['GlycanComposition'] == g]['GlycanTypeCategory'].iloc[0] == glycan_type
            )

            if type_count > 0:
                # Add light background
                rect = Rectangle((x_pos - 0.5, -0.5), type_count, len(peptide_order),
                                 linewidth=1, edgecolor='gray', facecolor='none',
                                 linestyle='--', alpha=0.3)
                ax_main.add_patch(rect)

                # Add glycan type label at top
                ax_main.text(x_pos + type_count / 2 - 0.5, len(peptide_order) + 0.5,
                             glycan_type, ha='center', va='bottom', fontsize=10,
                             fontweight='bold', color=glycan_type_colors[glycan_type])

                x_pos += type_count

        # Set axis properties
        ax_main.set_xlim(-0.5, len(glycan_order) - 0.5)
        ax_main.set_ylim(-0.5, len(peptide_order) - 0.5)

        # Set x-axis labels (glycan compositions)
        ax_main.set_xticks(range(len(glycan_order)))
        ax_main.set_xticklabels(glycan_order, rotation=90, fontsize=8)

        # Set y-axis labels (peptides)
        ax_main.set_yticks(range(len(peptide_order)))
        ax_main.set_yticklabels(peptide_order, fontsize=8)

        ax_main.set_xlabel('Glycan Composition (grouped by type)', fontsize=12)
        ax_main.set_ylabel('Peptide (sorted by VIP score)', fontsize=12)

        # Invert y-axis to have highest VIP at top
        ax_main.invert_yaxis()

        # Grid
        ax_main.grid(True, alpha=0.2, linestyle=':', linewidth=0.5)
        ax_main.set_axisbelow(True)

        # Create legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=glycan_type_colors[gt], label=gt, alpha=0.7)
            for gt in glycan_type_order
        ]
        legend_elements.append(Patch(facecolor='gray', label='Intensity: transparency', alpha=0.5))

        ax_main.legend(handles=legend_elements, loc='upper left',
                       bbox_to_anchor=(1.02, 1), frameon=True, fontsize=10,
                       title='Glycan Type')

        # Title
        fig.suptitle(f'Glycopeptide Dot Heatmap - Sample {sample_name}\n'
                     f'(Top {len(peptide_order)} peptides by VIP score)',
                     fontsize=14, fontweight='bold', y=0.995)

        # Save plot
        output_file = self.output_dir / f'glycopeptide_dot_heatmap_{sample_name}.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved glycopeptide dot heatmap to {output_file}")

        # Save trace data
        trace_data = df_plot[['Peptide', 'GlycanComposition', 'GlycanTypeCategory',
                              'PeptideVIP', sample_name]].copy()
        save_trace_data(trace_data, self.output_dir,
                        f'glycopeptide_dot_heatmap_{sample_name}_data.csv')

        plt.close()
