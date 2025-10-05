"""
Glycopeptide Comparison Heatmap Plot Module for pGlyco Auto Combine
Compares Cancer vs Normal groups with improved visualization
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
import logging
from utils import replace_empty_with_zero, save_trace_data

logger = logging.getLogger(__name__)


class GlycopeptideComparisonHeatmapMixin:
    """Mixin class for glycopeptide comparison heatmap (Cancer vs Normal)"""

    def plot_glycopeptide_comparison_heatmap(
        self,
        df: pd.DataFrame,
        vip_scores: pd.DataFrame,
        figsize: tuple = (24, 16),
        max_peptides: int = 50,
        max_glycans_per_type: int = 15
    ):
        """
        Create dot-based heatmap comparing Cancer vs Normal groups

        Layout (top to bottom):
        1. Top Panel: Aggregated intensity line plots (Cancer vs Normal)
        2. Middle Panel: Gradient colored bar showing glycan type regions
        3. Main Panel: Dot heatmap with side-by-side comparison
        4. Bottom: Glycan composition labels (rotated 90° vertically)

        Features:
        - Y-axis: Peptides sorted by VIP score (descending)
        - X-axis: Glycan compositions grouped by type (HM, F, S, SF, C/H)
        - Dots: Side-by-side - Circle (○ left) = Cancer, Square (□ right) = Normal
        - Color: By glycan type with transparency by intensity
        - Type labels: Clean labels (HM, F, S, SF, C/H) on colored bar (middle)
        - Composition labels: Detailed labels (e.g., H(5)N(4)A(2)) below heatmap (rotated)

        Args:
            df: Annotated DataFrame with all samples
            vip_scores: VIP scores DataFrame from PLS-DA
            figsize: Figure size
            max_peptides: Maximum number of peptides to show
            max_glycans_per_type: Maximum glycans per type to show
        """
        logger.info("Creating glycopeptide comparison heatmap (Cancer vs Normal)...")

        # Define glycan type order and colors
        glycan_type_order = ['HM', 'F', 'S', 'SF', 'C/H']
        glycan_type_colors = {
            'HM': '#00CC00',      # Green
            'F': '#FF0000',       # Red
            'S': '#FF69B4',       # Pink
            'SF': '#FFA500',      # Orange
            'C/H': '#0000FF'      # Blue
        }

        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]

        # Get top peptides by VIP score
        top_peptides = vip_scores.nlargest(max_peptides * 5, 'VIP_Score')['Peptide'].unique()[:max_peptides]

        # Filter data to top peptides
        df_filtered = df[df['Peptide'].isin(top_peptides)].copy()

        # Add VIP scores to filtered data
        peptide_vip_map = vip_scores.groupby('Peptide')['VIP_Score'].mean().to_dict()
        df_filtered['PeptideVIP'] = df_filtered['Peptide'].map(peptide_vip_map)

        # Calculate aggregated intensities for Cancer and Normal
        df_filtered['Cancer_Mean'] = replace_empty_with_zero(df_filtered[cancer_samples]).mean(axis=1)
        df_filtered['Normal_Mean'] = replace_empty_with_zero(df_filtered[normal_samples]).mean(axis=1)

        # Select glycans by type (limit per type)
        selected_glycans = []
        for glycan_type in glycan_type_order:
            type_glycans = df_filtered[df_filtered['GlycanTypeCategory'] == glycan_type]

            if len(type_glycans) > 0:
                type_glycans = type_glycans.copy()
                # Sort by total intensity across both groups
                type_glycans['Total_Intensity'] = type_glycans['Cancer_Mean'] + type_glycans['Normal_Mean']

                # Get top N glycans of this type
                top_type_glycans = type_glycans.nlargest(
                    min(max_glycans_per_type, len(type_glycans)),
                    'Total_Intensity'
                )
                selected_glycans.append(top_type_glycans)

        if len(selected_glycans) == 0:
            logger.warning("No glycopeptides found for comparison")
            return

        df_plot = pd.concat(selected_glycans, ignore_index=False)

        # Create peptide order (by VIP score, descending)
        peptide_order = df_plot.groupby('Peptide')['PeptideVIP'].mean().sort_values(ascending=False).index.tolist()

        # Create glycan order (grouped by type)
        glycan_order = []
        glycan_type_positions = {}  # Track position of each glycan type
        current_pos = 0

        # Helper function for natural/numeric sorting of glycan compositions
        def glycan_sort_key(glycan_comp):
            """Extract numbers from glycan composition for natural sorting"""
            import re
            # Extract all monosaccharide types and their counts
            # Pattern: Letter(number)
            parts = re.findall(r'([A-Z]+)\((\d+)\)', glycan_comp)
            # Create a tuple of (monosaccharide, count as int) for sorting
            # Sort by: H, N, A, F, G order (typical glycan structure)
            monosaccharide_order = {'H': 0, 'N': 1, 'A': 2, 'F': 3, 'G': 4}
            sort_tuple = []
            for mono, count in parts:
                order = monosaccharide_order.get(mono, 99)
                sort_tuple.append((order, int(count)))
            return tuple(sort_tuple)

        for glycan_type in glycan_type_order:
            type_glycans = df_plot[df_plot['GlycanTypeCategory'] == glycan_type]['GlycanComposition'].unique()
            if len(type_glycans) > 0:
                # Sort glycan compositions by numeric values (natural sorting)
                type_glycans_sorted = sorted(type_glycans, key=glycan_sort_key)

                glycan_type_positions[glycan_type] = {
                    'start': current_pos,
                    'end': current_pos + len(type_glycans_sorted)
                }
                glycan_order.extend(type_glycans_sorted)
                current_pos += len(type_glycans_sorted)

        # Prepare indices
        peptide_to_idx = {p: i for i, p in enumerate(peptide_order)}
        glycan_to_idx = {g: i for i, g in enumerate(glycan_order)}

        # Create figure with three panels
        # Layout: Line plot (top) → Color bar (middle) → Heatmap (bottom)
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(3, 1, height_ratios=[1, 0.3, 5.5], hspace=0.05)

        ax_top = fig.add_subplot(gs[0])      # Line plot
        ax_colorbar = fig.add_subplot(gs[1])  # Color bar (was bottom, now middle)
        ax_main = fig.add_subplot(gs[2])     # Main heatmap

        # === TOP PANEL: Aggregated intensity comparison ===
        cancer_intensities = []
        normal_intensities = []

        for glycan in glycan_order:
            glycan_data = df_plot[df_plot['GlycanComposition'] == glycan]
            cancer_int = glycan_data['Cancer_Mean'].sum()
            normal_int = glycan_data['Normal_Mean'].sum()
            cancer_intensities.append(cancer_int)
            normal_intensities.append(normal_int)

        x_pos = np.arange(len(glycan_order))
        ax_top.plot(x_pos, cancer_intensities, color='#E74C3C', linewidth=2.5,
                   marker='o', markersize=6, label='Cancer', alpha=0.85,
                   markeredgecolor='white', markeredgewidth=0.5)
        ax_top.plot(x_pos, normal_intensities, color='#3498DB', linewidth=2.5,
                   marker='s', markersize=6, label='Normal', alpha=0.85,
                   markeredgecolor='white', markeredgewidth=0.5)

        ax_top.set_xlim(-0.5, len(glycan_order) - 0.5)
        ax_top.set_ylabel('Aggregated\nIntensity', fontsize=12, fontweight='bold')
        ax_top.set_xticks([])
        ax_top.grid(axis='y', alpha=0.4, linestyle='--', linewidth=0.8)
        ax_top.legend(loc='upper right', fontsize=11, framealpha=0.9, edgecolor='#333')
        ax_top.spines['bottom'].set_visible(False)
        ax_top.spines['top'].set_visible(False)
        ax_top.spines['right'].set_visible(False)

        # === MIDDLE PANEL: Solid colored bar for glycan types ===
        ax_colorbar.set_xlim(-0.5, len(glycan_order) - 0.5)
        ax_colorbar.set_ylim(0, 1)

        # Draw solid colored rectangles for each glycan type group
        for glycan_type, pos_info in glycan_type_positions.items():
            width = pos_info['end'] - pos_info['start']
            rect = Rectangle((pos_info['start'] - 0.5, 0), width, 1,
                            facecolor=glycan_type_colors[glycan_type],
                            edgecolor='#333', linewidth=1.5, alpha=0.85)
            ax_colorbar.add_patch(rect)

            # Add glycan type label in center - larger and more visible
            center_x = (pos_info['start'] + pos_info['end']) / 2 - 0.5
            ax_colorbar.text(center_x, 0.5, glycan_type, ha='center', va='center',
                          fontsize=16, fontweight='bold', color='white',
                          bbox=dict(boxstyle='round,pad=0.4', facecolor='black', alpha=0.6,
                                   edgecolor='white', linewidth=1.5))

        # Hide axes
        ax_colorbar.set_xticks([])
        ax_colorbar.set_yticks([])
        ax_colorbar.spines['top'].set_visible(False)
        ax_colorbar.spines['bottom'].set_visible(False)
        ax_colorbar.spines['left'].set_visible(False)
        ax_colorbar.spines['right'].set_visible(False)

        # === MAIN PANEL: Dot heatmap ===

        # Calculate max intensity PER GLYCAN TYPE for relative transparency
        max_intensity_by_type = {}
        for glycan_type in glycan_type_order:
            type_data = df_plot[df_plot['GlycanTypeCategory'] == glycan_type]
            if len(type_data) > 0:
                type_intensities = np.concatenate([
                    type_data['Cancer_Mean'].values,
                    type_data['Normal_Mean'].values
                ])
                max_intensity_by_type[glycan_type] = max(type_intensities.max(), 1)
            else:
                max_intensity_by_type[glycan_type] = 1

        logger.info("Max intensities by glycan type (for transparency):")
        for gtype, max_int in max_intensity_by_type.items():
            logger.info(f"  {gtype}: {max_int:.2e}")

        # Plot dots for both Cancer and Normal
        for idx, row in df_plot.iterrows():
            peptide = row['Peptide']
            glycan = row['GlycanComposition']
            glycan_type = row['GlycanTypeCategory']

            if peptide not in peptide_to_idx or glycan not in glycan_to_idx:
                continue

            y_pos = peptide_to_idx[peptide]
            x_pos = glycan_to_idx[glycan]

            color = glycan_type_colors.get(glycan_type, '#CCCCCC')

            # Get max intensity for THIS glycan type (relative transparency)
            type_max_intensity = max_intensity_by_type[glycan_type]

            # Symbol visualization: × for Cancer (red), + for Normal (blue)
            # Symbols placed directly on grid intersections

            # Cancer symbol (× cross) - Red
            cancer_intensity = row['Cancer_Mean']
            if cancer_intensity > 0:
                alpha = min(0.3 + (cancer_intensity / type_max_intensity) * 0.7, 1.0)
                ax_main.scatter(x_pos, y_pos, s=400, c='#E74C3C', alpha=alpha,
                              marker='x', linewidths=3.0, zorder=4,
                              label='_nolegend_')

            # Normal symbol (+ plus) - Blue
            normal_intensity = row['Normal_Mean']
            if normal_intensity > 0:
                alpha = min(0.3 + (normal_intensity / type_max_intensity) * 0.7, 1.0)
                ax_main.scatter(x_pos, y_pos, s=400, c='#3498DB', alpha=alpha,
                              marker='+', linewidths=3.0, zorder=3,
                              label='_nolegend_')

        # Add light vertical separators between glycan types
        for glycan_type, pos_info in glycan_type_positions.items():
            if pos_info['end'] < len(glycan_order):  # Don't draw after last group
                ax_main.axvline(pos_info['end'] - 0.5, color='gray',
                              linestyle='--', linewidth=1.5, alpha=0.4, zorder=1)

        # Set axis properties
        ax_main.set_xlim(-0.5, len(glycan_order) - 0.5)
        ax_main.set_ylim(-0.5, len(peptide_order) - 0.5)

        # Add x-axis glycan composition labels (rotated 45 degrees clockwise)
        ax_main.set_xticks(range(len(glycan_order)))
        ax_main.set_xticklabels(glycan_order, rotation=45, fontsize=11, ha='right')
        ax_main.set_xlabel('Glycan Composition', fontsize=12, fontweight='bold', labelpad=10)

        # Set y-axis labels (peptides) - larger font
        ax_main.set_yticks(range(len(peptide_order)))
        ax_main.set_yticklabels(peptide_order, fontsize=10)
        ax_main.set_ylabel('Peptide (sorted by VIP score)', fontsize=13, fontweight='bold')

        # Invert y-axis to have highest VIP at top
        ax_main.invert_yaxis()

        # Grid - enhanced visibility
        ax_main.grid(True, alpha=0.5, linestyle='-', linewidth=1.0, color='#BBBBBB', zorder=0)
        ax_main.set_axisbelow(True)

        # Add minor grid for better cell separation
        ax_main.set_xticks([i - 0.5 for i in range(1, len(glycan_order))], minor=True)
        ax_main.set_yticks([i - 0.5 for i in range(1, len(peptide_order))], minor=True)
        ax_main.grid(which='minor', alpha=0.3, linestyle='-', linewidth=0.5, color='#DDDDDD', zorder=0)

        # === LEGEND ===
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D

        legend_elements = []

        # Glycan type colors
        for gt in glycan_type_order:
            if gt in glycan_type_positions:
                legend_elements.append(Patch(facecolor=glycan_type_colors[gt],
                                            label=f'{gt}', alpha=0.85, edgecolor='#333'))

        # Group indicators - explain symbol visualization (larger, more visible)
        legend_elements.append(Line2D([0], [0], marker='x', color='w',
                                     markerfacecolor='#E74C3C', markersize=14,
                                     markeredgewidth=3.0,
                                     label='× Cancer', markeredgecolor='#E74C3C'))
        legend_elements.append(Line2D([0], [0], marker='+', color='w',
                                     markerfacecolor='#3498DB', markersize=14,
                                     markeredgewidth=3.0,
                                     label='+ Normal (Control)', markeredgecolor='#3498DB'))
        legend_elements.append(Patch(facecolor='gray', label='Darkness = Intensity (relative to type)',
                                    alpha=0.5, edgecolor='#333', linewidth=1))

        ax_main.legend(handles=legend_elements, loc='upper left',
                      bbox_to_anchor=(1.01, 1), frameon=True, fontsize=12,
                      title='Legend', title_fontsize=14, framealpha=0.95,
                      edgecolor='#333', fancybox=True, shadow=True)

        # === TITLE ===
        fig.suptitle('Glycopeptide Comparison Heatmap: Cancer vs Normal\n'
                    f'(Top {len(peptide_order)} peptides by VIP score)',
                    fontsize=18, fontweight='bold', y=0.995)

        # Save plot
        output_file = self.output_dir / 'glycopeptide_comparison_heatmap.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved glycopeptide comparison heatmap to {output_file}")

        # === SAVE COMPREHENSIVE TRACE DATA ===

        # NOTE: The Cancer_Mean and Normal_Mean are AGGREGATED values
        # (averaged across samples) used for plotting.
        # Individual sample columns contain ORIGINAL raw intensities from integrated.csv

        # Prepare enhanced trace data with all individual sample values
        trace_data = df_plot.copy()

        # Add plot position information
        trace_data['Plot_X_Position'] = trace_data['GlycanComposition'].map(glycan_to_idx)
        trace_data['Plot_Y_Position'] = trace_data['Peptide'].map(peptide_to_idx)

        # Add statistical information
        cancer_samples = [col for col in df_plot.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df_plot.columns if col.startswith('N') and col[1:].isdigit()]

        # Calculate individual sample statistics
        cancer_values = replace_empty_with_zero(trace_data[cancer_samples])
        normal_values = replace_empty_with_zero(trace_data[normal_samples])

        trace_data['Cancer_StdDev'] = cancer_values.std(axis=1)
        trace_data['Normal_StdDev'] = normal_values.std(axis=1)
        trace_data['Cancer_SampleCount'] = (cancer_values > 0).sum(axis=1)
        trace_data['Normal_SampleCount'] = (normal_values > 0).sum(axis=1)
        trace_data['Cancer_Min'] = cancer_values.min(axis=1)
        trace_data['Cancer_Max'] = cancer_values.max(axis=1)
        trace_data['Normal_Min'] = normal_values.min(axis=1)
        trace_data['Normal_Max'] = normal_values.max(axis=1)

        # Calculate fold change
        trace_data['Fold_Change'] = trace_data.apply(
            lambda row: row['Cancer_Mean'] / row['Normal_Mean'] if row['Normal_Mean'] > 0 else np.inf,
            axis=1
        )
        trace_data['Log2_Fold_Change'] = trace_data.apply(
            lambda row: np.log2(row['Cancer_Mean'] / row['Normal_Mean'])
                        if (row['Cancer_Mean'] > 0 and row['Normal_Mean'] > 0)
                        else np.nan,
            axis=1
        )

        # Add intensity transparency values used in plot (relative to glycan type max)
        def calculate_alpha(row, group_col):
            """Calculate alpha transparency relative to glycan type max"""
            intensity = row[group_col]
            if intensity <= 0:
                return 0
            glycan_type = row['GlycanTypeCategory']
            type_max = max_intensity_by_type.get(glycan_type, 1)
            return min(0.3 + (intensity / type_max) * 0.7, 1.0)

        trace_data['Cancer_Alpha'] = trace_data.apply(lambda row: calculate_alpha(row, 'Cancer_Mean'), axis=1)
        trace_data['Normal_Alpha'] = trace_data.apply(lambda row: calculate_alpha(row, 'Normal_Mean'), axis=1)

        # Add flags for presence in plot
        trace_data['Cancer_Dot_Plotted'] = trace_data['Cancer_Mean'] > 0
        trace_data['Normal_Dot_Plotted'] = trace_data['Normal_Mean'] > 0

        # Reorder columns for readability
        key_columns = [
            'Peptide', 'GlycanComposition', 'GlycanTypeCategory',
            'Plot_X_Position', 'Plot_Y_Position',
            'PeptideVIP',
            'Cancer_Mean', 'Cancer_StdDev', 'Cancer_SampleCount', 'Cancer_Min', 'Cancer_Max',
            'Normal_Mean', 'Normal_StdDev', 'Normal_SampleCount', 'Normal_Min', 'Normal_Max',
            'Fold_Change', 'Log2_Fold_Change',
            'Cancer_Alpha', 'Normal_Alpha',
            'Cancer_Dot_Plotted', 'Normal_Dot_Plotted'
        ]

        # All individual sample columns
        sample_columns = cancer_samples + normal_samples

        # Final column order
        final_columns = key_columns + sample_columns

        # Save comprehensive trace data
        save_trace_data(trace_data[final_columns], self.output_dir,
                       'glycopeptide_comparison_heatmap_data.csv')

        # Also save a summary file for quick review
        summary_data = trace_data[key_columns].copy()
        save_trace_data(summary_data, self.output_dir,
                       'glycopeptide_comparison_heatmap_summary.csv')

        plt.close()

        logger.info(f"Comparison heatmap complete:")
        logger.info(f"  - Cancer samples: {len(cancer_samples)}")
        logger.info(f"  - Normal samples: {len(normal_samples)}")
        logger.info(f"  - Peptides shown: {len(peptide_order)}")
        logger.info(f"  - Glycans shown: {len(glycan_order)}")
        for gt, pos in glycan_type_positions.items():
            logger.info(f"    - {gt}: {pos['end'] - pos['start']} glycans")
