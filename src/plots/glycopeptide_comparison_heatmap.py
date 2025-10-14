"""
Glycopeptide Comparison Heatmap Plot Module for pGlyco Auto Combine
Compares Cancer vs Normal groups with improved visualization

UPDATED: Now uses centralized data preparation for consistency
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import logging
from ..utils import save_trace_data, calculate_fold_change, get_sample_columns
from ..data_preparation import (
    DataPreparationConfig,
    prepare_visualization_data,
    calculate_group_statistics_standardized
)
from .heatmap_helpers import (
    glycan_sort_key,
    get_symbol_info,
    create_peptide_order,
    create_glycan_order,
    create_index_mappings
)
from .plot_config import (
    EXTENDED_CATEGORY_COLORS, save_publication_figure, DPI_COMPLEX, COLOR_CANCER, COLOR_NORMAL,
    TITLE_SIZE, AXIS_LABEL_SIZE, TICK_LABEL_SIZE, LEGEND_SIZE, ANNOTATION_SIZE,
    PLOT_LINE_LINEWIDTH, PLOT_LINE_LINEWIDTH_THICK,
    LINE_MEDIUM, LINE_MEDIUM_THIN, LINE_NORMAL,
    GRID_LINEWIDTH, GRID_LINEWIDTH_THIN, GRID_LINEWIDTH_THICK,
    EDGE_LINEWIDTH_THICK,
    ALPHA_VERY_LIGHT, ALPHA_LIGHT, ALPHA_MEDIUM_LIGHT, ALPHA_MEDIUM, ALPHA_MEDIUM_HIGH,
    ALPHA_MODERATE, ALPHA_HIGH, ALPHA_VERY_HIGH, ALPHA_NEAR_OPAQUE, ALPHA_ALMOST_OPAQUE,
    POINT_ALPHA, OVERLAY_ALPHA, FRAME_ALPHA, ANNOTATION_ALPHA,
    SCATTER_SIZE_SMALL, SCATTER_SIZE_LARGE, MARKER_SIZE_SMALL,
    LINE_MARKER_SIZE, LEGEND_MARKER_SIZE, LEGEND_MARKER_SIZE_LARGE,
    EDGE_COLOR_WHITE, FRAME_EDGE_COLOR, SEPARATOR_EDGE_COLOR, MARKER_EDGE_COLOR_LIGHT,
    # Marker style constants (Phase 10.3.9)
    MARKER_CANCER, MARKER_NORMAL, MARKER_ANNOTATION, MARKER_ANNOTATION_ALT,
    # Linestyle constants (Phase 10.3.8)
    THRESHOLD_LINESTYLE, GRID_LINESTYLE_MAJOR,
    # Zorder constants (Phase 10.3.7)
    ZORDER_BACKGROUND, ZORDER_GRID, ZORDER_SEPARATOR,
    ZORDER_DATA_LOW, ZORDER_DATA_HIGH,
    ZORDER_THRESHOLD, ZORDER_ANNOTATION,
    ZORDER_OVERLAY, ZORDER_EFFECT,
    ZORDER_TOP, ZORDER_ABSOLUTE_TOP
)

logger = logging.getLogger(__name__)


class GlycopeptideComparisonHeatmapMixin:
    """Mixin class for glycopeptide comparison heatmap (Cancer vs Normal)"""

    def plot_glycopeptide_comparison_heatmap(
        self,
        df: pd.DataFrame,
        vip_scores: pd.DataFrame,
        config: DataPreparationConfig = None,
        figsize: tuple = (24, 16),
        max_peptides: int = 50,
        max_glycans_per_type: int = 15
    ):
        """
        Create dot-based heatmap comparing Cancer vs Normal groups

        Layout (top to bottom):
        1. Top Panel: Average intensity line plots (Cancer vs Normal)
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

        # Use default config if not provided
        if config is None:
            config = DataPreparationConfig(
                min_detection_pct=0.30,  # STANDARDIZED: 30% detection (same as other visualizations)
                min_samples=5,
                missing_data_method='skipna'
            )

        # Define glycan type order and use user-specified color rule from plot_config
        glycan_type_order = ['HM', 'F', 'S', 'SF', 'C/H']
        glycan_type_colors = EXTENDED_CATEGORY_COLORS  # User-specified color rule

        # STANDARDIZED DATA PREPARATION (eliminates double-filtering)
        # Data is already filtered by DataPipeline in main.py - no need to refilter
        df_with_vip = prepare_visualization_data(
            df=df,
            config=config,
            vip_scores=vip_scores,
            merge_method='left',  # FIXED: Keep all glycopeptides (no pre-filtering via inner join)
            apply_detection_filter=False,  # Data already filtered in main.py
            log_prefix="[Comparison Heatmap] "
        )

        if len(df_with_vip) == 0:
            logger.error("No glycopeptides pass detection filter!")
            return

        # Select top glycopeptides by VIP score, then get first N unique peptides
        df_sorted = df_with_vip.sort_values('VIP_Score', ascending=False)

        # Get first max_peptides unique peptides (peptide count will be approximately max_peptides)
        peptides_seen = []
        selected_rows = []

        for idx, row in df_sorted.iterrows():
            peptide = row['Peptide']

            # Track unique peptides
            if peptide not in peptides_seen:
                peptides_seen.append(peptide)

            # Include all glycopeptides for the top max_peptides unique peptides
            if len(peptides_seen) <= max_peptides:
                selected_rows.append(idx)
            elif peptide in peptides_seen[:max_peptides]:
                # Include remaining glycoforms of already-selected peptides
                selected_rows.append(idx)

        df_filtered = df_with_vip.loc[selected_rows].copy()

        # Select top glycans by type (limit per type to avoid overcrowding)
        selected_glycans = []
        for glycan_type in glycan_type_order:
            type_glycans = df_filtered[df_filtered['GlycanTypeCategory'] == glycan_type]

            if len(type_glycans) > 0:
                type_glycans = type_glycans.copy()
                # Sort by VIP score (preserve ranking by discriminative power)
                type_glycans_sorted = type_glycans.sort_values('VIP_Score', ascending=False)

                # Get top N glycans of this type
                top_type_glycans = type_glycans_sorted.head(
                    min(max_glycans_per_type, len(type_glycans))
                )
                selected_glycans.append(top_type_glycans)

        if len(selected_glycans) == 0:
            logger.warning("No glycopeptides found for comparison")
            return

        df_plot = pd.concat(selected_glycans, ignore_index=False)

        # Create peptide order (by max VIP score per peptide, descending)
        peptide_order = create_peptide_order(df_plot)

        # Create glycan order (grouped by type with natural sorting)
        glycan_order, glycan_type_positions = create_glycan_order(
            df_plot, glycan_type_order, include_positions=True
        )

        # Prepare indices
        peptide_to_idx, glycan_to_idx = create_index_mappings(peptide_order, glycan_order)

        # Create figure with three panels
        # Layout: Line plot (top) → Color bar (middle) → Heatmap (bottom)
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(3, 1, height_ratios=[1, 0.3, 5.5], hspace=0.05)

        ax_top = fig.add_subplot(gs[0])      # Line plot
        ax_colorbar = fig.add_subplot(gs[1])  # Color bar (was bottom, now middle)
        ax_main = fig.add_subplot(gs[2])     # Main heatmap

        # === TOP PANEL: Average intensity comparison ===
        cancer_intensities = []
        normal_intensities = []

        for glycan in glycan_order:
            glycan_data = df_plot[df_plot['GlycanComposition'] == glycan]
            cancer_int = glycan_data['Cancer_Mean'].mean()
            normal_int = glycan_data['Normal_Mean'].mean()
            cancer_intensities.append(cancer_int)
            normal_intensities.append(normal_int)

        x_pos = np.arange(len(glycan_order))
        ax_top.plot(x_pos, cancer_intensities, color=COLOR_CANCER, linewidth=PLOT_LINE_LINEWIDTH_THICK,
                    marker=MARKER_CANCER, markersize=LINE_MARKER_SIZE, label='Cancer', alpha=ALPHA_NEAR_OPAQUE,
                    markeredgecolor=MARKER_EDGE_COLOR_LIGHT, markeredgewidth=GRID_LINEWIDTH)
        ax_top.plot(x_pos, normal_intensities, color=COLOR_NORMAL, linewidth=PLOT_LINE_LINEWIDTH_THICK,
                    marker=MARKER_NORMAL, markersize=LINE_MARKER_SIZE, label='Normal', alpha=ALPHA_NEAR_OPAQUE,
                    markeredgecolor=MARKER_EDGE_COLOR_LIGHT, markeredgewidth=GRID_LINEWIDTH)

        ax_top.set_xlim(-0.5, len(glycan_order) - 0.5)
        ax_top.set_ylabel('Average\nIntensity', fontsize=AXIS_LABEL_SIZE, fontweight='bold')

        # Add vertical grid lines at glycan positions for PERFECT alignment with heatmap below
        ax_top.set_xticks(range(len(glycan_order)), minor=False)
        ax_top.set_xticklabels([])  # Hide labels but keep ticks for grid
        ax_top.grid(axis='y', alpha=ALPHA_MEDIUM, linestyle=GRID_LINESTYLE_MAJOR, linewidth=LINE_MEDIUM_THIN)
        ax_top.grid(axis='x', alpha=ALPHA_LIGHT, linestyle=GRID_LINESTYLE_MAJOR, linewidth=GRID_LINEWIDTH_THICK, color='#BBBBBB')

        ax_top.legend(loc='upper right', fontsize=ANNOTATION_SIZE, framealpha=FRAME_ALPHA, edgecolor=FRAME_EDGE_COLOR)
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
                             edgecolor=FRAME_EDGE_COLOR, linewidth=EDGE_LINEWIDTH_THICK, alpha=ALPHA_NEAR_OPAQUE)
            ax_colorbar.add_patch(rect)

            # Add glycan type label in center - larger and more visible
            center_x = (pos_info['start'] + pos_info['end']) / 2 - 0.5
            ax_colorbar.text(center_x, 0.5, glycan_type, ha='center', va='center',
                             fontsize=TITLE_SIZE, fontweight='bold', color='white',
                             bbox=dict(boxstyle='round,pad=0.4', facecolor='black', alpha=OVERLAY_ALPHA,
                                       edgecolor=EDGE_COLOR_WHITE, linewidth=EDGE_LINEWIDTH_THICK))

        # Hide axes
        ax_colorbar.set_xticks([])
        ax_colorbar.set_yticks([])
        ax_colorbar.spines['top'].set_visible(False)
        ax_colorbar.spines['bottom'].set_visible(False)
        ax_colorbar.spines['left'].set_visible(False)
        ax_colorbar.spines['right'].set_visible(False)

        # === MAIN PANEL: Symbol-based heatmap with qualitative highlighting ===
        logger.info("Creating symbol markers (× for Cancer, + for Normal) with bold highlighting for qualitative differences...")

        # Plot symbols for each glycopeptide
        # Design: × (Cancer) and + (Normal) with bold linewidth when qualitatively different
        for idx, row in df_plot.iterrows():
            peptide = row['Peptide']
            glycan = row['GlycanComposition']
            glycan_type = row['GlycanTypeCategory']

            if peptide not in peptide_to_idx or glycan not in glycan_to_idx:
                continue

            y_pos = peptide_to_idx[peptide]
            x_pos = glycan_to_idx[glycan]

            # Check presence in each group
            cancer_intensity = row['Cancer_Mean']
            normal_intensity = row['Normal_Mean']
            has_cancer = not np.isnan(cancer_intensity) and cancer_intensity > 0
            has_normal = not np.isnan(normal_intensity) and normal_intensity > 0

            # Qualitative difference: only one group present
            is_qualitatively_different = (has_cancer and not has_normal) or (has_normal and not has_cancer)

            # Set linewidth based on qualitative difference
            if is_qualitatively_different:
                linewidth = 5.0  # Bold for qualitative difference (NOTE: Using literal 5.0 - custom heavy emphasis for standard heatmap)
            else:
                linewidth = PLOT_LINE_LINEWIDTH_THICK  # Normal for both groups present

            # Plot Cancer marker (× symbol)
            if has_cancer:
                ax_main.scatter(x_pos, y_pos, s=SCATTER_SIZE_LARGE,
                               marker=MARKER_ANNOTATION,
                               c=glycan_type_colors[glycan_type],
                               linewidths=linewidth,
                               alpha=POINT_ALPHA,
                               zorder=ZORDER_DATA_HIGH,
                               label='_nolegend_')

            # Plot Normal marker (+ symbol)
            if has_normal:
                ax_main.scatter(x_pos, y_pos, s=SCATTER_SIZE_LARGE,
                               marker=MARKER_ANNOTATION_ALT,
                               c=glycan_type_colors[glycan_type],
                               linewidths=linewidth,
                               alpha=POINT_ALPHA,
                               zorder=ZORDER_DATA_LOW,
                               label='_nolegend_')

        # Add light vertical separators between glycan types
        for glycan_type, pos_info in glycan_type_positions.items():
            if pos_info['end'] < len(glycan_order):  # Don't draw after last group
                ax_main.axvline(pos_info['end'] - 0.5, color='gray',
                                linestyle=THRESHOLD_LINESTYLE, linewidth=EDGE_LINEWIDTH_THICK, alpha=ALPHA_MEDIUM, zorder=ZORDER_SEPARATOR)

        # Set axis properties with equal aspect ratio for square cells
        ax_main.set_aspect('equal', adjustable='box')
        ax_main.set_xlim(-0.5, len(glycan_order) - 0.5)
        ax_main.set_ylim(-0.5, len(peptide_order) - 0.5)

        # Add x-axis glycan composition labels (rotated 45 degrees clockwise)
        ax_main.set_xticks(range(len(glycan_order)))
        ax_main.set_xticklabels(glycan_order, rotation=45, fontsize=TICK_LABEL_SIZE, ha='right')
        ax_main.set_xlabel('Glycan Composition', fontsize=AXIS_LABEL_SIZE, fontweight='bold', labelpad=10)

        # Set y-axis labels (peptides) - larger font
        ax_main.set_yticks(range(len(peptide_order)))
        ax_main.set_yticklabels(peptide_order, fontsize=TICK_LABEL_SIZE)
        ax_main.set_ylabel('Peptide (sorted by VIP score)', fontsize=AXIS_LABEL_SIZE, fontweight='bold')

        # Invert y-axis to have highest VIP at top
        ax_main.invert_yaxis()

        # Grid - enhanced visibility
        ax_main.grid(True, alpha=ALPHA_MEDIUM_HIGH, linestyle=GRID_LINESTYLE_MAJOR, linewidth=LINE_NORMAL, color='#BBBBBB', zorder=ZORDER_BACKGROUND)
        ax_main.set_axisbelow(True)

        # Add minor grid for better cell separation
        ax_main.set_xticks([i - 0.5 for i in range(1, len(glycan_order))], minor=True)
        ax_main.set_yticks([i - 0.5 for i in range(1, len(peptide_order))], minor=True)
        ax_main.grid(which='minor', alpha=ALPHA_MEDIUM_LIGHT, linestyle=GRID_LINESTYLE_MAJOR, linewidth=GRID_LINEWIDTH, color='#DDDDDD', zorder=ZORDER_BACKGROUND)

        # === LEGEND ===
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D

        legend_elements = []

        # Glycan type colors
        for gt in glycan_type_order:
            if gt in glycan_type_positions:
                legend_elements.append(Patch(facecolor=glycan_type_colors[gt],
                                             label=f'{gt}', alpha=ALPHA_NEAR_OPAQUE, edgecolor=FRAME_EDGE_COLOR))

        # Symbol indicators - linewidth shows qualitative difference
        legend_elements.append(Line2D([0], [0], marker=MARKER_ANNOTATION, color='w',
                                      markerfacecolor='gray', markersize=LEGEND_MARKER_SIZE_LARGE,
                                      markeredgewidth=PLOT_LINE_LINEWIDTH_THICK, markeredgecolor=SEPARATOR_EDGE_COLOR,
                                      label='× Cancer (both groups)'))
        legend_elements.append(Line2D([0], [0], marker=MARKER_ANNOTATION_ALT, color='w',
                                      markerfacecolor='gray', markersize=LEGEND_MARKER_SIZE_LARGE,
                                      markeredgewidth=PLOT_LINE_LINEWIDTH_THICK, markeredgecolor=SEPARATOR_EDGE_COLOR,
                                      label='+ Normal (both groups)'))
        legend_elements.append(Line2D([0], [0], marker=MARKER_ANNOTATION, color='w',
                                      markerfacecolor='gray', markersize=LEGEND_MARKER_SIZE_LARGE,
                                      markeredgewidth=5.0, markeredgecolor=SEPARATOR_EDGE_COLOR,
                                      label='× Cancer only (bold)'))
        legend_elements.append(Line2D([0], [0], marker=MARKER_ANNOTATION_ALT, color='w',
                                      markerfacecolor='gray', markersize=LEGEND_MARKER_SIZE_LARGE,
                                      markeredgewidth=5.0, markeredgecolor=SEPARATOR_EDGE_COLOR,
                                      label='+ Normal only (bold)'))

        ax_main.legend(handles=legend_elements, loc='upper left',
                       bbox_to_anchor=(1.01, 1), frameon=True, fontsize=LEGEND_SIZE,
                       title='Legend', title_fontsize=AXIS_LABEL_SIZE, framealpha=ALPHA_ALMOST_OPAQUE,
                       edgecolor=FRAME_EDGE_COLOR, fancybox=True, shadow=True)

        # === TITLE ===
        fig.suptitle('Glycopeptide Comparison Heatmap: Cancer vs Normal\n'
                     f'(Top {len(peptide_order)} peptides by VIP score)',
                     fontsize=TITLE_SIZE, fontweight='bold', y=0.995)

        # Save plot
        output_file = self.output_dir / 'glycopeptide_comparison_heatmap.png'
        save_publication_figure(fig, output_file, dpi=DPI_COMPLEX)
        logger.info(f"Saved glycopeptide comparison heatmap to {output_file} (optimized, {DPI_COMPLEX} DPI)")

        # === SAVE COMPREHENSIVE TRACE DATA ===

        # NOTE: The Cancer_Mean and Normal_Mean are MEAN values
        # (averaged across samples within each glycopeptide).
        # The top panel shows the AVERAGE of these means per glycan composition.
        # Individual sample columns contain ORIGINAL raw intensities from integrated.csv

        # Prepare enhanced trace data with all individual sample values
        trace_data = df_plot.copy()

        # Add plot position information
        trace_data['Plot_X_Position'] = trace_data['GlycanComposition'].map(glycan_to_idx)
        trace_data['Plot_Y_Position'] = trace_data['Peptide'].map(peptide_to_idx)

        # Add statistical information
        cancer_samples, normal_samples = get_sample_columns(df_plot)

        # Calculate individual sample statistics using proper missing data handling
        # SCIENTIFIC VALIDITY: Uses skipna=True to avoid zero-bias
        cancer_trace_stats = calculate_group_statistics_standardized(trace_data, cancer_samples, method='skipna')
        normal_trace_stats = calculate_group_statistics_standardized(trace_data, normal_samples, method='skipna')

        trace_data['Cancer_StdDev'] = cancer_trace_stats['std']
        trace_data['Normal_StdDev'] = normal_trace_stats['std']
        trace_data['Cancer_SampleCount'] = cancer_trace_stats['count']
        trace_data['Normal_SampleCount'] = normal_trace_stats['count']
        trace_data['Cancer_Detection_Pct'] = cancer_trace_stats['count'] / len(cancer_samples)
        trace_data['Normal_Detection_Pct'] = normal_trace_stats['count'] / len(normal_samples)
        trace_data['Cancer_Min'] = cancer_trace_stats['min']
        trace_data['Cancer_Max'] = cancer_trace_stats['max']
        trace_data['Normal_Min'] = normal_trace_stats['min']
        trace_data['Normal_Max'] = normal_trace_stats['max']

        # Calculate fold change using scientifically improved method
        # SCIENTIFIC VALIDITY: Log2 fold change with pseudocount handles zeros and is symmetric
        trace_data['Fold_Change'] = trace_data.apply(
            lambda row: calculate_fold_change(row['Cancer_Mean'], row['Normal_Mean'], log_scale=False),
            axis=1
        )
        trace_data['Log2_Fold_Change'] = trace_data.apply(
            lambda row: calculate_fold_change(row['Cancer_Mean'], row['Normal_Mean'], log_scale=True),
            axis=1
        )

        # Add symbol information based on group presence
        # NOTE: Using literal 5.0 for linewidth_bold - custom heavy emphasis for standard heatmap
        symbol_info = trace_data.apply(
            lambda row: get_symbol_info(row, linewidth_bold=5.0, linewidth_normal=PLOT_LINE_LINEWIDTH_THICK),
            axis=1
        )
        trace_data['Cancer_Symbol'] = [x[0] for x in symbol_info]
        trace_data['Normal_Symbol'] = [x[1] for x in symbol_info]
        trace_data['Symbol_Linewidth'] = [x[2] for x in symbol_info]
        trace_data['Group_Presence'] = [x[3] for x in symbol_info]

        # Add uniform alpha for all symbols
        trace_data['Symbol_Alpha'] = POINT_ALPHA

        # Add flags for presence in plot
        trace_data['Cancer_Plotted'] = trace_data['Cancer_Symbol'].notna()
        trace_data['Normal_Plotted'] = trace_data['Normal_Symbol'].notna()

        # Reorder columns for readability
        key_columns = [
            'Peptide', 'GlycanComposition', 'GlycanTypeCategory',
            'Plot_X_Position', 'Plot_Y_Position',
            'VIP_Score',  # Individual glycopeptide VIP score
            'Cancer_Mean', 'Cancer_StdDev', 'Cancer_SampleCount', 'Cancer_Detection_Pct', 'Cancer_Min', 'Cancer_Max',
            'Normal_Mean', 'Normal_StdDev', 'Normal_SampleCount', 'Normal_Detection_Pct', 'Normal_Min', 'Normal_Max',
            'Fold_Change', 'Log2_Fold_Change',
            'Cancer_Symbol', 'Normal_Symbol', 'Symbol_Linewidth', 'Symbol_Alpha', 'Group_Presence',
            'Cancer_Plotted', 'Normal_Plotted'
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

        logger.info("Comparison heatmap complete:")
        logger.info(f"  - Cancer samples: {len(cancer_samples)}")
        logger.info(f"  - Normal samples: {len(normal_samples)}")
        logger.info(f"  - Peptides shown: {len(peptide_order)}")
        logger.info(f"  - Glycans shown: {len(glycan_order)}")
        for gt, pos in glycan_type_positions.items():
            logger.info(f"    - {gt}: {pos['end'] - pos['start']} glycans")

    def plot_glycopeptide_comparison_heatmap_full(
        self,
        df: pd.DataFrame,
        vip_scores: pd.DataFrame,
        config: DataPreparationConfig = None
    ):
        """
        Create FULL-SCALE dot-based heatmap comparing Cancer vs Normal groups
        Shows ALL glycopeptides without VIP filtering (complete landscape view)

        Layout (top to bottom):
        1. Top Panel: Average intensity line plots (Cancer vs Normal)
        2. Middle Panel: Gradient colored bar showing glycan type regions
        3. Main Panel: Dot heatmap with side-by-side comparison (ALL glycopeptides)
        4. Bottom: Glycan composition labels (rotated 45°)

        Features:
        - Y-axis: ALL peptides sorted by VIP score (descending)
        - X-axis: ALL glycan compositions grouped by type (HM, F, S, SF, C/H)
        - Dots: Side-by-side - × (Cancer), + (Normal)
        - Bold symbols (thick linewidth) when qualitatively different
        - Color: By glycan type with uniform transparency
        - Complete landscape: 374 peptides × 238 glycans = 2,314 glycopeptides

        Args:
            df: Annotated DataFrame with all samples
            vip_scores: VIP scores DataFrame from PLS-DA
            config: Data preparation configuration
        """
        logger.info("Creating FULL-SCALE glycopeptide comparison heatmap (ALL glycopeptides)...")

        # Use default config if not provided
        if config is None:
            config = DataPreparationConfig(
                min_detection_pct=0.30,
                min_samples=5,
                missing_data_method='skipna'
            )

        # Define glycan type order and use user-specified color rule
        glycan_type_order = ['HM', 'F', 'S', 'SF', 'C/H']
        glycan_type_colors = EXTENDED_CATEGORY_COLORS

        # STANDARDIZED DATA PREPARATION (no filtering)
        df_with_vip = prepare_visualization_data(
            df=df,
            config=config,
            vip_scores=vip_scores,
            merge_method='left',
            apply_detection_filter=False,  # Data already filtered in main.py
            log_prefix="[Full Comparison Heatmap] "
        )

        if len(df_with_vip) == 0:
            logger.error("No glycopeptides available!")
            return

        logger.info(f"Processing {len(df_with_vip)} total glycopeptides (complete dataset)")

        # Use ALL glycopeptides without filtering
        df_plot = df_with_vip.copy()

        # Create peptide order (by max VIP score per peptide, descending)
        peptide_order = create_peptide_order(df_plot)

        # Create glycan order (grouped by type with natural sorting)
        glycan_order, glycan_type_positions = create_glycan_order(
            df_plot, glycan_type_order, include_positions=True
        )

        # Prepare indices
        peptide_to_idx, glycan_to_idx = create_index_mappings(peptide_order, glycan_order)

        # Dynamic figure sizing based on data dimensions
        # Base size: 0.2 inches per peptide/glycan
        n_peptides = len(peptide_order)
        n_glycans = len(glycan_order)

        fig_height = max(30, n_peptides * 0.2 + 5)  # Min 30", +5" for top panels/title
        fig_width = max(30, n_glycans * 0.23 + 10)  # Min 30", +10" for labels/legend

        logger.info(f"Figure dimensions: {fig_width:.1f}\" × {fig_height:.1f}\" ({n_peptides} peptides × {n_glycans} glycans)")

        # Create figure with three panels
        fig = plt.figure(figsize=(fig_width, fig_height))
        gs = fig.add_gridspec(3, 1, height_ratios=[1, 0.3, 5.5], hspace=0.05)

        ax_top = fig.add_subplot(gs[0])
        ax_colorbar = fig.add_subplot(gs[1])
        ax_main = fig.add_subplot(gs[2])

        # === TOP PANEL: Average intensity comparison ===
        cancer_intensities = []
        normal_intensities = []

        for glycan in glycan_order:
            glycan_data = df_plot[df_plot['GlycanComposition'] == glycan]
            cancer_int = glycan_data['Cancer_Mean'].mean()
            normal_int = glycan_data['Normal_Mean'].mean()
            cancer_intensities.append(cancer_int)
            normal_intensities.append(normal_int)

        x_pos = np.arange(len(glycan_order))
        ax_top.plot(x_pos, cancer_intensities, color=COLOR_CANCER, linewidth=PLOT_LINE_LINEWIDTH,
                    marker=MARKER_CANCER, markersize=MARKER_SIZE_SMALL, label='Cancer', alpha=ALPHA_NEAR_OPAQUE,
                    markeredgecolor=MARKER_EDGE_COLOR_LIGHT, markeredgewidth=GRID_LINEWIDTH_THIN)
        ax_top.plot(x_pos, normal_intensities, color=COLOR_NORMAL, linewidth=PLOT_LINE_LINEWIDTH,
                    marker=MARKER_NORMAL, markersize=MARKER_SIZE_SMALL, label='Normal', alpha=ALPHA_NEAR_OPAQUE,
                    markeredgecolor=MARKER_EDGE_COLOR_LIGHT, markeredgewidth=GRID_LINEWIDTH_THIN)

        ax_top.set_xlim(-0.5, len(glycan_order) - 0.5)
        ax_top.set_ylabel('Average\nIntensity', fontsize=ANNOTATION_SIZE, fontweight='bold')

        # Add vertical grid lines at glycan positions for PERFECT alignment with heatmap below
        ax_top.set_xticks(range(len(glycan_order)), minor=False)
        ax_top.set_xticklabels([])  # Hide labels but keep ticks for grid
        ax_top.grid(axis='y', alpha=ALPHA_MEDIUM, linestyle=GRID_LINESTYLE_MAJOR, linewidth=GRID_LINEWIDTH_THICK)
        ax_top.grid(axis='x', alpha=ALPHA_LIGHT, linestyle=GRID_LINESTYLE_MAJOR, linewidth=GRID_LINEWIDTH, color='#BBBBBB')

        ax_top.legend(loc='upper right', fontsize=ANNOTATION_SIZE, framealpha=FRAME_ALPHA, edgecolor=FRAME_EDGE_COLOR)
        ax_top.spines['bottom'].set_visible(False)
        ax_top.spines['top'].set_visible(False)
        ax_top.spines['right'].set_visible(False)

        # === MIDDLE PANEL: Colored bar for glycan types ===
        ax_colorbar.set_xlim(-0.5, len(glycan_order) - 0.5)
        ax_colorbar.set_ylim(0, 1)

        for glycan_type, pos_info in glycan_type_positions.items():
            width = pos_info['end'] - pos_info['start']
            rect = Rectangle((pos_info['start'] - 0.5, 0), width, 1,
                             facecolor=glycan_type_colors[glycan_type],
                             edgecolor=FRAME_EDGE_COLOR, linewidth=LINE_MEDIUM, alpha=ALPHA_NEAR_OPAQUE)
            ax_colorbar.add_patch(rect)

            # Add glycan type label in center
            center_x = (pos_info['start'] + pos_info['end']) / 2 - 0.5
            ax_colorbar.text(center_x, 0.5, glycan_type, ha='center', va='center',
                             fontsize=AXIS_LABEL_SIZE, fontweight='bold', color='white',
                             bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=OVERLAY_ALPHA,
                                       edgecolor=EDGE_COLOR_WHITE, linewidth=LINE_NORMAL))

        ax_colorbar.set_xticks([])
        ax_colorbar.set_yticks([])
        ax_colorbar.spines['top'].set_visible(False)
        ax_colorbar.spines['bottom'].set_visible(False)
        ax_colorbar.spines['left'].set_visible(False)
        ax_colorbar.spines['right'].set_visible(False)

        # === MAIN PANEL: Symbol-based heatmap ===
        logger.info("Plotting symbols for all glycopeptides...")

        # Plot symbols for each glycopeptide
        for idx, row in df_plot.iterrows():
            peptide = row['Peptide']
            glycan = row['GlycanComposition']
            glycan_type = row['GlycanTypeCategory']

            if peptide not in peptide_to_idx or glycan not in glycan_to_idx:
                continue

            y_pos = peptide_to_idx[peptide]
            x_pos = glycan_to_idx[glycan]

            # Check presence in each group
            cancer_intensity = row['Cancer_Mean']
            normal_intensity = row['Normal_Mean']
            has_cancer = not np.isnan(cancer_intensity) and cancer_intensity > 0
            has_normal = not np.isnan(normal_intensity) and normal_intensity > 0

            # Qualitative difference: only one group present
            is_qualitatively_different = (has_cancer and not has_normal) or (has_normal and not has_cancer)

            # Set linewidth based on qualitative difference (reduced for full scale)
            if is_qualitatively_different:
                linewidth = 3.0  # Bold for qualitative difference (NOTE: Using literal 3.0 - custom full-scale emphasis)
            else:
                linewidth = EDGE_LINEWIDTH_THICK  # Normal for both groups present (reduced from 2.5)

            # Plot Cancer marker (× symbol) - reduced size
            if has_cancer:
                ax_main.scatter(x_pos, y_pos, s=SCATTER_SIZE_SMALL,  # Reduced from 400
                               marker=MARKER_ANNOTATION,
                               c=glycan_type_colors[glycan_type],
                               linewidths=linewidth,
                               alpha=POINT_ALPHA,
                               zorder=ZORDER_DATA_HIGH,
                               label='_nolegend_')

            # Plot Normal marker (+ symbol) - reduced size
            if has_normal:
                ax_main.scatter(x_pos, y_pos, s=SCATTER_SIZE_SMALL,  # Reduced from 400
                               marker=MARKER_ANNOTATION_ALT,
                               c=glycan_type_colors[glycan_type],
                               linewidths=linewidth,
                               alpha=POINT_ALPHA,
                               zorder=ZORDER_DATA_LOW,
                               label='_nolegend_')

        # Add light vertical separators between glycan types
        for glycan_type, pos_info in glycan_type_positions.items():
            if pos_info['end'] < len(glycan_order):
                ax_main.axvline(pos_info['end'] - 0.5, color='gray',
                                linestyle=THRESHOLD_LINESTYLE, linewidth=LINE_NORMAL, alpha=ALPHA_MEDIUM, zorder=ZORDER_SEPARATOR)

        # Set axis properties with equal aspect ratio
        ax_main.set_aspect('equal', adjustable='box')
        ax_main.set_xlim(-0.5, len(glycan_order) - 0.5)
        ax_main.set_ylim(-0.5, len(peptide_order) - 0.5)

        # Add x-axis glycan composition labels (rotated 45°) - smaller font
        ax_main.set_xticks(range(len(glycan_order)))
        ax_main.set_xticklabels(glycan_order, rotation=45, fontsize=ANNOTATION_SIZE, ha='right')
        ax_main.set_xlabel('Glycan Composition', fontsize=AXIS_LABEL_SIZE, fontweight='bold', labelpad=10)

        # Set y-axis labels (peptides) - smaller font
        ax_main.set_yticks(range(len(peptide_order)))
        ax_main.set_yticklabels(peptide_order, fontsize=ANNOTATION_SIZE)
        ax_main.set_ylabel('Peptide (sorted by VIP score)', fontsize=AXIS_LABEL_SIZE, fontweight='bold')

        # Invert y-axis to have highest VIP at top
        ax_main.invert_yaxis()

        # Grid - lighter for full scale
        ax_main.grid(True, alpha=ALPHA_MEDIUM_LIGHT, linestyle=GRID_LINESTYLE_MAJOR, linewidth=GRID_LINEWIDTH, color='#CCCCCC', zorder=ZORDER_BACKGROUND)
        ax_main.set_axisbelow(True)

        # Add minor grid
        ax_main.set_xticks([i - 0.5 for i in range(1, len(glycan_order))], minor=True)
        ax_main.set_yticks([i - 0.5 for i in range(1, len(peptide_order))], minor=True)
        ax_main.grid(which='minor', alpha=ALPHA_VERY_LIGHT, linestyle=GRID_LINESTYLE_MAJOR, linewidth=GRID_LINEWIDTH_THIN, color='#DDDDDD', zorder=ZORDER_BACKGROUND)

        # === LEGEND ===
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D

        legend_elements = []

        # Glycan type colors
        for gt in glycan_type_order:
            if gt in glycan_type_positions:
                legend_elements.append(Patch(facecolor=glycan_type_colors[gt],
                                             label=f'{gt}', alpha=ALPHA_NEAR_OPAQUE, edgecolor=FRAME_EDGE_COLOR))

        # Symbol indicators
        legend_elements.append(Line2D([0], [0], marker=MARKER_ANNOTATION, color='w',
                                      markerfacecolor='gray', markersize=LEGEND_MARKER_SIZE,
                                      markeredgewidth=EDGE_LINEWIDTH_THICK, markeredgecolor=SEPARATOR_EDGE_COLOR,
                                      label='× Cancer (both groups)'))
        legend_elements.append(Line2D([0], [0], marker=MARKER_ANNOTATION_ALT, color='w',
                                      markerfacecolor='gray', markersize=LEGEND_MARKER_SIZE,
                                      markeredgewidth=EDGE_LINEWIDTH_THICK, markeredgecolor=SEPARATOR_EDGE_COLOR,
                                      label='+ Normal (both groups)'))
        legend_elements.append(Line2D([0], [0], marker=MARKER_ANNOTATION, color='w',
                                      markerfacecolor='gray', markersize=LEGEND_MARKER_SIZE,
                                      markeredgewidth=3.0, markeredgecolor=SEPARATOR_EDGE_COLOR,
                                      label='× Cancer only (bold)'))
        legend_elements.append(Line2D([0], [0], marker=MARKER_ANNOTATION_ALT, color='w',
                                      markerfacecolor='gray', markersize=LEGEND_MARKER_SIZE,
                                      markeredgewidth=3.0, markeredgecolor=SEPARATOR_EDGE_COLOR,
                                      label='+ Normal only (bold)'))

        ax_main.legend(handles=legend_elements, loc='upper left',
                       bbox_to_anchor=(1.01, 1), frameon=True, fontsize=LEGEND_SIZE,
                       title='Legend', title_fontsize=AXIS_LABEL_SIZE, framealpha=ALPHA_ALMOST_OPAQUE,
                       edgecolor=FRAME_EDGE_COLOR, fancybox=True, shadow=True)

        # === TITLE ===
        fig.suptitle('Full-Scale Glycopeptide Comparison Heatmap: Cancer vs Normal\n'
                     f'(ALL {len(peptide_order)} peptides × {len(glycan_order)} glycans = {len(df_plot)} glycopeptides)',
                     fontsize=TITLE_SIZE, fontweight='bold', y=0.997)

        # Save plot
        output_file = self.output_dir / 'glycopeptide_comparison_heatmap_full.png'
        save_publication_figure(fig, output_file, dpi=DPI_COMPLEX)
        logger.info(f"Saved full-scale heatmap to {output_file} (optimized, {DPI_COMPLEX} DPI)")

        # === SAVE COMPREHENSIVE TRACE DATA ===
        trace_data = df_plot.copy()

        # Add plot position information
        trace_data['Plot_X_Position'] = trace_data['GlycanComposition'].map(glycan_to_idx)
        trace_data['Plot_Y_Position'] = trace_data['Peptide'].map(peptide_to_idx)

        # Add statistical information
        cancer_samples, normal_samples = get_sample_columns(df_plot)

        cancer_trace_stats = calculate_group_statistics_standardized(trace_data, cancer_samples, method='skipna')
        normal_trace_stats = calculate_group_statistics_standardized(trace_data, normal_samples, method='skipna')

        trace_data['Cancer_StdDev'] = cancer_trace_stats['std']
        trace_data['Normal_StdDev'] = normal_trace_stats['std']
        trace_data['Cancer_SampleCount'] = cancer_trace_stats['count']
        trace_data['Normal_SampleCount'] = normal_trace_stats['count']
        trace_data['Cancer_Detection_Pct'] = cancer_trace_stats['count'] / len(cancer_samples)
        trace_data['Normal_Detection_Pct'] = normal_trace_stats['count'] / len(normal_samples)
        trace_data['Cancer_Min'] = cancer_trace_stats['min']
        trace_data['Cancer_Max'] = cancer_trace_stats['max']
        trace_data['Normal_Min'] = normal_trace_stats['min']
        trace_data['Normal_Max'] = normal_trace_stats['max']

        # Calculate fold change
        trace_data['Fold_Change'] = trace_data.apply(
            lambda row: calculate_fold_change(row['Cancer_Mean'], row['Normal_Mean'], log_scale=False),
            axis=1
        )
        trace_data['Log2_Fold_Change'] = trace_data.apply(
            lambda row: calculate_fold_change(row['Cancer_Mean'], row['Normal_Mean'], log_scale=True),
            axis=1
        )

        # Add symbol information
        # NOTE: Using literal 3.0 for linewidth_bold - custom full-scale emphasis for trace data
        symbol_info = trace_data.apply(
            lambda row: get_symbol_info(row, linewidth_bold=3.0, linewidth_normal=EDGE_LINEWIDTH_THICK),
            axis=1
        )
        trace_data['Cancer_Symbol'] = [x[0] for x in symbol_info]
        trace_data['Normal_Symbol'] = [x[1] for x in symbol_info]
        trace_data['Symbol_Linewidth'] = [x[2] for x in symbol_info]
        trace_data['Group_Presence'] = [x[3] for x in symbol_info]
        trace_data['Symbol_Alpha'] = POINT_ALPHA
        trace_data['Cancer_Plotted'] = trace_data['Cancer_Symbol'].notna()
        trace_data['Normal_Plotted'] = trace_data['Normal_Symbol'].notna()

        # Reorder columns
        key_columns = [
            'Peptide', 'GlycanComposition', 'GlycanTypeCategory',
            'Plot_X_Position', 'Plot_Y_Position',
            'VIP_Score',
            'Cancer_Mean', 'Cancer_StdDev', 'Cancer_SampleCount', 'Cancer_Detection_Pct', 'Cancer_Min', 'Cancer_Max',
            'Normal_Mean', 'Normal_StdDev', 'Normal_SampleCount', 'Normal_Detection_Pct', 'Normal_Min', 'Normal_Max',
            'Fold_Change', 'Log2_Fold_Change',
            'Cancer_Symbol', 'Normal_Symbol', 'Symbol_Linewidth', 'Symbol_Alpha', 'Group_Presence',
            'Cancer_Plotted', 'Normal_Plotted'
        ]

        sample_columns = cancer_samples + normal_samples
        final_columns = key_columns + sample_columns

        # Save comprehensive trace data
        save_trace_data(trace_data[final_columns], self.output_dir,
                        'glycopeptide_comparison_heatmap_full_data.csv')

        # Save summary
        summary_data = trace_data[key_columns].copy()
        save_trace_data(summary_data, self.output_dir,
                        'glycopeptide_comparison_heatmap_full_summary.csv')

        plt.close()

        logger.info("Full-scale comparison heatmap complete:")
        logger.info(f"  - Cancer samples: {len(cancer_samples)}")
        logger.info(f"  - Normal samples: {len(normal_samples)}")
        logger.info(f"  - Peptides shown: {len(peptide_order)} (ALL unique peptides)")
        logger.info(f"  - Glycans shown: {len(glycan_order)} (ALL unique glycans)")
        logger.info(f"  - Total glycopeptides: {len(df_plot)}")
        for gt, pos in glycan_type_positions.items():
            logger.info(f"    - {gt}: {pos['end'] - pos['start']} glycans")

    def plot_glycopeptide_comparison_heatmap_by_type(
        self,
        df: pd.DataFrame,
        vip_scores: pd.DataFrame,
        glycan_type: str,
        config: DataPreparationConfig = None
    ):
        """
        Create glycopeptide comparison heatmap for a SPECIFIC glycan type
        Shows ALL glycopeptides of the specified type (complete landscape view)

        Layout (top to bottom):
        1. Top Panel: Average intensity line plots (Cancer vs Normal)
        2. Main Panel: Dot heatmap with side-by-side comparison
        3. Bottom: Glycan composition labels (rotated 45°)

        Features:
        - Filtered to single glycan type (HM, F, S, SF, or C/H)
        - Y-axis: ALL peptides with this glycan type, sorted by VIP score
        - X-axis: ALL glycan compositions of this type, naturally sorted
        - Symbols: × (Cancer), + (Normal) colored by glycan type
        - Bold symbols (thick linewidth) when qualitatively different
        - NO color bar (redundant since all glycans are same type)

        Args:
            df: Annotated DataFrame with all samples
            vip_scores: VIP scores DataFrame from PLS-DA
            glycan_type: Glycan type to filter ('HM', 'F', 'S', 'SF', or 'C/H')
            config: Data preparation configuration
        """
        logger.info(f"Creating glycopeptide comparison heatmap for {glycan_type} type...")

        # Use default config if not provided
        if config is None:
            config = DataPreparationConfig(
                min_detection_pct=0.30,
                min_samples=5,
                missing_data_method='skipna'
            )

        # Get glycan type color
        glycan_type_color = EXTENDED_CATEGORY_COLORS[glycan_type]

        # STANDARDIZED DATA PREPARATION
        df_with_vip = prepare_visualization_data(
            df=df,
            config=config,
            vip_scores=vip_scores,
            merge_method='left',
            apply_detection_filter=False,
            log_prefix=f"[{glycan_type} Heatmap] "
        )

        if len(df_with_vip) == 0:
            logger.error("No glycopeptides available!")
            return

        # Filter to specific glycan type
        df_filtered = df_with_vip[df_with_vip['GlycanTypeCategory'] == glycan_type].copy()

        if len(df_filtered) == 0:
            logger.warning(f"No glycopeptides found for type {glycan_type}")
            return

        logger.info(f"Processing {len(df_filtered)} glycopeptides of type {glycan_type}")

        # Use ALL glycopeptides of this type
        df_plot = df_filtered.copy()

        # Create peptide order (by max VIP score per peptide, descending)
        peptide_order = create_peptide_order(df_plot)

        # Create glycan order (naturally sorted, no type grouping needed)
        glycan_compositions = df_plot['GlycanComposition'].unique()
        glycan_order = sorted(glycan_compositions, key=glycan_sort_key)

        # Prepare indices
        peptide_to_idx, glycan_to_idx = create_index_mappings(peptide_order, glycan_order)

        # Dynamic figure sizing
        n_peptides = len(peptide_order)
        n_glycans = len(glycan_order)

        fig_height = max(20, n_peptides * 0.2 + 5)
        fig_width = max(20, n_glycans * 0.23 + 10)

        logger.info(f"Figure dimensions: {fig_width:.1f}\" × {fig_height:.1f}\" ({n_peptides} peptides × {n_glycans} glycans)")

        # Create figure with TWO panels (no color bar)
        fig = plt.figure(figsize=(fig_width, fig_height))
        gs = fig.add_gridspec(2, 1, height_ratios=[1, 6], hspace=0.08)

        ax_top = fig.add_subplot(gs[0])
        ax_main = fig.add_subplot(gs[1])

        # === TOP PANEL: Average intensity comparison ===
        cancer_intensities = []
        normal_intensities = []

        for glycan in glycan_order:
            glycan_data = df_plot[df_plot['GlycanComposition'] == glycan]
            cancer_int = glycan_data['Cancer_Mean'].mean()
            normal_int = glycan_data['Normal_Mean'].mean()
            cancer_intensities.append(cancer_int)
            normal_intensities.append(normal_int)

        x_pos = np.arange(len(glycan_order))
        ax_top.plot(x_pos, cancer_intensities, color=COLOR_CANCER, linewidth=PLOT_LINE_LINEWIDTH,
                    marker=MARKER_CANCER, markersize=MARKER_SIZE_SMALL, label='Cancer', alpha=ALPHA_NEAR_OPAQUE,
                    markeredgecolor=MARKER_EDGE_COLOR_LIGHT, markeredgewidth=GRID_LINEWIDTH_THIN)
        ax_top.plot(x_pos, normal_intensities, color=COLOR_NORMAL, linewidth=PLOT_LINE_LINEWIDTH,
                    marker=MARKER_NORMAL, markersize=MARKER_SIZE_SMALL, label='Normal', alpha=ALPHA_NEAR_OPAQUE,
                    markeredgecolor=MARKER_EDGE_COLOR_LIGHT, markeredgewidth=GRID_LINEWIDTH_THIN)

        ax_top.set_xlim(-0.5, len(glycan_order) - 0.5)
        ax_top.set_ylabel('Average\nIntensity', fontsize=ANNOTATION_SIZE, fontweight='bold')

        # Add vertical grid lines at glycan positions for PERFECT alignment with heatmap below
        ax_top.set_xticks(range(len(glycan_order)), minor=False)
        ax_top.set_xticklabels([])  # Hide labels but keep ticks for grid
        ax_top.grid(axis='y', alpha=ALPHA_MEDIUM, linestyle=GRID_LINESTYLE_MAJOR, linewidth=GRID_LINEWIDTH_THICK)
        ax_top.grid(axis='x', alpha=ALPHA_LIGHT, linestyle=GRID_LINESTYLE_MAJOR, linewidth=GRID_LINEWIDTH, color='#BBBBBB')

        ax_top.legend(loc='upper right', fontsize=ANNOTATION_SIZE, framealpha=FRAME_ALPHA, edgecolor=FRAME_EDGE_COLOR)
        ax_top.spines['bottom'].set_visible(False)
        ax_top.spines['top'].set_visible(False)
        ax_top.spines['right'].set_visible(False)

        # === MAIN PANEL: Symbol-based heatmap ===
        logger.info(f"Plotting symbols (all colored {glycan_type_color})...")

        # Plot symbols for each glycopeptide - ALL SAME COLOR
        for idx, row in df_plot.iterrows():
            peptide = row['Peptide']
            glycan = row['GlycanComposition']

            if peptide not in peptide_to_idx or glycan not in glycan_to_idx:
                continue

            y_pos = peptide_to_idx[peptide]
            x_pos = glycan_to_idx[glycan]

            # Check presence in each group
            cancer_intensity = row['Cancer_Mean']
            normal_intensity = row['Normal_Mean']
            has_cancer = not np.isnan(cancer_intensity) and cancer_intensity > 0
            has_normal = not np.isnan(normal_intensity) and normal_intensity > 0

            # Qualitative difference: only one group present
            is_qualitatively_different = (has_cancer and not has_normal) or (has_normal and not has_cancer)

            # Set linewidth based on qualitative difference
            if is_qualitatively_different:
                linewidth = 3.0  # Bold (NOTE: Using literal 3.0 - custom by-type emphasis)
            else:
                linewidth = EDGE_LINEWIDTH_THICK  # Normal

            # Plot Cancer marker (× symbol) - colored by glycan type
            if has_cancer:
                ax_main.scatter(x_pos, y_pos, s=SCATTER_SIZE_SMALL,
                               marker=MARKER_ANNOTATION,
                               c=glycan_type_color,  # All same color
                               linewidths=linewidth,
                               alpha=POINT_ALPHA,
                               zorder=ZORDER_DATA_HIGH,
                               label='_nolegend_')

            # Plot Normal marker (+ symbol) - colored by glycan type
            if has_normal:
                ax_main.scatter(x_pos, y_pos, s=SCATTER_SIZE_SMALL,
                               marker=MARKER_ANNOTATION_ALT,
                               c=glycan_type_color,  # All same color
                               linewidths=linewidth,
                               alpha=POINT_ALPHA,
                               zorder=ZORDER_DATA_LOW,
                               label='_nolegend_')

        # Set axis properties with equal aspect ratio
        ax_main.set_aspect('equal', adjustable='box')
        ax_main.set_xlim(-0.5, len(glycan_order) - 0.5)
        ax_main.set_ylim(-0.5, len(peptide_order) - 0.5)

        # Add x-axis glycan composition labels (rotated 45°)
        ax_main.set_xticks(range(len(glycan_order)))
        ax_main.set_xticklabels(glycan_order, rotation=45, fontsize=ANNOTATION_SIZE, ha='right')
        ax_main.set_xlabel('Glycan Composition', fontsize=AXIS_LABEL_SIZE, fontweight='bold', labelpad=10)

        # Set y-axis labels (peptides)
        ax_main.set_yticks(range(len(peptide_order)))
        ax_main.set_yticklabels(peptide_order, fontsize=ANNOTATION_SIZE)
        ax_main.set_ylabel('Peptide (sorted by VIP score)', fontsize=AXIS_LABEL_SIZE, fontweight='bold')

        # Invert y-axis to have highest VIP at top
        ax_main.invert_yaxis()

        # Grid
        ax_main.grid(True, alpha=ALPHA_MEDIUM_LIGHT, linestyle=GRID_LINESTYLE_MAJOR, linewidth=GRID_LINEWIDTH, color='#CCCCCC', zorder=ZORDER_BACKGROUND)
        ax_main.set_axisbelow(True)

        # Add minor grid
        ax_main.set_xticks([i - 0.5 for i in range(1, len(glycan_order))], minor=True)
        ax_main.set_yticks([i - 0.5 for i in range(1, len(peptide_order))], minor=True)
        ax_main.grid(which='minor', alpha=ALPHA_VERY_LIGHT, linestyle=GRID_LINESTYLE_MAJOR, linewidth=GRID_LINEWIDTH_THIN, color='#DDDDDD', zorder=ZORDER_BACKGROUND)

        # === LEGEND ===
        from matplotlib.lines import Line2D

        legend_elements = []

        # Symbol indicators with glycan type color
        legend_elements.append(Line2D([0], [0], marker=MARKER_ANNOTATION, color='w',
                                      markerfacecolor=glycan_type_color, markersize=LEGEND_MARKER_SIZE,
                                      markeredgewidth=EDGE_LINEWIDTH_THICK, markeredgecolor=glycan_type_color,
                                      label='× Cancer (both groups)'))
        legend_elements.append(Line2D([0], [0], marker=MARKER_ANNOTATION_ALT, color='w',
                                      markerfacecolor=glycan_type_color, markersize=LEGEND_MARKER_SIZE,
                                      markeredgewidth=EDGE_LINEWIDTH_THICK, markeredgecolor=glycan_type_color,
                                      label='+ Normal (both groups)'))
        legend_elements.append(Line2D([0], [0], marker=MARKER_ANNOTATION, color='w',
                                      markerfacecolor=glycan_type_color, markersize=LEGEND_MARKER_SIZE,
                                      markeredgewidth=3.0, markeredgecolor=glycan_type_color,
                                      label='× Cancer only (bold)'))
        legend_elements.append(Line2D([0], [0], marker=MARKER_ANNOTATION_ALT, color='w',
                                      markerfacecolor=glycan_type_color, markersize=LEGEND_MARKER_SIZE,
                                      markeredgewidth=3.0, markeredgecolor=glycan_type_color,
                                      label='+ Normal only (bold)'))

        ax_main.legend(handles=legend_elements, loc='upper left',
                       bbox_to_anchor=(1.01, 1), frameon=True, fontsize=LEGEND_SIZE,
                       title=f'{glycan_type} Type', title_fontsize=AXIS_LABEL_SIZE, framealpha=ALPHA_ALMOST_OPAQUE,
                       edgecolor=FRAME_EDGE_COLOR, fancybox=True, shadow=True)

        # === TITLE ===
        glycan_type_names = {
            'HM': 'High-Mannose',
            'F': 'Fucosylated',
            'S': 'Sialylated',
            'SF': 'Sialo-Fucosylated',
            'C/H': 'Complex/Hybrid'
        }
        fig.suptitle(f'Glycopeptide Comparison: {glycan_type_names[glycan_type]} Type\n'
                     f'Cancer vs Normal ({len(peptide_order)} peptides × {len(glycan_order)} glycans = {len(df_plot)} glycopeptides)',
                     fontsize=TITLE_SIZE, fontweight='bold', y=0.997)

        # Save plot with glycan-type-specific filename
        glycan_type_filename = glycan_type.replace('/', '_')  # C/H → C_H
        output_file = self.output_dir / f'glycopeptide_comparison_heatmap_{glycan_type_filename}.png'
        save_publication_figure(fig, output_file, dpi=DPI_COMPLEX)
        logger.info(f"Saved {glycan_type} heatmap to {output_file} ({DPI_COMPLEX} DPI)")

        # === SAVE TRACE DATA ===
        trace_data = df_plot.copy()

        # Add plot position information
        trace_data['Plot_X_Position'] = trace_data['GlycanComposition'].map(glycan_to_idx)
        trace_data['Plot_Y_Position'] = trace_data['Peptide'].map(peptide_to_idx)

        # Add statistical information
        cancer_samples, normal_samples = get_sample_columns(df_plot)

        cancer_trace_stats = calculate_group_statistics_standardized(trace_data, cancer_samples, method='skipna')
        normal_trace_stats = calculate_group_statistics_standardized(trace_data, normal_samples, method='skipna')

        trace_data['Cancer_StdDev'] = cancer_trace_stats['std']
        trace_data['Normal_StdDev'] = normal_trace_stats['std']
        trace_data['Cancer_SampleCount'] = cancer_trace_stats['count']
        trace_data['Normal_SampleCount'] = normal_trace_stats['count']
        trace_data['Cancer_Detection_Pct'] = cancer_trace_stats['count'] / len(cancer_samples)
        trace_data['Normal_Detection_Pct'] = normal_trace_stats['count'] / len(normal_samples)
        trace_data['Cancer_Min'] = cancer_trace_stats['min']
        trace_data['Cancer_Max'] = cancer_trace_stats['max']
        trace_data['Normal_Min'] = normal_trace_stats['min']
        trace_data['Normal_Max'] = normal_trace_stats['max']

        # Calculate fold change
        trace_data['Fold_Change'] = trace_data.apply(
            lambda row: calculate_fold_change(row['Cancer_Mean'], row['Normal_Mean'], log_scale=False),
            axis=1
        )
        trace_data['Log2_Fold_Change'] = trace_data.apply(
            lambda row: calculate_fold_change(row['Cancer_Mean'], row['Normal_Mean'], log_scale=True),
            axis=1
        )

        # Add symbol information
        # NOTE: Using literal 3.0 for linewidth_bold - custom by-type emphasis for trace data
        symbol_info = trace_data.apply(
            lambda row: get_symbol_info(row, linewidth_bold=3.0, linewidth_normal=EDGE_LINEWIDTH_THICK),
            axis=1
        )
        trace_data['Cancer_Symbol'] = [x[0] for x in symbol_info]
        trace_data['Normal_Symbol'] = [x[1] for x in symbol_info]
        trace_data['Symbol_Linewidth'] = [x[2] for x in symbol_info]
        trace_data['Group_Presence'] = [x[3] for x in symbol_info]
        trace_data['Symbol_Color'] = glycan_type_color
        trace_data['Symbol_Alpha'] = POINT_ALPHA
        trace_data['Cancer_Plotted'] = trace_data['Cancer_Symbol'].notna()
        trace_data['Normal_Plotted'] = trace_data['Normal_Symbol'].notna()

        # Reorder columns
        key_columns = [
            'Peptide', 'GlycanComposition', 'GlycanTypeCategory',
            'Plot_X_Position', 'Plot_Y_Position',
            'VIP_Score',
            'Cancer_Mean', 'Cancer_StdDev', 'Cancer_SampleCount', 'Cancer_Detection_Pct', 'Cancer_Min', 'Cancer_Max',
            'Normal_Mean', 'Normal_StdDev', 'Normal_SampleCount', 'Normal_Detection_Pct', 'Normal_Min', 'Normal_Max',
            'Fold_Change', 'Log2_Fold_Change',
            'Cancer_Symbol', 'Normal_Symbol', 'Symbol_Linewidth', 'Symbol_Color', 'Symbol_Alpha', 'Group_Presence',
            'Cancer_Plotted', 'Normal_Plotted'
        ]

        sample_columns = cancer_samples + normal_samples
        final_columns = key_columns + sample_columns

        # Save trace data
        save_trace_data(trace_data[final_columns], self.output_dir,
                        f'glycopeptide_comparison_heatmap_{glycan_type_filename}_data.csv')

        # Save summary
        summary_data = trace_data[key_columns].copy()
        save_trace_data(summary_data, self.output_dir,
                        f'glycopeptide_comparison_heatmap_{glycan_type_filename}_summary.csv')

        plt.close()

        logger.info(f"{glycan_type} heatmap complete:")
        logger.info(f"  - Peptides: {len(peptide_order)}")
        logger.info(f"  - Glycans: {len(glycan_order)}")
        logger.info(f"  - Total glycopeptides: {len(df_plot)}")
