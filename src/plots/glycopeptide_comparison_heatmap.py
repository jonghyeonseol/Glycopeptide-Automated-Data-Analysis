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

    def _plot_top_panel(
        self,
        ax,
        glycan_order: list,
        df_plot: pd.DataFrame,
        marker_size: str = 'large'
    ):
        """
        Plot average intensity comparison panel (Cancer vs Normal)

        Args:
            ax: Matplotlib axes object for top panel
            glycan_order: Ordered list of glycan compositions
            df_plot: DataFrame with plot data
            marker_size: 'large' for standard heatmap, 'small' for full-scale/by-type
        """
        # Calculate average intensities per glycan
        cancer_intensities = []
        normal_intensities = []

        for glycan in glycan_order:
            glycan_data = df_plot[df_plot['GlycanComposition'] == glycan]
            cancer_int = glycan_data['Cancer_Mean'].mean()
            normal_int = glycan_data['Normal_Mean'].mean()
            cancer_intensities.append(cancer_int)
            normal_intensities.append(normal_int)

        # Set styling based on marker size
        if marker_size == 'large':
            linewidth = PLOT_LINE_LINEWIDTH_THICK
            markersize = LINE_MARKER_SIZE
            ylabel_fontsize = AXIS_LABEL_SIZE
            legend_fontsize = ANNOTATION_SIZE
            markeredgewidth = GRID_LINEWIDTH
            y_grid_linewidth = LINE_MEDIUM_THIN
            x_grid_linewidth = GRID_LINEWIDTH_THICK
        else:  # 'small'
            linewidth = PLOT_LINE_LINEWIDTH
            markersize = MARKER_SIZE_SMALL
            ylabel_fontsize = ANNOTATION_SIZE
            legend_fontsize = ANNOTATION_SIZE
            markeredgewidth = GRID_LINEWIDTH_THIN
            y_grid_linewidth = GRID_LINEWIDTH_THICK
            x_grid_linewidth = GRID_LINEWIDTH

        # Plot lines
        x_pos = np.arange(len(glycan_order))
        ax.plot(x_pos, cancer_intensities, color=COLOR_CANCER, linewidth=linewidth,
                marker=MARKER_CANCER, markersize=markersize, label='Cancer', alpha=ALPHA_NEAR_OPAQUE,
                markeredgecolor=MARKER_EDGE_COLOR_LIGHT, markeredgewidth=markeredgewidth)
        ax.plot(x_pos, normal_intensities, color=COLOR_NORMAL, linewidth=linewidth,
                marker=MARKER_NORMAL, markersize=markersize, label='Normal', alpha=ALPHA_NEAR_OPAQUE,
                markeredgecolor=MARKER_EDGE_COLOR_LIGHT, markeredgewidth=markeredgewidth)

        # Set axis properties
        ax.set_xlim(-0.5, len(glycan_order) - 0.5)
        ax.set_ylabel('Average\nIntensity', fontsize=ylabel_fontsize, fontweight='bold')

        # Add vertical grid lines at glycan positions for PERFECT alignment with heatmap below
        ax.set_xticks(range(len(glycan_order)), minor=False)
        ax.set_xticklabels([])  # Hide labels but keep ticks for grid
        ax.grid(axis='y', alpha=ALPHA_MEDIUM, linestyle=GRID_LINESTYLE_MAJOR, linewidth=y_grid_linewidth)
        ax.grid(axis='x', alpha=ALPHA_LIGHT, linestyle=GRID_LINESTYLE_MAJOR, linewidth=x_grid_linewidth, color='#BBBBBB')

        ax.legend(loc='upper right', fontsize=legend_fontsize, framealpha=FRAME_ALPHA, edgecolor=FRAME_EDGE_COLOR)
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    def _plot_colorbar_panel(
        self,
        ax,
        glycan_order: list,
        glycan_type_positions: dict,
        glycan_type_colors: dict,
        label_size: str = 'large'
    ):
        """
        Plot glycan type color bar panel

        Args:
            ax: Matplotlib axes object for colorbar panel
            glycan_order: Ordered list of glycan compositions
            glycan_type_positions: Dict mapping glycan_type -> {'start': int, 'end': int}
            glycan_type_colors: Dict mapping glycan_type -> color
            label_size: 'large' for standard heatmap, 'small' for full-scale
        """
        ax.set_xlim(-0.5, len(glycan_order) - 0.5)
        ax.set_ylim(0, 1)

        # Set styling based on label size
        if label_size == 'large':
            edge_linewidth = EDGE_LINEWIDTH_THICK
            label_fontsize = TITLE_SIZE
            label_padding = 0.4
            bbox_linewidth = EDGE_LINEWIDTH_THICK
        else:  # 'small'
            edge_linewidth = LINE_MEDIUM
            label_fontsize = AXIS_LABEL_SIZE
            label_padding = 0.3
            bbox_linewidth = LINE_NORMAL

        # Draw solid colored rectangles for each glycan type group
        for glycan_type, pos_info in glycan_type_positions.items():
            width = pos_info['end'] - pos_info['start']
            rect = Rectangle((pos_info['start'] - 0.5, 0), width, 1,
                             facecolor=glycan_type_colors[glycan_type],
                             edgecolor=FRAME_EDGE_COLOR, linewidth=edge_linewidth, alpha=ALPHA_NEAR_OPAQUE)
            ax.add_patch(rect)

            # Add glycan type label in center
            center_x = (pos_info['start'] + pos_info['end']) / 2 - 0.5
            ax.text(center_x, 0.5, glycan_type, ha='center', va='center',
                    fontsize=label_fontsize, fontweight='bold', color='white',
                    bbox=dict(boxstyle=f'round,pad={label_padding}', facecolor='black', alpha=OVERLAY_ALPHA,
                              edgecolor=EDGE_COLOR_WHITE, linewidth=bbox_linewidth))

        # Hide axes
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)

    def _plot_symbol_heatmap(
        self,
        ax,
        df_plot: pd.DataFrame,
        peptide_to_idx: dict,
        glycan_to_idx: dict,
        peptide_order: list,
        glycan_order: list,
        marker_size: str = 'large',
        linewidth_bold: float = 5.0,
        linewidth_normal: float = None,
        glycan_type_colors: dict = None,
        fixed_color: str = None,
        glycan_type_positions: dict = None,
        separator_linewidth: float = None
    ):
        """
        Plot symbol-based heatmap with Cancer (×) and Normal (+) markers

        Args:
            ax: Matplotlib axes object for heatmap
            df_plot: DataFrame with plot data
            peptide_to_idx: Dict mapping peptide -> y position
            glycan_to_idx: Dict mapping glycan composition -> x position
            peptide_order: Ordered list of peptides
            glycan_order: Ordered list of glycan compositions
            marker_size: 'large' for standard heatmap, 'small' for full-scale/by-type
            linewidth_bold: Linewidth for qualitatively different samples
            linewidth_normal: Linewidth for samples present in both groups
            glycan_type_colors: Dict mapping glycan_type -> color (for multi-type heatmaps)
            fixed_color: Single color to use for all symbols (for type-specific heatmaps)
            glycan_type_positions: Dict of glycan type positions for vertical separators
            separator_linewidth: Linewidth for vertical separators (None = use default based on marker_size)
        """
        if linewidth_normal is None:
            linewidth_normal = PLOT_LINE_LINEWIDTH_THICK

        # Set scatter size based on marker size
        if marker_size == 'large':
            scatter_size = SCATTER_SIZE_LARGE
            default_separator_linewidth = EDGE_LINEWIDTH_THICK
        else:  # 'small'
            scatter_size = SCATTER_SIZE_SMALL
            default_separator_linewidth = LINE_NORMAL

        if separator_linewidth is None:
            separator_linewidth = default_separator_linewidth

        # Plot symbols for each glycopeptide
        for idx, row in df_plot.iterrows():
            peptide = row['Peptide']
            glycan = row['GlycanComposition']

            if peptide not in peptide_to_idx or glycan not in glycan_to_idx:
                continue

            y_pos = peptide_to_idx[peptide]
            x_pos = glycan_to_idx[glycan]

            # Determine color
            if fixed_color:
                color = fixed_color
            else:
                glycan_type = row['GlycanTypeCategory']
                color = glycan_type_colors[glycan_type]

            # Check presence in each group
            cancer_intensity = row['Cancer_Mean']
            normal_intensity = row['Normal_Mean']
            has_cancer = not np.isnan(cancer_intensity) and cancer_intensity > 0
            has_normal = not np.isnan(normal_intensity) and normal_intensity > 0

            # Qualitative difference: only one group present
            is_qualitatively_different = (has_cancer and not has_normal) or (has_normal and not has_cancer)

            # Set linewidth based on qualitative difference
            if is_qualitatively_different:
                linewidth = linewidth_bold
            else:
                linewidth = linewidth_normal

            # Plot Cancer marker (× symbol)
            if has_cancer:
                ax.scatter(x_pos, y_pos, s=scatter_size,
                          marker=MARKER_ANNOTATION,
                          c=color,
                          linewidths=linewidth,
                          alpha=POINT_ALPHA,
                          zorder=ZORDER_DATA_HIGH,
                          label='_nolegend_')

            # Plot Normal marker (+ symbol)
            if has_normal:
                ax.scatter(x_pos, y_pos, s=scatter_size,
                          marker=MARKER_ANNOTATION_ALT,
                          c=color,
                          linewidths=linewidth,
                          alpha=POINT_ALPHA,
                          zorder=ZORDER_DATA_LOW,
                          label='_nolegend_')

        # Add vertical separators between glycan types (if positions provided)
        if glycan_type_positions:
            for glycan_type, pos_info in glycan_type_positions.items():
                if pos_info['end'] < len(glycan_order):  # Don't draw after last group
                    ax.axvline(pos_info['end'] - 0.5, color='gray',
                              linestyle=THRESHOLD_LINESTYLE, linewidth=separator_linewidth,
                              alpha=ALPHA_MEDIUM, zorder=ZORDER_SEPARATOR)

        # Set axis properties with equal aspect ratio
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(-0.5, len(glycan_order) - 0.5)
        ax.set_ylim(-0.5, len(peptide_order) - 0.5)

        # Add x-axis glycan composition labels (rotated 45 degrees)
        ax.set_xticks(range(len(glycan_order)))
        if marker_size == 'large':
            ax.set_xticklabels(glycan_order, rotation=45, fontsize=TICK_LABEL_SIZE, ha='right')
        else:  # 'small'
            ax.set_xticklabels(glycan_order, rotation=45, fontsize=ANNOTATION_SIZE, ha='right')
        ax.set_xlabel('Glycan Composition', fontsize=AXIS_LABEL_SIZE, fontweight='bold', labelpad=10)

        # Set y-axis labels (peptides)
        ax.set_yticks(range(len(peptide_order)))
        if marker_size == 'large':
            ax.set_yticklabels(peptide_order, fontsize=TICK_LABEL_SIZE)
        else:  # 'small'
            ax.set_yticklabels(peptide_order, fontsize=ANNOTATION_SIZE)
        ax.set_ylabel('Peptide (sorted by VIP score)', fontsize=AXIS_LABEL_SIZE, fontweight='bold')

        # Invert y-axis to have highest VIP at top
        ax.invert_yaxis()

    def _prepare_and_save_trace_data(
        self,
        df_plot: pd.DataFrame,
        peptide_to_idx: dict,
        glycan_to_idx: dict,
        filename_base: str,
        linewidth_bold: float = 5.0,
        linewidth_normal: float = None,
        include_symbol_color: bool = False,
        symbol_color: str = None
    ):
        """
        Prepare and save comprehensive trace data with plot positions and statistics

        Args:
            df_plot: DataFrame with plot data
            peptide_to_idx: Dict mapping peptide -> y position
            glycan_to_idx: Dict mapping glycan composition -> x position
            filename_base: Base filename for output (e.g., 'glycopeptide_comparison_heatmap')
            linewidth_bold: Linewidth for qualitatively different samples
            linewidth_normal: Linewidth for samples present in both groups (defaults to PLOT_LINE_LINEWIDTH_THICK)
            include_symbol_color: If True, add Symbol_Color column
            symbol_color: Color to use if include_symbol_color=True
        """
        if linewidth_normal is None:
            linewidth_normal = PLOT_LINE_LINEWIDTH_THICK

        # Prepare enhanced trace data
        trace_data = df_plot.copy()

        # Add plot position information
        trace_data['Plot_X_Position'] = trace_data['GlycanComposition'].map(glycan_to_idx)
        trace_data['Plot_Y_Position'] = trace_data['Peptide'].map(peptide_to_idx)

        # Add statistical information
        cancer_samples, normal_samples = get_sample_columns(df_plot)

        # Calculate individual sample statistics using proper missing data handling
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
        symbol_info = trace_data.apply(
            lambda row: get_symbol_info(row, linewidth_bold=linewidth_bold, linewidth_normal=linewidth_normal),
            axis=1
        )
        trace_data['Cancer_Symbol'] = [x[0] for x in symbol_info]
        trace_data['Normal_Symbol'] = [x[1] for x in symbol_info]
        trace_data['Symbol_Linewidth'] = [x[2] for x in symbol_info]
        trace_data['Group_Presence'] = [x[3] for x in symbol_info]

        if include_symbol_color:
            trace_data['Symbol_Color'] = symbol_color

        trace_data['Symbol_Alpha'] = POINT_ALPHA
        trace_data['Cancer_Plotted'] = trace_data['Cancer_Symbol'].notna()
        trace_data['Normal_Plotted'] = trace_data['Normal_Symbol'].notna()

        # Reorder columns for readability
        key_columns = [
            'Peptide', 'GlycanComposition', 'GlycanTypeCategory',
            'Plot_X_Position', 'Plot_Y_Position',
            'VIP_Score',
            'Cancer_Mean', 'Cancer_StdDev', 'Cancer_SampleCount', 'Cancer_Detection_Pct', 'Cancer_Min', 'Cancer_Max',
            'Normal_Mean', 'Normal_StdDev', 'Normal_SampleCount', 'Normal_Detection_Pct', 'Normal_Min', 'Normal_Max',
            'Fold_Change', 'Log2_Fold_Change',
            'Cancer_Symbol', 'Normal_Symbol', 'Symbol_Linewidth'
        ]

        if include_symbol_color:
            key_columns.append('Symbol_Color')

        key_columns.extend(['Symbol_Alpha', 'Group_Presence', 'Cancer_Plotted', 'Normal_Plotted'])

        # All individual sample columns
        sample_columns = cancer_samples + normal_samples
        final_columns = key_columns + sample_columns

        # Save comprehensive trace data
        save_trace_data(trace_data[final_columns], self.output_dir, f'{filename_base}_data.csv')

        # Save summary
        summary_data = trace_data[key_columns].copy()
        save_trace_data(summary_data, self.output_dir, f'{filename_base}_summary.csv')

    def _select_top_glycopeptides(
        self,
        df_with_vip: pd.DataFrame,
        max_peptides: int,
        max_glycans_per_type: int,
        glycan_type_order: list
    ) -> pd.DataFrame:
        """
        Select top N peptides and top N glycans per type for standard heatmap variant

        Args:
            df_with_vip: DataFrame with VIP scores merged
            max_peptides: Maximum number of unique peptides to include
            max_glycans_per_type: Maximum number of glycans per type to include
            glycan_type_order: Ordered list of glycan types

        Returns:
            DataFrame with selected glycopeptides
        """
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
            return pd.DataFrame()

        return pd.concat(selected_glycans, ignore_index=False)

    def _plot_comparison_heatmap_base(
        self,
        df: pd.DataFrame,
        vip_scores: pd.DataFrame,
        heatmap_variant: str,
        config: DataPreparationConfig = None,
        figsize: tuple = None,
        max_peptides: int = 50,
        max_glycans_per_type: int = 15,
        glycan_type: str = None
    ):
        """
        Unified base method for all three heatmap variants

        Args:
            df: Annotated DataFrame with all samples
            vip_scores: VIP scores DataFrame from PLS-DA
            heatmap_variant: 'standard', 'full', or 'by_type'
            config: Data preparation configuration
            figsize: Figure size (used for 'standard' variant)
            max_peptides: Max peptides for 'standard' variant
            max_glycans_per_type: Max glycans per type for 'standard' variant
            glycan_type: Glycan type filter for 'by_type' variant
        """
        # Variant-specific logging
        variant_names = {
            'standard': 'glycopeptide comparison heatmap (Cancer vs Normal)',
            'full': 'FULL-SCALE glycopeptide comparison heatmap (ALL glycopeptides)',
            'by_type': f'glycopeptide comparison heatmap for {glycan_type} type'
        }
        logger.info(f"Creating {variant_names[heatmap_variant]}...")

        # === CONFIGURATION SETUP ===
        if config is None:
            config = DataPreparationConfig(
                min_detection_pct=0.30,
                min_samples=5,
                missing_data_method='skipna'
            )

        # Define glycan type order and colors
        glycan_type_order = ['HM', 'F', 'S', 'SF', 'C/H']
        glycan_type_colors = EXTENDED_CATEGORY_COLORS

        # Variant-specific styling parameters
        if heatmap_variant == 'standard':
            marker_size = 'large'
            linewidth_bold = 5.0
            linewidth_normal = PLOT_LINE_LINEWIDTH_THICK
            label_size = 'large'
            legend_marker_size = LEGEND_MARKER_SIZE_LARGE
            grid_alpha_major = ALPHA_MEDIUM_HIGH
            grid_color = '#BBBBBB'
        else:  # 'full' or 'by_type'
            marker_size = 'small'
            linewidth_bold = 3.0
            linewidth_normal = EDGE_LINEWIDTH_THICK
            label_size = 'small'
            legend_marker_size = LEGEND_MARKER_SIZE
            grid_alpha_major = ALPHA_MEDIUM_LIGHT
            grid_color = '#CCCCCC'

        # === DATA PREPARATION ===
        log_prefix = {
            'standard': "[Comparison Heatmap] ",
            'full': "[Full Comparison Heatmap] ",
            'by_type': f"[{glycan_type} Heatmap] "
        }[heatmap_variant]

        df_with_vip = prepare_visualization_data(
            df=df,
            config=config,
            vip_scores=vip_scores,
            merge_method='left',
            apply_detection_filter=False,
            log_prefix=log_prefix
        )

        if len(df_with_vip) == 0:
            logger.error("No glycopeptides available!")
            return

        # === DATA SELECTION (Conditional based on variant) ===
        if heatmap_variant == 'standard':
            # Select top N peptides and glycans
            df_plot = self._select_top_glycopeptides(
                df_with_vip, max_peptides, max_glycans_per_type, glycan_type_order
            )
            if len(df_plot) == 0:
                logger.warning("No glycopeptides found for comparison")
                return

        elif heatmap_variant == 'full':
            # Use ALL glycopeptides
            logger.info(f"Processing {len(df_with_vip)} total glycopeptides (complete dataset)")
            df_plot = df_with_vip.copy()

        elif heatmap_variant == 'by_type':
            # Filter to specific glycan type
            df_filtered = df_with_vip[df_with_vip['GlycanTypeCategory'] == glycan_type].copy()
            if len(df_filtered) == 0:
                logger.warning(f"No glycopeptides found for type {glycan_type}")
                return
            logger.info(f"Processing {len(df_filtered)} glycopeptides of type {glycan_type}")
            df_plot = df_filtered.copy()

        # === CREATE ORDERINGS ===
        peptide_order = create_peptide_order(df_plot)

        # Glycan ordering depends on variant
        if heatmap_variant in ['standard', 'full']:
            glycan_order, glycan_type_positions = create_glycan_order(
                df_plot, glycan_type_order, include_positions=True
            )
        else:  # 'by_type'
            glycan_compositions = df_plot['GlycanComposition'].unique()
            glycan_order = sorted(glycan_compositions, key=glycan_sort_key)
            glycan_type_positions = None

        peptide_to_idx, glycan_to_idx = create_index_mappings(peptide_order, glycan_order)

        # === FIGURE SIZING (Conditional based on variant) ===
        n_peptides = len(peptide_order)
        n_glycans = len(glycan_order)

        if heatmap_variant == 'standard':
            fig_width, fig_height = figsize
        elif heatmap_variant == 'full':
            fig_height = max(30, n_peptides * 0.2 + 5)
            fig_width = max(30, n_glycans * 0.23 + 10)
            logger.info(f"Figure dimensions: {fig_width:.1f}\" × {fig_height:.1f}\" ({n_peptides} peptides × {n_glycans} glycans)")
        elif heatmap_variant == 'by_type':
            fig_height = max(20, n_peptides * 0.2 + 5)
            fig_width = max(20, n_glycans * 0.23 + 10)
            logger.info(f"Figure dimensions: {fig_width:.1f}\" × {fig_height:.1f}\" ({n_peptides} peptides × {n_glycans} glycans)")

        # === FIGURE CREATION (Conditional panel layout) ===
        fig = plt.figure(figsize=(fig_width, fig_height))

        if heatmap_variant in ['standard', 'full']:
            # 3-panel layout with colorbar
            gs = fig.add_gridspec(3, 1, height_ratios=[1, 0.3, 5.5], hspace=0.05)
            ax_top = fig.add_subplot(gs[0])
            ax_colorbar = fig.add_subplot(gs[1])
            ax_main = fig.add_subplot(gs[2])
        else:  # 'by_type'
            # 2-panel layout without colorbar
            gs = fig.add_gridspec(2, 1, height_ratios=[1, 6], hspace=0.08)
            ax_top = fig.add_subplot(gs[0])
            ax_main = fig.add_subplot(gs[1])
            ax_colorbar = None

        # === PLOT TOP PANEL ===
        self._plot_top_panel(ax_top, glycan_order, df_plot, marker_size=marker_size)

        # === PLOT COLORBAR PANEL (if applicable) ===
        if ax_colorbar is not None:
            self._plot_colorbar_panel(
                ax_colorbar, glycan_order, glycan_type_positions,
                glycan_type_colors, label_size=label_size
            )

        # === PLOT MAIN HEATMAP ===
        if heatmap_variant == 'standard':
            logger.info("Creating symbol markers (× for Cancer, + for Normal) with bold highlighting for qualitative differences...")
        else:
            if heatmap_variant == 'by_type':
                glycan_type_color = glycan_type_colors[glycan_type]
                logger.info(f"Plotting symbols (all colored {glycan_type_color})...")
            else:
                logger.info("Plotting symbols for all glycopeptides...")

        # Determine color scheme
        if heatmap_variant == 'by_type':
            color_param = {'fixed_color': glycan_type_colors[glycan_type]}
        else:
            color_param = {'glycan_type_colors': glycan_type_colors, 'glycan_type_positions': glycan_type_positions}

        self._plot_symbol_heatmap(
            ax_main, df_plot, peptide_to_idx, glycan_to_idx,
            peptide_order, glycan_order,
            marker_size=marker_size,
            linewidth_bold=linewidth_bold,
            linewidth_normal=linewidth_normal,
            **color_param
        )

        # === GRID SETUP ===
        ax_main.grid(True, alpha=grid_alpha_major, linestyle=GRID_LINESTYLE_MAJOR,
                    linewidth=LINE_NORMAL if heatmap_variant == 'standard' else GRID_LINEWIDTH,
                    color=grid_color, zorder=ZORDER_BACKGROUND)
        ax_main.set_axisbelow(True)

        # Minor grid
        ax_main.set_xticks([i - 0.5 for i in range(1, len(glycan_order))], minor=True)
        ax_main.set_yticks([i - 0.5 for i in range(1, len(peptide_order))], minor=True)
        minor_grid_alpha = ALPHA_MEDIUM_LIGHT if heatmap_variant == 'standard' else ALPHA_VERY_LIGHT
        ax_main.grid(which='minor', alpha=minor_grid_alpha, linestyle=GRID_LINESTYLE_MAJOR,
                    linewidth=GRID_LINEWIDTH if heatmap_variant == 'standard' else GRID_LINEWIDTH_THIN,
                    color='#DDDDDD', zorder=ZORDER_BACKGROUND)

        # === LEGEND (Conditional) ===
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D

        legend_elements = []

        if heatmap_variant in ['standard', 'full']:
            # Multi-type legend with glycan type patches
            for gt in glycan_type_order:
                if gt in glycan_type_positions:
                    legend_elements.append(Patch(
                        facecolor=glycan_type_colors[gt],
                        label=f'{gt}',
                        alpha=ALPHA_NEAR_OPAQUE,
                        edgecolor=FRAME_EDGE_COLOR
                    ))

            # Symbol indicators
            legend_elements.extend([
                Line2D([0], [0], marker=MARKER_ANNOTATION, color='w',
                      markerfacecolor='gray', markersize=legend_marker_size,
                      markeredgewidth=linewidth_normal, markeredgecolor=SEPARATOR_EDGE_COLOR,
                      label='× Cancer (both groups)'),
                Line2D([0], [0], marker=MARKER_ANNOTATION_ALT, color='w',
                      markerfacecolor='gray', markersize=legend_marker_size,
                      markeredgewidth=linewidth_normal, markeredgecolor=SEPARATOR_EDGE_COLOR,
                      label='+ Normal (both groups)'),
                Line2D([0], [0], marker=MARKER_ANNOTATION, color='w',
                      markerfacecolor='gray', markersize=legend_marker_size,
                      markeredgewidth=linewidth_bold, markeredgecolor=SEPARATOR_EDGE_COLOR,
                      label='× Cancer only (bold)'),
                Line2D([0], [0], marker=MARKER_ANNOTATION_ALT, color='w',
                      markerfacecolor='gray', markersize=legend_marker_size,
                      markeredgewidth=linewidth_bold, markeredgecolor=SEPARATOR_EDGE_COLOR,
                      label='+ Normal only (bold)')
            ])

            legend_title = 'Legend'

        else:  # 'by_type'
            # Single-type legend with colored symbols
            glycan_type_color = glycan_type_colors[glycan_type]
            legend_elements.extend([
                Line2D([0], [0], marker=MARKER_ANNOTATION, color='w',
                      markerfacecolor=glycan_type_color, markersize=legend_marker_size,
                      markeredgewidth=linewidth_normal, markeredgecolor=glycan_type_color,
                      label='× Cancer (both groups)'),
                Line2D([0], [0], marker=MARKER_ANNOTATION_ALT, color='w',
                      markerfacecolor=glycan_type_color, markersize=legend_marker_size,
                      markeredgewidth=linewidth_normal, markeredgecolor=glycan_type_color,
                      label='+ Normal (both groups)'),
                Line2D([0], [0], marker=MARKER_ANNOTATION, color='w',
                      markerfacecolor=glycan_type_color, markersize=legend_marker_size,
                      markeredgewidth=linewidth_bold, markeredgecolor=glycan_type_color,
                      label='× Cancer only (bold)'),
                Line2D([0], [0], marker=MARKER_ANNOTATION_ALT, color='w',
                      markerfacecolor=glycan_type_color, markersize=legend_marker_size,
                      markeredgewidth=linewidth_bold, markeredgecolor=glycan_type_color,
                      label='+ Normal only (bold)')
            ])

            legend_title = f'{glycan_type} Type'

        ax_main.legend(handles=legend_elements, loc='upper left',
                      bbox_to_anchor=(1.01, 1), frameon=True, fontsize=LEGEND_SIZE,
                      title=legend_title, title_fontsize=AXIS_LABEL_SIZE,
                      framealpha=ALPHA_ALMOST_OPAQUE, edgecolor=FRAME_EDGE_COLOR,
                      fancybox=True, shadow=True)

        # === TITLE (Conditional) ===
        if heatmap_variant == 'standard':
            title = f'Glycopeptide Comparison Heatmap: Cancer vs Normal\n(Top {len(peptide_order)} peptides by VIP score)'
            y_pos = 0.995
        elif heatmap_variant == 'full':
            title = f'Full-Scale Glycopeptide Comparison Heatmap: Cancer vs Normal\n(ALL {len(peptide_order)} peptides × {len(glycan_order)} glycans = {len(df_plot)} glycopeptides)'
            y_pos = 0.997
        elif heatmap_variant == 'by_type':
            glycan_type_names = {
                'HM': 'High-Mannose',
                'F': 'Fucosylated',
                'S': 'Sialylated',
                'SF': 'Sialo-Fucosylated',
                'C/H': 'Complex/Hybrid'
            }
            title = f'Glycopeptide Comparison: {glycan_type_names[glycan_type]} Type\nCancer vs Normal ({len(peptide_order)} peptides × {len(glycan_order)} glycans = {len(df_plot)} glycopeptides)'
            y_pos = 0.997

        fig.suptitle(title, fontsize=TITLE_SIZE, fontweight='bold', y=y_pos)

        # === SAVE FIGURE ===
        if heatmap_variant == 'standard':
            output_file = self.output_dir / 'glycopeptide_comparison_heatmap.png'
            filename_base = 'glycopeptide_comparison_heatmap'
        elif heatmap_variant == 'full':
            output_file = self.output_dir / 'glycopeptide_comparison_heatmap_full.png'
            filename_base = 'glycopeptide_comparison_heatmap_full'
        elif heatmap_variant == 'by_type':
            glycan_type_filename = glycan_type.replace('/', '_')
            output_file = self.output_dir / f'glycopeptide_comparison_heatmap_{glycan_type_filename}.png'
            filename_base = f'glycopeptide_comparison_heatmap_{glycan_type_filename}'

        save_publication_figure(fig, output_file, dpi=DPI_COMPLEX)
        logger.info(f"Saved heatmap to {output_file} (optimized, {DPI_COMPLEX} DPI)")

        # === SAVE TRACE DATA ===
        if heatmap_variant == 'by_type':
            self._prepare_and_save_trace_data(
                df_plot, peptide_to_idx, glycan_to_idx,
                filename_base,
                linewidth_bold=linewidth_bold,
                linewidth_normal=linewidth_normal,
                include_symbol_color=True,
                symbol_color=glycan_type_colors[glycan_type]
            )
        else:
            self._prepare_and_save_trace_data(
                df_plot, peptide_to_idx, glycan_to_idx,
                filename_base,
                linewidth_bold=linewidth_bold,
                linewidth_normal=linewidth_normal
            )

        plt.close()

        # === LOGGING ===
        cancer_samples, normal_samples = get_sample_columns(df_plot)

        if heatmap_variant == 'standard':
            logger.info("Comparison heatmap complete:")
            logger.info(f"  - Cancer samples: {len(cancer_samples)}")
            logger.info(f"  - Normal samples: {len(normal_samples)}")
            logger.info(f"  - Peptides shown: {len(peptide_order)}")
            logger.info(f"  - Glycans shown: {len(glycan_order)}")
            for gt, pos in glycan_type_positions.items():
                logger.info(f"    - {gt}: {pos['end'] - pos['start']} glycans")

        elif heatmap_variant == 'full':
            logger.info("Full-scale comparison heatmap complete:")
            logger.info(f"  - Cancer samples: {len(cancer_samples)}")
            logger.info(f"  - Normal samples: {len(normal_samples)}")
            logger.info(f"  - Peptides shown: {len(peptide_order)} (ALL unique peptides)")
            logger.info(f"  - Glycans shown: {len(glycan_order)} (ALL unique glycans)")
            logger.info(f"  - Total glycopeptides: {len(df_plot)}")
            for gt, pos in glycan_type_positions.items():
                logger.info(f"    - {gt}: {pos['end'] - pos['start']} glycans")

        elif heatmap_variant == 'by_type':
            logger.info(f"{glycan_type} heatmap complete:")
            logger.info(f"  - Peptides: {len(peptide_order)}")
            logger.info(f"  - Glycans: {len(glycan_order)}")
            logger.info(f"  - Total glycopeptides: {len(df_plot)}")

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
        4. Bottom: Glycan composition labels (rotated 90 degrees vertically)

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
        return self._plot_comparison_heatmap_base(
            df=df,
            vip_scores=vip_scores,
            heatmap_variant='standard',
            config=config,
            figsize=figsize,
            max_peptides=max_peptides,
            max_glycans_per_type=max_glycans_per_type
        )

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
        4. Bottom: Glycan composition labels (rotated 45 degrees)

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
        return self._plot_comparison_heatmap_base(
            df=df,
            vip_scores=vip_scores,
            heatmap_variant='full',
            config=config
        )

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
        3. Bottom: Glycan composition labels (rotated 45 degrees)

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
        return self._plot_comparison_heatmap_base(
            df=df,
            vip_scores=vip_scores,
            heatmap_variant='by_type',
            config=config,
            glycan_type=glycan_type
        )
