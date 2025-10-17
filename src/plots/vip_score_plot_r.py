"""
VIP Score Plot Module using Seaborn for pGlyco Auto Combine
Handles VIP score visualizations with Python graphics

UPDATED: Switched from R/ggplot2 to Python/Seaborn for better design control
"""

import pandas as pd
import logging
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
from typing import Optional
from ..utils import get_sample_columns, save_trace_data
from ..data_preparation import (
    DataPreparationConfig,
    calculate_group_statistics_standardized
)
from .plot_config import (
    # Font sizes
    VIP_FEATURE_NAME_SIZE,
    VIP_SUBTITLE_FONTSIZE,
    VIP_LEGEND_FONTSIZE,
    AXIS_LABEL_SIZE,
    TITLE_SIZE,
    # Figure dimensions
    VIP_FIGURE_WIDTH,
    VIP_FIGURE_HEIGHT,
    VIP_GROUPED_HEIGHT_MULTIPLIER,
    # VIP styling
    VIP_SIGNIFICANCE_THRESHOLD,
    VIP_THRESHOLD_LINE_COLOR,
    VIP_THRESHOLD_LINE_WIDTH_SEABORN,
    VIP_THRESHOLD_LINE_ALPHA,
    VIP_THRESHOLD_TEXT_Y_OFFSET,
    VIP_THRESHOLD_TEXT_Y_OFFSET_GROUPED,
    # Bar styling
    VIP_BAR_HEIGHT,
    VIP_BAR_HEIGHT_GROUPED,
    VIP_BAR_ALPHA,
    VIP_BAR_EDGE_LINEWIDTH,
    # Grid styling
    VIP_GRID_ALPHA,
    VIP_GRID_LINEWIDTH,
    # Layout
    VIP_XLIM_EXPANSION,
    VIP_GROUP_SPACING,
    VIP_SPINE_LINEWIDTH,
    # Colors
    COLOR_CANCER,
    COLOR_NORMAL,
    COLOR_GRAY,
    # DPI and saving
    DPI_SUPPLEMENTARY,
    save_publication_figure
)

logger = logging.getLogger(__name__)

# Message templates (Phase 2: Eliminate magic string duplication)
VIP_SUBTITLE_TEXT = 'Bar color indicates higher intensity group (Red: Cancer, Blue: Normal)'
VIP_THRESHOLD_LABEL = 'VIP = 1.0'
VIP_AXIS_LABEL = 'VIP Score'

# Seaborn style configuration
sns.set_style("whitegrid", {
    'grid.linestyle': '-',
    'grid.linewidth': 0.3,
    'grid.color': '#E5E5E5',
    'axes.edgecolor': 'black',
    'axes.linewidth': 0.8
})


class VIPDataStrategy:
    """
    Strategy patterns for VIP data filtering and aggregation

    Extracts nested functions from plot methods to improve testability
    and eliminate the nested function anti-pattern.

    Phase 2 Refactoring: Consolidates 7 nested functions across 3 methods.
    """

    @staticmethod
    def mask_glycopeptide(df: pd.DataFrame, peptide: str, glycan: str) -> pd.Series:
        """
        Filter DataFrame by exact peptide + glycan combination

        Args:
            df: DataFrame to filter
            peptide: Peptide sequence
            glycan: Glycan composition string

        Returns:
            Boolean mask for matching rows
        """
        return (df['Peptide'] == peptide) & (df['GlycanComposition'] == glycan)

    @staticmethod
    def mask_glycan_composition(df: pd.DataFrame, glycan: str) -> pd.Series:
        """
        Filter DataFrame by glycan composition

        Args:
            df: DataFrame to filter
            glycan: Glycan composition string

        Returns:
            Boolean mask for matching rows
        """
        return df['GlycanComposition'] == glycan

    @staticmethod
    def mask_peptide(df: pd.DataFrame, peptide: str) -> pd.Series:
        """
        Filter DataFrame by peptide sequence

        Args:
            df: DataFrame to filter
            peptide: Peptide sequence

        Returns:
            Boolean mask for matching rows
        """
        return df['Peptide'] == peptide

    @staticmethod
    def aggregate_mean(stats: dict) -> float:
        """
        Extract mean value from statistics dictionary

        Used for individual glycopeptides where we want the mean intensity.

        Args:
            stats: Statistics dictionary from calculate_group_statistics_standardized

        Returns:
            Mean value, or 0 if all NaN
        """
        return stats['mean'].iloc[0] if not stats['mean'].isna().all() else 0

    @staticmethod
    def aggregate_sum(stats: dict) -> float:
        """
        Sum all values from statistics dictionary

        Used for aggregating across multiple glycoforms or peptides.

        Args:
            stats: Statistics dictionary from calculate_group_statistics_standardized

        Returns:
            Sum of all values
        """
        return stats['sum'].sum()


class VIPScorePlotRMixin:
    """
    Mixin class for Seaborn-based VIP score plots

    NOTE: Method names retain '_r' suffix for backwards compatibility with main.py.
    This is a legacy artifact from the R/ggplot2 implementation. The code is now
    pure Python/Seaborn, but changing method names would be a breaking API change.

    Future consideration (v4.0.0): Remove '_r' suffix in a major version update.
    """

    @staticmethod
    def _assign_bar_colors_vectorized(fold_change_df: pd.DataFrame) -> pd.Series:
        """
        Vectorized color assignment based on fold change direction

        Performance: ~3-5x faster than loop-based approach

        Args:
            fold_change_df: DataFrame with 'Cancer' and 'Normal' columns

        Returns:
            Series of color codes (Red=Cancer enriched, Blue=Normal enriched, Gray=Equal)
        """
        colors = pd.Series(COLOR_GRAY, index=fold_change_df.index)  # Default gray

        cancer_higher = fold_change_df['Cancer'] > fold_change_df['Normal']
        normal_higher = fold_change_df['Normal'] > fold_change_df['Cancer']

        colors[cancer_higher] = COLOR_CANCER  # Red
        colors[normal_higher] = COLOR_NORMAL  # Blue

        return colors

    @staticmethod
    def _add_vip_legend(ax, fontsize: int = None):
        """
        Add standard Cancer/Normal enrichment legend

        Args:
            ax: Matplotlib axes object
            fontsize: Font size for legend (defaults to VIP_LEGEND_FONTSIZE)
        """
        if fontsize is None:
            fontsize = VIP_LEGEND_FONTSIZE

        cancer_patch = mpatches.Patch(color=COLOR_CANCER, label='Cancer enriched',
                                      alpha=VIP_BAR_ALPHA)
        normal_patch = mpatches.Patch(color=COLOR_NORMAL, label='Normal enriched',
                                      alpha=VIP_BAR_ALPHA)
        ax.legend(handles=[cancer_patch, normal_patch], loc='lower right',
                  frameon=True, fancybox=False, shadow=False, fontsize=fontsize)

    @staticmethod
    def _validate_vip_input(vip_df: pd.DataFrame, required_cols: list, top_n: int, plot_type: str) -> Optional[int]:
        """
        Centralized input validation for VIP plot methods

        Eliminates code duplication across all plot_vip_scores_* methods.

        Args:
            vip_df: VIP scores DataFrame to validate
            required_cols: List of required column names
            top_n: Number of top features to show
            plot_type: Type of plot for logging (e.g., 'glycopeptide', 'peptide')

        Returns:
            Adjusted top_n value, or None if validation fails (empty DataFrame)

        Raises:
            ValueError: If required columns are missing or top_n is invalid
        """
        # Check for empty DataFrame
        if vip_df.empty:
            logger.warning(f"VIP DataFrame is empty, skipping {plot_type} VIP plot")
            return None

        # Validate required columns
        if not all(col in vip_df.columns for col in required_cols):
            raise ValueError(f"VIP DataFrame missing required columns: {required_cols}")

        # Validate top_n parameter
        if top_n < 1:
            raise ValueError(f"top_n must be >= 1, got {top_n}")

        # Adjust top_n if it exceeds available features
        if top_n > len(vip_df):
            logger.warning(f"top_n ({top_n}) exceeds available features ({len(vip_df)}), showing all")
            return len(vip_df)

        return top_n

    @staticmethod
    def _calculate_fold_change_color(cancer_mean: float, normal_mean: float) -> str:
        """
        Calculate bar color based on fold change direction

        Helper method to reduce code duplication in grouped plotting.

        Args:
            cancer_mean: Mean intensity in cancer group
            normal_mean: Mean intensity in normal group

        Returns:
            Color code (COLOR_CANCER, COLOR_NORMAL, or COLOR_GRAY)
        """
        if cancer_mean > normal_mean:
            return COLOR_CANCER  # Red = Cancer enriched
        elif normal_mean > cancer_mean:
            return COLOR_NORMAL  # Blue = Normal enriched
        else:
            return COLOR_GRAY  # Gray = Equal

    def _get_vip_output_path(self, plot_type: str, extension: str = 'png'):
        """
        Generate standardized output file path for VIP plots

        Phase 3 Refactoring: Centralizes file path generation.

        Args:
            plot_type: Type of VIP plot ('glycopeptide', 'glycan_composition', 'peptide', 'peptide_grouped')
            extension: File extension (default: 'png')

        Returns:
            Path object for output file
        """
        filename = f'vip_score_{plot_type}_r.{extension}'
        return self.output_dir / filename

    def _get_vip_trace_filename(self, plot_type: str) -> str:
        """
        Generate standardized trace data filename for VIP plots

        Phase 3 Refactoring: Centralizes trace data file naming.

        Args:
            plot_type: Type of VIP plot ('glycopeptide', 'glycan_composition', 'peptide', 'peptide_grouped')

        Returns:
            Filename string for trace data CSV
        """
        return f'vip_score_{plot_type}_data.csv'

    def _build_glycoform_plot_rows(
        self,
        df: pd.DataFrame,
        vip_df: pd.DataFrame,
        top_peptides: list,
        cancer_samples: list,
        normal_samples: list
    ) -> list:
        """
        Build plot row data for grouped peptide VIP visualization

        Extracts complex glycoform processing logic from plot_vip_scores_peptide_grouped_r
        to improve maintainability and testability.

        Args:
            df: Annotated DataFrame with intensity data
            vip_df: DataFrame with VIP scores
            top_peptides: List of top peptide names
            cancer_samples: List of cancer sample column names
            normal_samples: List of normal sample column names

        Returns:
            List of dictionaries containing plot row data
        """
        plot_rows = []
        y_pos = 0
        config = DataPreparationConfig(missing_data_method='skipna')

        for peptide in reversed(top_peptides):  # Reverse for top-to-bottom plotting
            # Get all glycoforms for this peptide, sorted by VIP score descending
            peptide_glycoforms = vip_df[vip_df['Peptide'] == peptide].sort_values(
                'VIP_Score', ascending=False).reset_index(drop=True)

            # Process each glycoform
            for idx, row in peptide_glycoforms.iterrows():
                glycan_comp = row['GlycanComposition']
                vip_score = row['VIP_Score']

                # Calculate statistics using standardized methods
                mask = (df['Peptide'] == peptide) & (df['GlycanComposition'] == glycan_comp)
                if mask.sum() > 0:
                    glycopeptide_row = df[mask]

                    cancer_stats = calculate_group_statistics_standardized(
                        glycopeptide_row, cancer_samples, method=config.missing_data_method
                    )
                    normal_stats = calculate_group_statistics_standardized(
                        glycopeptide_row, normal_samples, method=config.missing_data_method
                    )

                    cancer_mean = cancer_stats['mean'].iloc[0] if not cancer_stats['mean'].isna().all() else 0
                    normal_mean = normal_stats['mean'].iloc[0] if not normal_stats['mean'].isna().all() else 0
                else:
                    cancer_mean = 0
                    normal_mean = 0

                # Determine bar color using helper method
                bar_color = self._calculate_fold_change_color(cancer_mean, normal_mean)

                plot_rows.append({
                    'Peptide': peptide,
                    'GlycanComposition': glycan_comp,
                    'VIP_Score': vip_score,
                    'y_pos': y_pos,
                    'BarColor': bar_color,
                    'is_first': idx == 0,
                    'is_last': idx == len(peptide_glycoforms) - 1,
                    'group_size': len(peptide_glycoforms)
                })
                y_pos += 1

            # Add spacing between peptide groups
            y_pos += VIP_GROUP_SPACING

        return plot_rows

    def _create_vip_plot_seaborn(self, vip_data: pd.DataFrame, fold_change_direction: pd.DataFrame,
                                  title: str, output_file: str) -> None:
        """
        Create clean VIP score plot using Python/Seaborn (Publication-quality design)

        Args:
            vip_data: DataFrame with columns [Feature, VIP_Score]
            fold_change_direction: DataFrame with columns [Feature, Cancer, Normal] for color determination
            title: Plot title
            output_file: Output file path
        """
        # Prepare plot data with fold change direction for color coding
        plot_data = vip_data.copy()

        # Vectorized color assignment (3-5x faster than loop)
        plot_data['BarColor'] = self._assign_bar_colors_vectorized(fold_change_direction)

        # Sort by VIP score for top-to-bottom display
        plot_data = plot_data.sort_values('VIP_Score', ascending=True).reset_index(drop=True)

        # Create figure
        fig, ax = plt.subplots(figsize=(VIP_FIGURE_WIDTH, VIP_FIGURE_HEIGHT))

        # Plot horizontal bars
        y_positions = np.arange(len(plot_data))
        bars = ax.barh(y_positions, plot_data['VIP_Score'],
                       height=VIP_BAR_HEIGHT, color=plot_data['BarColor'],
                       alpha=VIP_BAR_ALPHA, edgecolor='white', linewidth=VIP_BAR_EDGE_LINEWIDTH)

        # VIP significance threshold line
        ax.axvline(VIP_SIGNIFICANCE_THRESHOLD, color=VIP_THRESHOLD_LINE_COLOR,
                   linestyle='--', linewidth=VIP_THRESHOLD_LINE_WIDTH_SEABORN,
                   alpha=VIP_THRESHOLD_LINE_ALPHA, zorder=0)
        ax.text(VIP_SIGNIFICANCE_THRESHOLD, VIP_THRESHOLD_TEXT_Y_OFFSET, VIP_THRESHOLD_LABEL,
                ha='left', va='top', fontsize=9, color=VIP_THRESHOLD_LINE_COLOR,
                style='italic')

        # Set labels
        ax.set_yticks(y_positions)
        ax.set_yticklabels(plot_data['Feature'], fontsize=VIP_FEATURE_NAME_SIZE)
        ax.set_xlabel(VIP_AXIS_LABEL, fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax.set_title(title, fontsize=TITLE_SIZE, fontweight='bold', pad=15)

        # Add subtitle
        ax.text(0.5, 1.02, VIP_SUBTITLE_TEXT,
                ha='center', va='bottom', transform=ax.transAxes,
                fontsize=VIP_SUBTITLE_FONTSIZE, color='#666666', style='italic')

        # Clean up spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(VIP_SPINE_LINEWIDTH)
        ax.spines['bottom'].set_linewidth(VIP_SPINE_LINEWIDTH)

        # Grid styling (horizontal only)
        ax.grid(axis='x', alpha=VIP_GRID_ALPHA, linestyle='-', linewidth=VIP_GRID_LINEWIDTH)
        ax.set_axisbelow(True)

        # Add legend using helper method
        self._add_vip_legend(ax)

        # Set x-axis limits
        max_vip = plot_data['VIP_Score'].max()
        ax.set_xlim(0, max_vip * VIP_XLIM_EXPANSION)

        # Tight layout
        plt.tight_layout()

        # Save figure using standardized function
        save_publication_figure(fig, output_file, dpi=DPI_SUPPLEMENTARY)
        plt.close(fig)

        logger.info(f"Saved Seaborn-based VIP score plot to {output_file} (optimized, {DPI_SUPPLEMENTARY} DPI)")

    def _prepare_vip_heatmap_data_generic(
        self,
        df: pd.DataFrame,
        features_df: pd.DataFrame,
        mask_fn,
        aggregation_fn
    ) -> tuple:
        """
        Generic heatmap data preparation for VIP score plots.

        This helper consolidates the common heatmap data preparation pattern used by:
        - plot_vip_scores_glycopeptide_r()
        - plot_vip_scores_glycan_composition_r()
        - plot_vip_scores_peptide_r()

        Uses Strategy Pattern for flexible filtering and aggregation.

        Args:
            df: Annotated DataFrame
            features_df: DataFrame with 'Feature' and 'VIP_Score' columns
            mask_fn: Callable(df, row) -> boolean mask for filtering data
            aggregation_fn: Callable(stats_dict, row) -> aggregated value

        Returns:
            Tuple of (heatmap_df, vip_plot_data)
        """
        # Get sample columns and config
        cancer_samples, normal_samples = get_sample_columns(df)
        config = DataPreparationConfig(missing_data_method='skipna')

        # Prepare heatmap data using provided strategies
        heatmap_data = []

        for _, row in features_df.iterrows():
            mask = mask_fn(df, row)

            if mask.sum() > 0:
                filtered_rows = df[mask]

                # Use standardized statistics calculation
                cancer_stats = calculate_group_statistics_standardized(
                    filtered_rows, cancer_samples, method=config.missing_data_method
                )
                normal_stats = calculate_group_statistics_standardized(
                    filtered_rows, normal_samples, method=config.missing_data_method
                )

                # Apply aggregation strategy
                cancer_value = aggregation_fn(cancer_stats, row)
                normal_value = aggregation_fn(normal_stats, row)

                heatmap_data.append({
                    'Feature': row['Feature'],
                    'Cancer': cancer_value,
                    'Normal': normal_value
                })
            else:
                heatmap_data.append({
                    'Feature': row['Feature'],
                    'Cancer': 0,
                    'Normal': 0
                })

        heatmap_df = pd.DataFrame(heatmap_data)
        vip_plot_data = features_df[['Feature', 'VIP_Score']].copy()

        return heatmap_df, vip_plot_data

    def plot_vip_scores_glycopeptide_r(self, df: pd.DataFrame, vip_df: pd.DataFrame, top_n: int = 10) -> None:
        """
        Plot top VIP scores by glycopeptide using Python/Seaborn

        Args:
            df: Annotated DataFrame
            vip_df: DataFrame with VIP scores by glycopeptide
            top_n: Number of top glycopeptides to show (must be >= 1)

        Raises:
            ValueError: If vip_df is empty or missing required columns
        """
        try:
            # Input validation using centralized method
            required_cols = ['Peptide', 'GlycanComposition', 'VIP_Score']
            top_n = self._validate_vip_input(vip_df, required_cols, top_n, 'glycopeptide')
            if top_n is None:
                return

            top_n_data = vip_df.head(top_n).copy()

            # Create highly readable labels
            def format_feature_label(peptide, glycan):
                """Create clear, readable labels - NO truncation for maximum clarity"""
                # Format: "PEPTIDE | H(5)N(4)A(1)"
                return f"{peptide} | {glycan}"

            top_n_data['Feature'] = top_n_data.apply(
                lambda row: format_feature_label(row['Peptide'], row['GlycanComposition']), axis=1
            )

            # Prepare strategies using VIPDataStrategy class
            strategy = VIPDataStrategy()

            def mask_fn(df_inner, row):
                return strategy.mask_glycopeptide(df_inner, row['Peptide'], row['GlycanComposition'])

            def agg_fn(stats, row):
                return strategy.aggregate_mean(stats)

            # Prepare heatmap data using unified helper with strategy pattern
            heatmap_df, vip_plot_data = self._prepare_vip_heatmap_data_generic(
                df=df,
                features_df=top_n_data,
                mask_fn=mask_fn,
                aggregation_fn=agg_fn
            )

            # Generate output paths using centralized helper
            output_file = self._get_vip_output_path('glycopeptide')
            trace_filename = self._get_vip_trace_filename('glycopeptide')

            # Note: heatmap_df now used only for fold change direction (bar coloring)
            self._create_vip_plot_seaborn(vip_plot_data, heatmap_df,
                                           f'Top {top_n} Glycopeptides by VIP Score',
                                           str(output_file))

            # Save trace data for reproducibility
            save_trace_data(vip_plot_data, self.output_dir, trace_filename)

        except Exception as e:
            logger.error(f"Failed to create glycopeptide VIP plot: {e}")
            raise

    def plot_vip_scores_glycan_composition_r(self, df: pd.DataFrame, vip_df: pd.DataFrame, top_n: int = 10) -> None:
        """
        Plot VIP scores by GlycanComposition using Python/Seaborn

        Args:
            df: Annotated DataFrame
            vip_df: DataFrame with all VIP scores by glycopeptide
            top_n: Number of top GlycanCompositions to show (must be >= 1)

        Raises:
            ValueError: If vip_df is empty or missing required columns
        """
        try:
            # Input validation using centralized method
            required_cols = ['GlycanComposition', 'VIP_Score']
            top_n = self._validate_vip_input(vip_df, required_cols, top_n, 'glycan composition')
            if top_n is None:
                return

            # Group by GlycanComposition and get the max VIP score for each
            glycan_vip = vip_df.groupby('GlycanComposition')['VIP_Score'].max().nlargest(top_n).reset_index()
            glycan_vip['Feature'] = glycan_vip['GlycanComposition']

            # Prepare strategies using VIPDataStrategy class
            strategy = VIPDataStrategy()

            def mask_fn(df_inner, row):
                return strategy.mask_glycan_composition(df_inner, row['GlycanComposition'])

            def agg_fn(stats, row):
                return strategy.aggregate_sum(stats)

            # Prepare heatmap data using unified helper with strategy pattern
            heatmap_df, vip_plot_data = self._prepare_vip_heatmap_data_generic(
                df=df,
                features_df=glycan_vip,
                mask_fn=mask_fn,
                aggregation_fn=agg_fn
            )

            # Generate output paths using centralized helper
            output_file = self._get_vip_output_path('glycan_composition')
            trace_filename = self._get_vip_trace_filename('glycan_composition')

            # Note: heatmap_df now used only for fold change direction (bar coloring)
            self._create_vip_plot_seaborn(vip_plot_data, heatmap_df,
                                           f'Top {top_n} Glycan Compositions by VIP Score',
                                           str(output_file))

            # Save trace data for reproducibility
            save_trace_data(vip_plot_data, self.output_dir, trace_filename)

        except Exception as e:
            logger.error(f"Failed to create glycan composition VIP plot: {e}")
            raise

    def plot_vip_scores_peptide_r(self, df: pd.DataFrame, vip_df: pd.DataFrame, top_n: int = 10) -> None:
        """
        Plot VIP scores by Peptide using Python/Seaborn

        Args:
            df: Annotated DataFrame
            vip_df: DataFrame with all VIP scores
            top_n: Number of top peptides to show (must be >= 1)

        Raises:
            ValueError: If vip_df is empty or missing required columns
        """
        try:
            # Input validation using centralized method
            required_cols = ['Peptide', 'VIP_Score']
            top_n = self._validate_vip_input(vip_df, required_cols, top_n, 'peptide')
            if top_n is None:
                return

            # Group by Peptide and get the max VIP score for each
            peptide_vip = vip_df.groupby('Peptide')['VIP_Score'].max().nlargest(top_n).reset_index()
            peptide_vip['Feature'] = peptide_vip['Peptide']

            # Prepare strategies using VIPDataStrategy class
            strategy = VIPDataStrategy()

            def mask_fn(df_inner, row):
                return strategy.mask_peptide(df_inner, row['Peptide'])

            def agg_fn(stats, row):
                return strategy.aggregate_sum(stats)

            # Prepare heatmap data using unified helper with strategy pattern
            heatmap_df, vip_plot_data = self._prepare_vip_heatmap_data_generic(
                df=df,
                features_df=peptide_vip,
                mask_fn=mask_fn,
                aggregation_fn=agg_fn
            )

            # Generate output paths using centralized helper
            output_file = self._get_vip_output_path('peptide')
            trace_filename = self._get_vip_trace_filename('peptide')

            # Note: heatmap_df now used only for fold change direction (bar coloring)
            self._create_vip_plot_seaborn(vip_plot_data, heatmap_df,
                                           f'Top {top_n} Peptides by VIP Score',
                                           str(output_file))

            # Save trace data for reproducibility
            save_trace_data(vip_plot_data, self.output_dir, trace_filename)

        except Exception as e:
            logger.error(f"Failed to create peptide VIP plot: {e}")
            raise

    def plot_vip_scores_peptide_grouped_r(self, df: pd.DataFrame, vip_df: pd.DataFrame, top_n: int = 10) -> None:
        """
        Plot VIP scores with peptide grouping (bracket notation for multiple glycoforms) using Python/Seaborn

        Args:
            df: Annotated DataFrame
            vip_df: DataFrame with all VIP scores by glycopeptide
            top_n: Number of top peptides to show (must be >= 1)

        Raises:
            ValueError: If vip_df is empty or missing required columns
        """
        try:
            # Input validation using centralized method
            required_cols = ['Peptide', 'GlycanComposition', 'VIP_Score']
            top_n = self._validate_vip_input(vip_df, required_cols, top_n, 'grouped peptide')
            if top_n is None:
                return

            # Get top peptides by max VIP score
            top_peptides = vip_df.groupby('Peptide')['VIP_Score'].max().nlargest(top_n).index.tolist()

            # Get sample columns (C1-C24, N1-N24)
            cancer_samples, normal_samples = get_sample_columns(df)

            # Build plot data using extracted helper method
            plot_rows = self._build_glycoform_plot_rows(
                df, vip_df, top_peptides, cancer_samples, normal_samples
            )
            plot_data = pd.DataFrame(plot_rows)

            # Generate output paths using centralized helper
            output_file = self._get_vip_output_path('peptide_grouped')
            trace_filename = self._get_vip_trace_filename('peptide_grouped')

            # Create Seaborn plot - unified design with other VIP plots
            self._create_grouped_vip_plot_seaborn(plot_data, output_file, top_n)

            # Save trace data for reproducibility
            save_trace_data(plot_data[['Peptide', 'GlycanComposition', 'VIP_Score', 'BarColor']],
                            self.output_dir, trace_filename)

        except Exception as e:
            logger.error(f"Failed to create grouped peptide VIP plot: {e}")
            raise

    def _create_grouped_vip_plot_seaborn(self, plot_data: pd.DataFrame, output_file: str, top_n: int) -> None:
        """
        Create grouped VIP plot showing glycoforms per peptide using Python/Seaborn

        Args:
            plot_data: DataFrame with peptide/glycan/VIP/color data
            output_file: Output file path
            top_n: Number of top peptides shown
        """
        # Sort by VIP score for bottom-to-top display
        plot_data = plot_data.sort_values('VIP_Score', ascending=True).reset_index(drop=True)

        # Create labels combining peptide and glycan for grouped entries
        labels = []
        for _, row in plot_data.iterrows():
            if row['group_size'] > 1:
                labels.append(f"   {row['GlycanComposition']}")  # Indent glycans
            else:
                labels.append(row['Peptide'])
        plot_data['DisplayLabel'] = labels

        # Create figure (using height multiplier constant)
        fig, ax = plt.subplots(figsize=(VIP_FIGURE_WIDTH, VIP_FIGURE_HEIGHT * VIP_GROUPED_HEIGHT_MULTIPLIER))

        # Plot horizontal bars (using constants)
        y_positions = np.arange(len(plot_data))
        bars = ax.barh(y_positions, plot_data['VIP_Score'],
                       height=VIP_BAR_HEIGHT_GROUPED, color=plot_data['BarColor'],
                       alpha=VIP_BAR_ALPHA, edgecolor='white', linewidth=VIP_BAR_EDGE_LINEWIDTH)

        # VIP significance threshold line (using constants)
        ax.axvline(VIP_SIGNIFICANCE_THRESHOLD, color=VIP_THRESHOLD_LINE_COLOR,
                   linestyle='--', linewidth=VIP_THRESHOLD_LINE_WIDTH_SEABORN,
                   alpha=VIP_THRESHOLD_LINE_ALPHA, zorder=0)
        ax.text(VIP_SIGNIFICANCE_THRESHOLD, VIP_THRESHOLD_TEXT_Y_OFFSET_GROUPED, VIP_THRESHOLD_LABEL,
                ha='left', va='top', fontsize=9, color=VIP_THRESHOLD_LINE_COLOR,
                style='italic')

        # Set labels (using constants)
        ax.set_yticks(y_positions)
        ax.set_yticklabels(plot_data['DisplayLabel'], fontsize=VIP_FEATURE_NAME_SIZE * 0.9)

        # Add peptide names for grouped entries (bold)
        for peptide in plot_data['Peptide'].unique():
            peptide_data = plot_data[plot_data['Peptide'] == peptide]
            if len(peptide_data) > 1:
                # Add peptide label on left for grouped glycoforms
                first_idx = peptide_data.index[0]
                y_pos = y_positions[plot_data.index.get_loc(first_idx)]
                ax.text(-0.15, y_pos, peptide,
                        ha='right', va='center', fontsize=VIP_FEATURE_NAME_SIZE,
                        fontweight='bold', transform=ax.get_yaxis_transform())

        ax.set_xlabel(VIP_AXIS_LABEL, fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax.set_title(f'Top {top_n} Peptides by VIP Score (Grouped by Glycoforms)',
                     fontsize=TITLE_SIZE, fontweight='bold', pad=15)

        # Add subtitle (using constant)
        ax.text(0.5, 1.02, VIP_SUBTITLE_TEXT,
                ha='center', va='bottom', transform=ax.transAxes,
                fontsize=VIP_SUBTITLE_FONTSIZE, color='#666666', style='italic')

        # Clean up spines (using constant)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(VIP_SPINE_LINEWIDTH)
        ax.spines['bottom'].set_linewidth(VIP_SPINE_LINEWIDTH)

        # Grid styling (using constants)
        ax.grid(axis='x', alpha=VIP_GRID_ALPHA, linestyle='-', linewidth=VIP_GRID_LINEWIDTH)
        ax.set_axisbelow(True)

        # Add legend using helper method
        self._add_vip_legend(ax)

        # Set x-axis limits (using constant)
        max_vip = plot_data['VIP_Score'].max()
        ax.set_xlim(0, max_vip * VIP_XLIM_EXPANSION)

        # Tight layout
        plt.tight_layout()

        # Save figure using standardized function
        save_publication_figure(fig, output_file, dpi=DPI_SUPPLEMENTARY)
        plt.close(fig)

        logger.info(f"Saved Seaborn-based grouped VIP plot to {output_file} (optimized, {DPI_SUPPLEMENTARY} DPI)")


