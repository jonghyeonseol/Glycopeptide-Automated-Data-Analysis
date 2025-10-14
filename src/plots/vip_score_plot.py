"""
VIP Score Plot Module for pGlyco Auto Combine
Handles VIP score visualizations with heatmap

Dependencies:
    External:
        - pandas: Data manipulation
        - matplotlib: Plotting backend
        - seaborn: Statistical visualization (heatmap)

    Internal:
        - src.utils: get_sample_columns
        - src.data_preparation: DataPreparationConfig, calculate_group_statistics_standardized
        - src.plots.plot_config: save_publication_figure

UPDATED: Now uses centralized data preparation for consistency
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from matplotlib.gridspec import GridSpec
from ..utils import get_sample_columns
from ..data_preparation import (
    DataPreparationConfig,
    calculate_group_statistics_standardized
)
from .plot_config import (
    save_publication_figure, DPI_MAIN,
    TITLE_SIZE, AXIS_LABEL_SIZE, TICK_LABEL_SIZE,
    ALPHA_MEDIUM_LIGHT,  # Alpha constants
    EDGE_COLOR_BLACK,# Edge color standardization
    HEATMAP_CMAP_FOLDCHANGE,  # Colormap constants (Phase 10.3.11)
    # Zorder constants (Phase 10.3.7)
    ZORDER_BACKGROUND, ZORDER_GRID, ZORDER_SEPARATOR,
    ZORDER_DATA_LOW, ZORDER_DATA_HIGH,
    ZORDER_THRESHOLD, ZORDER_ANNOTATION,
    ZORDER_OVERLAY, ZORDER_EFFECT,
    ZORDER_TOP, ZORDER_ABSOLUTE_TOP
)

logger = logging.getLogger(__name__)


class VIPScorePlotMixin:
    """Mixin class for VIP score-related plots"""

    def _plot_vip_scores_base(
        self, df: pd.DataFrame, display_data: pd.DataFrame,
        mask_fn, aggregation_method: str,
        ylabel: str, title: str, output_suffix: str,
        label_column: str, figsize: tuple = (10, 6)
    ):
        """
        Unified base method for VIP score visualization with heatmap

        Eliminates ~235 lines of code duplication across 3 VIP score methods.

        Args:
            df: Annotated DataFrame with all glycopeptides
            display_data: DataFrame with VIP_Score column and data to display (top N items)
            mask_fn: Function(df, row) -> boolean mask for selecting data
            aggregation_method: 'mean' or 'sum' - how to aggregate cancer/normal values
            ylabel: Y-axis label for VIP score plot
            title: Plot title
            output_suffix: Suffix for output filename (e.g., 'glycopeptide', 'glycan_composition')
            label_column: Column name in display_data to use for y-axis labels
            figsize: Figure size tuple

        Returns:
            None (saves figure to output directory)

        Pattern Used:
            Strategy Pattern - mask_fn and aggregation_method parameterize the differing behavior
        """
        # Get sample columns
        cancer_samples, normal_samples = get_sample_columns(df)

        # STANDARDIZED: Prepare heatmap data using centralized statistics
        config = DataPreparationConfig(missing_data_method='skipna')
        heatmap_data = []

        for _, row in display_data.iterrows():
            mask = mask_fn(df, row)

            if mask.sum() > 0:
                matched_rows = df[mask]

                # Use standardized statistics calculation
                cancer_stats = calculate_group_statistics_standardized(
                    matched_rows, cancer_samples, method=config.missing_data_method
                )
                normal_stats = calculate_group_statistics_standardized(
                    matched_rows, normal_samples, method=config.missing_data_method
                )

                # Extract values based on aggregation method
                if aggregation_method == 'mean':
                    cancer_value = cancer_stats['mean'].iloc[0] if not cancer_stats['mean'].isna().all() else 0
                    normal_value = normal_stats['mean'].iloc[0] if not normal_stats['mean'].isna().all() else 0
                elif aggregation_method == 'sum':
                    cancer_value = cancer_stats['sum'].sum()
                    normal_value = normal_stats['sum'].sum()
                else:
                    raise ValueError(f"Unknown aggregation method: {aggregation_method}")

                heatmap_data.append([cancer_value, normal_value])
            else:
                heatmap_data.append([0, 0])

        heatmap_df = pd.DataFrame(heatmap_data, columns=['Cancer', 'Normal'])

        # Binary classification: 1.0 for higher value, 0.0 for lower value
        heatmap_normalized = pd.DataFrame(0.0, index=heatmap_df.index, columns=heatmap_df.columns)
        for idx in heatmap_df.index:
            if heatmap_df.loc[idx, 'Cancer'] > heatmap_df.loc[idx, 'Normal']:
                heatmap_normalized.loc[idx, 'Cancer'] = 1.0
                heatmap_normalized.loc[idx, 'Normal'] = 0.0
            else:
                heatmap_normalized.loc[idx, 'Cancer'] = 0.0
                heatmap_normalized.loc[idx, 'Normal'] = 1.0

        # Create figure with GridSpec
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(1, 2, width_ratios=[3, 1], wspace=0.05, figure=fig)

        ax_vip = fig.add_subplot(gs[0])
        ax_heatmap = fig.add_subplot(gs[1], sharey=ax_vip)

        # Plot VIP scores (left)
        ax_vip.scatter(display_data['VIP_Score'], range(len(display_data)),
                       c='#555555', s=100, edgecolors=EDGE_COLOR_BLACK, linewidths=1.5, zorder=ZORDER_DATA_LOW)
        ax_vip.set_yticks(range(len(display_data)))
        ax_vip.set_yticklabels(display_data[label_column], fontsize=TICK_LABEL_SIZE)
        ax_vip.set_xlabel('VIP Score', fontsize=AXIS_LABEL_SIZE)
        ax_vip.set_ylabel(ylabel, fontsize=AXIS_LABEL_SIZE)
        ax_vip.set_title(title, fontsize=TITLE_SIZE, fontweight='bold', pad=20)
        ax_vip.invert_yaxis()
        ax_vip.grid(axis='x', alpha=ALPHA_MEDIUM_LIGHT)

        # Plot heatmap (right) with square cells
        sns.heatmap(heatmap_normalized, ax=ax_heatmap, cmap=HEATMAP_CMAP_FOLDCHANGE,
                    cbar_kws={'label': 'Relative Intensity', 'shrink': 0.5,
                              'ticks': [0, 1]},
                    linewidths=0.5, linecolor='white',
                    vmin=0, vmax=1, yticklabels=False, square=True)

        # Customize colorbar labels
        cbar = ax_heatmap.collections[0].colorbar
        cbar.set_ticks([0, 1])
        cbar.set_ticklabels(['Low', 'High'])

        ax_heatmap.set_xticks([0.5, 1.5])
        ax_heatmap.set_xticklabels(['Cancer', 'Normal'], rotation=0, fontsize=TICK_LABEL_SIZE)
        ax_heatmap.set_xlabel('')
        ax_heatmap.tick_params(left=False)

        plt.tight_layout()

        output_file = self.output_dir / f'vip_score_{output_suffix}.png'
        save_publication_figure(fig, output_file, dpi=DPI_MAIN)
        logger.info(f"Saved VIP score {output_suffix} plot to {output_file} (optimized, {DPI_MAIN} DPI)")

        plt.close()

    def plot_vip_scores_glycopeptide(
        self, df: pd.DataFrame, vip_df: pd.DataFrame,
        figsize: tuple = (10, 6), top_n: int = 10
    ):
        """
        Plot top VIP scores by glycopeptide with heatmap showing Cancer/Normal intensity

        Args:
            df: Annotated DataFrame
            vip_df: DataFrame with VIP scores by glycopeptide
            figsize: Figure size
            top_n: Number of top glycopeptides to show
        """
        # Prepare display data
        top_n_data = vip_df.head(top_n).copy()
        top_n_data['Label'] = top_n_data['Peptide'] + ' | ' + top_n_data['GlycanComposition']

        # Define mask function for glycopeptide matching
        def mask_fn(df, row):
            return (df['Peptide'] == row['Peptide']) & (df['GlycanComposition'] == row['GlycanComposition'])

        # Call unified base method
        self._plot_vip_scores_base(
            df=df,
            display_data=top_n_data,
            mask_fn=mask_fn,
            aggregation_method='mean',
            ylabel='Glycopeptide',
            title=f'Top {top_n} Glycopeptides by VIP Score',
            output_suffix='glycopeptide',
            label_column='Label',
            figsize=figsize
        )

    def plot_vip_scores_glycan_composition(
        self, df: pd.DataFrame, vip_df: pd.DataFrame,
        figsize: tuple = (10, 6), top_n: int = 10
    ):
        """
        Plot VIP scores by GlycanComposition with heatmap

        Args:
            df: Annotated DataFrame
            vip_df: DataFrame with all VIP scores by glycopeptide
            figsize: Figure size
            top_n: Number of top GlycanCompositions to show
        """
        # Prepare display data: group by GlycanComposition and get max VIP score
        glycan_vip = vip_df.groupby('GlycanComposition')['VIP_Score'].max().nlargest(top_n).reset_index()

        # Define mask function for glycan composition matching
        def mask_fn(df, row):
            return df['GlycanComposition'] == row['GlycanComposition']

        # Call unified base method
        self._plot_vip_scores_base(
            df=df,
            display_data=glycan_vip,
            mask_fn=mask_fn,
            aggregation_method='sum',
            ylabel='Glycan Composition',
            title=f'Top {top_n} Glycan Compositions by VIP Score',
            output_suffix='glycan_composition',
            label_column='GlycanComposition',
            figsize=figsize
        )

    def plot_vip_scores_peptide(
        self, df: pd.DataFrame, vip_df: pd.DataFrame,
        figsize: tuple = (10, 6), top_n: int = 10
    ):
        """
        Plot VIP scores by Peptide with heatmap

        Args:
            df: Annotated DataFrame
            vip_df: DataFrame with all VIP scores
            figsize: Figure size
            top_n: Number of top peptides to show
        """
        # Prepare display data: group by Peptide and get max VIP score
        peptide_vip = vip_df.groupby('Peptide')['VIP_Score'].max().nlargest(top_n).reset_index()

        # Define mask function for peptide matching
        def mask_fn(df, row):
            return df['Peptide'] == row['Peptide']

        # Call unified base method
        self._plot_vip_scores_base(
            df=df,
            display_data=peptide_vip,
            mask_fn=mask_fn,
            aggregation_method='sum',
            ylabel='Peptide',
            title=f'Top {top_n} Peptides by VIP Score',
            output_suffix='peptide',
            label_column='Peptide',
            figsize=figsize
        )
