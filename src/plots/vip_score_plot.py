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
        top_n_data = vip_df.head(top_n).copy()

        # Create labels
        top_n_data['Label'] = top_n_data['Peptide'] + ' | ' + top_n_data['GlycanComposition']

        # Get sample columns
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        cancer_samples + normal_samples

        # STANDARDIZED: Prepare heatmap data using centralized mean calculation
        config = DataPreparationConfig(missing_data_method='skipna')
        heatmap_data = []

        for _, row in top_n_data.iterrows():
            mask = (df['Peptide'] == row['Peptide']) & (df['GlycanComposition'] == row['GlycanComposition'])
            if mask.sum() > 0:
                glycopeptide_row = df[mask]

                # Use standardized statistics calculation
                cancer_stats = calculate_group_statistics_standardized(
                    glycopeptide_row, cancer_samples, method=config.missing_data_method
                )
                normal_stats = calculate_group_statistics_standardized(
                    glycopeptide_row, normal_samples, method=config.missing_data_method
                )

                cancer_mean = cancer_stats['mean'].iloc[0] if not cancer_stats['mean'].isna().all() else 0
                normal_mean = normal_stats['mean'].iloc[0] if not normal_stats['mean'].isna().all() else 0

                heatmap_data.append([cancer_mean, normal_mean])
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
        ax_vip.scatter(top_n_data['VIP_Score'], range(len(top_n_data)),
                       c='#555555', s=100, edgecolors=EDGE_COLOR_BLACK, linewidths=1.5, zorder=ZORDER_DATA_LOW)
        ax_vip.set_yticks(range(len(top_n_data)))
        ax_vip.set_yticklabels(top_n_data['Label'], fontsize=TICK_LABEL_SIZE)
        ax_vip.set_xlabel('VIP Score', fontsize=AXIS_LABEL_SIZE)
        ax_vip.set_ylabel('Glycopeptide', fontsize=AXIS_LABEL_SIZE)
        ax_vip.set_title(f'Top {top_n} Glycopeptides by VIP Score', fontsize=TITLE_SIZE, fontweight='bold', pad=20)
        ax_vip.invert_yaxis()
        ax_vip.grid(axis='x', alpha=ALPHA_MEDIUM_LIGHT)

        # Plot heatmap (right) with square cells
        sns.heatmap(heatmap_normalized, ax=ax_heatmap, cmap='RdBu_r',
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

        output_file = self.output_dir / 'vip_score_glycopeptide.png'
        save_publication_figure(fig, output_file, dpi=DPI_MAIN)
        logger.info(f"Saved VIP score glycopeptide plot to {output_file} (optimized, {DPI_MAIN} DPI)")

        plt.close()

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
        # Group by GlycanComposition and get the max VIP score for each
        glycan_vip = vip_df.groupby('GlycanComposition')['VIP_Score'].max().nlargest(top_n).reset_index()

        # Get sample columns
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        cancer_samples + normal_samples

        # STANDARDIZED: Prepare heatmap data using centralized statistics
        config = DataPreparationConfig(missing_data_method='skipna')
        heatmap_data = []

        for _, row in glycan_vip.iterrows():
            glycan_comp = row['GlycanComposition']
            mask = df['GlycanComposition'] == glycan_comp

            if mask.sum() > 0:
                glycan_rows = df[mask]

                # Use standardized statistics calculation
                cancer_stats = calculate_group_statistics_standardized(
                    glycan_rows, cancer_samples, method=config.missing_data_method
                )
                normal_stats = calculate_group_statistics_standardized(
                    glycan_rows, normal_samples, method=config.missing_data_method
                )

                # Sum across all peptides with this glycan
                cancer_total = cancer_stats['sum'].sum()
                normal_total = normal_stats['sum'].sum()

                heatmap_data.append([cancer_total, normal_total])
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
        ax_vip.scatter(glycan_vip['VIP_Score'], range(len(glycan_vip)),
                       c='#555555', s=100, edgecolors=EDGE_COLOR_BLACK, linewidths=1.5, zorder=ZORDER_DATA_LOW)
        ax_vip.set_yticks(range(len(glycan_vip)))
        ax_vip.set_yticklabels(glycan_vip['GlycanComposition'], fontsize=TICK_LABEL_SIZE)
        ax_vip.set_xlabel('VIP Score', fontsize=AXIS_LABEL_SIZE)
        ax_vip.set_ylabel('Glycan Composition', fontsize=AXIS_LABEL_SIZE)
        ax_vip.set_title(f'Top {top_n} Glycan Compositions by VIP Score', fontsize=TITLE_SIZE, fontweight='bold', pad=20)
        ax_vip.invert_yaxis()
        ax_vip.grid(axis='x', alpha=ALPHA_MEDIUM_LIGHT)

        # Plot heatmap (right) with square cells
        sns.heatmap(heatmap_normalized, ax=ax_heatmap, cmap='RdBu_r',
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

        output_file = self.output_dir / 'vip_score_glycan_composition.png'
        save_publication_figure(fig, output_file, dpi=DPI_MAIN)
        logger.info(f"Saved VIP score glycan composition plot to {output_file} (optimized, {DPI_MAIN} DPI)")

        plt.close()

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
        # Group by Peptide and get the max VIP score for each
        peptide_vip = vip_df.groupby('Peptide')['VIP_Score'].max().nlargest(top_n).reset_index()

        # Get sample columns
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        cancer_samples + normal_samples

        # STANDARDIZED: Prepare heatmap data using centralized statistics
        config = DataPreparationConfig(missing_data_method='skipna')
        heatmap_data = []

        for _, row in peptide_vip.iterrows():
            peptide = row['Peptide']
            mask = df['Peptide'] == peptide

            if mask.sum() > 0:
                peptide_rows = df[mask]

                # Use standardized statistics calculation
                cancer_stats = calculate_group_statistics_standardized(
                    peptide_rows, cancer_samples, method=config.missing_data_method
                )
                normal_stats = calculate_group_statistics_standardized(
                    peptide_rows, normal_samples, method=config.missing_data_method
                )

                # Sum across all glycoforms of this peptide
                cancer_total = cancer_stats['sum'].sum()
                normal_total = normal_stats['sum'].sum()

                heatmap_data.append([cancer_total, normal_total])
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
        ax_vip.scatter(peptide_vip['VIP_Score'], range(len(peptide_vip)),
                       c='#555555', s=100, edgecolors=EDGE_COLOR_BLACK, linewidths=1.5, zorder=ZORDER_DATA_LOW)
        ax_vip.set_yticks(range(len(peptide_vip)))
        ax_vip.set_yticklabels(peptide_vip['Peptide'], fontsize=TICK_LABEL_SIZE)
        ax_vip.set_xlabel('VIP Score', fontsize=AXIS_LABEL_SIZE)
        ax_vip.set_ylabel('Peptide', fontsize=AXIS_LABEL_SIZE)
        ax_vip.set_title(f'Top {top_n} Peptides by VIP Score', fontsize=TITLE_SIZE, fontweight='bold', pad=20)
        ax_vip.invert_yaxis()
        ax_vip.grid(axis='x', alpha=ALPHA_MEDIUM_LIGHT)

        # Plot heatmap (right) with square cells
        sns.heatmap(heatmap_normalized, ax=ax_heatmap, cmap='RdBu_r',
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

        output_file = self.output_dir / 'vip_score_peptide.png'
        save_publication_figure(fig, output_file, dpi=DPI_MAIN)
        logger.info(f"Saved VIP score peptide plot to {output_file} (optimized, {DPI_MAIN} DPI)")

        plt.close()
