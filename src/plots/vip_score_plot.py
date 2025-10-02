"""
VIP Score Plot Module for pGlyco Auto Combine
Handles VIP score visualizations with heatmap
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from matplotlib.gridspec import GridSpec

logger = logging.getLogger(__name__)


class VIPScorePlotMixin:
    """Mixin class for VIP score-related plots"""

    def plot_vip_scores_glycopeptide(self, df: pd.DataFrame, vip_df: pd.DataFrame, figsize: tuple = (14, 10), top_n: int = 30):
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
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']
        sample_cols = [col for col in df.columns if col not in metadata_cols]
        cancer_samples = [col for col in sample_cols if col.startswith('C')]
        normal_samples = [col for col in sample_cols if col.startswith('N')]

        # Prepare heatmap data (Cancer mean, Normal mean for each glycopeptide)
        heatmap_data = []
        for _, row in top_n_data.iterrows():
            mask = (df['Peptide'] == row['Peptide']) & (df['GlycanComposition'] == row['GlycanComposition'])
            if mask.sum() > 0:
                glycopeptide_row = df[mask].iloc[0]

                cancer_mean = glycopeptide_row[cancer_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0).mean()
                normal_mean = glycopeptide_row[normal_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0).mean()

                heatmap_data.append([cancer_mean, normal_mean])
            else:
                heatmap_data.append([0, 0])

        heatmap_df = pd.DataFrame(heatmap_data, columns=['Cancer', 'Normal'])

        # Normalize each row (feature) to 0-1 for heatmap
        heatmap_normalized = heatmap_df.div(heatmap_df.max(axis=1), axis=0).fillna(0)

        # Create figure with GridSpec
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(1, 2, width_ratios=[3, 1], wspace=0.05, figure=fig)

        ax_vip = fig.add_subplot(gs[0])
        ax_heatmap = fig.add_subplot(gs[1], sharey=ax_vip)

        # Plot VIP scores (left)
        ax_vip.scatter(top_n_data['VIP_Score'], range(len(top_n_data)),
                      c='#555555', s=100, edgecolors='black', linewidths=1.5, zorder=3)
        ax_vip.set_yticks(range(len(top_n_data)))
        ax_vip.set_yticklabels(top_n_data['Label'], fontsize=8)
        ax_vip.set_xlabel('VIP Score', fontsize=12)
        ax_vip.set_ylabel('Glycopeptide', fontsize=12)
        ax_vip.set_title(f'Top {top_n} Glycopeptides by VIP Score', fontsize=14, fontweight='bold', pad=20)
        ax_vip.invert_yaxis()
        ax_vip.grid(axis='x', alpha=0.3)

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
        ax_heatmap.set_xticklabels(['Cancer', 'Normal'], rotation=0, fontsize=10)
        ax_heatmap.set_xlabel('')
        ax_heatmap.tick_params(left=False)

        plt.tight_layout()

        output_file = self.output_dir / 'vip_score_glycopeptide.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved VIP score glycopeptide plot to {output_file}")

        plt.close()

    def plot_vip_scores_glycan_composition(self, df: pd.DataFrame, vip_df: pd.DataFrame, figsize: tuple = (14, 10), top_n: int = 30):
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
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']
        sample_cols = [col for col in df.columns if col not in metadata_cols]
        cancer_samples = [col for col in sample_cols if col.startswith('C')]
        normal_samples = [col for col in sample_cols if col.startswith('N')]

        # Prepare heatmap data
        heatmap_data = []
        for _, row in glycan_vip.iterrows():
            glycan_comp = row['GlycanComposition']
            mask = df['GlycanComposition'] == glycan_comp

            if mask.sum() > 0:
                glycan_rows = df[mask]

                # Calculate total intensity across all peptides with this glycan
                cancer_total = 0
                normal_total = 0
                for _, gly_row in glycan_rows.iterrows():
                    cancer_total += gly_row[cancer_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0).sum()
                    normal_total += gly_row[normal_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0).sum()

                heatmap_data.append([cancer_total, normal_total])
            else:
                heatmap_data.append([0, 0])

        heatmap_df = pd.DataFrame(heatmap_data, columns=['Cancer', 'Normal'])

        # Normalize each row to 0-1
        heatmap_normalized = heatmap_df.div(heatmap_df.max(axis=1), axis=0).fillna(0)

        # Create figure with GridSpec
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(1, 2, width_ratios=[3, 1], wspace=0.05, figure=fig)

        ax_vip = fig.add_subplot(gs[0])
        ax_heatmap = fig.add_subplot(gs[1], sharey=ax_vip)

        # Plot VIP scores (left)
        ax_vip.scatter(glycan_vip['VIP_Score'], range(len(glycan_vip)),
                      c='#555555', s=100, edgecolors='black', linewidths=1.5, zorder=3)
        ax_vip.set_yticks(range(len(glycan_vip)))
        ax_vip.set_yticklabels(glycan_vip['GlycanComposition'], fontsize=9)
        ax_vip.set_xlabel('VIP Score', fontsize=12)
        ax_vip.set_ylabel('Glycan Composition', fontsize=12)
        ax_vip.set_title(f'Top {top_n} Glycan Compositions by VIP Score', fontsize=14, fontweight='bold', pad=20)
        ax_vip.invert_yaxis()
        ax_vip.grid(axis='x', alpha=0.3)

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
        ax_heatmap.set_xticklabels(['Cancer', 'Normal'], rotation=0, fontsize=10)
        ax_heatmap.set_xlabel('')
        ax_heatmap.tick_params(left=False)

        plt.tight_layout()

        output_file = self.output_dir / 'vip_score_glycan_composition.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved VIP score glycan composition plot to {output_file}")

        plt.close()

    def plot_vip_scores_peptide(self, df: pd.DataFrame, vip_df: pd.DataFrame, figsize: tuple = (14, 10), top_n: int = 30):
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
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']
        sample_cols = [col for col in df.columns if col not in metadata_cols]
        cancer_samples = [col for col in sample_cols if col.startswith('C')]
        normal_samples = [col for col in sample_cols if col.startswith('N')]

        # Prepare heatmap data
        heatmap_data = []
        for _, row in peptide_vip.iterrows():
            peptide = row['Peptide']
            mask = df['Peptide'] == peptide

            if mask.sum() > 0:
                peptide_rows = df[mask]

                # Calculate total intensity across all glycoforms of this peptide
                cancer_total = 0
                normal_total = 0
                for _, pep_row in peptide_rows.iterrows():
                    cancer_total += pep_row[cancer_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0).sum()
                    normal_total += pep_row[normal_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0).sum()

                heatmap_data.append([cancer_total, normal_total])
            else:
                heatmap_data.append([0, 0])

        heatmap_df = pd.DataFrame(heatmap_data, columns=['Cancer', 'Normal'])

        # Normalize each row to 0-1
        heatmap_normalized = heatmap_df.div(heatmap_df.max(axis=1), axis=0).fillna(0)

        # Create figure with GridSpec
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(1, 2, width_ratios=[3, 1], wspace=0.05, figure=fig)

        ax_vip = fig.add_subplot(gs[0])
        ax_heatmap = fig.add_subplot(gs[1], sharey=ax_vip)

        # Plot VIP scores (left)
        ax_vip.scatter(peptide_vip['VIP_Score'], range(len(peptide_vip)),
                      c='#555555', s=100, edgecolors='black', linewidths=1.5, zorder=3)
        ax_vip.set_yticks(range(len(peptide_vip)))
        ax_vip.set_yticklabels(peptide_vip['Peptide'], fontsize=8)
        ax_vip.set_xlabel('VIP Score', fontsize=12)
        ax_vip.set_ylabel('Peptide', fontsize=12)
        ax_vip.set_title(f'Top {top_n} Peptides by VIP Score', fontsize=14, fontweight='bold', pad=20)
        ax_vip.invert_yaxis()
        ax_vip.grid(axis='x', alpha=0.3)

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
        ax_heatmap.set_xticklabels(['Cancer', 'Normal'], rotation=0, fontsize=10)
        ax_heatmap.set_xlabel('')
        ax_heatmap.tick_params(left=False)

        plt.tight_layout()

        output_file = self.output_dir / 'vip_score_peptide.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved VIP score peptide plot to {output_file}")

        plt.close()
