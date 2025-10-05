"""
Volcano Plot Module for pGlyco Auto Combine
Handles differential expression visualization
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from scipy import stats
from adjustText import adjust_text
from ..utils import replace_empty_with_zero, save_trace_data

logger = logging.getLogger(__name__)


class VolcanoPlotMixin:
    """Mixin class for Volcano plot visualization"""

    def plot_volcano(self, df: pd.DataFrame, vip_df: pd.DataFrame,
                     fdr_threshold: float = 0.05, fc_threshold: float = 1.5,
                     figsize: tuple = (12, 10)):
        """
        Create volcano plot showing log2(fold change) vs -log10(FDR)

        Args:
            df: Annotated DataFrame with intensity data
            vip_df: DataFrame with VIP scores
            fdr_threshold: FDR significance threshold (default 0.05)
            fc_threshold: Fold change threshold (default 1.5)
            figsize: Figure size (width, height)
        """
        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]

        # Prepare data for volcano plot
        volcano_data = []

        for idx, row in df.iterrows():
            peptide = row['Peptide']
            glycan_comp = row['GlycanComposition']
            glycopeptide = f"{peptide}_{glycan_comp}"

            # Get intensity values
            cancer_values = replace_empty_with_zero(row[cancer_samples]).values.astype(float)
            normal_values = replace_empty_with_zero(row[normal_samples]).values.astype(float)

            # Filter out zeros for statistics
            cancer_nonzero = cancer_values[cancer_values > 0]
            normal_nonzero = normal_values[normal_values > 0]

            # Skip if insufficient data
            if len(cancer_nonzero) < 3 or len(normal_nonzero) < 3:
                continue

            # Calculate means
            cancer_mean = np.mean(cancer_nonzero)
            normal_mean = np.mean(normal_nonzero)

            # Calculate fold change (avoid division by zero)
            if normal_mean > 0:
                fold_change = cancer_mean / normal_mean
                log2_fc = np.log2(fold_change)
            else:
                log2_fc = np.nan

            # Perform Wilcoxon rank-sum test (Mann-Whitney U)
            try:
                statistic, p_value = stats.mannwhitneyu(cancer_nonzero, normal_nonzero, alternative='two-sided')
            except:
                p_value = np.nan

            # Get VIP score if available
            vip_score = 0
            vip_match = vip_df[(vip_df['Peptide'] == peptide) &
                              (vip_df['GlycanComposition'] == glycan_comp)]
            if not vip_match.empty:
                vip_score = vip_match['VIP_Score'].values[0]

            # Get glycan type category if available
            glycan_type = row.get('GlycanTypeCategory', 'C/H')

            volcano_data.append({
                'Glycopeptide': glycopeptide,
                'Peptide': peptide,
                'GlycanComposition': glycan_comp,
                'GlycanTypeCategory': glycan_type,
                'Log2FC': log2_fc,
                'P_value': p_value,
                'VIP_Score': vip_score,
                'Cancer_Mean': cancer_mean,
                'Normal_Mean': normal_mean
            })

        # Create DataFrame
        volcano_df = pd.DataFrame(volcano_data)

        # Remove rows with NaN
        volcano_df = volcano_df.dropna(subset=['Log2FC', 'P_value'])

        # Calculate -log10(p-value)
        volcano_df['-Log10P'] = -np.log10(volcano_df['P_value'])

        # FDR correction using Benjamini-Hochberg
        from statsmodels.stats.multitest import multipletests
        _, fdr_values, _, _ = multipletests(volcano_df['P_value'].values, method='fdr_bh')
        volcano_df['FDR'] = fdr_values
        volcano_df['-Log10FDR'] = -np.log10(fdr_values)

        # Classification
        volcano_df['Regulation'] = 'Non-significant'

        # Up-regulated: log2FC > threshold AND FDR < threshold
        up_mask = (volcano_df['Log2FC'] > np.log2(fc_threshold)) & (volcano_df['FDR'] < fdr_threshold)
        volcano_df.loc[up_mask, 'Regulation'] = 'Up in Cancer'

        # Down-regulated: log2FC < -threshold AND FDR < threshold
        down_mask = (volcano_df['Log2FC'] < -np.log2(fc_threshold)) & (volcano_df['FDR'] < fdr_threshold)
        volcano_df.loc[down_mask, 'Regulation'] = 'Down in Cancer'

        # Create plot
        fig, ax = plt.subplots(figsize=figsize)

        # Define colors
        colors = {
            'Up in Cancer': '#E74C3C',
            'Down in Cancer': '#3498DB',
            'Non-significant': '#95A5A6'
        }

        # Plot points by regulation status
        for regulation, color in colors.items():
            mask = volcano_df['Regulation'] == regulation
            subset = volcano_df[mask]

            # Size by VIP score (scale between 20 and 200)
            sizes = 20 + (subset['VIP_Score'] / subset['VIP_Score'].max() * 180) if subset['VIP_Score'].max() > 0 else 50

            ax.scatter(subset['Log2FC'], subset['-Log10FDR'],
                      c=color, s=sizes, alpha=0.6,
                      edgecolors='black', linewidth=0.5,
                      label=f"{regulation} (n={len(subset)})", zorder=3)

        # Add threshold lines
        ax.axhline(-np.log10(fdr_threshold), color='gray', linestyle='--', linewidth=1, alpha=0.5, zorder=1)
        ax.axvline(np.log2(fc_threshold), color='gray', linestyle='--', linewidth=1, alpha=0.5, zorder=1)
        ax.axvline(-np.log2(fc_threshold), color='gray', linestyle='--', linewidth=1, alpha=0.5, zorder=1)

        # Define glycan type colors
        glycan_type_colors = {
            'HM': '#00CC00',   # Green
            'F': '#FF0000',    # Red
            'S': '#FF69B4',    # Pink
            'SF': '#FFA500',   # Orange
            'C/H': '#0000FF'   # Blue
        }

        # Find significant increases and decreases (p < 0.05 AND |FC| > 2)
        # Rank by combined score: |log2FC| * -log10(p-value) for better prioritization

        # Increased group: high log2FC and low p-value
        increased_candidates = volcano_df[
            (volcano_df['P_value'] < 0.05) & (volcano_df['Log2FC'] > 1)
        ].copy()
        if len(increased_candidates) > 0:
            increased_candidates['Score'] = increased_candidates['Log2FC'] * increased_candidates['-Log10P']
            significant_increases = increased_candidates.nlargest(3, 'Score')
        else:
            significant_increases = pd.DataFrame()

        # Decreased group: low (negative) log2FC and low p-value
        decreased_candidates = volcano_df[
            (volcano_df['P_value'] < 0.05) & (volcano_df['Log2FC'] < -1)
        ].copy()
        if len(decreased_candidates) > 0:
            decreased_candidates['Score'] = abs(decreased_candidates['Log2FC']) * decreased_candidates['-Log10P']
            significant_decreases = decreased_candidates.nlargest(3, 'Score')
        else:
            significant_decreases = pd.DataFrame()

        # Combine for annotation
        to_annotate = pd.concat([significant_increases, significant_decreases])

        texts = []
        for _, row in to_annotate.iterrows():
            # Create label in format "PEPTIDE_H(5)N(4)A(2)"
            label = f"{row['Peptide']}_{row['GlycanComposition']}"

            # Get color based on glycan type
            glycan_type = row.get('GlycanTypeCategory', 'C/H')
            text_color = glycan_type_colors.get(glycan_type, '#000000')

            texts.append(ax.text(row['Log2FC'], row['-Log10FDR'], label,
                               fontsize=7, ha='center', va='bottom',
                               color=text_color, fontweight='bold',
                               bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                       edgecolor=text_color, alpha=0.8, linewidth=1.5)))

        # Adjust text to avoid overlap
        if len(texts) > 0:
            try:
                adjust_text(texts, arrowprops=dict(arrowstyle='->', color='gray', lw=0.8, alpha=0.7))
            except:
                pass  # If adjust_text fails, just show overlapping labels

        # Labels and title
        ax.set_xlabel(f'Log2 Fold Change (Cancer / Normal)', fontsize=12, fontweight='bold')
        ax.set_ylabel(f'-Log10 FDR', fontsize=12, fontweight='bold')
        ax.set_title('Volcano Plot: Differential Glycopeptide Expression\nCancer vs Normal',
                    fontsize=14, fontweight='bold')

        # Legend
        ax.legend(loc='upper right', frameon=True, fontsize=10)

        # Grid
        ax.grid(True, alpha=0.3, zorder=0)

        # Add statistics text
        n_up = len(volcano_df[volcano_df['Regulation'] == 'Up in Cancer'])
        n_down = len(volcano_df[volcano_df['Regulation'] == 'Down in Cancer'])
        n_ns = len(volcano_df[volcano_df['Regulation'] == 'Non-significant'])

        stats_text = f"FC threshold: {fc_threshold}x | FDR < {fdr_threshold}\n"
        stats_text += f"Up: {n_up} | Down: {n_down} | NS: {n_ns}"

        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
               fontsize=9, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'volcano_plot.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved volcano plot to {output_file}")

        # Save trace data
        save_trace_data(volcano_df, self.output_dir, 'volcano_plot_data.csv')

        plt.close()

        return volcano_df
