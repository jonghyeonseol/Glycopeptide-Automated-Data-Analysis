"""
VIP Score Plot Module for pGlyco Auto Combine
Handles VIP score visualizations
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class VIPScorePlotMixin:
    """Mixin class for VIP score-related plots"""

    def plot_vip_scores_glycopeptide(self, df: pd.DataFrame, vip_df: pd.DataFrame, figsize: tuple = (10, 8)):
        """
        Plot top 10 VIP scores by glycopeptide (dot plot with Cancer/Normal coloring)

        Args:
            df: Annotated DataFrame
            vip_df: DataFrame with VIP scores by glycopeptide
            figsize: Figure size
        """
        top_10 = vip_df.head(10).copy()

        # Create labels
        top_10['Label'] = top_10['Peptide'] + ' | ' + top_10['GlycanComposition']

        # Calculate Cancer vs Normal mean intensity for color determination
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']
        sample_cols = [col for col in df.columns if col not in metadata_cols]
        cancer_samples = [col for col in sample_cols if col.startswith('C')]
        normal_samples = [col for col in sample_cols if col.startswith('N')]

        colors = []
        for _, row in top_10.iterrows():
            # Find this glycopeptide in the dataframe
            mask = (df['Peptide'] == row['Peptide']) & (df['GlycanComposition'] == row['GlycanComposition'])
            if mask.sum() > 0:
                glycopeptide_row = df[mask].iloc[0]

                # Calculate means
                cancer_mean = glycopeptide_row[cancer_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0).mean()
                normal_mean = glycopeptide_row[normal_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0).mean()

                # Determine color
                if cancer_mean > normal_mean:
                    colors.append('#E74C3C')  # Red - High in Cancer
                else:
                    colors.append('#3498DB')  # Blue - High in Normal
            else:
                colors.append('#CCCCCC')  # Gray - not found

        fig, ax = plt.subplots(figsize=figsize)

        # Plot dots
        ax.scatter(top_10['VIP_Score'], range(len(top_10)), c=colors, s=150, edgecolors='black', linewidths=1.5, zorder=3)

        ax.set_yticks(range(len(top_10)))
        ax.set_yticklabels(top_10['Label'], fontsize=9)
        ax.set_xlabel('VIP Score', fontsize=12)
        ax.set_ylabel('Glycopeptide', fontsize=12)
        ax.set_title('Top 10 Glycopeptides by VIP Score', fontsize=14, fontweight='bold')
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.3)

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#E74C3C', edgecolor='black', label='High in Cancer'),
            Patch(facecolor='#3498DB', edgecolor='black', label='High in Normal')
        ]
        ax.legend(handles=legend_elements, loc='lower right', frameon=True, fontsize=10)

        plt.tight_layout()

        output_file = self.output_dir / 'vip_score_glycopeptide.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved VIP score glycopeptide plot to {output_file}")

        plt.close()

    def plot_vip_scores_glycan_type(self, df: pd.DataFrame, vip_df: pd.DataFrame,
                                     classification_col: str = 'SecondaryClassification',
                                     figsize: tuple = (12, 8)):
        """
        Plot VIP scores by glycan type (dot plot showing top 10 glycopeptides per type)

        Args:
            df: Annotated DataFrame
            vip_df: DataFrame with all VIP scores by glycopeptide
            classification_col: Column to use for classification
            figsize: Figure size
        """
        # Merge VIP scores with classification info
        vip_with_class = vip_df.merge(
            df[['Peptide', 'GlycanComposition', classification_col]],
            on=['Peptide', 'GlycanComposition'],
            how='left'
        )

        # Get sample columns for Cancer/Normal comparison
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']
        sample_cols = [col for col in df.columns if col not in metadata_cols]
        cancer_samples = [col for col in sample_cols if col.startswith('C')]
        normal_samples = [col for col in sample_cols if col.startswith('N')]

        # Get top glycan types by mean VIP score (limit to top 6 for readability)
        top_types = vip_with_class.groupby(classification_col)['VIP_Score'].mean().nlargest(6).index.tolist()

        # Exclude Truncated and Outlier if present
        top_types = [t for t in top_types if t not in ['Truncated', 'Outlier', 'Unknown']][:5]

        # Collect data for plotting
        plot_data = []
        for glycan_type in top_types:
            # Get top 10 glycopeptides for this type
            type_data = vip_with_class[vip_with_class[classification_col] == glycan_type].nlargest(10, 'VIP_Score')

            for _, row in type_data.iterrows():
                # Find intensity in original df for color determination
                mask = (df['Peptide'] == row['Peptide']) & (df['GlycanComposition'] == row['GlycanComposition'])
                if mask.sum() > 0:
                    glycopeptide_row = df[mask].iloc[0]
                    cancer_mean = glycopeptide_row[cancer_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0).mean()
                    normal_mean = glycopeptide_row[normal_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0).mean()

                    color = '#E74C3C' if cancer_mean > normal_mean else '#3498DB'
                else:
                    color = '#CCCCCC'

                plot_data.append({
                    'GlycanType': glycan_type,
                    'VIP_Score': row['VIP_Score'],
                    'Color': color
                })

        plot_df = pd.DataFrame(plot_data)

        # Create plot
        fig, ax = plt.subplots(figsize=figsize)

        # Plot each glycan type
        y_positions = []
        y_labels = []
        for i, glycan_type in enumerate(top_types):
            type_df = plot_df[plot_df['GlycanType'] == glycan_type]
            y_pos = [i] * len(type_df)

            # Add jitter to y positions for visibility
            y_pos_jittered = [y + np.random.uniform(-0.2, 0.2) for y in y_pos]

            ax.scatter(type_df['VIP_Score'], y_pos_jittered, c=type_df['Color'],
                      s=80, edgecolors='black', linewidths=1, alpha=0.7, zorder=3)

            y_positions.append(i)
            y_labels.append(glycan_type)

        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels, fontsize=11)
        ax.set_xlabel('VIP Score', fontsize=12)
        ax.set_ylabel('Glycan Type', fontsize=12)
        ax.set_title('VIP Scores by Glycan Type (Top 10 per type)', fontsize=14, fontweight='bold')
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.3)

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#E74C3C', edgecolor='black', label='High in Cancer'),
            Patch(facecolor='#3498DB', edgecolor='black', label='High in Normal')
        ]
        ax.legend(handles=legend_elements, loc='lower right', frameon=True, fontsize=10)

        plt.tight_layout()

        output_file = self.output_dir / 'vip_score_glycan_type.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved VIP score glycan type plot to {output_file}")

        plt.close()

    def plot_vip_scores_peptide_top10(self, df: pd.DataFrame, vip_df: pd.DataFrame, figsize: tuple = (12, 8)):
        """
        Plot VIP scores for top 10 peptides (dot plot showing top 10 glycopeptides per peptide)

        Args:
            df: Annotated DataFrame
            vip_df: DataFrame with all VIP scores
            figsize: Figure size
        """
        # Get top 10 peptides by mean VIP score
        peptide_vip = vip_df.groupby('Peptide')['VIP_Score'].mean().nlargest(10).index.tolist()

        # Get sample columns for Cancer/Normal comparison
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']
        sample_cols = [col for col in df.columns if col not in metadata_cols]
        cancer_samples = [col for col in sample_cols if col.startswith('C')]
        normal_samples = [col for col in sample_cols if col.startswith('N')]

        # Collect data for plotting
        plot_data = []
        for peptide in peptide_vip:
            # Get top 10 glycopeptides for this peptide
            peptide_data = vip_df[vip_df['Peptide'] == peptide].nlargest(10, 'VIP_Score')

            for _, row in peptide_data.iterrows():
                # Find intensity in original df for color determination
                mask = (df['Peptide'] == row['Peptide']) & (df['GlycanComposition'] == row['GlycanComposition'])
                if mask.sum() > 0:
                    glycopeptide_row = df[mask].iloc[0]
                    cancer_mean = glycopeptide_row[cancer_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0).mean()
                    normal_mean = glycopeptide_row[normal_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0).mean()

                    color = '#E74C3C' if cancer_mean > normal_mean else '#3498DB'
                else:
                    color = '#CCCCCC'

                plot_data.append({
                    'Peptide': peptide,
                    'VIP_Score': row['VIP_Score'],
                    'Color': color
                })

        plot_df = pd.DataFrame(plot_data)

        # Create plot
        fig, ax = plt.subplots(figsize=figsize)

        # Plot each peptide
        y_positions = []
        y_labels = []
        for i, peptide in enumerate(peptide_vip):
            peptide_df = plot_df[plot_df['Peptide'] == peptide]
            y_pos = [i] * len(peptide_df)

            # Add jitter to y positions for visibility
            y_pos_jittered = [y + np.random.uniform(-0.2, 0.2) for y in y_pos]

            ax.scatter(peptide_df['VIP_Score'], y_pos_jittered, c=peptide_df['Color'],
                      s=80, edgecolors='black', linewidths=1, alpha=0.7, zorder=3)

            y_positions.append(i)
            y_labels.append(peptide)

        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels, fontsize=9)
        ax.set_xlabel('VIP Score', fontsize=12)
        ax.set_ylabel('Peptide', fontsize=12)
        ax.set_title('VIP Scores by Peptide (Top 10 peptides, Top 10 glycoforms each)', fontsize=14, fontweight='bold')
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.3)

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#E74C3C', edgecolor='black', label='High in Cancer'),
            Patch(facecolor='#3498DB', edgecolor='black', label='High in Normal')
        ]
        ax.legend(handles=legend_elements, loc='lower right', frameon=True, fontsize=10)

        plt.tight_layout()

        output_file = self.output_dir / 'vip_score_peptide_top10.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved VIP score peptide (top 10) plot to {output_file}")

        plt.close()

    def plot_vip_scores_peptide_all(self, df: pd.DataFrame, vip_df: pd.DataFrame, figsize: tuple = (14, 10)):
        """
        Plot VIP scores for all peptides (dot plot showing all glycopeptide combinations)

        Args:
            df: Annotated DataFrame
            vip_df: DataFrame with all VIP scores
            figsize: Figure size
        """
        # Get all peptides sorted by mean VIP score
        peptide_mean_vip = vip_df.groupby('Peptide')['VIP_Score'].mean().sort_values(ascending=False)
        all_peptides = peptide_mean_vip.head(20).index.tolist()  # Limit to top 20 for readability

        # Get sample columns for Cancer/Normal comparison
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']
        sample_cols = [col for col in df.columns if col not in metadata_cols]
        cancer_samples = [col for col in sample_cols if col.startswith('C')]
        normal_samples = [col for col in sample_cols if col.startswith('N')]

        # Collect data for plotting
        plot_data = []
        for peptide in all_peptides:
            # Get ALL glycopeptides for this peptide
            peptide_data = vip_df[vip_df['Peptide'] == peptide]

            for _, row in peptide_data.iterrows():
                # Find intensity in original df for color determination
                mask = (df['Peptide'] == row['Peptide']) & (df['GlycanComposition'] == row['GlycanComposition'])
                if mask.sum() > 0:
                    glycopeptide_row = df[mask].iloc[0]
                    cancer_mean = glycopeptide_row[cancer_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0).mean()
                    normal_mean = glycopeptide_row[normal_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0).mean()

                    color = '#E74C3C' if cancer_mean > normal_mean else '#3498DB'
                else:
                    color = '#CCCCCC'

                plot_data.append({
                    'Peptide': peptide,
                    'VIP_Score': row['VIP_Score'],
                    'Color': color
                })

        plot_df = pd.DataFrame(plot_data)

        # Create plot
        fig, ax = plt.subplots(figsize=figsize)

        # Plot each peptide
        y_positions = []
        y_labels = []
        for i, peptide in enumerate(all_peptides):
            peptide_df = plot_df[plot_df['Peptide'] == peptide]
            y_pos = [i] * len(peptide_df)

            # Add jitter to y positions for visibility
            y_pos_jittered = [y + np.random.uniform(-0.3, 0.3) for y in y_pos]

            ax.scatter(peptide_df['VIP_Score'], y_pos_jittered, c=peptide_df['Color'],
                      s=60, edgecolors='black', linewidths=0.8, alpha=0.6, zorder=3)

            y_positions.append(i)
            y_labels.append(peptide)

        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels, fontsize=8)
        ax.set_xlabel('VIP Score', fontsize=12)
        ax.set_ylabel('Peptide', fontsize=12)
        ax.set_title('VIP Scores by Peptide (All glycoforms, Top 20 peptides)', fontsize=14, fontweight='bold')
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.3)

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#E74C3C', edgecolor='black', label='High in Cancer'),
            Patch(facecolor='#3498DB', edgecolor='black', label='High in Normal')
        ]
        ax.legend(handles=legend_elements, loc='lower right', frameon=True, fontsize=10)

        plt.tight_layout()

        output_file = self.output_dir / 'vip_score_peptide_all.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved VIP score peptide (all) plot to {output_file}")

        plt.close()
