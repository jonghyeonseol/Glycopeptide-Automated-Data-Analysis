"""
Visualizer Module for pGlyco Auto Combine
Handles data visualization (PCA, boxplot, heatmap)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class GlycanVisualizer:
    """Create visualizations for glycoproteomics data"""

    def __init__(self, output_dir: str, dpi: int = 300, colors: dict = None):
        """
        Initialize GlycanVisualizer

        Args:
            output_dir: Directory to save plots
            dpi: Resolution for saved figures
            colors: Color scheme for glycan types
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.dpi = dpi

        # Default color scheme
        self.colors = colors or {
            'Non': '#CCCCCC',
            'Sialylated': '#E74C3C',
            'Fucosylated': '#3498DB',
            'Both': '#9B59B6'
        }

        # Set style
        sns.set_style("whitegrid")
        plt.rcParams['font.size'] = 10
        plt.rcParams['axes.labelsize'] = 12
        plt.rcParams['axes.titlesize'] = 14

    def plot_pca(self, pca_results: dict, figsize: tuple = (10, 8)):
        """
        Create PCA scatter plot

        Args:
            pca_results: Dictionary containing PCA results from analyzer
            figsize: Figure size (width, height)
        """
        pca_df = pca_results['pca_df']
        explained_variance = pca_results['explained_variance']

        fig, ax = plt.subplots(figsize=figsize)

        # Plot by group
        groups = pca_df['Group'].unique()
        group_colors = {'Cancer': '#E74C3C', 'Normal': '#3498DB'}
        group_markers = {'Cancer': 'o', 'Normal': 's'}

        for group in groups:
            mask = pca_df['Group'] == group
            ax.scatter(
                pca_df.loc[mask, 'PC1'],
                pca_df.loc[mask, 'PC2'],
                label=group,
                color=group_colors[group],
                marker=group_markers[group],
                s=100,
                alpha=0.7,
                edgecolors='black',
                linewidths=1
            )

        # Add labels for samples
        for idx, row in pca_df.iterrows():
            ax.annotate(
                idx,
                (row['PC1'], row['PC2']),
                xytext=(5, 5),
                textcoords='offset points',
                fontsize=8,
                alpha=0.7
            )

        ax.set_xlabel(f'PC1 ({explained_variance[0]*100:.2f}% variance)')
        ax.set_ylabel(f'PC2 ({explained_variance[1]*100:.2f}% variance)')
        ax.set_title('PCA: Cancer vs Normal Samples')
        ax.legend(loc='best', frameon=True)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        output_file = self.output_dir / 'pca_plot.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved PCA plot to {output_file}")

        plt.close()

    def plot_pca_by_glycan_type(self, df: pd.DataFrame, pca_results: dict, figsize: tuple = (10, 8)):
        """
        Create PCA biplot colored by glycan type contribution

        Args:
            df: Annotated DataFrame
            pca_results: Dictionary containing PCA results
            figsize: Figure size
        """
        pca_df = pca_results['pca_df']
        explained_variance = pca_results['explained_variance']

        fig, ax = plt.subplots(figsize=figsize)

        # Plot samples by group
        groups = pca_df['Group'].unique()
        group_colors = {'Cancer': '#E74C3C', 'Normal': '#3498DB'}
        group_markers = {'Cancer': 'o', 'Normal': 's'}

        for group in groups:
            mask = pca_df['Group'] == group
            ax.scatter(
                pca_df.loc[mask, 'PC1'],
                pca_df.loc[mask, 'PC2'],
                label=group,
                color=group_colors[group],
                marker=group_markers[group],
                s=100,
                alpha=0.7,
                edgecolors='black',
                linewidths=1
            )

        ax.set_xlabel(f'PC1 ({explained_variance[0]*100:.2f}% variance)')
        ax.set_ylabel(f'PC2 ({explained_variance[1]*100:.2f}% variance)')
        ax.set_title('PCA: Sample Distribution')
        ax.legend(loc='best', frameon=True)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        output_file = self.output_dir / 'pca_samples.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved PCA plot to {output_file}")

        plt.close()

    def plot_boxplot(self, boxplot_data: pd.DataFrame, figsize: tuple = (12, 6)):
        """
        Create boxplot comparing glycan types between groups

        Args:
            boxplot_data: Long-format DataFrame from analyzer
            figsize: Figure size
        """
        fig, ax = plt.subplots(figsize=figsize)

        # Create boxplot
        sns.boxplot(
            data=boxplot_data,
            x='GlycanType',
            y='Intensity',
            hue='Group',
            palette={'Cancer': '#E74C3C', 'Normal': '#3498DB'},
            ax=ax
        )

        ax.set_xlabel('Glycan Type')
        ax.set_ylabel('Log2(Intensity + 1)')
        ax.set_title('Glycan Intensity Distribution by Type and Group')
        ax.legend(title='Group', loc='best')

        plt.tight_layout()

        output_file = self.output_dir / 'boxplot_glycan_types.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved boxplot to {output_file}")

        plt.close()

    def plot_glycan_type_distribution(self, df: pd.DataFrame, figsize: tuple = (10, 6)):
        """
        Create bar plot showing distribution of glycan types

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        fig, ax = plt.subplots(figsize=figsize)

        # Count glycan types
        type_counts = df['GlycanType'].value_counts()

        # Create bar plot with custom colors
        bars = ax.bar(
            range(len(type_counts)),
            type_counts.values,
            color=[self.colors.get(t, '#CCCCCC') for t in type_counts.index]
        )

        ax.set_xticks(range(len(type_counts)))
        ax.set_xticklabels(type_counts.index, rotation=0)
        ax.set_xlabel('Glycan Type')
        ax.set_ylabel('Count')
        ax.set_title('Distribution of Glycan Types')

        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width()/2.,
                height,
                f'{int(height)}',
                ha='center',
                va='bottom'
            )

        plt.tight_layout()

        output_file = self.output_dir / 'glycan_type_distribution.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved glycan type distribution to {output_file}")

        plt.close()

    def plot_heatmap(self, df: pd.DataFrame, figsize: tuple = (14, 10), top_n: int = 50):
        """
        Create heatmap of top glycopeptides

        Args:
            df: Annotated DataFrame
            figsize: Figure size
            top_n: Number of top glycopeptides to show
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType']

        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Get intensity matrix
        intensity_matrix = df[sample_cols].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

        # Calculate mean intensity across all samples
        df['MeanIntensity'] = intensity_matrix.mean(axis=1)

        # Select top N glycopeptides
        top_glycopeptides = df.nlargest(top_n, 'MeanIntensity')

        # Prepare data for heatmap
        heatmap_data = top_glycopeptides[sample_cols].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

        # Log transform
        heatmap_data = np.log2(heatmap_data + 1)

        # Create row labels (Peptide_GlycanComposition_GlycanType)
        row_labels = [
            f"{row['Peptide']}_{row['GlycanComposition']}_{row['GlycanType']}"
            for _, row in top_glycopeptides.iterrows()
        ]

        heatmap_data.index = row_labels

        # Create heatmap
        fig, ax = plt.subplots(figsize=figsize)

        sns.heatmap(
            heatmap_data,
            cmap='YlOrRd',
            cbar_kws={'label': 'Log2(Intensity + 1)'},
            xticklabels=True,
            yticklabels=True,
            ax=ax,
            linewidths=0.5,
            linecolor='lightgray'
        )

        ax.set_xlabel('Sample')
        ax.set_ylabel('Glycopeptide')
        ax.set_title(f'Top {top_n} Glycopeptides Heatmap')

        plt.tight_layout()

        output_file = self.output_dir / 'heatmap_top_glycopeptides.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved heatmap to {output_file}")

        plt.close()

    def plot_all(self, df: pd.DataFrame, pca_results: dict, boxplot_data: pd.DataFrame):
        """
        Generate all plots

        Args:
            df: Annotated DataFrame
            pca_results: PCA results from analyzer
            boxplot_data: Boxplot data from analyzer
        """
        logger.info("Generating all visualizations...")

        self.plot_pca(pca_results)
        self.plot_pca_by_glycan_type(df, pca_results)
        self.plot_boxplot(boxplot_data)
        self.plot_glycan_type_distribution(df)
        self.plot_heatmap(df)

        logger.info(f"All visualizations saved to {self.output_dir}")


if __name__ == "__main__":
    # Test
    pass
