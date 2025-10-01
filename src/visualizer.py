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
from scipy import stats
from matplotlib.patches import Ellipse
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform

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

    def _draw_confidence_ellipse(self, ax, x, y, color, alpha=0.2, n_std=1.96):
        """
        Draw 95% confidence ellipse for a group of points

        Args:
            ax: Matplotlib axis
            x, y: Data points
            color: Ellipse color
            alpha: Transparency
            n_std: Number of standard deviations (1.96 for 95% confidence)
        """
        if len(x) < 3:
            return

        # Calculate covariance matrix
        cov = np.cov(x, y)

        # Calculate eigenvalues and eigenvectors
        eigenvalues, eigenvectors = np.linalg.eig(cov)

        # Get the largest eigenvalue and its corresponding eigenvector
        order = eigenvalues.argsort()[::-1]
        eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]

        # Calculate ellipse angle
        angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))

        # Width and height are "full" widths, not radius
        width, height = 2 * n_std * np.sqrt(eigenvalues)

        # Draw ellipse
        ellipse = Ellipse(
            xy=(np.mean(x), np.mean(y)),
            width=width,
            height=height,
            angle=angle,
            facecolor=color,
            edgecolor=color,
            alpha=alpha,
            linewidth=2,
            linestyle='--'
        )
        ax.add_patch(ellipse)

    def plot_pca(self, pca_results: dict, figsize: tuple = (10, 8)):
        """
        Create PCA scatter plot with 95% confidence ellipses

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
            pc1_vals = pca_df.loc[mask, 'PC1'].values
            pc2_vals = pca_df.loc[mask, 'PC2'].values

            # Draw 95% confidence ellipse
            self._draw_confidence_ellipse(
                ax,
                pc1_vals,
                pc2_vals,
                group_colors[group],
                alpha=0.2,
                n_std=1.96  # 95% confidence
            )

            # Plot scatter points
            ax.scatter(
                pc1_vals,
                pc2_vals,
                label=group,
                color=group_colors[group],
                marker=group_markers[group],
                s=100,
                alpha=0.7,
                edgecolors='black',
                linewidths=1,
                zorder=3
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
        ax.set_title('PCA: Cancer vs Normal Samples (95% Confidence Ellipses)')
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
        Create boxplot comparing glycan types between groups with statistical significance

        Args:
            boxplot_data: Long-format DataFrame from analyzer
            figsize: Figure size
        """
        fig, ax = plt.subplots(figsize=figsize)

        # Define fixed order for glycan types
        glycan_order = ['Non', 'Sialylated', 'Fucosylated', 'Both']

        # Filter to only include existing glycan types
        existing_types = [gt for gt in glycan_order if gt in boxplot_data['GlycanType'].unique()]

        # Create boxplot with ordered categories
        sns.boxplot(
            data=boxplot_data,
            x='GlycanType',
            y='Intensity',
            hue='Group',
            order=existing_types,
            palette={'Cancer': '#E74C3C', 'Normal': '#3498DB'},
            ax=ax
        )

        # Perform statistical tests and add significance markers
        y_max = boxplot_data['Intensity'].max()
        y_range = boxplot_data['Intensity'].max() - boxplot_data['Intensity'].min()

        for i, glycan_type in enumerate(existing_types):
            # Get data for each group
            cancer_data = boxplot_data[
                (boxplot_data['GlycanType'] == glycan_type) &
                (boxplot_data['Group'] == 'Cancer')
            ]['Intensity'].values

            normal_data = boxplot_data[
                (boxplot_data['GlycanType'] == glycan_type) &
                (boxplot_data['Group'] == 'Normal')
            ]['Intensity'].values

            # Skip if either group has insufficient data
            if len(cancer_data) < 3 or len(normal_data) < 3:
                continue

            # Perform Mann-Whitney U test (non-parametric)
            try:
                statistic, p_value = stats.mannwhitneyu(cancer_data, normal_data, alternative='two-sided')

                # Determine significance level
                if p_value < 0.001:
                    sig_marker = '***'
                elif p_value < 0.01:
                    sig_marker = '**'
                elif p_value < 0.05:
                    sig_marker = '*'
                else:
                    sig_marker = 'ns'

                # Add significance marker if significant
                if sig_marker != 'ns':
                    # Position marker above the boxes
                    y_position = y_max + y_range * 0.05 * (i % 2 + 1)
                    ax.text(
                        i,
                        y_position,
                        sig_marker,
                        ha='center',
                        va='bottom',
                        fontsize=14,
                        fontweight='bold'
                    )

                    logger.info(f"{glycan_type}: p={p_value:.4f} ({sig_marker})")

            except Exception as e:
                logger.warning(f"Statistical test failed for {glycan_type}: {str(e)}")

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

    def plot_heatmap(self, df: pd.DataFrame, figsize: tuple = (16, 12), top_n: int = 50):
        """
        Create clustered heatmap of top glycopeptides with hierarchical clustering

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
        df_copy = df.copy()
        df_copy['MeanIntensity'] = intensity_matrix.mean(axis=1)

        # Select top N glycopeptides
        top_glycopeptides = df_copy.nlargest(top_n, 'MeanIntensity')

        # Prepare data for heatmap
        heatmap_data = top_glycopeptides[sample_cols].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

        # Log transform
        heatmap_data = np.log2(heatmap_data + 1)

        # Transpose for sample clustering (samples as rows)
        heatmap_data_t = heatmap_data.T

        # Create row labels (Peptide_GlycanComposition_GlycanType)
        row_labels = [
            f"{row['Peptide']}_{row['GlycanComposition']}_{row['GlycanType']}"
            for _, row in top_glycopeptides.iterrows()
        ]

        heatmap_data.index = row_labels

        # Create color annotation for sample groups
        sample_colors = []
        for col in heatmap_data.columns:
            if col.startswith('C'):
                sample_colors.append('#E74C3C')  # Cancer - red
            else:
                sample_colors.append('#3498DB')  # Normal - blue

        # Create clustermap with hierarchical clustering on both axes
        g = sns.clustermap(
            heatmap_data,
            cmap='YlOrRd',
            figsize=figsize,
            dendrogram_ratio=0.15,
            cbar_pos=(0.02, 0.83, 0.03, 0.15),
            cbar_kws={'label': 'Log2(Intensity + 1)'},
            col_cluster=True,  # Cluster samples
            row_cluster=True,  # Cluster glycopeptides
            xticklabels=True,
            yticklabels=True,
            linewidths=0.5,
            linecolor='lightgray',
            col_colors=sample_colors,
            method='average',  # Linkage method
            metric='euclidean'  # Distance metric
        )

        # Adjust labels
        g.ax_heatmap.set_xlabel('Sample', fontsize=12)
        g.ax_heatmap.set_ylabel('Glycopeptide', fontsize=12)
        g.fig.suptitle(f'Top {top_n} Glycopeptides Heatmap with Hierarchical Clustering',
                       fontsize=14, y=0.98)

        # Rotate sample labels
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

        # Add legend for sample colors
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#E74C3C', label='Cancer'),
            Patch(facecolor='#3498DB', label='Normal')
        ]
        g.ax_heatmap.legend(
            handles=legend_elements,
            loc='upper left',
            bbox_to_anchor=(1.2, 1),
            frameon=True
        )

        output_file = self.output_dir / 'heatmap_top_glycopeptides.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved clustered heatmap to {output_file}")

        plt.close()

    def plot_heatmap_full_profile(self, df: pd.DataFrame, figsize: tuple = (18, 14)):
        """
        Create clustered heatmap of ALL glycopeptides (full glycan profile per sample)

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType']

        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Get intensity matrix for ALL glycopeptides
        intensity_matrix = df[sample_cols].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

        # Filter out glycopeptides with all zeros (not detected in any sample)
        row_sums = intensity_matrix.sum(axis=1)
        df_filtered = df[row_sums > 0].copy()
        heatmap_data = intensity_matrix[row_sums > 0].copy()

        logger.info(f"Full profile heatmap: {len(heatmap_data)} glycopeptides across {len(sample_cols)} samples")

        # Log transform
        heatmap_data = np.log2(heatmap_data + 1)

        # Create row labels (Peptide_GlycanComposition_GlycanType)
        row_labels = [
            f"{row['Peptide']}_{row['GlycanComposition']}_{row['GlycanType']}"
            for _, row in df_filtered.iterrows()
        ]

        heatmap_data.index = row_labels

        # Create color annotation for sample groups
        sample_colors = []
        for col in heatmap_data.columns:
            if col.startswith('C'):
                sample_colors.append('#E74C3C')  # Cancer - red
            else:
                sample_colors.append('#3498DB')  # Normal - blue

        # Create clustermap with hierarchical clustering on both axes
        g = sns.clustermap(
            heatmap_data,
            cmap='YlOrRd',
            figsize=figsize,
            dendrogram_ratio=0.1,
            cbar_pos=(0.02, 0.85, 0.03, 0.12),
            cbar_kws={'label': 'Log2(Intensity + 1)'},
            col_cluster=True,  # Cluster samples
            row_cluster=True,  # Cluster glycopeptides
            xticklabels=True,
            yticklabels=False,  # Too many glycopeptides to show labels
            linewidths=0,
            col_colors=sample_colors,
            method='average',  # Linkage method
            metric='euclidean'  # Distance metric
        )

        # Adjust labels
        g.ax_heatmap.set_xlabel('Sample', fontsize=12)
        g.ax_heatmap.set_ylabel('Glycopeptide', fontsize=12)
        g.fig.suptitle(f'Full Glycan Profile Heatmap ({len(heatmap_data)} glycopeptides)',
                       fontsize=14, y=0.98)

        # Rotate sample labels
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

        # Add legend for sample colors
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#E74C3C', label='Cancer'),
            Patch(facecolor='#3498DB', label='Normal')
        ]
        g.ax_heatmap.legend(
            handles=legend_elements,
            loc='upper left',
            bbox_to_anchor=(1.15, 1),
            frameon=True
        )

        output_file = self.output_dir / 'heatmap_full_glycan_profile.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved full glycan profile heatmap to {output_file}")

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
        self.plot_heatmap_full_profile(df)

        logger.info(f"All visualizations saved to {self.output_dir}")


if __name__ == "__main__":
    # Test
    pass
