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
            'Sialylated': "#FF00EA",
            'Fucosylated': "#DB3434",
            'Both': "#FFB700"
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
            hue_order=['Normal', 'Cancer'],
            palette={'Cancer': '#E74C3C', 'Normal': '#3498DB'},
            width=0.6,
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
                    # Position marker inside the plot area, near the top
                    y_position = y_max - y_range * 0.02 - (y_range * 0.03 * (i % 2))
                    ax.text(
                        i,
                        y_position,
                        sig_marker,
                        ha='center',
                        va='top',
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

    def plot_boxplot_extended(self, boxplot_data: pd.DataFrame, figsize: tuple = (14, 6)):
        """
        Create extended boxplot comparing 5 glycan categories between groups

        Args:
            boxplot_data: Long-format DataFrame from analyzer (extended)
            figsize: Figure size
        """
        fig, ax = plt.subplots(figsize=figsize)

        # Define fixed order for extended categories
        category_order = ['HM', 'C/H', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        # Filter to only include existing categories
        existing_categories = [cat for cat in category_order if cat in boxplot_data['ExtendedCategory'].unique()]

        # Create boxplot with ordered categories
        sns.boxplot(
            data=boxplot_data,
            x='ExtendedCategory',
            y='Intensity',
            hue='Group',
            order=existing_categories,
            hue_order=['Normal', 'Cancer'],
            palette={'Cancer': '#E74C3C', 'Normal': '#3498DB'},
            width=0.6,
            ax=ax
        )

        # Perform statistical tests and add significance markers
        y_max = boxplot_data['Intensity'].max()
        y_range = boxplot_data['Intensity'].max() - boxplot_data['Intensity'].min()

        for i, category in enumerate(existing_categories):
            # Get data for each group
            cancer_data = boxplot_data[
                (boxplot_data['ExtendedCategory'] == category) &
                (boxplot_data['Group'] == 'Cancer')
            ]['Intensity'].values

            normal_data = boxplot_data[
                (boxplot_data['ExtendedCategory'] == category) &
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
                    # Position marker inside the plot area, near the top
                    y_position = y_max - y_range * 0.02 - (y_range * 0.03 * (i % 2))
                    ax.text(
                        i,
                        y_position,
                        sig_marker,
                        ha='center',
                        va='top',
                        fontsize=14,
                        fontweight='bold'
                    )

                    logger.info(f"{category}: p={p_value:.4f} ({sig_marker})")

            except Exception as e:
                logger.warning(f"Statistical test failed for {category}: {str(e)}")

        ax.set_xlabel('Glycan Category', fontsize=12)
        ax.set_ylabel('Log2(Intensity + 1)', fontsize=12)
        ax.set_title('Extended Glycan Category Distribution (HM, C/H, Fucosylated, Sialylated, Sialofucosylated)', fontsize=14)
        ax.legend(title='Group', loc='upper left', bbox_to_anchor=(1, 1), frameon=True)

        plt.tight_layout()

        output_file = self.output_dir / 'boxplot_extended_categories.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved extended boxplot to {output_file}")

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
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid']

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
        
        # Position title at the top with more space to avoid overlap
        g.fig.suptitle(f'Top {top_n} Glycopeptides Heatmap with Hierarchical Clustering',
                       fontsize=14, y=1.00, fontweight='bold')

        # Rotate sample labels
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

        # Add legend for sample colors - position to avoid clustering line
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#E74C3C', label='Cancer'),
            Patch(facecolor='#3498DB', label='Normal')
        ]
        g.ax_heatmap.legend(
            handles=legend_elements,
            loc='lower left',
            bbox_to_anchor=(1.05, 0.5),
            frameon=True,
            fontsize=10
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
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid']

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
        
        # Position title at the top with more space to avoid overlap
        g.fig.suptitle(f'Full Glycan Profile Heatmap ({len(heatmap_data)} glycopeptides)',
                       fontsize=14, y=1.00, fontweight='bold')

        # Rotate sample labels
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

        # Add legend for sample colors - position to avoid clustering line
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#E74C3C', label='Cancer'),
            Patch(facecolor='#3498DB', label='Normal')
        ]
        g.ax_heatmap.legend(
            handles=legend_elements,
            loc='lower left',
            bbox_to_anchor=(1.05, 0.5),
            frameon=True,
            fontsize=10
        )

        output_file = self.output_dir / 'heatmap_full_glycan_profile.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved full glycan profile heatmap to {output_file}")

        plt.close()

    def plot_histogram_by_sample(self, df: pd.DataFrame, figsize: tuple = (20, 12)):
        """
        Create histogram showing glycan type intensities per sample (original data)

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid']

        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Define glycan categories to display
        # Priority order: High mannose, C/H, Sialylated, Fucosylated, Both
        categories = []

        # Calculate total intensity per sample per category
        data_for_plot = []

        for sample in sample_cols:
            sample_data = {'Sample': sample}

            # Get intensity matrix
            intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

            # High mannose (priority 1)
            high_mannose_mask = df['IsHighMannose']
            sample_data['High mannose'] = intensity_col[high_mannose_mask].sum()

            # C/H (priority 2)
            ch_mask = df['IsComplexHybrid']
            sample_data['C/H'] = intensity_col[ch_mask].sum()

            # Both (Sialylated AND Fucosylated) (priority 3)
            both_mask = df['IsSialylated'] & df['IsFucosylated']
            sample_data['Both'] = intensity_col[both_mask].sum()

            # Sialylated only (priority 4)
            sia_only_mask = df['IsSialylated'] & ~df['IsFucosylated'] & ~df['IsHighMannose'] & ~df['IsComplexHybrid']
            sample_data['Sialylated'] = intensity_col[sia_only_mask].sum()

            # Fucosylated only (priority 5)
            fuc_only_mask = df['IsFucosylated'] & ~df['IsSialylated'] & ~df['IsHighMannose'] & ~df['IsComplexHybrid']
            sample_data['Fucosylated'] = intensity_col[fuc_only_mask].sum()

            data_for_plot.append(sample_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot)
        plot_df = plot_df.set_index('Sample')

        # Reorder columns
        column_order = ['High mannose', 'C/H', 'Fucosylated', 'Sialylated', 'Both']
        plot_df = plot_df[column_order]

        # Create stacked bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Define colors for each category
        colors = {
            'High mannose': '#2ECC71',  # Green
            'C/H': "#1500FF",           # Blue
            'Fucosylated': "#FF0000",   # Red
            'Sialylated': "#FF00FB",    # Pink
            'Both': "#FF9D00"           # Orange
        }

        plot_df.plot(
            kind='bar',
            stacked=True,
            ax=ax,
            color=[colors[col] for col in column_order],
            width=0.8,
            edgecolor='black',
            linewidth=0.5
        )

        ax.set_xlabel('Sample', fontsize=12)
        ax.set_ylabel('Total Signal Intensity', fontsize=12)
        ax.set_title('Glycan Type Distribution by Sample (Original Data)', fontsize=14)
        ax.legend(title='Glycan Type', loc='upper left', bbox_to_anchor=(1, 1), frameon=True)

        # Rotate x-axis labels
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=90, ha='right')

        # Use scientific notation for y-axis
        ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

        plt.tight_layout()

        output_file = self.output_dir / 'histogram_glycan_types_by_sample.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved histogram to {output_file}")

        plt.close()

    def plot_histogram_normalized(self, df: pd.DataFrame, figsize: tuple = (20, 12)):
        """
        Create histogram showing glycan type intensities per sample (TIC normalized data)

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid']

        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Get intensity matrix
        intensity_matrix = df[sample_cols].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

        # Apply TIC normalization (same as analyzer.py)
        sample_sums = intensity_matrix.sum(axis=0)
        median_sum = sample_sums.median()
        sample_sums_safe = sample_sums.replace(0, 1)
        intensity_normalized = intensity_matrix / sample_sums_safe * median_sum

        # Calculate total intensity per sample per category (normalized)
        data_for_plot = []

        for i, sample in enumerate(sample_cols):
            sample_data = {'Sample': sample}

            # Get normalized intensity column
            intensity_col = intensity_normalized.iloc[:, i]

            # High mannose (priority 1)
            high_mannose_mask = df['IsHighMannose']
            sample_data['High mannose'] = intensity_col[high_mannose_mask].sum()

            # C/H (priority 2)
            ch_mask = df['IsComplexHybrid']
            sample_data['C/H'] = intensity_col[ch_mask].sum()

            # Both (Sialylated AND Fucosylated) (priority 3)
            both_mask = df['IsSialylated'] & df['IsFucosylated']
            sample_data['Both'] = intensity_col[both_mask].sum()

            # Sialylated only (priority 4)
            sia_only_mask = df['IsSialylated'] & ~df['IsFucosylated'] & ~df['IsHighMannose'] & ~df['IsComplexHybrid']
            sample_data['Sialylated'] = intensity_col[sia_only_mask].sum()

            # Fucosylated only (priority 5)
            fuc_only_mask = df['IsFucosylated'] & ~df['IsSialylated'] & ~df['IsHighMannose'] & ~df['IsComplexHybrid']
            sample_data['Fucosylated'] = intensity_col[fuc_only_mask].sum()

            data_for_plot.append(sample_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot)
        plot_df = plot_df.set_index('Sample')

        # Reorder columns
        column_order = ['High mannose', 'C/H', 'Fucosylated', 'Sialylated', 'Both']
        plot_df = plot_df[column_order]

        # Create stacked bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Define colors for each category
        colors = {
            'High mannose': '#2ECC71',  # Green
            'C/H': "#1500FF",           # Blue
            'Fucosylated': "#FF0000",   # Red
            'Sialylated': "#FF00FB",    # Pink
            'Both': "#FF9D00"           # Orange
        }

        plot_df.plot(
            kind='bar',
            stacked=True,
            ax=ax,
            color=[colors[col] for col in column_order],
            width=0.8,
            edgecolor='black',
            linewidth=0.5
        )

        ax.set_xlabel('Sample', fontsize=12)
        ax.set_ylabel('Total Signal Intensity (TIC Normalized)', fontsize=12)
        ax.set_title('Glycan Type Distribution by Sample (TIC Normalized Data)', fontsize=14, fontweight='bold')
        ax.legend(title='Glycan Type', loc='upper left', bbox_to_anchor=(1, 1), frameon=True)

        # Rotate x-axis labels
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=90, ha='right')

        # Use scientific notation for y-axis
        ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

        plt.tight_layout()

        output_file = self.output_dir / 'histogram_glycan_types_by_sample_normalized.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved normalized histogram to {output_file}")

        plt.close()

    def plot_histogram_primary_classification(self, df: pd.DataFrame, normalization: str = 'raw',
                                               figsize: tuple = (20, 10)):
        """
        Create histogram for primary classification (Truncated, High Mannose, ComplexHybrid)
        across all samples

        Args:
            df: Annotated DataFrame
            normalization: 'raw' (normalize raw data then sum) or 'aggregated' (sum then normalize)
            figsize: Figure size
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']

        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Get intensity matrix
        intensity_matrix = df[sample_cols].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

        # Primary classification categories (exclude Outlier for visualization)
        primary_categories = ['Truncated', 'High Mannose', 'ComplexHybrid']

        data_for_plot = []

        for sample in sample_cols:
            sample_data = {'Sample': sample}

            if normalization == 'raw':
                # Normalize raw data using min-max scaling
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

                # Min-max normalization
                min_val = intensity_col.min()
                max_val = intensity_col.max()
                if max_val > min_val:
                    intensity_col = (intensity_col - min_val) / (max_val - min_val)
                else:
                    intensity_col = intensity_col * 0  # All zeros if no variation
            else:
                # Use raw intensities for aggregation
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

            # Sum by primary classification
            for category in primary_categories:
                mask = df['PrimaryClassification'] == category
                sample_data[category] = intensity_col[mask].sum()

            data_for_plot.append(sample_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot)
        plot_df = plot_df.set_index('Sample')

        # Apply min-max normalization after aggregation if needed
        if normalization == 'aggregated':
            for col in primary_categories:
                min_val = plot_df[col].min()
                max_val = plot_df[col].max()
                if max_val > min_val:
                    plot_df[col] = (plot_df[col] - min_val) / (max_val - min_val)

        # Create grouped bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Define colors
        colors_primary = {
            'Truncated': '#CCCCCC',
            'High Mannose': '#2ECC71',
            'ComplexHybrid': '#3498DB'
        }

        plot_df[primary_categories].plot(
            kind='bar',
            ax=ax,
            color=[colors_primary[cat] for cat in primary_categories],
            width=0.8,
            edgecolor='black',
            linewidth=0.5
        )

        norm_text = 'Raw Normalized then Summed' if normalization == 'raw' else 'Summed then Normalized'
        ax.set_xlabel('Sample', fontsize=12)
        ax.set_ylabel('Normalized Intensity Sum', fontsize=12)
        ax.set_title(f'Primary Classification Distribution by Sample ({norm_text})', fontsize=14)
        ax.legend(title='Primary Classification', loc='upper left', bbox_to_anchor=(1, 1), frameon=True)

        plt.setp(ax.xaxis.get_majorticklabels(), rotation=90, ha='right')
        plt.tight_layout()

        output_file = self.output_dir / f'histogram_primary_{normalization}_normalized.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved primary classification histogram to {output_file}")

        plt.close()

    def plot_histogram_secondary_classification(self, df: pd.DataFrame, normalization: str = 'raw',
                                                 figsize: tuple = (20, 10)):
        """
        Create histogram for secondary classification across all samples
        (excludes Truncated and Outlier)

        Args:
            df: Annotated DataFrame
            normalization: 'raw' (normalize raw data then sum) or 'aggregated' (sum then normalize)
            figsize: Figure size
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']

        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Secondary classification categories (5 categories, exclude Truncated and Outlier)
        secondary_categories = ['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        data_for_plot = []

        for sample in sample_cols:
            sample_data = {'Sample': sample}

            if normalization == 'raw':
                # Normalize raw data using min-max scaling
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

                # Min-max normalization
                min_val = intensity_col.min()
                max_val = intensity_col.max()
                if max_val > min_val:
                    intensity_col = (intensity_col - min_val) / (max_val - min_val)
                else:
                    intensity_col = intensity_col * 0
            else:
                # Use raw intensities for aggregation
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

            # Sum by secondary classification
            for category in secondary_categories:
                mask = df['SecondaryClassification'] == category
                sample_data[category] = intensity_col[mask].sum()

            data_for_plot.append(sample_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot)
        plot_df = plot_df.set_index('Sample')

        # Apply min-max normalization after aggregation if needed
        if normalization == 'aggregated':
            for col in secondary_categories:
                min_val = plot_df[col].min()
                max_val = plot_df[col].max()
                if max_val > min_val:
                    plot_df[col] = (plot_df[col] - min_val) / (max_val - min_val)

        # Create grouped bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Define colors
        colors_secondary = {
            'High Mannose': '#2ECC71',
            'Complex/Hybrid': '#3498DB',
            'Fucosylated': '#E74C3C',
            'Sialylated': '#9B59B6',
            'Sialofucosylated': '#F39C12'
        }

        plot_df[secondary_categories].plot(
            kind='bar',
            ax=ax,
            color=[colors_secondary[cat] for cat in secondary_categories],
            width=0.8,
            edgecolor='black',
            linewidth=0.5
        )

        norm_text = 'Raw Normalized then Summed' if normalization == 'raw' else 'Summed then Normalized'
        ax.set_xlabel('Sample', fontsize=12)
        ax.set_ylabel('Normalized Intensity Sum', fontsize=12)
        ax.set_title(f'Secondary Classification Distribution by Sample ({norm_text})', fontsize=14)
        ax.legend(title='Secondary Classification', loc='upper left', bbox_to_anchor=(1, 1), frameon=True)

        plt.setp(ax.xaxis.get_majorticklabels(), rotation=90, ha='right')
        plt.tight_layout()

        output_file = self.output_dir / f'histogram_secondary_{normalization}_normalized.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved secondary classification histogram to {output_file}")

        plt.close()

    def plot_histogram_cancer_vs_normal_primary(self, df: pd.DataFrame, figsize: tuple = (12, 8)):
        """
        Create histogram comparing Cancer vs Normal for primary classification
        (Aggregated with Log2 scaling)

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']

        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Separate Cancer and Normal samples
        cancer_samples = [col for col in sample_cols if col.startswith('C')]
        normal_samples = [col for col in sample_cols if col.startswith('N')]

        # Primary classification categories (exclude Outlier)
        primary_categories = ['Truncated', 'High Mannose', 'ComplexHybrid']

        # Calculate sums for each group
        data_for_plot = []

        for group_name, group_samples in [('Cancer', cancer_samples), ('Normal', normal_samples)]:
            group_data = {'Group': group_name}

            # Get intensity matrix for this group
            intensity_matrix = df[group_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

            # Sum across all samples in the group, then by category
            for category in primary_categories:
                mask = df['PrimaryClassification'] == category
                group_data[category] = intensity_matrix[mask].sum().sum()

            data_for_plot.append(group_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot)
        plot_df = plot_df.set_index('Group')

        # Apply Log2 transformation (no min-max normalization)
        for col in primary_categories:
            plot_df[col] = np.log2(plot_df[col] + 1)

        # Create grouped bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Transpose for proper grouping (categories as x-axis)
        plot_df_t = plot_df.T

        # Define colors
        group_colors = {'Cancer': '#E74C3C', 'Normal': '#3498DB'}

        plot_df_t.plot(
            kind='bar',
            ax=ax,
            color=[group_colors[g] for g in plot_df_t.columns],
            width=0.7,
            edgecolor='black',
            linewidth=1
        )

        ax.set_xlabel('Primary Classification', fontsize=12)
        ax.set_ylabel('Log2(Total Intensity + 1)', fontsize=12)
        ax.set_title('Primary Classification: Cancer vs Normal (Log2 Scaled)', fontsize=14)
        ax.legend(title='Group', loc='upper right', frameon=True)

        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
        plt.tight_layout()

        output_file = self.output_dir / 'histogram_primary_cancer_vs_normal.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved primary Cancer vs Normal histogram to {output_file}")

        plt.close()

    def plot_histogram_cancer_vs_normal_secondary(self, df: pd.DataFrame, figsize: tuple = (12, 8)):
        """
        Create histogram comparing Cancer vs Normal for secondary classification
        (Aggregated with Log2 scaling)

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']

        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Separate Cancer and Normal samples
        cancer_samples = [col for col in sample_cols if col.startswith('C')]
        normal_samples = [col for col in sample_cols if col.startswith('N')]

        # Secondary classification categories (5 categories)
        secondary_categories = ['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        # Calculate sums for each group
        data_for_plot = []

        for group_name, group_samples in [('Cancer', cancer_samples), ('Normal', normal_samples)]:
            group_data = {'Group': group_name}

            # Get intensity matrix for this group
            intensity_matrix = df[group_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

            # Sum across all samples in the group, then by category
            for category in secondary_categories:
                mask = df['SecondaryClassification'] == category
                group_data[category] = intensity_matrix[mask].sum().sum()

            data_for_plot.append(group_data)

        # Create DataFrame
        plot_df = pd.DataFrame(data_for_plot)
        plot_df = plot_df.set_index('Group')

        # Apply Log2 transformation (no min-max normalization)
        for col in secondary_categories:
            plot_df[col] = np.log2(plot_df[col] + 1)

        # Create grouped bar plot
        fig, ax = plt.subplots(figsize=figsize)

        # Transpose for proper grouping (categories as x-axis)
        plot_df_t = plot_df.T

        # Define colors
        group_colors = {'Cancer': '#E74C3C', 'Normal': '#3498DB'}

        plot_df_t.plot(
            kind='bar',
            ax=ax,
            color=[group_colors[g] for g in plot_df_t.columns],
            width=0.7,
            edgecolor='black',
            linewidth=1
        )

        ax.set_xlabel('Secondary Classification', fontsize=12)
        ax.set_ylabel('Log2(Total Intensity + 1)', fontsize=12)
        ax.set_title('Secondary Classification: Cancer vs Normal (Log2 Scaled)', fontsize=14)
        ax.legend(title='Group', loc='upper right', frameon=True)

        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
        plt.tight_layout()

        output_file = self.output_dir / 'histogram_secondary_cancer_vs_normal.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved secondary Cancer vs Normal histogram to {output_file}")

        plt.close()

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

    def plot_boxplot_primary_classification(self, df: pd.DataFrame, normalization: str = 'raw', figsize: tuple = (12, 8)):
        """
        Create boxplot for primary classification

        Args:
            df: Annotated DataFrame
            normalization: 'raw' or 'aggregated'
            figsize: Figure size
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']
        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Primary classification categories
        primary_categories = ['Truncated', 'High Mannose', 'ComplexHybrid']

        # Prepare long format data
        data_for_plot = []

        for sample in sample_cols:
            if normalization == 'raw':
                # Min-max normalize per sample
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)
                min_val, max_val = intensity_col.min(), intensity_col.max()
                if max_val > min_val:
                    intensity_col = (intensity_col - min_val) / (max_val - min_val)
            else:
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

            group = 'Cancer' if sample.startswith('C') else 'Normal'

            for idx, row in df.iterrows():
                if row['PrimaryClassification'] in primary_categories:
                    data_for_plot.append({
                        'Sample': sample,
                        'Group': group,
                        'Classification': row['PrimaryClassification'],
                        'Intensity': intensity_col.iloc[idx] if normalization == 'raw' else intensity_col.iloc[idx]
                    })

        plot_df = pd.DataFrame(data_for_plot)

        # Apply aggregated normalization if needed
        if normalization == 'aggregated':
            plot_df['Intensity'] = np.log2(plot_df['Intensity'] + 1)

        fig, ax = plt.subplots(figsize=figsize)

        sns.boxplot(data=plot_df, x='Classification', y='Intensity', hue='Group',
                   palette={'Cancer': '#E74C3C', 'Normal': '#3498DB'}, ax=ax)

        norm_text = 'Raw Normalized' if normalization == 'raw' else 'Log2 Scaled'
        ax.set_xlabel('Primary Classification', fontsize=12)
        ax.set_ylabel(f'Intensity ({norm_text})', fontsize=12)
        ax.set_title(f'Primary Classification Distribution ({norm_text})', fontsize=14)
        ax.legend(title='Group', loc='upper right')

        plt.tight_layout()

        output_file = self.output_dir / f'boxplot_primary_{normalization}_normalized.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved primary boxplot to {output_file}")

        plt.close()

    def plot_boxplot_secondary_classification(self, df: pd.DataFrame, normalization: str = 'raw', figsize: tuple = (14, 8)):
        """
        Create boxplot for secondary classification

        Args:
            df: Annotated DataFrame
            normalization: 'raw' or 'aggregated'
            figsize: Figure size
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']
        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Secondary classification categories
        secondary_categories = ['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        # Prepare long format data
        data_for_plot = []

        for sample in sample_cols:
            if normalization == 'raw':
                # Min-max normalize per sample
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)
                min_val, max_val = intensity_col.min(), intensity_col.max()
                if max_val > min_val:
                    intensity_col = (intensity_col - min_val) / (max_val - min_val)
            else:
                intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

            group = 'Cancer' if sample.startswith('C') else 'Normal'

            for idx, row in df.iterrows():
                if row['SecondaryClassification'] in secondary_categories:
                    data_for_plot.append({
                        'Sample': sample,
                        'Group': group,
                        'Classification': row['SecondaryClassification'],
                        'Intensity': intensity_col.iloc[idx] if normalization == 'raw' else intensity_col.iloc[idx]
                    })

        plot_df = pd.DataFrame(data_for_plot)

        # Apply aggregated normalization if needed
        if normalization == 'aggregated':
            plot_df['Intensity'] = np.log2(plot_df['Intensity'] + 1)

        fig, ax = plt.subplots(figsize=figsize)

        sns.boxplot(data=plot_df, x='Classification', y='Intensity', hue='Group',
                   order=secondary_categories,
                   palette={'Cancer': '#E74C3C', 'Normal': '#3498DB'}, ax=ax)

        norm_text = 'Raw Normalized' if normalization == 'raw' else 'Log2 Scaled'
        ax.set_xlabel('Secondary Classification', fontsize=12)
        ax.set_ylabel(f'Intensity ({norm_text})', fontsize=12)
        ax.set_title(f'Secondary Classification Distribution ({norm_text})', fontsize=14)
        ax.legend(title='Group', loc='upper right')
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

        plt.tight_layout()

        output_file = self.output_dir / f'boxplot_secondary_{normalization}_normalized.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved secondary boxplot to {output_file}")

        plt.close()

    def plot_boxplot_cancer_vs_normal_primary(self, df: pd.DataFrame, figsize: tuple = (10, 6)):
        """
        Create boxplot for Cancer vs Normal comparison (Primary classification)

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']
        sample_cols = [col for col in df.columns if col not in metadata_cols]

        primary_categories = ['Truncated', 'High Mannose', 'ComplexHybrid']

        # Prepare data
        data_for_plot = []
        for sample in sample_cols:
            intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)
            group = 'Cancer' if sample.startswith('C') else 'Normal'

            for idx, row in df.iterrows():
                if row['PrimaryClassification'] in primary_categories:
                    data_for_plot.append({
                        'Group': group,
                        'Classification': row['PrimaryClassification'],
                        'Intensity': np.log2(intensity_col.iloc[idx] + 1)
                    })

        plot_df = pd.DataFrame(data_for_plot)

        fig, ax = plt.subplots(figsize=figsize)

        sns.boxplot(data=plot_df, x='Classification', y='Intensity', hue='Group',
                   palette={'Cancer': '#E74C3C', 'Normal': '#3498DB'}, ax=ax)

        ax.set_xlabel('Primary Classification', fontsize=12)
        ax.set_ylabel('Log2(Intensity + 1)', fontsize=12)
        ax.set_title('Primary Classification: Cancer vs Normal', fontsize=14)
        ax.legend(title='Group')

        plt.tight_layout()

        output_file = self.output_dir / 'boxplot_primary_cancer_vs_normal.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved primary Cancer vs Normal boxplot to {output_file}")

        plt.close()

    def plot_boxplot_cancer_vs_normal_secondary(self, df: pd.DataFrame, figsize: tuple = (12, 6)):
        """
        Create boxplot for Cancer vs Normal comparison (Secondary classification)

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid', 'N_count',
                        'PrimaryClassification', 'SecondaryClassification']
        sample_cols = [col for col in df.columns if col not in metadata_cols]

        secondary_categories = ['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        # Prepare data
        data_for_plot = []
        for sample in sample_cols:
            intensity_col = df[sample].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)
            group = 'Cancer' if sample.startswith('C') else 'Normal'

            for idx, row in df.iterrows():
                if row['SecondaryClassification'] in secondary_categories:
                    data_for_plot.append({
                        'Group': group,
                        'Classification': row['SecondaryClassification'],
                        'Intensity': np.log2(intensity_col.iloc[idx] + 1)
                    })

        plot_df = pd.DataFrame(data_for_plot)

        fig, ax = plt.subplots(figsize=figsize)

        sns.boxplot(data=plot_df, x='Classification', y='Intensity', hue='Group',
                   order=secondary_categories,
                   palette={'Cancer': '#E74C3C', 'Normal': '#3498DB'}, ax=ax)

        ax.set_xlabel('Secondary Classification', fontsize=12)
        ax.set_ylabel('Log2(Intensity + 1)', fontsize=12)
        ax.set_title('Secondary Classification: Cancer vs Normal', fontsize=14)
        ax.legend(title='Group')
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

        plt.tight_layout()

        output_file = self.output_dir / 'boxplot_secondary_cancer_vs_normal.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved secondary Cancer vs Normal boxplot to {output_file}")

        plt.close()

    def plot_all(self, df: pd.DataFrame, pca_results: dict, boxplot_data: pd.DataFrame, boxplot_data_extended: pd.DataFrame = None):
        """
        Generate all plots

        Args:
            df: Annotated DataFrame
            pca_results: PCA results from analyzer
            boxplot_data: Boxplot data from analyzer
            boxplot_data_extended: Extended boxplot data from analyzer (optional)
        """
        logger.info("Generating all visualizations...")

        self.plot_pca(pca_results)
        self.plot_pca_by_glycan_type(df, pca_results)
        self.plot_boxplot(boxplot_data)

        # Plot extended boxplot if data is provided
        if boxplot_data_extended is not None:
            self.plot_boxplot_extended(boxplot_data_extended)

        self.plot_glycan_type_distribution(df)
        self.plot_heatmap(df)
        self.plot_heatmap_full_profile(df)
        self.plot_histogram_by_sample(df)
        self.plot_histogram_normalized(df)  # Add normalized histogram

        logger.info(f"All visualizations saved to {self.output_dir}")


if __name__ == "__main__":
    # Test
    pass
