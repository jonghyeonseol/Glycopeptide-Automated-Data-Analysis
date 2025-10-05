"""
PCA Plot Module for pGlyco Auto Combine
Handles PCA visualization
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from matplotlib.patches import Ellipse
from adjustText import adjust_text
from ..utils import save_trace_data
from .plot_config import (
    PCA_FIGSIZE, PCA_LABEL_FONTSIZE, PCA_POINT_SIZE,
    PCA_POINT_LINEWIDTH, PCA_POINT_ALPHA, GROUP_PALETTE,
    apply_standard_axis_style, apply_standard_legend
)

logger = logging.getLogger(__name__)


class PCAPlotMixin:
    """Mixin class for PCA-related plots"""

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
        fig, ax = plt.subplots(figsize=figsize)

        # Extract data
        pca_df = pca_results['pca_df']
        explained_var = pca_results['explained_variance']

        # Separate Cancer and Normal samples
        cancer_mask = pca_df['Group'] == 'Cancer'
        normal_mask = pca_df['Group'] == 'Normal'

        # Use standardized colors from plot_config
        cancer_color = GROUP_PALETTE['Cancer']
        normal_color = GROUP_PALETTE['Normal']

        # Plot Cancer samples (Prism style: larger, bolder points)
        cancer_scores = pca_df[cancer_mask]
        ax.scatter(cancer_scores['PC1'], cancer_scores['PC2'],
                  c=cancer_color, s=PCA_POINT_SIZE, alpha=PCA_POINT_ALPHA,
                  edgecolors='black', linewidth=PCA_POINT_LINEWIDTH,
                  label='Cancer', zorder=3)

        # Plot Normal samples (Prism style: larger, bolder points)
        normal_scores = pca_df[normal_mask]
        ax.scatter(normal_scores['PC1'], normal_scores['PC2'],
                  c=normal_color, s=PCA_POINT_SIZE, alpha=PCA_POINT_ALPHA,
                  edgecolors='black', linewidth=PCA_POINT_LINEWIDTH,
                  label='Normal', zorder=3)

        # Draw 95% confidence ellipses
        self._draw_confidence_ellipse(ax, cancer_scores['PC1'].values, cancer_scores['PC2'].values,
                                     color=cancer_color, alpha=0.15)
        self._draw_confidence_ellipse(ax, normal_scores['PC1'].values, normal_scores['PC2'].values,
                                     color=normal_color, alpha=0.15)

        # Add sample labels (for all points)
        texts = []
        for idx, row in pca_df.iterrows():
            # Color label by group
            label_color = cancer_color if row['Group'] == 'Cancer' else normal_color

            # Add text annotation using standardized font size
            texts.append(ax.text(row['PC1'], row['PC2'], idx,
                               fontsize=PCA_LABEL_FONTSIZE, ha='center', va='center',
                               color=label_color, fontweight='bold',
                               bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                       edgecolor=label_color, alpha=0.7, linewidth=1),
                               zorder=5))

        # Adjust text to avoid overlap
        try:
            adjust_text(texts,
                       arrowprops=dict(arrowstyle='->', color='gray', lw=0.5, alpha=0.6),
                       expand_points=(1.5, 1.5),
                       force_text=(0.5, 0.5))
        except:
            logger.warning("adjustText optimization failed, showing overlapping labels")

        # Apply standardized styling
        apply_standard_axis_style(
            ax,
            xlabel=f'PC1 ({explained_var[0]*100:.2f}%)',
            ylabel=f'PC2 ({explained_var[1]*100:.2f}%)',
            title='PCA: Cancer vs Normal Samples (with Sample Labels)',
            grid=True
        )

        # Apply standardized legend (positioned outside plot area)
        apply_standard_legend(ax)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'pca_plot.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved PCA plot to {output_file}")

        # Save trace data
        trace_data = pca_df.copy()
        trace_data['PC1_variance'] = explained_var[0]
        trace_data['PC2_variance'] = explained_var[1]
        save_trace_data(trace_data, self.output_dir, 'pca_plot_data.csv')

        plt.close()

    def plot_pca_by_glycan_type(self, df: pd.DataFrame, pca_results: dict, figsize: tuple = (12, 8)):
        """
        Create PCA plot colored by sample names (alternative view)

        Args:
            df: Annotated DataFrame
            pca_results: PCA results from analyzer
            figsize: Figure size
        """
        fig, ax = plt.subplots(figsize=figsize)

        # Extract data
        pca_df = pca_results['pca_df']
        explained_var = pca_results['explained_variance']

        # Plot with sample labels
        for idx, row in pca_df.iterrows():
            color = '#E74C3C' if row['Group'] == 'Cancer' else '#3498DB'
            marker = 'o' if row['Group'] == 'Cancer' else 's'

            ax.scatter(row['PC1'], row['PC2'], c=color, s=100,
                      alpha=0.6, edgecolors='black', linewidth=1,
                      marker=marker, zorder=3)

            # Add sample label
            ax.annotate(idx, (row['PC1'], row['PC2']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.7)

        ax.set_xlabel(f'PC1 ({explained_var[0]*100:.2f}%)', fontsize=12)
        ax.set_ylabel(f'PC2 ({explained_var[1]*100:.2f}%)', fontsize=12)
        ax.set_title('PCA with Sample Labels', fontsize=14, fontweight='bold')

        # Custom legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#E74C3C', edgecolor='black', label='Cancer'),
            Patch(facecolor='#3498DB', edgecolor='black', label='Normal')
        ]
        ax.legend(handles=legend_elements, loc='best', frameon=True)

        ax.grid(True, alpha=0.3)
        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'pca_samples.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved PCA plot to {output_file}")

        # Save trace data
        trace_data = pca_df.copy()
        trace_data['PC1_variance'] = explained_var[0]
        trace_data['PC2_variance'] = explained_var[1]
        save_trace_data(trace_data, self.output_dir, 'pca_samples_data.csv')

        plt.close()
