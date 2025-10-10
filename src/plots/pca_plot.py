"""
PCA Plot Module for pGlyco Auto Combine
Handles PCA visualization
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
from matplotlib.patches import Ellipse
from adjustText import adjust_text
from ..utils import save_trace_data
from .plot_config import (
    PCA_LABEL_FONTSIZE, PCA_POINT_SIZE,
    PCA_POINT_LINEWIDTH, PCA_POINT_ALPHA,
    apply_standard_axis_style, apply_standard_legend,
    add_sample_size_annotation,  # Phase 2.2 enhancement
    get_group_style,  # Phase 3 enhancement
    save_publication_figure, PCA_DPI,
    create_fancy_bbox, apply_publication_theme,  # ✨ Enhanced styling
    GRADIENT_ALPHA_START, GRADIENT_ALPHA_END  # ✨ Gradient fills
)

logger = logging.getLogger(__name__)


class PCAPlotMixin:
    """Mixin class for PCA-related plots"""

    def _draw_confidence_ellipse(self, ax, x, y, color, alpha=0.2, n_std=1.96):
        """
        Draw 95% confidence ellipse with ENHANCED gradient fill

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

        # ✨ ENHANCED: Draw multiple ellipses with gradient alpha (center→edge fade)
        n_gradient_levels = 5
        alphas = np.linspace(GRADIENT_ALPHA_START * 0.4, GRADIENT_ALPHA_END * 0.3, n_gradient_levels)

        for i, alpha_val in enumerate(alphas):
            scale = 1.0 - (i * 0.15)  # Shrink towards center
            ellipse = Ellipse(
                xy=(np.mean(x), np.mean(y)),
                width=width * scale,
                height=height * scale,
                angle=angle,
                facecolor=color,
                edgecolor='none',  # No edge for gradient layers
                alpha=alpha_val,
                linewidth=0,
                zorder=1
            )
            ax.add_patch(ellipse)

        # ✨ Draw outer boundary with dashed line
        ellipse_boundary = Ellipse(
            xy=(np.mean(x), np.mean(y)),
            width=width,
            height=height,
            angle=angle,
            facecolor='none',
            edgecolor=color,
            alpha=0.6,
            linewidth=2,
            linestyle='--',
            zorder=2
        )
        ax.add_patch(ellipse_boundary)

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

        # Phase 3: Use color + shape encoding for colorblind accessibility
        cancer_color, cancer_marker = get_group_style('Cancer')
        normal_color, normal_marker = get_group_style('Normal')

        # Plot Cancer samples (Circle markers - colorblind-safe)
        cancer_scores = pca_df[cancer_mask]
        ax.scatter(cancer_scores['PC1'], cancer_scores['PC2'],
                   c=cancer_color, s=PCA_POINT_SIZE, alpha=PCA_POINT_ALPHA,
                   edgecolors='black', linewidth=PCA_POINT_LINEWIDTH,
                   marker=cancer_marker,  # Phase 3: Circle for Cancer
                   label='Cancer (○)', zorder=3)

        # Plot Normal samples (Square markers - colorblind-safe)
        normal_scores = pca_df[normal_mask]
        ax.scatter(normal_scores['PC1'], normal_scores['PC2'],
                   c=normal_color, s=PCA_POINT_SIZE, alpha=PCA_POINT_ALPHA,
                   edgecolors='black', linewidth=PCA_POINT_LINEWIDTH,
                   marker=normal_marker,  # Phase 3: Square for Normal
                   label='Normal (□)', zorder=3)

        # ✨ Draw ENHANCED 95% confidence ellipses with gradient fill
        self._draw_confidence_ellipse(ax, cancer_scores['PC1'].values, cancer_scores['PC2'].values,
                                      color=cancer_color, alpha=0.15)
        self._draw_confidence_ellipse(ax, normal_scores['PC1'].values, normal_scores['PC2'].values,
                                      color=normal_color, alpha=0.15)

        # Add sample labels (for all points)
        texts = []
        for idx, row in pca_df.iterrows():
            # Color label by group
            label_color = cancer_color if row['Group'] == 'Cancer' else normal_color

            # ✨ ENHANCED: Add text annotation with fancy bbox
            bbox_props = create_fancy_bbox(facecolor='white', edgecolor=label_color,
                                           alpha=0.85, linewidth=1.2)
            texts.append(ax.text(row['PC1'], row['PC2'], idx,
                                 fontsize=PCA_LABEL_FONTSIZE, ha='center', va='center',
                                 color=label_color, fontweight='bold',
                                 bbox=bbox_props,
                                 zorder=5))

        # Adjust text to avoid overlap
        try:
            adjust_text(texts,
                        arrowprops=dict(arrowstyle='->', color='gray', lw=0.5, alpha=0.6),
                        expand_points=(1.5, 1.5),
                        force_text=(0.5, 0.5))
        except Exception:
            logger.warning("adjustText optimization failed, showing overlapping labels")

        # Apply standardized styling
        apply_standard_axis_style(
            ax,
            xlabel =f'PC1 ({explained_var[0] * 100:.2f}%)',
            ylabel =f'PC2 ({explained_var[1] * 100:.2f}%)',
            title='PCA: Cancer vs Normal Samples (with Sample Labels)',
            grid=True
        )

        # Apply standardized legend (positioned outside plot area)
        apply_standard_legend(ax)

        # Add sample size annotation (Phase 2.2 enhancement)
        n_cancer = cancer_mask.sum()
        n_normal = normal_mask.sum()
        add_sample_size_annotation(ax, n_cancer=n_cancer, n_normal=n_normal,
                                   location='lower right', fontsize=10)

        plt.tight_layout()

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        # Save plot
        output_file = self.output_dir / 'pca_plot.png'
        save_publication_figure(fig, output_file, dpi=PCA_DPI)
        logger.info(f"✨ Saved ENHANCED PCA plot to {output_file} ({PCA_DPI} DPI)")

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

        # Phase 3: Plot with color + shape encoding
        for idx, row in pca_df.iterrows():
            color, marker = get_group_style(row['Group'])

            ax.scatter(row['PC1'], row['PC2'], c=color, s=100,
                       alpha=0.6, edgecolors='black', linewidth=1,
                       marker=marker, zorder=3)  # Phase 3: Colorblind-safe markers

            # Add sample label
            ax.annotate(idx, (row['PC1'], row['PC2']),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=8, alpha=0.7)

        ax.set_xlabel(f'PC1 ({explained_var[0] * 100:.2f}%)', fontsize=12)
        ax.set_ylabel(f'PC2 ({explained_var[1] * 100:.2f}%)', fontsize=12)
        ax.set_title('PCA with Sample Labels', fontsize=14, fontweight='bold')

        # Phase 3: Custom legend with shape markers
        from matplotlib.lines import Line2D
        cancer_color, cancer_marker = get_group_style('Cancer')
        normal_color, normal_marker = get_group_style('Normal')

        legend_elements = [
            Line2D([0], [0], marker=cancer_marker, color='w',
                   markerfacecolor=cancer_color, markeredgecolor='black',
                   markersize=10, label='Cancer (○)'),
            Line2D([0], [0], marker=normal_marker, color='w',
                   markerfacecolor=normal_color, markeredgecolor='black',
                   markersize=10, label='Normal (□)')
        ]
        ax.legend(handles=legend_elements, loc='best', frameon=True)

        # Add sample size annotation (Phase 2.2 enhancement)
        n_cancer = (pca_df['Group'] == 'Cancer').sum()
        n_normal = (pca_df['Group'] == 'Normal').sum()
        add_sample_size_annotation(ax, n_cancer=n_cancer, n_normal=n_normal,
                                   location='lower right', fontsize=10)

        ax.grid(True, alpha=0.3)
        plt.tight_layout()

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        # Save plot
        output_file = self.output_dir / 'pca_samples.png'
        save_publication_figure(fig, output_file, dpi=PCA_DPI)
        logger.info(f"✨ Saved ENHANCED PCA plot to {output_file} ({PCA_DPI} DPI)")

        # Save trace data
        trace_data = pca_df.copy()
        trace_data['PC1_variance'] = explained_var[0]
        trace_data['PC2_variance'] = explained_var[1]
        save_trace_data(trace_data, self.output_dir, 'pca_samples_data.csv')

        plt.close()
