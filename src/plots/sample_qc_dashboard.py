"""
Sample QC Dashboard Module for pGlyco Auto Combine
Creates comprehensive per-sample quality control visualization

Dependencies:
    External:
        - pandas: Data manipulation
        - numpy: Numerical computations
        - matplotlib: Plotting backend
        - scipy.stats: Statistical functions (chi2 distribution for outlier detection)

    Internal:
        - src.utils: save_trace_data, get_sample_columns
        - src.plots.plot_config: QC_* constants, save_publication_figure

Phase 4.3: Critical for identifying problematic samples
Shows quality metrics, outliers, and technical variability per sample
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import logging

from ..utils import save_trace_data, get_sample_columns
from .plot_config import (
    QC_FIGSIZE, QC_THRESHOLD_DETECTION_RATE, QC_THRESHOLD_CV_GOOD,
    QC_THRESHOLD_CV_ACCEPTABLE, QC_THRESHOLD_OUTLIER_P,
    QC_LABEL_FONTSIZE, QC_TITLE_FONTSIZE, QC_LEGEND_FONTSIZE, QC_XTICKLABEL_FONTSIZE,
    QC_DPI, save_publication_figure, COLOR_CANCER, COLOR_NORMAL,
    EDGE_LINEWIDTH_THIN, PLOT_LINE_LINEWIDTH, LINE_BOLD,
    LINE_ALPHA, ALPHA_MEDIUM_LIGHT,
    EDGE_COLOR_BLACK,  # Edge color standardization
    # Linestyle constants (Phase 10.3.8)
    THRESHOLD_LINESTYLE
)

logger = logging.getLogger(__name__)


class SampleQCDashboardMixin:
    """Mixin class for per-sample QC dashboard"""

    def plot_sample_qc_dashboard(self, df: pd.DataFrame, figsize: tuple = None):
        """
        Create comprehensive per-sample QC dashboard

        CRITICAL FOR DATA QUALITY:
        - Panel 1: Total intensity per sample (TIC)
        - Panel 2: Detection count (# glycopeptides detected)
        - Panel 3: Detection rate (% of all glycopeptides)
        - Panel 4: Median intensity (central tendency)
        - Panel 5: CV distribution (technical variability)
        - Panel 6: Outlier detection (Mahalanobis distance)

        Args:
            df: Annotated DataFrame with intensity data
            figsize: Figure size (width, height), defaults to QC_FIGSIZE
        """
        logger.info("Creating per-sample QC dashboard...")

        # Use centralized figsize
        if figsize is None:
            figsize = QC_FIGSIZE

        # Get sample columns using centralized function
        cancer_samples, normal_samples = get_sample_columns(df)
        all_samples = cancer_samples + normal_samples
        logger.debug(f"  Analyzing {len(all_samples)} samples ({len(cancer_samples)} Cancer, {len(normal_samples)} Normal)")

        # Extract intensity matrix
        intensity_matrix = df[all_samples].copy()
        intensity_matrix = intensity_matrix.replace('', np.nan).astype(float)

        # Calculate QC metrics for each sample
        qc_metrics = self._calculate_qc_metrics(intensity_matrix, all_samples, cancer_samples)
        logger.debug(f"  QC metrics calculated for {len(qc_metrics)} samples")

        # Create 2x3 subplot layout
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        fig.suptitle('Per-Sample QC Dashboard: Data Quality Assessment',
                     fontsize=QC_TITLE_FONTSIZE + 5, fontweight='bold', y=0.995)

        # Sample colors (Cancer=red, Normal=blue)
        sample_colors = [COLOR_CANCER if s in cancer_samples else COLOR_NORMAL for s in all_samples]

        # =====================================================================
        # PANEL 1: Total Intensity (TIC) per sample
        # =====================================================================
        ax1 = axes[0, 0]
        logger.debug("  Panel 1: Total intensity (TIC)...")

        _ = ax1.bar(range(len(all_samples)), qc_metrics['total_intensity'],
                    color=sample_colors, edgecolor=EDGE_COLOR_BLACK, linewidth=EDGE_LINEWIDTH_THIN)

        # Add group means as horizontal lines
        cancer_mean_tic = qc_metrics.loc[qc_metrics['group'] == 'Cancer', 'total_intensity'].mean()
        normal_mean_tic = qc_metrics.loc[qc_metrics['group'] == 'Normal', 'total_intensity'].mean()

        ax1.axhline(cancer_mean_tic, color=COLOR_CANCER, linestyle=THRESHOLD_LINESTYLE, linewidth=PLOT_LINE_LINEWIDTH,
                    label=f'Cancer mean: {cancer_mean_tic:.2e}', alpha=LINE_ALPHA)
        ax1.axhline(normal_mean_tic, color=COLOR_NORMAL, linestyle=THRESHOLD_LINESTYLE, linewidth=PLOT_LINE_LINEWIDTH,
                    label=f'Normal mean: {normal_mean_tic:.2e}', alpha=LINE_ALPHA)

        ax1.set_xlabel('Sample', fontsize=QC_LABEL_FONTSIZE, fontweight='bold')
        ax1.set_ylabel('Total Intensity (TIC)', fontsize=QC_LABEL_FONTSIZE, fontweight='bold')
        ax1.set_title('Total Ion Current (TIC)', fontsize=QC_TITLE_FONTSIZE, fontweight='bold')
        ax1.set_xticks(range(len(all_samples)))
        ax1.set_xticklabels(all_samples, rotation=90, fontsize=QC_XTICKLABEL_FONTSIZE)
        ax1.legend(loc='upper right', fontsize=QC_LEGEND_FONTSIZE)
        ax1.grid(axis='y', alpha=ALPHA_MEDIUM_LIGHT)
        ax1.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

        # =====================================================================
        # PANEL 2: Detection Count
        # =====================================================================
        ax2 = axes[0, 1]
        logger.debug("  Panel 2: Detection count...")

        _ = ax2.bar(range(len(all_samples)), qc_metrics['detection_count'],
                        color=sample_colors, edgecolor=EDGE_COLOR_BLACK, linewidth=EDGE_LINEWIDTH_THIN)

        # Add group means
        cancer_mean_det = qc_metrics.loc[qc_metrics['group'] == 'Cancer', 'detection_count'].mean()
        normal_mean_det = qc_metrics.loc[qc_metrics['group'] == 'Normal', 'detection_count'].mean()

        ax2.axhline(cancer_mean_det, color=COLOR_CANCER, linestyle=THRESHOLD_LINESTYLE, linewidth=PLOT_LINE_LINEWIDTH,
                    label=f'Cancer mean: {cancer_mean_det:.0f}', alpha=LINE_ALPHA)
        ax2.axhline(normal_mean_det, color=COLOR_NORMAL, linestyle=THRESHOLD_LINESTYLE, linewidth=PLOT_LINE_LINEWIDTH,
                    label=f'Normal mean: {normal_mean_det:.0f}', alpha=LINE_ALPHA)

        ax2.set_xlabel('Sample', fontsize=QC_LABEL_FONTSIZE, fontweight='bold')
        ax2.set_ylabel('Detected Glycopeptides', fontsize=QC_LABEL_FONTSIZE, fontweight='bold')
        ax2.set_title('Detection Count', fontsize=QC_TITLE_FONTSIZE, fontweight='bold')
        ax2.set_xticks(range(len(all_samples)))
        ax2.set_xticklabels(all_samples, rotation=90, fontsize=QC_XTICKLABEL_FONTSIZE)
        ax2.legend(loc='upper right', fontsize=QC_LEGEND_FONTSIZE)
        ax2.grid(axis='y', alpha=ALPHA_MEDIUM_LIGHT)

        # =====================================================================
        # PANEL 3: Detection Rate (%) with centralized threshold
        # =====================================================================
        ax3 = axes[0, 2]
        logger.debug("  Panel 3: Detection rate...")

        _ = ax3.bar(range(len(all_samples)), qc_metrics['detection_rate_pct'],
                        color=sample_colors, edgecolor=EDGE_COLOR_BLACK, linewidth=EDGE_LINEWIDTH_THIN)

        # Add threshold line using centralized constant
        ax3.axhline(QC_THRESHOLD_DETECTION_RATE, color='red', linestyle=THRESHOLD_LINESTYLE, linewidth=PLOT_LINE_LINEWIDTH,
                    label=f'{QC_THRESHOLD_DETECTION_RATE}% threshold', alpha=LINE_ALPHA)

        ax3.set_xlabel('Sample', fontsize=QC_LABEL_FONTSIZE, fontweight='bold')
        ax3.set_ylabel('Detection Rate (%)', fontsize=QC_LABEL_FONTSIZE, fontweight='bold')
        ax3.set_title('Detection Completeness', fontsize=QC_TITLE_FONTSIZE, fontweight='bold')
        ax3.set_xticks(range(len(all_samples)))
        ax3.set_xticklabels(all_samples, rotation=90, fontsize=QC_XTICKLABEL_FONTSIZE)
        ax3.set_ylim(0, 100)
        ax3.legend(loc='upper right', fontsize=QC_LEGEND_FONTSIZE)
        ax3.grid(axis='y', alpha=ALPHA_MEDIUM_LIGHT)

        # =====================================================================
        # PANEL 4: Median Intensity
        # =====================================================================
        ax4 = axes[1, 0]
        logger.debug("  Panel 4: Median intensity...")

        _ = ax4.bar(range(len(all_samples)), qc_metrics['median_intensity'],
                        color=sample_colors, edgecolor=EDGE_COLOR_BLACK, linewidth=EDGE_LINEWIDTH_THIN)

        # Add group means
        cancer_mean_med = qc_metrics.loc[qc_metrics['group'] == 'Cancer', 'median_intensity'].mean()
        normal_mean_med = qc_metrics.loc[qc_metrics['group'] == 'Normal', 'median_intensity'].mean()

        ax4.axhline(cancer_mean_med, color=COLOR_CANCER, linestyle=THRESHOLD_LINESTYLE, linewidth=PLOT_LINE_LINEWIDTH,
                    label=f'Cancer mean: {cancer_mean_med:.2e}', alpha=LINE_ALPHA)
        ax4.axhline(normal_mean_med, color=COLOR_NORMAL, linestyle=THRESHOLD_LINESTYLE, linewidth=PLOT_LINE_LINEWIDTH,
                    label=f'Normal mean: {normal_mean_med:.2e}', alpha=LINE_ALPHA)

        ax4.set_xlabel('Sample', fontsize=QC_LABEL_FONTSIZE, fontweight='bold')
        ax4.set_ylabel('Median Intensity', fontsize=QC_LABEL_FONTSIZE, fontweight='bold')
        ax4.set_title('Central Tendency', fontsize=QC_TITLE_FONTSIZE, fontweight='bold')
        ax4.set_xticks(range(len(all_samples)))
        ax4.set_xticklabels(all_samples, rotation=90, fontsize=QC_XTICKLABEL_FONTSIZE)
        ax4.legend(loc='upper right', fontsize=QC_LEGEND_FONTSIZE)
        ax4.grid(axis='y', alpha=ALPHA_MEDIUM_LIGHT)
        ax4.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

        # =====================================================================
        # PANEL 5: CV Distribution with centralized thresholds
        # =====================================================================
        ax5 = axes[1, 1]
        logger.debug("  Panel 5: CV distribution...")

        _ = ax5.bar(range(len(all_samples)), qc_metrics['cv_median_pct'],
                        color=sample_colors, edgecolor=EDGE_COLOR_BLACK, linewidth=EDGE_LINEWIDTH_THIN)

        # Add threshold lines using centralized constants
        ax5.axhline(QC_THRESHOLD_CV_GOOD, color='green', linestyle=THRESHOLD_LINESTYLE, linewidth=PLOT_LINE_LINEWIDTH,
                    label=f'Good QC (<{QC_THRESHOLD_CV_GOOD}%)', alpha=LINE_ALPHA)
        ax5.axhline(QC_THRESHOLD_CV_ACCEPTABLE, color='orange', linestyle=THRESHOLD_LINESTYLE, linewidth=PLOT_LINE_LINEWIDTH,
                    label=f'Acceptable (<{QC_THRESHOLD_CV_ACCEPTABLE}%)', alpha=LINE_ALPHA)

        ax5.set_xlabel('Sample', fontsize=QC_LABEL_FONTSIZE, fontweight='bold')
        ax5.set_ylabel('Median CV (%)', fontsize=QC_LABEL_FONTSIZE, fontweight='bold')
        ax5.set_title('Technical Variability', fontsize=QC_TITLE_FONTSIZE, fontweight='bold')
        ax5.set_xticks(range(len(all_samples)))
        ax5.set_xticklabels(all_samples, rotation=90, fontsize=QC_XTICKLABEL_FONTSIZE)
        ax5.legend(loc='upper right', fontsize=QC_LEGEND_FONTSIZE)
        ax5.grid(axis='y', alpha=ALPHA_MEDIUM_LIGHT)

        # =====================================================================
        # PANEL 6: Outlier Detection with dynamic df calculation
        # =====================================================================
        ax6 = axes[1, 2]
        logger.debug("  Panel 6: Outlier detection (Mahalanobis distance)...")

        # Calculate Mahalanobis distance for outlier detection
        outlier_scores = self._calculate_mahalanobis_distance(qc_metrics)
        qc_metrics['mahalanobis_distance'] = outlier_scores

        bars6 = ax6.bar(range(len(all_samples)), outlier_scores,
                        color=sample_colors, edgecolor=EDGE_COLOR_BLACK, linewidth=EDGE_LINEWIDTH_THIN)

        # ========================================
        # DYNAMIC DF CALCULATION (NEW)
        # ========================================
        # Calculate df from the actual metric columns used in Mahalanobis calculation
        metric_cols = ['total_intensity', 'detection_count', 'detection_rate_pct',
                       'median_intensity', 'cv_median_pct']
        n_metrics = len(metric_cols)

        # Use centralized confidence level
        threshold = stats.chi2.ppf(QC_THRESHOLD_OUTLIER_P, df=n_metrics)
        logger.debug(f"  Outlier threshold: {threshold:.2f} (chi2, p={QC_THRESHOLD_OUTLIER_P}, df={n_metrics})")

        ax6.axhline(threshold, color='red', linestyle=THRESHOLD_LINESTYLE, linewidth=PLOT_LINE_LINEWIDTH,
                    label=f'Outlier threshold (p<{1-QC_THRESHOLD_OUTLIER_P:.2f})', alpha=LINE_ALPHA)

        ax6.set_xlabel('Sample', fontsize=QC_LABEL_FONTSIZE, fontweight='bold')
        ax6.set_ylabel('Mahalanobis Distance', fontsize=QC_LABEL_FONTSIZE, fontweight='bold')
        ax6.set_title('Outlier Detection', fontsize=QC_TITLE_FONTSIZE, fontweight='bold')
        ax6.set_xticks(range(len(all_samples)))
        ax6.set_xticklabels(all_samples, rotation=90, fontsize=QC_XTICKLABEL_FONTSIZE)
        ax6.legend(loc='upper right', fontsize=QC_LEGEND_FONTSIZE)
        ax6.grid(axis='y', alpha=ALPHA_MEDIUM_LIGHT)

        # Highlight outliers
        outliers = np.where(outlier_scores > threshold)[0]
        if len(outliers) > 0:
            logger.debug(f"  Found {len(outliers)} outlier(s)")
            for outlier_idx in outliers:
                bars6[outlier_idx].set_edgecolor('red')
                bars6[outlier_idx].set_linewidth(LINE_BOLD)
        else:
            logger.debug("  No outliers detected")

        plt.tight_layout()

        # ========================================
        # SAVE with standardized function
        # ========================================
        output_file = self.output_dir / 'sample_qc_dashboard.png'
        save_publication_figure(fig, output_file, dpi=QC_DPI)
        logger.info(f"Saved sample QC dashboard to {output_file} (optimized, {QC_DPI} DPI)")

        # ========================================
        # TRACE DATA for reproducibility
        # ========================================
        save_trace_data(qc_metrics, self.output_dir, 'sample_qc_metrics.csv')

        plt.close()

        # Summary report
        n_outliers = len(outliers)
        outlier_samples = [all_samples[i] for i in outliers]

        logger.info("✓ Per-sample QC analysis complete")
        logger.info(f"  Samples analyzed: {len(all_samples)}")
        logger.info(f"  Outliers detected: {n_outliers}")
        if n_outliers > 0:
            logger.info(f"  Outlier samples: {', '.join(outlier_samples)}")
        else:
            logger.info("  No outliers detected - all samples pass QC")

    def _calculate_qc_metrics(self, intensity_matrix: pd.DataFrame,
                              all_samples: list, cancer_samples: list) -> pd.DataFrame:
        """
        Calculate QC metrics for each sample

        Metrics calculated:
            - total_intensity: Sum of all intensities (TIC)
            - detection_count: Number of detected glycopeptides
            - detection_rate_pct: Percentage of detected glycopeptides
            - median_intensity: Central tendency measure
            - cv_median_pct: Coefficient of variation (technical variability)

        Args:
            intensity_matrix: Intensity matrix (glycopeptides x samples)
            all_samples: List of sample names
            cancer_samples: List of cancer sample names (for group assignment)

        Returns:
            DataFrame with QC metrics per sample
        """
        qc_data = []
        cancer_samples_set = set(cancer_samples)  # For O(1) lookup

        for sample in all_samples:
            sample_intensities = intensity_matrix[sample]

            # Total intensity (TIC)
            total_intensity = sample_intensities.sum()

            # Detection count and rate
            detection_count = sample_intensities.notna().sum()
            detection_rate_pct = detection_count / len(intensity_matrix) * 100

            # Median intensity (central tendency)
            median_intensity = sample_intensities.median()

            # CV (coefficient of variation) - median across all glycopeptides
            # Note: CV = std / mean for each glycopeptide, then take median across all
            mean_int = sample_intensities.mean()
            std_int = sample_intensities.std()
            cv_median_pct = (std_int / mean_int * 100) if mean_int > 0 else 0

            # Group
            group = 'Cancer' if sample in cancer_samples_set else 'Normal'

            qc_data.append({
                'sample': sample,
                'group': group,
                'total_intensity': total_intensity,
                'detection_count': detection_count,
                'detection_rate_pct': detection_rate_pct,
                'median_intensity': median_intensity,
                'cv_median_pct': cv_median_pct
            })

        return pd.DataFrame(qc_data)

    def _calculate_mahalanobis_distance(self, qc_metrics: pd.DataFrame) -> np.ndarray:
        """
        Calculate Mahalanobis distance for outlier detection

        Mahalanobis distance accounts for correlations between metrics
        and is more robust than simple z-scores.

        Args:
            qc_metrics: QC metrics DataFrame

        Returns:
            Array of Mahalanobis distances (one per sample)
        """
        # Select numeric columns for multivariate analysis
        metric_cols = ['total_intensity', 'detection_count', 'detection_rate_pct',
                       'median_intensity', 'cv_median_pct']

        X = qc_metrics[metric_cols].values

        # Standardize (mean=0, std=1)
        X_centered = X - X.mean(axis=0)

        # Covariance matrix
        cov_matrix = np.cov(X_centered, rowvar=False)

        # Regularize if needed (avoid singular matrix)
        try:
            cov_inv = np.linalg.inv(cov_matrix)
        except np.linalg.LinAlgError:
            logger.warning("Singular covariance matrix - using pseudo-inverse")
            cov_inv = np.linalg.pinv(cov_matrix)

        # Mahalanobis distance for each sample
        mahalanobis_distances = []
        for i in range(len(X_centered)):
            diff = X_centered[i]
            distance = np.sqrt(diff @ cov_inv @ diff.T)
            mahalanobis_distances.append(distance)

        return np.array(mahalanobis_distances)
