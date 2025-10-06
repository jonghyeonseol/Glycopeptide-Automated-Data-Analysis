"""
Sample QC Dashboard Module for pGlyco Auto Combine
Creates comprehensive per-sample quality control visualization

Phase 4.3: Critical for identifying problematic samples
Shows quality metrics, outliers, and technical variability per sample
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import logging
from ..utils import save_trace_data, get_sample_columns

logger = logging.getLogger(__name__)


class SampleQCDashboardMixin:
    """Mixin class for per-sample QC dashboard"""

    def plot_sample_qc_dashboard(self, df: pd.DataFrame, figsize: tuple = (20, 12)):
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
            figsize: Figure size (width, height)
        """
        logger.info("Creating per-sample QC dashboard...")

        # Get sample columns
        cancer_samples, normal_samples = get_sample_columns(df)
        all_samples = cancer_samples + normal_samples

        # Extract intensity matrix
        intensity_matrix = df[all_samples].copy()
        intensity_matrix = intensity_matrix.replace('', np.nan).astype(float)

        # Calculate QC metrics for each sample
        qc_metrics = self._calculate_qc_metrics(intensity_matrix, all_samples)

        # Create 2x3 subplot layout
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        fig.suptitle('Per-Sample QC Dashboard: Data Quality Assessment',
                    fontsize=16, fontweight='bold', y=0.995)

        # Sample colors (Cancer vs Normal)
        sample_colors = ['#E74C3C' if s.startswith('C') else '#3498DB' for s in all_samples]
        sample_groups = ['Cancer' if s.startswith('C') else 'Normal' for s in all_samples]

        # =====================================================================
        # PANEL 1: Total Intensity (TIC) per sample
        # =====================================================================
        ax1 = axes[0, 0]
        logger.info("  Panel 1: Total intensity per sample...")

        bars1 = ax1.bar(range(len(all_samples)), qc_metrics['total_intensity'],
                       color=sample_colors, edgecolor='black', linewidth=0.5)

        # Add group means as horizontal lines
        cancer_mean_tic = qc_metrics.loc[qc_metrics['group'] == 'Cancer', 'total_intensity'].mean()
        normal_mean_tic = qc_metrics.loc[qc_metrics['group'] == 'Normal', 'total_intensity'].mean()

        ax1.axhline(cancer_mean_tic, color='#E74C3C', linestyle='--', linewidth=2,
                   label=f'Cancer mean: {cancer_mean_tic:.2e}', alpha=0.7)
        ax1.axhline(normal_mean_tic, color='#3498DB', linestyle='--', linewidth=2,
                   label=f'Normal mean: {normal_mean_tic:.2e}', alpha=0.7)

        ax1.set_xlabel('Sample', fontsize=10, fontweight='bold')
        ax1.set_ylabel('Total Intensity (TIC)', fontsize=10, fontweight='bold')
        ax1.set_title('Total Ion Current (TIC)', fontsize=11, fontweight='bold')
        ax1.set_xticks(range(len(all_samples)))
        ax1.set_xticklabels(all_samples, rotation=90, fontsize=7)
        ax1.legend(loc='upper right', fontsize=8)
        ax1.grid(axis='y', alpha=0.3)
        ax1.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

        # =====================================================================
        # PANEL 2: Detection Count
        # =====================================================================
        ax2 = axes[0, 1]
        logger.info("  Panel 2: Detection count per sample...")

        bars2 = ax2.bar(range(len(all_samples)), qc_metrics['detection_count'],
                       color=sample_colors, edgecolor='black', linewidth=0.5)

        # Add group means
        cancer_mean_det = qc_metrics.loc[qc_metrics['group'] == 'Cancer', 'detection_count'].mean()
        normal_mean_det = qc_metrics.loc[qc_metrics['group'] == 'Normal', 'detection_count'].mean()

        ax2.axhline(cancer_mean_det, color='#E74C3C', linestyle='--', linewidth=2,
                   label=f'Cancer mean: {cancer_mean_det:.0f}', alpha=0.7)
        ax2.axhline(normal_mean_det, color='#3498DB', linestyle='--', linewidth=2,
                   label=f'Normal mean: {normal_mean_det:.0f}', alpha=0.7)

        ax2.set_xlabel('Sample', fontsize=10, fontweight='bold')
        ax2.set_ylabel('Detected Glycopeptides', fontsize=10, fontweight='bold')
        ax2.set_title('Detection Count', fontsize=11, fontweight='bold')
        ax2.set_xticks(range(len(all_samples)))
        ax2.set_xticklabels(all_samples, rotation=90, fontsize=7)
        ax2.legend(loc='upper right', fontsize=8)
        ax2.grid(axis='y', alpha=0.3)

        # =====================================================================
        # PANEL 3: Detection Rate (%)
        # =====================================================================
        ax3 = axes[0, 2]
        logger.info("  Panel 3: Detection rate per sample...")

        bars3 = ax3.bar(range(len(all_samples)), qc_metrics['detection_rate_pct'],
                       color=sample_colors, edgecolor='black', linewidth=0.5)

        # Add threshold line at 50%
        ax3.axhline(50, color='red', linestyle='--', linewidth=2,
                   label='50% threshold', alpha=0.7)

        ax3.set_xlabel('Sample', fontsize=10, fontweight='bold')
        ax3.set_ylabel('Detection Rate (%)', fontsize=10, fontweight='bold')
        ax3.set_title('Detection Completeness', fontsize=11, fontweight='bold')
        ax3.set_xticks(range(len(all_samples)))
        ax3.set_xticklabels(all_samples, rotation=90, fontsize=7)
        ax3.set_ylim(0, 100)
        ax3.legend(loc='upper right', fontsize=8)
        ax3.grid(axis='y', alpha=0.3)

        # =====================================================================
        # PANEL 4: Median Intensity
        # =====================================================================
        ax4 = axes[1, 0]
        logger.info("  Panel 4: Median intensity per sample...")

        bars4 = ax4.bar(range(len(all_samples)), qc_metrics['median_intensity'],
                       color=sample_colors, edgecolor='black', linewidth=0.5)

        # Add group means
        cancer_mean_med = qc_metrics.loc[qc_metrics['group'] == 'Cancer', 'median_intensity'].mean()
        normal_mean_med = qc_metrics.loc[qc_metrics['group'] == 'Normal', 'median_intensity'].mean()

        ax4.axhline(cancer_mean_med, color='#E74C3C', linestyle='--', linewidth=2,
                   label=f'Cancer mean: {cancer_mean_med:.2e}', alpha=0.7)
        ax4.axhline(normal_mean_med, color='#3498DB', linestyle='--', linewidth=2,
                   label=f'Normal mean: {normal_mean_med:.2e}', alpha=0.7)

        ax4.set_xlabel('Sample', fontsize=10, fontweight='bold')
        ax4.set_ylabel('Median Intensity', fontsize=10, fontweight='bold')
        ax4.set_title('Central Tendency', fontsize=11, fontweight='bold')
        ax4.set_xticks(range(len(all_samples)))
        ax4.set_xticklabels(all_samples, rotation=90, fontsize=7)
        ax4.legend(loc='upper right', fontsize=8)
        ax4.grid(axis='y', alpha=0.3)
        ax4.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

        # =====================================================================
        # PANEL 5: CV Distribution (technical variability)
        # =====================================================================
        ax5 = axes[1, 1]
        logger.info("  Panel 5: CV distribution per sample...")

        bars5 = ax5.bar(range(len(all_samples)), qc_metrics['cv_median_pct'],
                       color=sample_colors, edgecolor='black', linewidth=0.5)

        # Add threshold lines
        ax5.axhline(20, color='green', linestyle='--', linewidth=2,
                   label='Good QC (<20%)', alpha=0.7)
        ax5.axhline(30, color='orange', linestyle='--', linewidth=2,
                   label='Acceptable (<30%)', alpha=0.7)

        ax5.set_xlabel('Sample', fontsize=10, fontweight='bold')
        ax5.set_ylabel('Median CV (%)', fontsize=10, fontweight='bold')
        ax5.set_title('Technical Variability', fontsize=11, fontweight='bold')
        ax5.set_xticks(range(len(all_samples)))
        ax5.set_xticklabels(all_samples, rotation=90, fontsize=7)
        ax5.legend(loc='upper right', fontsize=8)
        ax5.grid(axis='y', alpha=0.3)

        # =====================================================================
        # PANEL 6: Outlier Detection (Mahalanobis distance)
        # =====================================================================
        ax6 = axes[1, 2]
        logger.info("  Panel 6: Outlier detection (Mahalanobis distance)...")

        # Calculate Mahalanobis distance for outlier detection
        outlier_scores = self._calculate_mahalanobis_distance(qc_metrics)
        qc_metrics['mahalanobis_distance'] = outlier_scores

        bars6 = ax6.bar(range(len(all_samples)), outlier_scores,
                       color=sample_colors, edgecolor='black', linewidth=0.5)

        # Add threshold line at chi-square critical value (p=0.01, df=5)
        # 5 metrics: TIC, detection_count, detection_rate, median_intensity, cv_median
        threshold = stats.chi2.ppf(0.99, df=5)  # 99% confidence
        ax6.axhline(threshold, color='red', linestyle='--', linewidth=2,
                   label=f'Outlier threshold (p<0.01)', alpha=0.7)

        ax6.set_xlabel('Sample', fontsize=10, fontweight='bold')
        ax6.set_ylabel('Mahalanobis Distance', fontsize=10, fontweight='bold')
        ax6.set_title('Outlier Detection', fontsize=11, fontweight='bold')
        ax6.set_xticks(range(len(all_samples)))
        ax6.set_xticklabels(all_samples, rotation=90, fontsize=7)
        ax6.legend(loc='upper right', fontsize=8)
        ax6.grid(axis='y', alpha=0.3)

        # Highlight outliers
        outliers = np.where(outlier_scores > threshold)[0]
        if len(outliers) > 0:
            for outlier_idx in outliers:
                bars6[outlier_idx].set_edgecolor('red')
                bars6[outlier_idx].set_linewidth(3)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'sample_qc_dashboard.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved sample QC dashboard to {output_file}")

        # Save QC metrics as trace data
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
                             all_samples: list) -> pd.DataFrame:
        """
        Calculate QC metrics for each sample

        Args:
            intensity_matrix: Intensity matrix (glycopeptides x samples)
            all_samples: List of sample names

        Returns:
            DataFrame with QC metrics per sample
        """
        qc_data = []

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
            group = 'Cancer' if sample.startswith('C') else 'Normal'

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
