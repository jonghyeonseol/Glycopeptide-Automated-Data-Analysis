"""
Missing Data Visualization Module for pGlyco Auto Combine
Creates missingno-style matrix to show data completeness patterns

Dependencies:
    External:
        - pandas: Data manipulation
        - numpy: Numerical computations
        - matplotlib: Plotting backend, GridSpec layout

    Internal:
        - src.utils: save_trace_data, get_sample_columns
        - src.plots.plot_config: save_publication_figure

Phase 4.1: Critical for data integrity and transparency
Shows where data is missing and validates filtering decisions
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
from ..utils import save_trace_data, get_sample_columns
from .plot_config import (
    save_publication_figure, COLOR_CANCER, COLOR_NORMAL, DPI_MAIN,
    TITLE_SIZE, AXIS_LABEL_SIZE, TICK_LABEL_SIZE, LEGEND_SIZE, ANNOTATION_SIZE,
    LINE_THIN, LINE_THICK,  # Linewidth constants
    ALPHA_MEDIUM_LIGHT, ALPHA_VERY_HIGH, ALPHA_MEDIUM_HIGH,  # Alpha constants
    EDGE_COLOR_NONE,
    EDGE_COLOR_BLACK,# Edge color standardization,
    # Font family constants (Phase 10.3.10)
    FONT_DATA,
    # Linestyle constants (Phase 10.3.8)
    THRESHOLD_LINESTYLE,
    # Zorder constants (Phase 10.3.7)
    ZORDER_BACKGROUND, ZORDER_GRID, ZORDER_SEPARATOR,
    ZORDER_DATA_LOW, ZORDER_DATA_HIGH,
    ZORDER_THRESHOLD, ZORDER_ANNOTATION,
    ZORDER_OVERLAY, ZORDER_EFFECT,
    ZORDER_TOP, ZORDER_ABSOLUTE_TOP
)

logger = logging.getLogger(__name__)


class MissingDataPlotMixin:
    """Mixin class for missing data visualization"""

    def plot_missing_data_matrix(self, df: pd.DataFrame, figsize: tuple = (16, 10)):
        """
        Create missing data matrix visualization (missingno-style)

        Shows:
        - Which samples have missing data for which glycopeptides
        - Detection frequency per sample (bottom bar chart)
        - Detection frequency per glycopeptide (right bar chart)

        CRITICAL FOR DATA INTEGRITY:
        - Identifies systematic missingness patterns
        - Validates 30% detection filter
        - Shows if any samples are problematic (low detection)

        Args:
            df: Annotated DataFrame with intensity data
            figsize: Figure size (width, height)
        """
        logger.info("Creating missing data matrix visualization...")

        # Get sample columns
        cancer_samples, normal_samples = get_sample_columns(df)
        all_samples = cancer_samples + normal_samples

        # Extract intensity matrix
        intensity_matrix = df[all_samples].copy()

        # Convert empty strings to NaN for missing data
        intensity_matrix = intensity_matrix.replace('', np.nan)
        intensity_matrix = intensity_matrix.astype(float)

        # Create binary matrix (1 = detected, 0 = missing)
        detection_matrix = (~intensity_matrix.isna()).astype(int)

        # Calculate statistics
        n_glycopeptides = len(detection_matrix)
        n_samples = len(all_samples)
        total_possible = n_glycopeptides * n_samples
        total_detected = detection_matrix.sum().sum()
        overall_detection_rate = total_detected / total_possible * 100

        logger.info(f"Overall detection rate: {overall_detection_rate:.1f}%")
        logger.info(f"Total measurements: {total_detected:,} / {total_possible:,}")

        # Sample detection rates
        sample_detection = detection_matrix.sum(axis=0) / n_glycopeptides * 100

        # Glycopeptide detection rates
        glycopeptide_detection = detection_matrix.sum(axis=1) / n_samples * 100

        # Create figure with GridSpec for better layout
        from matplotlib.gridspec import GridSpec
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(3, 3, figure=fig, height_ratios=[0.15, 1, 0.15],
                      width_ratios=[1, 0.02, 0.15], hspace=0.05, wspace=0.05)

        # Main heatmap (center)
        ax_main = fig.add_subplot(gs[1, 0])

        # Sample detection bar chart (top)
        ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)

        # Glycopeptide detection bar chart (right)
        ax_right = fig.add_subplot(gs[1, 2], sharey=ax_main)

        # Colorbar axis
        ax_cbar = fig.add_subplot(gs[1, 1])

        # =====================================================================
        # MAIN HEATMAP - Subsample for visibility
        # =====================================================================
        # For large datasets, show every Nth glycopeptide for clarity
        if n_glycopeptides > 500:
            step = n_glycopeptides // 500
            display_matrix = detection_matrix.iloc[::step, :]
            display_detection = glycopeptide_detection.iloc[::step]
            logger.info(f"Subsampling: showing {len(display_matrix)} of {n_glycopeptides} glycopeptides")
        else:
            display_matrix = detection_matrix
            display_detection = glycopeptide_detection

        # Plot heatmap (white = missing, color = detected)
        im = ax_main.imshow(display_matrix.values, aspect='auto', cmap='RdYlGn',
                            interpolation='nearest', vmin=0, vmax=1)

        # Sample group colors (Cancer vs Normal)
        sample_colors = [COLOR_CANCER if s in cancer_samples else COLOR_NORMAL for s in all_samples]

        # Add sample group color bar at top
        for i, (sample, color) in enumerate(zip(all_samples, sample_colors)):
            ax_main.add_patch(plt.Rectangle((i - 0.5, -0.5), 1, 0.3,
                                            facecolor=color, edgecolor=EDGE_COLOR_NONE,
                                            clip_on=False, transform=ax_main.transData))

        # Labels and styling
        ax_main.set_xlabel('Sample', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax_main.set_ylabel('Glycopeptide (subsampled)', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax_main.set_xticks(range(len(all_samples)))
        ax_main.set_xticklabels(all_samples, rotation=90, fontsize=TICK_LABEL_SIZE)
        ax_main.set_yticks([])  # Too many to show
        ax_main.tick_params(axis='x', which='both', bottom=False, top=False)

        # =====================================================================
        # TOP BAR CHART - Sample detection rates
        # =====================================================================
        _ = ax_top.bar(range(len(all_samples)), sample_detection.values,
                       color=sample_colors, edgecolor=EDGE_COLOR_BLACK, linewidth=LINE_THIN)

        # Add threshold line at 30%
        ax_top.axhline(30, color='red', linestyle=THRESHOLD_LINESTYLE, linewidth=LINE_THICK,
                       label='30% filter threshold', zorder=ZORDER_THRESHOLD)

        ax_top.set_ylim(0, 100)
        ax_top.set_ylabel('Detection\nRate (%)', fontsize=ANNOTATION_SIZE, fontweight='bold')
        ax_top.set_title('Missing Data Matrix: Detection Completeness',
                         fontsize=TITLE_SIZE, fontweight='bold', pad=15)
        ax_top.legend(loc='upper right', fontsize=LEGEND_SIZE)
        ax_top.grid(axis='y', alpha=ALPHA_MEDIUM_LIGHT)
        plt.setp(ax_top.get_xticklabels(), visible=False)
        ax_top.tick_params(axis='x', which='both', bottom=False)

        # =====================================================================
        # RIGHT BAR CHART - Glycopeptide detection rates
        # =====================================================================
        ax_right.barh(range(len(display_detection)), display_detection.values,
                      color='#27AE60', edgecolor=EDGE_COLOR_BLACK, linewidth=LINE_THIN)

        # Add threshold line at 30%
        ax_right.axvline(30, color='red', linestyle=THRESHOLD_LINESTYLE, linewidth=LINE_THICK, zorder=ZORDER_THRESHOLD)

        ax_right.set_xlim(0, 100)
        ax_right.set_xlabel('Detection\nRate (%)', fontsize=ANNOTATION_SIZE, fontweight='bold')
        ax_right.grid(axis='x', alpha=ALPHA_MEDIUM_LIGHT)
        plt.setp(ax_right.get_yticklabels(), visible=False)
        ax_right.tick_params(axis='y', which='both', left=False)

        # =====================================================================
        # COLORBAR
        # =====================================================================
        cbar = plt.colorbar(im, cax=ax_cbar)
        cbar.set_ticks([0, 1])
        cbar.set_ticklabels(['Missing', 'Detected'])
        cbar.ax.tick_params(labelsize=10)

        # =====================================================================
        # STATISTICS TEXT BOX
        # =====================================================================
        stats_text = f"Overall Detection: {overall_detection_rate:.1f}%\n"
        stats_text += f"Total: {total_detected:,} / {total_possible:,}\n"
        stats_text += f"Glycopeptides: {n_glycopeptides}\n"
        stats_text += f"Samples: {n_samples} ({len(cancer_samples)}C, {len(normal_samples)}N)"

        ax_main.text(0.02, 0.98, stats_text, transform=ax_main.transAxes,
                     fontsize=ANNOTATION_SIZE, verticalalignment='top',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=ALPHA_VERY_HIGH),
                     family=FONT_DATA, zorder=ZORDER_TOP)

        # =====================================================================
        # DATA INTEGRITY NOTES
        # =====================================================================
        integrity_text = "DATA INTEGRITY VALIDATION:\n"

        # Check for problematic samples (< 10% detection)
        low_detection_samples = sample_detection[sample_detection < 10]
        if len(low_detection_samples) > 0:
            integrity_text += f"⚠ {len(low_detection_samples)} samples with <10% detection\n"
        else:
            integrity_text += "✓ All samples pass 10% detection threshold\n"

        # Check for Cancer vs Normal balance
        cancer_mean = sample_detection.loc[cancer_samples].mean()
        normal_mean = sample_detection.loc[normal_samples].mean()
        if abs(cancer_mean - normal_mean) > 10:
            integrity_text += f"⚠ Detection imbalance: C={cancer_mean:.1f}%, N={normal_mean:.1f}%\n"
        else:
            integrity_text += f"✓ Balanced detection: C={cancer_mean:.1f}%, N={normal_mean:.1f}%"

        fig.text(0.02, 0.02, integrity_text, fontsize=ANNOTATION_SIZE, family=FONT_DATA,
                 verticalalignment='bottom',
                 bbox=dict(boxstyle='round', facecolor='lightblue', alpha=ALPHA_MEDIUM_HIGH))

        # Save plot using standardized function
        output_file = self.output_dir / 'missing_data_matrix.png'
        save_publication_figure(fig, output_file, dpi=DPI_MAIN)
        logger.info(f"Saved missing data matrix to {output_file} (optimized, {DPI_MAIN} DPI)")

        # Save trace data - detection rates per sample
        trace_data = pd.DataFrame({
            'Sample': all_samples,
            'Group': ['Cancer' if s in cancer_samples else 'Normal' for s in all_samples],
            'Detection_Rate_Percent': sample_detection.values,
            'Detected_Count': detection_matrix.sum(axis=0).values,
            'Total_Possible': n_glycopeptides
        })
        save_trace_data(trace_data, self.output_dir, 'missing_data_sample_stats.csv')

        # Save glycopeptide-level detection rates
        glycopeptide_trace = pd.DataFrame({
            'Glycopeptide_Index': range(n_glycopeptides),
            'Detection_Rate_Percent': glycopeptide_detection.values,
            'Detected_Samples': detection_matrix.sum(axis=1).values,
            'Total_Samples': n_samples
        })
        save_trace_data(glycopeptide_trace, self.output_dir, 'missing_data_glycopeptide_stats.csv')

        plt.close()

        logger.info("✓ Missing data analysis complete - validates data integrity")
