"""
Enhanced Pie Chart Plot Module for pGlyco Auto Combine
Publication-quality pie charts with fold change visualization

DESIGN INSPIRATION:
- MetaboAnalyst: Fold change bar charts with statistical annotations
- GraphPad Prism: Clean, professional appearance with bold colors
- Perseus: Clear comparative visualization

NEW FEATURES:
- Side-by-side pie charts (Cancer vs Normal)
- Fold change bar chart panel showing enrichment
- Color-coded bars (red=Cancer enrichment, blue=Normal enrichment)
- Statistical significance markers
- Publication-ready 300 DPI output
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests  # FDR correction

from ..utils import replace_empty_with_zero, save_trace_data
from .plot_config import (
    LEGACY_GLYCAN_COLORS,
    EXTENDED_CATEGORY_COLORS,
    GROUP_PALETTE,
    DPI,
    TITLE_SIZE,
    AXIS_LABEL_SIZE,
    TICK_LABEL_SIZE,
    LEGEND_SIZE,
    ANNOTATION_SIZE,
    save_publication_figure, DPI_MAIN,
    apply_publication_theme,  # ✨ Enhanced styling
    DESIGN_SYSTEM_AVAILABLE
)

# Import premium design system if available
if DESIGN_SYSTEM_AVAILABLE:
    from .design_system import VisualEffects

logger = logging.getLogger(__name__)


class PieChartPlotMixin:
    """Mixin class for enhanced pie chart visualizations with fold change"""

    def plot_pie_chart_glycan_types(self, df: pd.DataFrame, figsize: tuple = (16, 12)):
        """
        Create publication-quality pie charts with fold change visualization
        Shows glycan type distribution (Cancer vs Normal) + fold change bar chart

        Design inspired by MetaboAnalyst and GraphPad Prism

        Args:
            df: Annotated DataFrame with intensity data
            figsize: Figure size (width, height) - increased for fold change panel
        """
        logger.info("Creating enhanced glycan type pie charts with fold change...")

        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]

        # Calculate total intensities per glycan type for each group
        glycan_types = ['Non', 'Sialylated', 'Fucosylated', 'Both']

        cancer_data = {}
        normal_data = {}

        for glycan_type in glycan_types:
            mask = df['GlycanType'] == glycan_type

            # Sum intensities across all samples in group
            cancer_total = replace_empty_with_zero(df[mask][cancer_samples]).sum().sum()
            normal_total = replace_empty_with_zero(df[mask][normal_samples]).sum().sum()

            cancer_data[glycan_type] = cancer_total
            normal_data[glycan_type] = normal_total

        # Create figure with GridSpec for pie charts + fold change panel
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(2, 2, height_ratios=[2, 1], hspace=0.3, wspace=0.2)

        ax1 = fig.add_subplot(gs[0, 0])  # Cancer pie
        ax2 = fig.add_subplot(gs[0, 1])  # Normal pie
        ax3 = fig.add_subplot(gs[1, :])  # Fold change bar chart (spans both columns)

        # Define colors using standard palette
        colors = [LEGACY_GLYCAN_COLORS.get(gt, '#CCCCCC') for gt in glycan_types]

        # ===== Cancer pie chart =====
        cancer_values = [cancer_data[gt] for gt in glycan_types]
        cancer_total_sum = sum(cancer_values)

        wedges1, texts1, autotexts1 = ax1.pie(
            cancer_values,
            labels=glycan_types,
            colors=colors,
            autopct=lambda pct: f'{pct:.1f}%' if pct > 1 else '',
            startangle=90,
            textprops={'fontsize': TICK_LABEL_SIZE, 'weight': 'bold'}
        )

        ax1.set_title('Cancer Group\nGlycan Type Distribution',
                      fontsize=TITLE_SIZE, weight='bold', pad=15)

        # ===== Normal pie chart =====
        normal_values = [normal_data[gt] for gt in glycan_types]
        normal_total_sum = sum(normal_values)

        wedges2, texts2, autotexts2 = ax2.pie(
            normal_values,
            labels=glycan_types,
            colors=colors,
            autopct=lambda pct: f'{pct:.1f}%' if pct > 1 else '',
            startangle=90,
            textprops={'fontsize': TICK_LABEL_SIZE, 'weight': 'bold'}
        )

        ax2.set_title('Normal Group\nGlycan Type Distribution',
                      fontsize=TITLE_SIZE, weight='bold', pad=15)

        # ===== Fold change bar chart (MetaboAnalyst style) =====
        fold_changes = []
        fc_colors = []
        p_values = []

        for glycan_type in glycan_types:
            cancer_val = cancer_data[glycan_type]
            normal_val = normal_data[glycan_type]

            # ========================================
            # PHASE 1.3 FIX: Correct fold change calculation
            # ========================================
            # Calculate fold change (Cancer / Normal)
            if normal_val > 0 and cancer_val > 0:
                fc = cancer_val / normal_val
            elif normal_val == 0 and cancer_val > 0:
                # Cancer detected, Normal not detected → Maximum fold change
                fc = 100.0  # Cap at 100-fold (scientifically reasonable)
                logger.debug(f"{glycan_type}: Cancer={cancer_val:.2e}, Normal=0 → FC=100 (capped)")
            elif normal_val > 0 and cancer_val == 0:
                # Normal detected, Cancer not detected → Minimum fold change
                fc = 0.01  # Cap at 0.01 (100-fold decrease)
                logger.debug(f"{glycan_type}: Cancer=0, Normal={normal_val:.2e} → FC=0.01 (capped)")
            else:
                # Both zero → No change
                fc = 1.0
                logger.debug(f"{glycan_type}: Both groups zero → FC=1.0")

            fold_changes.append(fc)

            # Color by enrichment direction (Prism style)
            if fc > 1.1:  # >10% enrichment in Cancer
                fc_colors.append(GROUP_PALETTE['Cancer'])  # Red
            elif fc < 0.9:  # >10% enrichment in Normal
                fc_colors.append(GROUP_PALETTE['Normal'])  # Blue
            else:
                fc_colors.append('#95A5A6')  # Gray - similar

            # ========================================
            # PHASE 2.3 FIX: Explicit handling of insufficient samples
            # ========================================
            # Calculate statistical significance (Mann-Whitney U test)
            mask = df['GlycanType'] == glycan_type
            if mask.sum() > 0:
                cancer_intensities = replace_empty_with_zero(df[mask][cancer_samples]).values.flatten()
                normal_intensities = replace_empty_with_zero(df[mask][normal_samples]).values.flatten()

                # Remove zeros for statistical test
                cancer_nonzero = cancer_intensities[cancer_intensities > 0]
                normal_nonzero = normal_intensities[normal_intensities > 0]

                if len(cancer_nonzero) >= 3 and len(normal_nonzero) >= 3:
                    try:
                        _, p_val = scipy_stats.mannwhitneyu(cancer_nonzero, normal_nonzero, alternative='two-sided')
                        if np.isnan(p_val) or np.isinf(p_val):
                            logger.warning(f"{glycan_type}: Mann-Whitney U test returned invalid p-value ({p_val}), using p=1.0")
                            p_values.append(1.0)
                        else:
                            p_values.append(p_val)
                    except Exception as e:
                        logger.warning(f"{glycan_type}: Mann-Whitney U test failed ({str(e)}), using p=1.0")
                        p_values.append(1.0)
                else:
                    # Insufficient samples for valid test
                    logger.debug(f"{glycan_type}: Insufficient samples (Cancer: {len(cancer_nonzero)}, Normal: {len(normal_nonzero)}), p=1.0 (non-significant)")
                    p_values.append(1.0)
            else:
                # No glycopeptides of this type
                logger.debug(f"{glycan_type}: No glycopeptides found, p=1.0")
                p_values.append(1.0)

        # ========================================
        # PHASE 1.4 FIX: Apply FDR correction for multiple testing
        # ========================================
        # Apply Benjamini-Hochberg FDR correction
        valid_indices = [i for i, p in enumerate(p_values) if p < 1.0]
        valid_p_values = [p_values[i] for i in valid_indices]

        fdr_values = [1.0] * len(p_values)  # Default to 1.0
        if len(valid_p_values) > 0:
            _, fdr_corrected, _, _ = multipletests(valid_p_values, method='fdr_bh')
            for i, idx in enumerate(valid_indices):
                fdr_values[idx] = fdr_corrected[i]

            # Log FDR correction results
            n_significant = sum(fdr < 0.05 for fdr in fdr_corrected)
            logger.info(f"FDR correction applied (glycan types): {len(valid_p_values)} tests, "
                       f"{n_significant} significant (FDR < 0.05)")
            for gt, raw_p, fdr in zip(glycan_types, p_values, fdr_values):
                if raw_p < 1.0:
                    logger.info(f"  {gt}: p={raw_p:.4f} → FDR={fdr:.4f}")

        # Plot bars
        x_pos = np.arange(len(glycan_types))
        bars = ax3.bar(x_pos, fold_changes, color=fc_colors,
                       edgecolor='black', linewidth=1.5, alpha=0.9)

        # Add fold change values on top of bars
        for i, (bar, fc, fdr) in enumerate(zip(bars, fold_changes, fdr_values)):
            height = bar.get_height()

            # Fold change text
            fc_text = f'{fc:.2f}x'
            y_pos = height + 0.05

            ax3.text(bar.get_x() + bar.get_width() / 2., y_pos,
                     fc_text, ha='center', va='bottom',
                     fontsize=ANNOTATION_SIZE, weight='bold')

            # Significance marker using FDR (not raw p-value) - Prism style
            if fdr < 0.001:
                sig = '***'
            elif fdr < 0.01:
                sig = '**'
            elif fdr < 0.05:
                sig = '*'
            else:
                sig = ''

            if sig:
                ax3.text(bar.get_x() + bar.get_width() / 2., y_pos + 0.15,
                         sig, ha='center', va='bottom',
                         fontsize=14, weight='bold', color='black')

        # Add reference line at FC=1.0 (no change)
        ax3.axhline(y=1.0, color='black', linestyle='--', linewidth=2, alpha=0.7, label='No change')

        # Style fold change axis
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels(glycan_types, fontsize=TICK_LABEL_SIZE, weight='bold')
        ax3.set_ylabel('Fold Change (Cancer / Normal)', fontsize=AXIS_LABEL_SIZE, weight='bold')
        ax3.set_title('Fold Change Analysis (Cancer vs Normal)',
                      fontsize=TITLE_SIZE, weight='bold', pad=15)

        # Add legend for colors
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=GROUP_PALETTE['Cancer'], edgecolor='black', label='Enriched in Cancer (>10%)'),
            Patch(facecolor=GROUP_PALETTE['Normal'], edgecolor='black', label='Enriched in Normal (>10%)'),
            Patch(facecolor='#95A5A6', edgecolor='black', label='Similar (<10% difference)')
        ]
        ax3.legend(handles=legend_elements, loc='upper right', fontsize=LEGEND_SIZE)

        # Apply Prism-style grid
        ax3.grid(True, alpha=0.15, linestyle='-', linewidth=0.5, axis='y', zorder=0)
        ax3.set_axisbelow(True)

        # Overall title
        fig.suptitle('Glycan Type Distribution & Comparative Analysis',
                     fontsize=TITLE_SIZE + 2, weight='bold', y=0.98, family='Inter')

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'pie_chart_glycan_types_enhanced.png'
        save_publication_figure(fig, output_file, dpi=DPI_MAIN)
        logger.info(f"Saved enhanced glycan type pie charts to {output_file} (optimized, {DPI_MAIN} DPI)")

        # Save trace data
        trace_data = pd.DataFrame({
            'GlycanType': glycan_types,
            'Cancer_Intensity': cancer_values,
            'Cancer_Percentage': [v / cancer_total_sum * 100 for v in cancer_values],
            'Normal_Intensity': normal_values,
            'Normal_Percentage': [v / normal_total_sum * 100 for v in normal_values],
            'Fold_Change': fold_changes,
            'P_Value': p_values
        })
        save_trace_data(trace_data, self.output_dir, 'pie_chart_glycan_types_enhanced_data.csv')

        plt.close()

    def plot_pie_chart_primary_classification(self, df: pd.DataFrame, figsize: tuple = (16, 12)):
        """
        Create enhanced pie charts for primary classification with fold change

        Args:
            df: Annotated DataFrame with intensity data
            figsize: Figure size (width, height)
        """
        logger.info("Creating enhanced primary classification pie charts with fold change...")

        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]

        # Primary classification categories
        primary_categories = ['Truncated', 'High Mannose', 'ComplexHybrid', 'Outlier']

        cancer_data = {}
        normal_data = {}

        for category in primary_categories:
            mask = df['PrimaryClassification'] == category

            # Sum intensities across all samples in group
            cancer_total = replace_empty_with_zero(df[mask][cancer_samples]).sum().sum()
            normal_total = replace_empty_with_zero(df[mask][normal_samples]).sum().sum()

            cancer_data[category] = cancer_total
            normal_data[category] = normal_total

        # Create figure with GridSpec
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(2, 2, height_ratios=[2, 1], hspace=0.3, wspace=0.2)

        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, :])

        # Define colors
        colors_primary = {
            'Truncated': '#CCCCCC',       # Gray
            'High Mannose': '#2ECC71',    # Green
            'ComplexHybrid': '#3498DB',   # Blue
            'Outlier': '#95A5A6'          # Light gray
        }
        colors = [colors_primary.get(cat, '#CCCCCC') for cat in primary_categories]

        # Cancer pie chart
        cancer_values = [cancer_data[cat] for cat in primary_categories]
        cancer_total_sum = sum(cancer_values)

        wedges1, texts1, autotexts1 = ax1.pie(
            cancer_values,
            labels=primary_categories,
            colors=colors,
            autopct=lambda pct: f'{pct:.1f}%' if pct > 1 else '',
            startangle=90,
            textprops={'fontsize': TICK_LABEL_SIZE, 'weight': 'bold'}
        )

        ax1.set_title('Cancer Group\nPrimary Classification',
                      fontsize=TITLE_SIZE, weight='bold', pad=15)

        # Normal pie chart
        normal_values = [normal_data[cat] for cat in primary_categories]
        normal_total_sum = sum(normal_values)

        wedges2, texts2, autotexts2 = ax2.pie(
            normal_values,
            labels=primary_categories,
            colors=colors,
            autopct=lambda pct: f'{pct:.1f}%' if pct > 1 else '',
            startangle=90,
            textprops={'fontsize': TICK_LABEL_SIZE, 'weight': 'bold'}
        )

        ax2.set_title('Normal Group\nPrimary Classification',
                      fontsize=TITLE_SIZE, weight='bold', pad=15)

        # Fold change bar chart
        fold_changes = []
        fc_colors = []

        for i, category in enumerate(primary_categories):
            cancer_val = cancer_data[category]
            normal_val = normal_data[category]

            # ========================================
            # PHASE 1.3 FIX: Correct fold change calculation
            # ========================================
            if normal_val > 0 and cancer_val > 0:
                fc = cancer_val / normal_val
            elif normal_val == 0 and cancer_val > 0:
                fc = 100.0  # Cap at 100-fold
                logger.debug(f"{category}: Cancer={cancer_val:.2e}, Normal=0 → FC=100 (capped)")
            elif normal_val > 0 and cancer_val == 0:
                fc = 0.01  # Cap at 0.01 (100-fold decrease)
                logger.debug(f"{category}: Cancer=0, Normal={normal_val:.2e} → FC=0.01 (capped)")
            else:
                fc = 1.0  # Both zero
                logger.debug(f"{category}: Both groups zero → FC=1.0")

            fold_changes.append(fc)

            # Color by enrichment
            if fc > 1.1:
                fc_colors.append(GROUP_PALETTE['Cancer'])
            elif fc < 0.9:
                fc_colors.append(GROUP_PALETTE['Normal'])
            else:
                fc_colors.append('#95A5A6')

        # Plot bars
        x_pos = np.arange(len(primary_categories))
        bars = ax3.bar(x_pos, fold_changes, color=fc_colors,
                       edgecolor='black', linewidth=1.5, alpha=0.9)

        # Add fold change values
        for bar, fc in zip(bars, fold_changes):
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width() / 2., height + 0.05,
                     f'{fc:.2f}x', ha='center', va='bottom',
                     fontsize=11, weight='bold')

        # Reference line
        ax3.axhline(y=1.0, color='black', linestyle='--', linewidth=2, alpha=0.7)

        # Style axis
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels(primary_categories, fontsize=TICK_LABEL_SIZE, weight='bold', rotation=15, ha='right')
        ax3.set_ylabel('Fold Change (Cancer / Normal)', fontsize=AXIS_LABEL_SIZE, weight='bold')
        ax3.set_title('Fold Change Analysis (Cancer vs Normal)',
                      fontsize=TITLE_SIZE, weight='bold', pad=15)
        ax3.grid(True, alpha=0.15, linestyle='-', linewidth=0.5, axis='y', zorder=0)
        ax3.set_axisbelow(True)

        # Overall title
        fig.suptitle('Primary Classification Distribution & Comparative Analysis',
                     fontsize=TITLE_SIZE + 2, weight='bold', y=0.98, family='Inter')

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'pie_chart_primary_classification_enhanced.png'
        save_publication_figure(fig, output_file, dpi=DPI_MAIN)
        logger.info(f"Saved enhanced primary classification pie charts to {output_file} (optimized, {DPI_MAIN} DPI)")

        # Save trace data
        trace_data = pd.DataFrame({
            'PrimaryClassification': primary_categories,
            'Cancer_Intensity': cancer_values,
            'Cancer_Percentage': [v / cancer_total_sum * 100 if cancer_total_sum > 0 else 0 for v in cancer_values],
            'Normal_Intensity': normal_values,
            'Normal_Percentage': [v / normal_total_sum * 100 if normal_total_sum > 0 else 0 for v in normal_values],
            'Fold_Change': fold_changes
        })
        save_trace_data(trace_data, self.output_dir, 'pie_chart_primary_classification_enhanced_data.csv')

        plt.close()

    def plot_pie_chart_secondary_classification(self, df: pd.DataFrame, figsize: tuple = (18, 12)):
        """
        Create enhanced pie charts for secondary classification with fold change

        Args:
            df: Annotated DataFrame with intensity data
            figsize: Figure size (width, height)
        """
        logger.info("Creating enhanced secondary classification pie charts with fold change...")

        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]

        # Secondary classification categories (5 categories)
        secondary_categories = ['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        cancer_data = {}
        normal_data = {}

        for category in secondary_categories:
            mask = df['SecondaryClassification'] == category

            # Sum intensities across all samples in group
            cancer_total = replace_empty_with_zero(df[mask][cancer_samples]).sum().sum()
            normal_total = replace_empty_with_zero(df[mask][normal_samples]).sum().sum()

            cancer_data[category] = cancer_total
            normal_data[category] = normal_total

        # Create figure with GridSpec
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(2, 2, height_ratios=[2, 1], hspace=0.3, wspace=0.2)

        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, :])

        # Define colors (using standardized palette)
        colors = [EXTENDED_CATEGORY_COLORS.get(cat, '#CCCCCC') for cat in secondary_categories]

        # Cancer pie chart
        cancer_values = [cancer_data[cat] for cat in secondary_categories]
        cancer_total_sum = sum(cancer_values)

        wedges1, texts1, autotexts1 = ax1.pie(
            cancer_values,
            labels=secondary_categories,
            colors=colors,
            autopct=lambda pct: f'{pct:.1f}%' if pct > 1 else '',
            startangle=90,
            textprops={'fontsize': 10, 'weight': 'bold'}
        )

        ax1.set_title('Cancer Group\nSecondary Classification',
                      fontsize=TITLE_SIZE, weight='bold', pad=15)

        # Normal pie chart
        normal_values = [normal_data[cat] for cat in secondary_categories]
        normal_total_sum = sum(normal_values)

        wedges2, texts2, autotexts2 = ax2.pie(
            normal_values,
            labels=secondary_categories,
            colors=colors,
            autopct=lambda pct: f'{pct:.1f}%' if pct > 1 else '',
            startangle=90,
            textprops={'fontsize': 10, 'weight': 'bold'}
        )

        ax2.set_title('Normal Group\nSecondary Classification',
                      fontsize=TITLE_SIZE, weight='bold', pad=15)

        # Fold change bar chart
        fold_changes = []
        fc_colors = []

        for category in secondary_categories:
            cancer_val = cancer_data[category]
            normal_val = normal_data[category]

            # ========================================
            # PHASE 1.3 FIX: Correct fold change calculation
            # ========================================
            if normal_val > 0 and cancer_val > 0:
                fc = cancer_val / normal_val
            elif normal_val == 0 and cancer_val > 0:
                fc = 100.0  # Cap at 100-fold
                logger.debug(f"{category}: Cancer={cancer_val:.2e}, Normal=0 → FC=100 (capped)")
            elif normal_val > 0 and cancer_val == 0:
                fc = 0.01  # Cap at 0.01 (100-fold decrease)
                logger.debug(f"{category}: Cancer=0, Normal={normal_val:.2e} → FC=0.01 (capped)")
            else:
                fc = 1.0  # Both zero
                logger.debug(f"{category}: Both groups zero → FC=1.0")

            fold_changes.append(fc)

            # Color by enrichment
            if fc > 1.1:
                fc_colors.append(GROUP_PALETTE['Cancer'])
            elif fc < 0.9:
                fc_colors.append(GROUP_PALETTE['Normal'])
            else:
                fc_colors.append('#95A5A6')

        # Plot bars
        x_pos = np.arange(len(secondary_categories))
        bars = ax3.bar(x_pos, fold_changes, color=fc_colors,
                       edgecolor='black', linewidth=1.5, alpha=0.9)

        # Add fold change values
        for bar, fc in zip(bars, fold_changes):
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width() / 2., height + 0.05,
                     f'{fc:.2f}x', ha='center', va='bottom',
                     fontsize=11, weight='bold')

        # Reference line
        ax3.axhline(y=1.0, color='black', linestyle='--', linewidth=2, alpha=0.7)

        # Style axis
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels(secondary_categories, fontsize=TICK_LABEL_SIZE, weight='bold', rotation=20, ha='right')
        ax3.set_ylabel('Fold Change (Cancer / Normal)', fontsize=AXIS_LABEL_SIZE, weight='bold')
        ax3.set_title('Fold Change Analysis (Cancer vs Normal)',
                      fontsize=TITLE_SIZE, weight='bold', pad=15)
        ax3.grid(True, alpha=0.15, linestyle='-', linewidth=0.5, axis='y', zorder=0)
        ax3.set_axisbelow(True)

        # Overall title
        fig.suptitle('Secondary Classification Distribution & Comparative Analysis',
                     fontsize=TITLE_SIZE + 2, weight='bold', y=0.98, family='Inter')

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'pie_chart_secondary_classification_enhanced.png'
        save_publication_figure(fig, output_file, dpi=DPI_MAIN)
        logger.info(f"Saved enhanced secondary classification pie charts to {output_file} (optimized, {DPI_MAIN} DPI)")

        # Save trace data
        trace_data = pd.DataFrame({
            'SecondaryClassification': secondary_categories,
            'Cancer_Intensity': cancer_values,
            'Cancer_Percentage': [v / cancer_total_sum * 100 if cancer_total_sum > 0 else 0 for v in cancer_values],
            'Normal_Intensity': normal_values,
            'Normal_Percentage': [v / normal_total_sum * 100 if normal_total_sum > 0 else 0 for v in normal_values],
            'Fold_Change': fold_changes
        })
        save_trace_data(trace_data, self.output_dir, 'pie_chart_secondary_classification_enhanced_data.csv')

        plt.close()

    def plot_pie_chart_significant_glycan_types(self, df: pd.DataFrame, vip_df: pd.DataFrame = None,
                                                 config=None, log2fc_threshold: float = 2.0,
                                                 fdr_threshold: float = 0.05, figsize: tuple = (16, 12)):
        """
        Create publication-quality pie charts showing ONLY highly significant glycan types

        Filters glycopeptides by |Log2 FC| >= 2 (4-fold change) AND FDR < 0.05
        Shows distribution of ONLY these highly significant features

        Design inspired by MetaboAnalyst and GraphPad Prism

        Args:
            df: Annotated DataFrame with intensity data
            vip_df: Optional VIP scores DataFrame
            config: DataPreparationConfig (uses default if None)
            log2fc_threshold: Minimum |Log2 FC| for significance (default: 2 = 4-fold change)
            fdr_threshold: FDR significance threshold (default: 0.05)
            figsize: Figure size (width, height)
        """
        from ..data_preparation import (
            DataPreparationConfig,
            prepare_visualization_data,
            calculate_statistical_significance
        )
        from ..utils import get_sample_columns

        logger.info(f"Creating pie charts for HIGHLY SIGNIFICANT glycan types (|Log2 FC| >= {log2fc_threshold})...")

        # Use default config if not provided
        if config is None:
            config = DataPreparationConfig(
                min_detection_pct=0.30,
                min_samples=5,
                missing_data_method='skipna'
            )

        # Prepare data with centralized statistics
        df_prep = prepare_visualization_data(
            df=df,
            config=config,
            vip_scores=vip_df,
            merge_method='left',
            apply_detection_filter=False,  # Already filtered in main.py
            log_prefix="[Significant Pie Chart] "
        )

        # Get sample columns
        cancer_samples, normal_samples = get_sample_columns(df)

        # Calculate statistical significance (FDR correction)
        df_prep = calculate_statistical_significance(
            df_prep=df_prep,
            cancer_samples=cancer_samples,
            normal_samples=normal_samples,
            method='mannwhitneyu',
            fdr_correction=True
        )

        # CRITICAL FILTER: Separate by direction of change
        # Upregulated: Log2 FC >= +2 (Cancer > Normal)
        # Downregulated: Log2 FC <= -2 (Cancer < Normal)
        upregulated_mask = (df_prep['Log2_Fold_Change'] >= log2fc_threshold) & (df_prep['FDR'] < fdr_threshold)
        downregulated_mask = (df_prep['Log2_Fold_Change'] <= -log2fc_threshold) & (df_prep['FDR'] < fdr_threshold)

        df_upregulated = df_prep[upregulated_mask].copy()
        df_downregulated = df_prep[downregulated_mask].copy()

        n_upregulated = len(df_upregulated)
        n_downregulated = len(df_downregulated)
        n_total = n_upregulated + n_downregulated

        logger.info(f"  Found {n_total} glycopeptides with |Log2 FC| >= {log2fc_threshold} and FDR < {fdr_threshold}")
        logger.info(f"    Upregulated in Cancer (Log2 FC >= +{log2fc_threshold}): {n_upregulated} features")
        logger.info(f"    Downregulated in Cancer (Log2 FC <= -{log2fc_threshold}): {n_downregulated} features")

        if n_total == 0:
            logger.warning("  No significant glycopeptides found with current thresholds!")
            logger.warning(f"  Try lowering log2fc_threshold (current: {log2fc_threshold}) or fdr_threshold (current: {fdr_threshold})")
            return

        # Aggregate by glycan type category (5 categories) for each direction
        glycan_types = ['HM', 'C/H', 'F', 'S', 'SF']

        upregulated_counts = {}
        downregulated_counts = {}
        upregulated_intensities = {}
        downregulated_intensities = {}

        for glycan_type in glycan_types:
            # Upregulated features
            up_mask = df_upregulated['GlycanTypeCategory'] == glycan_type
            upregulated_counts[glycan_type] = up_mask.sum()
            upregulated_intensities[glycan_type] = replace_empty_with_zero(
                df_upregulated[up_mask][cancer_samples + normal_samples]
            ).sum().sum() if up_mask.sum() > 0 else 0

            # Downregulated features
            down_mask = df_downregulated['GlycanTypeCategory'] == glycan_type
            downregulated_counts[glycan_type] = down_mask.sum()
            downregulated_intensities[glycan_type] = replace_empty_with_zero(
                df_downregulated[down_mask][cancer_samples + normal_samples]
            ).sum().sum() if down_mask.sum() > 0 else 0

        logger.info("  Glycan type distribution (by direction):")
        logger.info(f"    Upregulated: {dict(upregulated_counts)}")
        logger.info(f"    Downregulated: {dict(downregulated_counts)}")

        # Create figure with GridSpec for pie charts + fold change panel
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(2, 2, height_ratios=[2, 1], hspace=0.3, wspace=0.2)

        ax1 = fig.add_subplot(gs[0, 0])  # Upregulated pie
        ax2 = fig.add_subplot(gs[0, 1])  # Downregulated pie
        ax3 = fig.add_subplot(gs[1, :])  # Direction comparison bar chart

        # Define colors using user-requested 5-category palette
        colors = [EXTENDED_CATEGORY_COLORS.get(gt, '#CCCCCC') for gt in glycan_types]

        # ===== Upregulated pie chart (Cancer > Normal) =====
        upregulated_values = [upregulated_counts[gt] for gt in glycan_types]
        upregulated_total = sum(upregulated_values)

        if upregulated_total > 0:
            wedges1, texts1, autotexts1 = ax1.pie(
                upregulated_values,
                labels=glycan_types,
                colors=colors,
                autopct=lambda pct: f'{pct:.1f}%' if pct > 1 else '',
                startangle=90,
                textprops={'fontsize': TICK_LABEL_SIZE, 'weight': 'bold'}
            )
            ax1.set_title(f'Upregulated in Cancer\nLog2 FC ≥ +{log2fc_threshold}\n(n={n_upregulated} features)',
                          fontsize=TITLE_SIZE, weight='bold', pad=15)
        else:
            ax1.text(0.5, 0.5, 'No upregulated\nfeatures found',
                    ha='center', va='center', fontsize=14, transform=ax1.transAxes)
            ax1.set_title(f'Upregulated in Cancer\nLog2 FC ≥ +{log2fc_threshold}\n(n=0 features)',
                          fontsize=TITLE_SIZE, weight='bold', pad=15)

        # ===== Downregulated pie chart (Cancer < Normal) =====
        downregulated_values = [downregulated_counts[gt] for gt in glycan_types]
        downregulated_total = sum(downregulated_values)

        if downregulated_total > 0:
            wedges2, texts2, autotexts2 = ax2.pie(
                downregulated_values,
                labels=glycan_types,
                colors=colors,
                autopct=lambda pct: f'{pct:.1f}%' if pct > 1 else '',
                startangle=90,
                textprops={'fontsize': TICK_LABEL_SIZE, 'weight': 'bold'}
            )
            ax2.set_title(f'Downregulated in Cancer\nLog2 FC ≤ -{log2fc_threshold}\n(n={n_downregulated} features)',
                          fontsize=TITLE_SIZE, weight='bold', pad=15)
        else:
            ax2.text(0.5, 0.5, 'No downregulated\nfeatures found',
                    ha='center', va='center', fontsize=14, transform=ax2.transAxes)
            ax2.set_title(f'Downregulated in Cancer\nLog2 FC ≤ -{log2fc_threshold}\n(n=0 features)',
                          fontsize=TITLE_SIZE, weight='bold', pad=15)

        # ===== Directional comparison bar chart =====
        # Show feature counts by direction for each glycan type
        x_pos = np.arange(len(glycan_types))
        width = 0.35  # Width of bars

        # Create grouped bars
        up_values = [upregulated_counts[gt] for gt in glycan_types]
        down_values = [downregulated_counts[gt] for gt in glycan_types]

        bars1 = ax3.bar(x_pos - width/2, up_values, width,
                       label=f'Upregulated (n={n_upregulated})',
                       color=GROUP_PALETTE['Cancer'], edgecolor='black',
                       linewidth=1.5, alpha=0.9)
        bars2 = ax3.bar(x_pos + width/2, down_values, width,
                       label=f'Downregulated (n={n_downregulated})',
                       color=GROUP_PALETTE['Normal'], edgecolor='black',
                       linewidth=1.5, alpha=0.9)

        # Add count labels on bars
        for bar in bars1:
            height = bar.get_height()
            if height > 0:
                ax3.text(bar.get_x() + bar.get_width() / 2., height,
                        f'{int(height)}',
                        ha='center', va='bottom',
                        fontsize=ANNOTATION_SIZE, weight='bold')

        for bar in bars2:
            height = bar.get_height()
            if height > 0:
                ax3.text(bar.get_x() + bar.get_width() / 2., height,
                        f'{int(height)}',
                        ha='center', va='bottom',
                        fontsize=ANNOTATION_SIZE, weight='bold')

        # Style axis
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels(glycan_types, fontsize=TICK_LABEL_SIZE, weight='bold')
        ax3.set_ylabel('Number of Features', fontsize=AXIS_LABEL_SIZE, weight='bold')
        ax3.set_title(f'Directional Comparison by Glycan Type\n|Log2 FC| ≥ {log2fc_threshold}, FDR < {fdr_threshold}',
                      fontsize=TITLE_SIZE, weight='bold', pad=15)
        ax3.legend(loc='upper right', fontsize=LEGEND_SIZE)

        # Apply Prism-style grid
        ax3.grid(True, alpha=0.15, linestyle='-', linewidth=0.5, axis='y', zorder=0)
        ax3.set_axisbelow(True)

        # Overall title with emphasis on "DIRECTIONAL" (moved above figure)
        fig.suptitle(f'Directional Analysis: Highly Significant Glycan Types\n|Log2 FC| ≥ {log2fc_threshold} (≥4-fold change) | FDR < {fdr_threshold}',
                     fontsize=TITLE_SIZE + 2, weight='bold', y=1.01, family='Inter')

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'pie_chart_glycan_types_significant.png'
        save_publication_figure(fig, output_file, dpi=DPI_MAIN)
        logger.info(f"✨ Saved SIGNIFICANT glycan type pie charts to {output_file} (optimized, {DPI_MAIN} DPI)")

        # Save comprehensive trace data (directional)
        trace_data = pd.DataFrame({
            'GlycanType': glycan_types,
            'N_Upregulated': [upregulated_counts[gt] for gt in glycan_types],
            'N_Downregulated': [downregulated_counts[gt] for gt in glycan_types],
            'N_Total': [upregulated_counts[gt] + downregulated_counts[gt] for gt in glycan_types],
            'Upregulated_Percentage': [upregulated_counts[gt] / n_upregulated * 100 if n_upregulated > 0 else 0 for gt in glycan_types],
            'Downregulated_Percentage': [downregulated_counts[gt] / n_downregulated * 100 if n_downregulated > 0 else 0 for gt in glycan_types],
            'Log2FC_Threshold': log2fc_threshold,
            'FDR_Threshold': fdr_threshold
        })
        save_trace_data(trace_data, self.output_dir, 'pie_chart_glycan_types_significant_data.csv')

        # Save list of upregulated glycopeptides
        if n_upregulated > 0:
            upregulated_features = df_upregulated[['Peptide', 'GlycanComposition', 'GlycanTypeCategory',
                                                   'Log2_Fold_Change', 'FDR', 'Cancer_Mean', 'Normal_Mean']].copy()
            upregulated_features = upregulated_features.sort_values('FDR')
            save_trace_data(upregulated_features, self.output_dir, 'significant_glycopeptides_upregulated.csv')
            logger.info(f"  Saved {len(upregulated_features)} upregulated glycopeptides to significant_glycopeptides_upregulated.csv")

        # Save list of downregulated glycopeptides
        if n_downregulated > 0:
            downregulated_features = df_downregulated[['Peptide', 'GlycanComposition', 'GlycanTypeCategory',
                                                       'Log2_Fold_Change', 'FDR', 'Cancer_Mean', 'Normal_Mean']].copy()
            downregulated_features = downregulated_features.sort_values('FDR')
            save_trace_data(downregulated_features, self.output_dir, 'significant_glycopeptides_downregulated.csv')
            logger.info(f"  Saved {len(downregulated_features)} downregulated glycopeptides to significant_glycopeptides_downregulated.csv")

        # Also save combined list (for backward compatibility)
        all_significant = pd.concat([df_upregulated, df_downregulated]) if n_total > 0 else pd.DataFrame()
        if len(all_significant) > 0:
            all_significant = all_significant[['Peptide', 'GlycanComposition', 'GlycanTypeCategory',
                                              'Log2_Fold_Change', 'FDR', 'Cancer_Mean', 'Normal_Mean']].copy()
            all_significant = all_significant.sort_values('FDR')
            save_trace_data(all_significant, self.output_dir, 'significant_glycopeptides_list.csv')
            logger.info(f"  Saved combined list of {len(all_significant)} significant glycopeptides to significant_glycopeptides_list.csv")

        plt.close()
