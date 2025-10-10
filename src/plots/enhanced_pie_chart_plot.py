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

            # Calculate fold change (Cancer / Normal)
            if normal_val > 0:
                fc = cancer_val / normal_val
            else:
                fc = cancer_val if cancer_val > 0 else 1.0

            fold_changes.append(fc)

            # Color by enrichment direction (Prism style)
            if fc > 1.1:  # >10% enrichment in Cancer
                fc_colors.append(GROUP_PALETTE['Cancer'])  # Red
            elif fc < 0.9:  # >10% enrichment in Normal
                fc_colors.append(GROUP_PALETTE['Normal'])  # Blue
            else:
                fc_colors.append('#95A5A6')  # Gray - similar

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
                        p_values.append(p_val)
                    except Exception:
                        p_values.append(1.0)
                else:
                    p_values.append(1.0)
            else:
                p_values.append(1.0)

        # Plot bars
        x_pos = np.arange(len(glycan_types))
        bars = ax3.bar(x_pos, fold_changes, color=fc_colors,
                       edgecolor='black', linewidth=1.5, alpha=0.9)

        # Add fold change values on top of bars
        for i, (bar, fc, p_val) in enumerate(zip(bars, fold_changes, p_values)):
            height = bar.get_height()

            # Fold change text
            fc_text = f'{fc:.2f}x'
            y_pos = height + 0.05

            ax3.text(bar.get_x() + bar.get_width() / 2., y_pos,
                     fc_text, ha='center', va='bottom',
                     fontsize=ANNOTATION_SIZE, weight='bold')

            # Significance marker (Prism style)
            if p_val < 0.001:
                sig = '***'
            elif p_val < 0.01:
                sig = '**'
            elif p_val < 0.05:
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

            if normal_val > 0:
                fc = cancer_val / normal_val
            else:
                fc = cancer_val if cancer_val > 0 else 1.0

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

            if normal_val > 0:
                fc = cancer_val / normal_val
            else:
                fc = cancer_val if cancer_val > 0 else 1.0

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
