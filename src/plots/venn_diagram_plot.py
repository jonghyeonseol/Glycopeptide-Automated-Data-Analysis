"""
Venn Diagram Plot Module for pGlyco Auto Combine
Visualizes overlap between glycan type categories

Dependencies:
    External:
        - pandas: Data manipulation
        - matplotlib: Plotting backend
        - matplotlib_venn: Venn diagram visualization (venn3, venn3_circles)

    Internal:
        - src.utils: save_trace_data
        - src.plots.plot_config: VENN_* constants, save_publication_figure
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import logging

from ..utils import save_trace_data
from .plot_config import (
    VENN_FIGSIZE, VENN_COLORS, VENN_ALPHA, VENN_LINEWIDTH, VENN_LINESTYLE,
    VENN_LABEL_FONTSIZE, VENN_SUBSET_FONTSIZE, VENN_TITLE_FONTSIZE,
    VENN_STATS_FONTSIZE, VENN_SUBSET_LABEL_THRESHOLD, VENN_DPI,
    save_publication_figure
)

logger = logging.getLogger(__name__)


def format_venn_subset_label(value: int, use_compact: bool = False) -> str:
    """
    Format Venn diagram subset labels with smart formatting

    Prevents label overlap using compact notation for large numbers

    Args:
        value: Subset count value
        use_compact: If True, use K/M notation for large numbers

    Returns:
        Formatted label string

    Examples:
        >>> format_venn_subset_label(50)
        '50'
        >>> format_venn_subset_label(1500)
        '1,500'
        >>> format_venn_subset_label(1500, use_compact=True)
        '1.5K'
        >>> format_venn_subset_label(2500000, use_compact=True)
        '2.5M'
    """
    if use_compact:
        # Use K/M notation for compact display
        if value >= 1_000_000:
            return f'{value / 1_000_000:.1f}M'
        elif value >= 1_000:
            return f'{value / 1_000:.1f}K'
        else:
            return str(value)
    else:
        # Use comma separators for readability
        return f'{value:,}'


class VennDiagramPlotMixin:
    """Mixin class for Venn diagram visualization"""

    def plot_glycan_venn_diagram(self, df: pd.DataFrame, figsize: tuple = None):
        """
        Create Venn diagram showing overlap between glycan modification types

        IMPORTANT: Uses data copy to protect original DataFrame

        Args:
            df: Annotated DataFrame with glycan annotations
            figsize: Figure size (width, height), defaults to VENN_FIGSIZE
        """
        logger.info("Creating Venn diagram for glycan type overlap...")

        # Use centralized figsize
        if figsize is None:
            figsize = VENN_FIGSIZE

        # ========================================
        # DATA PROTECTION: Use copy to avoid modifying original DataFrame
        # ========================================
        df_copy = df.copy()

        # Check if internal columns exist, recreate if needed
        if 'IsSialylated' not in df_copy.columns or 'IsFucosylated' not in df_copy.columns:
            logger.debug("Recreating internal annotation columns...")
            df_copy['IsSialylated'] = df_copy['Sialylation'] == 'Sialylated'
            df_copy['IsFucosylated'] = df_copy['Fucosylation'] == 'Fucosylated'

        # Create sets for Venn diagram
        sialylated_set = set(df_copy[df_copy['IsSialylated']].index)
        fucosylated_set = set(df_copy[df_copy['IsFucosylated']].index)

        # High Mannose is defined by PrimaryClassification
        high_mannose_set = set(df_copy[df_copy['PrimaryClassification'] == 'High Mannose'].index)

        logger.debug(f"Set sizes: Sialylated={len(sialylated_set)}, "
                    f"Fucosylated={len(fucosylated_set)}, "
                    f"High Mannose={len(high_mannose_set)}")

        # Create figure
        fig, ax = plt.subplots(figsize=figsize)

        # ========================================
        # VENN DIAGRAM with centralized constants
        # ========================================
        venn = venn3(
            subsets=[sialylated_set, fucosylated_set, high_mannose_set],
            set_labels=('Sialylated', 'Fucosylated', 'High Mannose'),
            ax=ax,
            set_colors=VENN_COLORS,
            alpha=VENN_ALPHA
        )

        # Add circles for better visualization
        _ = venn3_circles(
            subsets=[sialylated_set, fucosylated_set, high_mannose_set],
            ax=ax,
            linewidth=VENN_LINEWIDTH,
            linestyle=VENN_LINESTYLE,
            color='gray'
        )

        # ========================================
        # ENHANCED LABELING with font size control
        # ========================================
        # Customize set labels (Sialylated, Fucosylated, High Mannose)
        for text in venn.set_labels:
            if text:
                text.set_fontsize(VENN_LABEL_FONTSIZE)
                text.set_fontweight('bold')

        # ========================================
        # ENHANCED: Adaptive subset label formatting
        # ========================================
        # Calculate max value for adaptive sizing
        subset_values = []
        for text in venn.subset_labels:
            if text and text.get_text():
                try:
                    val = int(text.get_text())
                    subset_values.append(val)
                except (ValueError, AttributeError):
                    pass

        max_value = max(subset_values) if subset_values else 0
        use_compact = max_value > 10000  # Use compact notation if any value > 10K

        # Customize subset labels with adaptive formatting
        for text in venn.subset_labels:
            if text and text.get_text():
                try:
                    current_value = int(text.get_text())

                    # Format label
                    formatted_label = format_venn_subset_label(current_value, use_compact=use_compact)
                    text.set_text(formatted_label)

                    # Adaptive font size: smaller for intersection regions (value < 10% of max)
                    if current_value < max_value * 0.1:
                        text.set_fontsize(VENN_SUBSET_FONTSIZE * 0.85)
                    else:
                        text.set_fontsize(VENN_SUBSET_FONTSIZE)

                    # Add subtle background box for readability
                    text.set_bbox(dict(boxstyle='round,pad=0.3',
                                      facecolor='white',
                                      alpha=0.7,
                                      edgecolor='none'))

                    logger.debug(f"  Subset label: {current_value} → '{formatted_label}'")

                except (ValueError, AttributeError) as e:
                    logger.warning(f"  Failed to parse subset label: {text.get_text()} ({e})")
                    pass

        # Title with centralized font size
        ax.set_title('Glycan Type Overlap Analysis\n(Venn Diagram)',
                     fontsize=VENN_TITLE_FONTSIZE, fontweight='bold', pad=20)

        # ========================================
        # MAIN EXTERNAL STATISTICS BOX (MetaboAnalyst style)
        # ========================================
        # Calculate statistics
        total_glycopeptides = len(df_copy)
        sial_only = len(sialylated_set - fucosylated_set - high_mannose_set)
        fuc_only = len(fucosylated_set - sialylated_set - high_mannose_set)
        hm_only = len(high_mannose_set - sialylated_set - fucosylated_set)
        sial_fuc = len(sialylated_set & fucosylated_set - high_mannose_set)
        all_three = len(sialylated_set & fucosylated_set & high_mannose_set)

        # Create statistics text
        stats_text = f"Total Glycopeptides: {total_glycopeptides}\n\n"
        stats_text += f"Sialylated only: {sial_only}\n"
        stats_text += f"Fucosylated only: {fuc_only}\n"
        stats_text += f"High Mannose only: {hm_only}\n"
        stats_text += f"Sialylated + Fucosylated: {sial_fuc}\n"
        stats_text += f"All three: {all_three}"

        # Position using normalized axes coordinates (0-1)
        # This ensures consistent positioning regardless of data scale
        ax.text(1.15, 0.5, stats_text,
                transform=ax.transAxes,  # Use axes coordinates, not data coordinates
                fontsize=VENN_STATS_FONTSIZE,
                verticalalignment='center',
                bbox=dict(boxstyle='round', facecolor='white',
                         alpha=0.8, edgecolor='gray', linewidth=1.5))

        logger.info(f"Venn diagram statistics: Total={total_glycopeptides}, "
                   f"Sial_only={sial_only}, Fuc_only={fuc_only}, HM_only={hm_only}, "
                   f"Sial+Fuc={sial_fuc}, All_three={all_three}")

        plt.tight_layout()

        # ========================================
        # SAVE with standardized function
        # ========================================
        output_file = self.output_dir / 'venn_diagram_glycan_types.png'
        save_publication_figure(fig, output_file, dpi=VENN_DPI)
        logger.info(f"Saved Venn diagram to {output_file} (optimized, {VENN_DPI} DPI)")

        # ========================================
        # TRACE DATA for reproducibility
        # ========================================
        venn_data = []
        for idx in df_copy.index:
            categories = []
            if idx in sialylated_set:
                categories.append('Sialylated')
            if idx in fucosylated_set:
                categories.append('Fucosylated')
            if idx in high_mannose_set:
                categories.append('High Mannose')

            venn_data.append({
                'Index': idx,
                'Peptide': df_copy.loc[idx, 'Peptide'],
                'GlycanComposition': df_copy.loc[idx, 'GlycanComposition'],
                'Categories': '|'.join(categories) if categories else 'None',
                'Is_Sialylated': idx in sialylated_set,
                'Is_Fucosylated': idx in fucosylated_set,
                'Is_HighMannose': idx in high_mannose_set
            })

        venn_df = pd.DataFrame(venn_data)
        save_trace_data(venn_df, self.output_dir, 'venn_diagram_glycan_types_data.csv')

        plt.close()

        logger.info("✓ Venn diagram creation complete")
        return venn_df
