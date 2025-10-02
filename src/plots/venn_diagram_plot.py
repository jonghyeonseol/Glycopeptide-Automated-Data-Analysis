"""
Venn Diagram Plot Module for pGlyco Auto Combine
Visualizes overlap between glycan type categories
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from pathlib import Path
import logging
from utils import save_trace_data

logger = logging.getLogger(__name__)


class VennDiagramPlotMixin:
    """Mixin class for Venn diagram visualization"""

    def plot_glycan_venn_diagram(self, df: pd.DataFrame, figsize: tuple = (12, 10)):
        """
        Create Venn diagram showing overlap between glycan modification types

        Args:
            df: Annotated DataFrame with glycan annotations
            figsize: Figure size (width, height)
        """
        # Check if internal columns exist
        if 'IsSialylated' not in df.columns or 'IsFucosylated' not in df.columns:
            logger.warning("Internal annotation columns not found, recreating...")
            # Recreate internal columns
            df['IsSialylated'] = df['Sialylation'] == 'Sialylated'
            df['IsFucosylated'] = df['Fucosylation'] == 'Fucosylated'

        # Create sets for Venn diagram
        sialylated_set = set(df[df['IsSialylated']].index)
        fucosylated_set = set(df[df['IsFucosylated']].index)

        # High Mannose is defined by PrimaryClassification
        high_mannose_set = set(df[df['PrimaryClassification'] == 'High Mannose'].index)

        # Create figure
        fig, ax = plt.subplots(figsize=figsize)

        # Create Venn diagram
        venn = venn3(
            subsets=[sialylated_set, fucosylated_set, high_mannose_set],
            set_labels=('Sialylated', 'Fucosylated', 'High Mannose'),
            ax=ax,
            set_colors=('#E74C3C', '#3498DB', '#8E44AD'),
            alpha=0.6
        )

        # Add circles for better visualization
        venn_circles = venn3_circles(
            subsets=[sialylated_set, fucosylated_set, high_mannose_set],
            ax=ax, linewidth=2, linestyle='--', color='gray'
        )

        # Customize labels
        for text in venn.set_labels:
            if text:
                text.set_fontsize(13)
                text.set_fontweight('bold')

        for text in venn.subset_labels:
            if text:
                text.set_fontsize(11)

        # Title
        ax.set_title('Glycan Type Overlap Analysis\n(Venn Diagram)',
                    fontsize=16, fontweight='bold', pad=20)

        # Add statistics text
        total_glycopeptides = len(df)
        sial_only = len(sialylated_set - fucosylated_set - high_mannose_set)
        fuc_only = len(fucosylated_set - sialylated_set - high_mannose_set)
        hm_only = len(high_mannose_set - sialylated_set - fucosylated_set)
        sial_fuc = len(sialylated_set & fucosylated_set - high_mannose_set)
        all_three = len(sialylated_set & fucosylated_set & high_mannose_set)

        stats_text = f"Total Glycopeptides: {total_glycopeptides}\n\n"
        stats_text += f"Sialylated only: {sial_only}\n"
        stats_text += f"Fucosylated only: {fuc_only}\n"
        stats_text += f"High Mannose only: {hm_only}\n"
        stats_text += f"Sialylated + Fucosylated: {sial_fuc}\n"
        stats_text += f"All three: {all_three}"

        ax.text(1.3, 0.5, stats_text, transform=ax.transData,
               fontsize=10, verticalalignment='center',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'venn_diagram_glycan_types.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved Venn diagram to {output_file}")

        # Prepare trace data
        venn_data = []
        for idx in df.index:
            categories = []
            if idx in sialylated_set:
                categories.append('Sialylated')
            if idx in fucosylated_set:
                categories.append('Fucosylated')
            if idx in high_mannose_set:
                categories.append('High Mannose')

            venn_data.append({
                'Index': idx,
                'Peptide': df.loc[idx, 'Peptide'],
                'GlycanComposition': df.loc[idx, 'GlycanComposition'],
                'Categories': '|'.join(categories) if categories else 'None',
                'Is_Sialylated': idx in sialylated_set,
                'Is_Fucosylated': idx in fucosylated_set,
                'Is_HighMannose': idx in high_mannose_set
            })

        venn_df = pd.DataFrame(venn_data)
        save_trace_data(venn_df, self.output_dir, 'venn_diagram_glycan_types_data.csv')

        plt.close()

        return venn_df
