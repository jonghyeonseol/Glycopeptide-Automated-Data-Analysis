"""
Distribution Plot Module for pGlyco Auto Combine
Handles distribution visualizations
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from ..utils import save_trace_data

logger = logging.getLogger(__name__)


class DistributionPlotMixin:
    """Mixin class for distribution-related plots"""

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

        # Save plot
        output_file = self.output_dir / 'glycan_type_distribution.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved glycan type distribution to {output_file}")

        # Save trace data
        trace_data = pd.DataFrame({'GlycanType': type_counts.index, 'Count': type_counts.values})
        save_trace_data(trace_data, self.output_dir, 'glycan_type_distribution_data.csv')

        plt.close()
