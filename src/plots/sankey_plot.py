"""
Sankey Diagram Plot Module for pGlyco Auto Combine
Visualizes glycan type flows through regulation and significance status

Author: pGlyco Auto Combine Pipeline
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import logging
from pathlib import Path
from ..utils import save_trace_data, calculate_fold_change
from ..data_preparation import (
    DataPreparationConfig,
    prepare_visualization_data,
    calculate_statistical_significance
)
from .plot_config import EXTENDED_CATEGORY_COLORS

logger = logging.getLogger(__name__)


class SankeyPlotMixin:
    """Mixin class for Sankey diagram visualizations"""

    def plot_glycan_type_sankey(
        self,
        df: pd.DataFrame,
        vip_scores: pd.DataFrame,
        config: DataPreparationConfig = None,
        log2fc_threshold: float = 1.0,
        fdr_threshold: float = 0.05
    ):
        """
        Create Sankey diagram showing glycan type flows through regulation and significance

        Flow structure:
        Glycan Type → Regulation Status → Significance Status

        Nodes:
        - Source: 5 Glycan Types (HM, F, S, SF, C/H)
        - Middle: 3 Regulation States (Upregulated, Downregulated, Unchanged)
        - Target: 2 Significance States (Significant, Non-significant)

        Args:
            df: Annotated DataFrame with all samples
            vip_scores: VIP scores DataFrame from PLS-DA
            config: Data preparation configuration
            log2fc_threshold: Log2 fold change threshold for regulation
            fdr_threshold: FDR threshold for significance
        """
        logger.info("Creating Sankey diagram for glycan type flows...")

        # Use default config if not provided
        if config is None:
            config = DataPreparationConfig(
                min_detection_pct=0.30,
                min_samples=5,
                missing_data_method='skipna'
            )

        # Prepare data with statistics
        df_prepared = prepare_visualization_data(
            df=df,
            config=config,
            vip_scores=vip_scores,
            merge_method='left',
            apply_detection_filter=False,
            log_prefix="[Sankey] "
        )

        if len(df_prepared) == 0:
            logger.error("No glycopeptides available!")
            return

        # Get sample lists
        cancer_samples = [col for col in df_prepared.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df_prepared.columns if col.startswith('N') and col[1:].isdigit()]

        # Calculate statistical significance
        logger.info("Calculating statistical significance for Sankey flows...")
        df_with_stats = calculate_statistical_significance(
            df_prep=df_prepared,
            cancer_samples=cancer_samples,
            normal_samples=normal_samples,
            method='mannwhitneyu',
            fdr_correction=True
        )

        # Classify by regulation status
        df_with_stats['Regulation'] = 'Unchanged'
        df_with_stats.loc[df_with_stats['Log2_Fold_Change'] >= log2fc_threshold, 'Regulation'] = 'Upregulated'
        df_with_stats.loc[df_with_stats['Log2_Fold_Change'] <= -log2fc_threshold, 'Regulation'] = 'Downregulated'

        # Classify by significance
        df_with_stats['Significance'] = 'Non-significant'
        df_with_stats.loc[df_with_stats['FDR'] < fdr_threshold, 'Significance'] = 'Significant'

        # Count flows
        flow_counts = df_with_stats.groupby(['GlycanTypeCategory', 'Regulation', 'Significance']).size().reset_index(name='Count')

        logger.info(f"Total flows: {len(flow_counts)}")
        logger.info(f"Total glycopeptides: {len(df_with_stats)}")

        # Define nodes
        glycan_types = ['HM', 'F', 'S', 'SF', 'C/H']
        regulation_states = ['Upregulated', 'Downregulated', 'Unchanged']
        significance_states = ['Significant', 'Non-significant']

        # Create node labels and indices
        node_labels = glycan_types + regulation_states + significance_states
        node_dict = {label: idx for idx, label in enumerate(node_labels)}

        # Define node colors
        node_colors = []
        # Glycan type colors
        for gt in glycan_types:
            node_colors.append(EXTENDED_CATEGORY_COLORS.get(gt, '#CCCCCC'))
        # Regulation colors
        node_colors.append('#E74C3C')  # Upregulated - Red
        node_colors.append('#3498DB')  # Downregulated - Blue
        node_colors.append('#95A5A6')  # Unchanged - Gray
        # Significance colors
        node_colors.append('#27AE60')  # Significant - Green
        node_colors.append('#ECF0F1')  # Non-significant - Light Gray

        # Build flows
        source_indices = []
        target_indices = []
        values = []
        link_colors = []

        # Flow 1: Glycan Type → Regulation
        for _, row in flow_counts.iterrows():
            glycan_type = row['GlycanTypeCategory']
            regulation = row['Regulation']
            significance = row['Significance']
            count = row['Count']

            if glycan_type in node_dict and regulation in node_dict:
                # Glycan Type → Regulation
                source_indices.append(node_dict[glycan_type])
                target_indices.append(node_dict[regulation])
                values.append(count)

                # Color link by glycan type with transparency
                base_color = EXTENDED_CATEGORY_COLORS.get(glycan_type, '#CCCCCC')
                # Convert hex to rgba with alpha=0.3
                rgba = f"rgba({int(base_color[1:3], 16)}, {int(base_color[3:5], 16)}, {int(base_color[5:7], 16)}, 0.3)"
                link_colors.append(rgba)

        # Flow 2: Regulation → Significance
        regulation_to_sig = df_with_stats.groupby(['Regulation', 'Significance']).size().reset_index(name='Count')
        for _, row in regulation_to_sig.iterrows():
            regulation = row['Regulation']
            significance = row['Significance']
            count = row['Count']

            if regulation in node_dict and significance in node_dict:
                source_indices.append(node_dict[regulation])
                target_indices.append(node_dict[significance])
                values.append(count)

                # Color link by regulation status with transparency
                if regulation == 'Upregulated':
                    link_colors.append('rgba(231, 76, 60, 0.3)')  # Red
                elif regulation == 'Downregulated':
                    link_colors.append('rgba(52, 152, 219, 0.3)')  # Blue
                else:
                    link_colors.append('rgba(149, 165, 166, 0.3)')  # Gray

        # Create Sankey diagram
        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=node_labels,
                color=node_colors,
                customdata=[f"{label}<br>{values[i] if i < len(values) else 0}" for i, label in enumerate(node_labels)],
                hovertemplate='%{label}<br>Count: %{value}<extra></extra>'
            ),
            link=dict(
                source=source_indices,
                target=target_indices,
                value=values,
                color=link_colors,
                hovertemplate='%{source.label} → %{target.label}<br>Count: %{value}<extra></extra>'
            )
        )])

        # Update layout
        fig.update_layout(
            title=dict(
                text=f"Glycan Type Flow Analysis<br><sub>Regulation (|Log2FC| ≥ {log2fc_threshold}) and Significance (FDR < {fdr_threshold})</sub>",
                x=0.5,
                xanchor='center',
                font=dict(size=18, family='Arial, sans-serif')
            ),
            font=dict(size=12, family='Arial, sans-serif'),
            height=800,
            width=1400,
            plot_bgcolor='white',
            paper_bgcolor='white',
            margin=dict(l=50, r=50, t=100, b=50)
        )

        # Save as PNG and HTML
        output_png = self.output_dir / 'sankey_glycan_type_flow.png'
        output_html = self.output_dir / 'sankey_glycan_type_flow.html'

        try:
            fig.write_image(str(output_png), width=1400, height=800, scale=2)
            logger.info(f"Saved Sankey diagram (PNG) to {output_png}")
        except Exception as e:
            logger.warning(f"Could not save PNG (requires kaleido): {e}")

        fig.write_html(str(output_html))
        logger.info(f"Saved Sankey diagram (HTML) to {output_html}")

        # Save trace data
        trace_data = df_with_stats[['Peptide', 'GlycanComposition', 'GlycanTypeCategory',
                                      'Log2_Fold_Change', 'FDR', 'P_Value',
                                      'Regulation', 'Significance',
                                      'Cancer_Mean', 'Normal_Mean', 'VIP_Score']].copy()

        save_trace_data(trace_data, self.output_dir, 'sankey_glycan_flow_data.csv')

        # Log summary statistics
        logger.info("Sankey diagram summary:")
        logger.info(f"  Glycan types: {len(glycan_types)}")
        logger.info(f"  Regulation states: {len(regulation_states)}")
        logger.info(f"  Significance states: {len(significance_states)}")
        logger.info(f"  Total flows: {len(source_indices)}")
        logger.info(f"  Total glycopeptides: {len(df_with_stats)}")

        # Log breakdown by category
        for gt in glycan_types:
            count = len(df_with_stats[df_with_stats['GlycanTypeCategory'] == gt])
            if count > 0:
                logger.info(f"    {gt}: {count} glycopeptides")

        return fig
