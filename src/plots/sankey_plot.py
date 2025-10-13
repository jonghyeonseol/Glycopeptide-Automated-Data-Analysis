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

    def plot_group_to_glycan_sankey(
        self,
        df: pd.DataFrame,
        vip_scores: pd.DataFrame,
        config: DataPreparationConfig = None,
        log2fc_threshold: float = 1.0,
        fdr_threshold: float = 0.05
    ):
        """
        Create Sankey diagram showing flows from Cancer/Normal groups to Glycan types

        Flow structure:
        Group (Cancer/Normal) → Glycan Type (HM, F, S, SF, C/H)

        Features:
        - Vertical deployment of both groups on left side
        - Glycan types on right side
        - Links colored by regulation status (upregulated/downregulated/unchanged)
        - Legend showing regulation categories

        Args:
            df: Annotated DataFrame with all samples
            vip_scores: VIP scores DataFrame from PLS-DA
            config: Data preparation configuration
            log2fc_threshold: Log2 fold change threshold for regulation (default: 1.0 = 2-fold)
            fdr_threshold: FDR threshold for significance (default: 0.05)
        """
        logger.info("Creating Group → Glycan Type Sankey diagram...")

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
            log_prefix="[Group Sankey] "
        )

        if len(df_prepared) == 0:
            logger.error("No glycopeptides available!")
            return

        # Get sample lists
        cancer_samples = [col for col in df_prepared.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df_prepared.columns if col.startswith('N') and col[1:].isdigit()]

        # Calculate statistical significance
        logger.info("Calculating statistical significance and regulation status...")
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

        # Determine group presence for each glycopeptide
        df_with_stats['Has_Cancer'] = df_with_stats['Cancer_Mean'].notna() & (df_with_stats['Cancer_Mean'] > 0)
        df_with_stats['Has_Normal'] = df_with_stats['Normal_Mean'].notna() & (df_with_stats['Normal_Mean'] > 0)

        # Define nodes
        groups = ['Cancer', 'Normal']
        glycan_types = ['HM', 'F', 'S', 'SF', 'C/H']

        # Create node labels and indices
        node_labels = groups + glycan_types
        node_dict = {label: idx for idx, label in enumerate(node_labels)}

        # Define node colors
        node_colors = []
        # Group colors
        node_colors.append('#E74C3C')  # Cancer - Red
        node_colors.append('#3498DB')  # Normal - Blue
        # Glycan type colors (from plot_config)
        for gt in glycan_types:
            node_colors.append(EXTENDED_CATEGORY_COLORS.get(gt, '#CCCCCC'))

        # Build flows from Cancer to Glycan Types
        cancer_flows = []
        for gt in glycan_types:
            gt_data = df_with_stats[(df_with_stats['GlycanTypeCategory'] == gt) & df_with_stats['Has_Cancer']]
            if len(gt_data) > 0:
                # Count by regulation status
                reg_counts = gt_data['Regulation'].value_counts()
                for reg_status, count in reg_counts.items():
                    cancer_flows.append({
                        'source': 'Cancer',
                        'target': gt,
                        'regulation': reg_status,
                        'count': count
                    })

        # Build flows from Normal to Glycan Types
        normal_flows = []
        for gt in glycan_types:
            gt_data = df_with_stats[(df_with_stats['GlycanTypeCategory'] == gt) & df_with_stats['Has_Normal']]
            if len(gt_data) > 0:
                # Count by regulation status
                reg_counts = gt_data['Regulation'].value_counts()
                for reg_status, count in reg_counts.items():
                    normal_flows.append({
                        'source': 'Normal',
                        'target': gt,
                        'regulation': reg_status,
                        'count': count
                    })

        # Combine all flows
        all_flows = cancer_flows + normal_flows

        # Aggregate flows by source-target pairs (sum across regulation statuses)
        # But keep track of dominant regulation for coloring
        flow_aggregated = {}
        flow_regulation = {}

        for flow in all_flows:
            key = (flow['source'], flow['target'])
            if key not in flow_aggregated:
                flow_aggregated[key] = 0
                flow_regulation[key] = {}

            flow_aggregated[key] += flow['count']
            reg = flow['regulation']
            if reg not in flow_regulation[key]:
                flow_regulation[key][reg] = 0
            flow_regulation[key][reg] += flow['count']

        # Build Sankey data structures
        source_indices = []
        target_indices = []
        values = []
        link_colors = []
        link_labels = []

        # Define regulation colors
        reg_colors = {
            'Upregulated': 'rgba(231, 76, 60, 0.6)',     # Red with transparency
            'Downregulated': 'rgba(52, 152, 219, 0.6)',  # Blue with transparency
            'Unchanged': 'rgba(149, 165, 166, 0.4)'      # Gray with transparency
        }

        for (source, target), total_count in flow_aggregated.items():
            source_idx = node_dict[source]
            target_idx = node_dict[target]

            source_indices.append(source_idx)
            target_indices.append(target_idx)
            values.append(total_count)

            # Determine dominant regulation status
            reg_breakdown = flow_regulation[(source, target)]
            dominant_reg = max(reg_breakdown, key=reg_breakdown.get)

            # Color by dominant regulation
            link_colors.append(reg_colors[dominant_reg])

            # Create hover label with regulation breakdown
            reg_text = '<br>'.join([f"{reg}: {cnt}" for reg, cnt in reg_breakdown.items()])
            link_labels.append(f"{source} → {target}<br>Total: {total_count}<br>{reg_text}")

        # Create Sankey diagram
        fig = go.Figure(data=[go.Sankey(
            arrangement='snap',  # Snap to grid for vertical alignment
            node=dict(
                pad=20,
                thickness=25,
                line=dict(color="black", width=1.0),
                label=node_labels,
                color=node_colors,
                x=[0.01, 0.01, 0.99, 0.99, 0.99, 0.99, 0.99],  # Position: Groups left, Glycan types right
                y=[0.2, 0.8, 0.1, 0.3, 0.5, 0.7, 0.9],  # Vertical spacing
                hovertemplate='%{label}<br>Total: %{value}<extra></extra>'
            ),
            link=dict(
                source=source_indices,
                target=target_indices,
                value=values,
                color=link_colors,
                customdata=link_labels,
                hovertemplate='%{customdata}<extra></extra>'
            )
        )])

        # Update layout with legend
        fig.update_layout(
            title=dict(
                text="Glycan Type Distribution: Cancer vs Normal Groups<br>"
                     f"<sub>Regulation threshold: |Log2FC| ≥ {log2fc_threshold} (2-fold), "
                     f"Significance: FDR < {fdr_threshold}</sub>",
                x=0.5,
                xanchor='center',
                font=dict(size=20, family='Arial, sans-serif', color='#2C3E50')
            ),
            font=dict(size=13, family='Arial, sans-serif'),
            height=1100,  # Increased from 1000 to accommodate legend below
            width=1200,
            plot_bgcolor='white',
            paper_bgcolor='#F8F9FA',
            margin=dict(l=80, r=80, t=150, b=200),  # Increased bottom margin to 200 to create space
            annotations=[
                # Add group labels at top
                dict(
                    x=0.01, y=1.03,
                    xref='paper', yref='paper',
                    text='<b>Groups</b>',
                    showarrow=False,
                    font=dict(size=14, family='Arial, sans-serif', color='#2C3E50'),
                    xanchor='center'
                ),
                dict(
                    x=0.99, y=1.03,
                    xref='paper', yref='paper',
                    text='<b>Glycan Types</b>',
                    showarrow=False,
                    font=dict(size=14, family='Arial, sans-serif', color='#2C3E50'),
                    xanchor='center'
                ),
                # Add legend box underneath main diagram
                # y=-0.12 places it well below the Sankey flows (which end at y=0)
                dict(
                    x=0.5, y=-0.12,
                    xref='paper', yref='paper',
                    text='<b>Link Colors (Regulation Status):</b><br>'
                         '<span style="color:#E74C3C">━━━</span> Upregulated (|Log2FC| ≥ 1.0) | '
                         '<span style="color:#3498DB">━━━</span> Downregulated (|Log2FC| ≤ -1.0) | '
                         '<span style="color:#95A5A6">━━━</span> Unchanged',
                    showarrow=False,
                    font=dict(size=13, family='Arial, sans-serif'),
                    align='center',
                    xanchor='center',
                    yanchor='top',
                    bgcolor='white',
                    bordercolor='#BDC3C7',
                    borderwidth=1.5,
                    borderpad=10
                )
            ]
        )

        # Save as PNG and HTML
        output_png = self.output_dir / 'sankey_group_to_glycan_type.png'
        output_html = self.output_dir / 'sankey_group_to_glycan_type.html'

        try:
            fig.write_image(str(output_png), width=1200, height=1100, scale=2)
            logger.info(f"Saved Group → Glycan Type Sankey (PNG) to {output_png}")
        except Exception as e:
            logger.warning(f"Could not save PNG (requires kaleido): {e}")

        fig.write_html(str(output_html))
        logger.info(f"Saved Group → Glycan Type Sankey (HTML) to {output_html}")

        # Save detailed trace data with regulation breakdown
        trace_data = df_with_stats[['Peptide', 'GlycanComposition', 'GlycanTypeCategory',
                                      'Has_Cancer', 'Has_Normal',
                                      'Cancer_Mean', 'Normal_Mean',
                                      'Log2_Fold_Change', 'FDR', 'P_Value',
                                      'Regulation', 'Significance',
                                      'VIP_Score']].copy()

        save_trace_data(trace_data, self.output_dir, 'sankey_group_to_glycan_data.csv')

        # Log summary statistics
        logger.info("Group → Glycan Type Sankey summary:")
        logger.info(f"  Total glycopeptides: {len(df_with_stats)}")
        logger.info(f"  In Cancer: {df_with_stats['Has_Cancer'].sum()}")
        logger.info(f"  In Normal: {df_with_stats['Has_Normal'].sum()}")
        logger.info(f"  Total flows: {len(source_indices)}")

        # Log regulation breakdown
        reg_summary = df_with_stats['Regulation'].value_counts()
        logger.info("  Regulation breakdown:")
        for reg, count in reg_summary.items():
            logger.info(f"    {reg}: {count} ({count/len(df_with_stats)*100:.1f}%)")

        # Log glycan type breakdown
        logger.info("  Glycan type distribution:")
        for gt in glycan_types:
            cancer_count = len(df_with_stats[(df_with_stats['GlycanTypeCategory'] == gt) & df_with_stats['Has_Cancer']])
            normal_count = len(df_with_stats[(df_with_stats['GlycanTypeCategory'] == gt) & df_with_stats['Has_Normal']])
            logger.info(f"    {gt}: Cancer={cancer_count}, Normal={normal_count}")

        return fig
