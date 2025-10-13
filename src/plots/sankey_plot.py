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
from ..utils import save_trace_data, calculate_fold_change, get_sample_columns
from ..data_preparation import (
    DataPreparationConfig,
    prepare_visualization_data,
    calculate_statistical_significance
)
from .plot_config import EXTENDED_CATEGORY_COLORS, COLOR_CANCER, COLOR_NORMAL

logger = logging.getLogger(__name__)

# ============================================================================
# Module-level Constants
# ============================================================================

# Regulation state names
REGULATION_UPREGULATED = 'Upregulated'
REGULATION_DOWNREGULATED = 'Downregulated'
REGULATION_UNCHANGED = 'Unchanged'

# Regulation node suffixes (for Group → Glycan Type diagram)
REGULATION_SUFFIX_UP = 'Up'
REGULATION_SUFFIX_DOWN = 'Down'

# Significance state names
SIGNIFICANCE_SIGNIFICANT = 'Significant'
SIGNIFICANCE_NON_SIGNIFICANT = 'Non-significant'

# Glycan type categories
GLYCAN_TYPES = ['HM', 'F', 'S', 'SF', 'C/H']

# Group names
GROUP_CANCER = 'Cancer'
GROUP_NORMAL = 'Normal'

# Color constants - Regulation states
REGULATION_COLOR_UP = COLOR_CANCER  # Red
REGULATION_COLOR_DOWN = COLOR_NORMAL  # Blue
REGULATION_COLOR_UNCHANGED = '#95A5A6'  # Gray

# Color constants - Groups
GROUP_COLOR_CANCER = COLOR_CANCER  # Red
GROUP_COLOR_NORMAL = COLOR_NORMAL  # Blue

# Color constants - Significance
SIGNIFICANCE_COLOR_SIGNIFICANT = '#27AE60'  # Green
SIGNIFICANCE_COLOR_NON_SIGNIFICANT = '#95A5A6'  # Darker Gray (improved contrast)

# Detection threshold for group presence
MIN_INTENSITY_THRESHOLD = 0  # Minimum intensity to consider glycopeptide as "present"

# Default transparency for links
LINK_ALPHA = 0.3  # 30% transparency
LINK_ALPHA_GROUP = 0.4  # 40% transparency for group-to-glycan links


class SankeyPlotMixin:
    """Mixin class for Sankey diagram visualizations"""

    # ========================================================================
    # Helper Methods (Private)
    # ========================================================================

    @staticmethod
    def _get_sample_lists(df: pd.DataFrame) -> tuple:
        """
        Extract cancer and normal sample column names from DataFrame.

        Uses centralized get_sample_columns() utility for consistent sample extraction.

        Args:
            df: DataFrame with sample columns

        Returns:
            Tuple of (cancer_samples, normal_samples) as lists
        """
        return get_sample_columns(df)

    @staticmethod
    def _validate_samples(cancer_samples: list, normal_samples: list) -> bool:
        """
        Validate that both cancer and normal samples are available.

        Args:
            cancer_samples: List of cancer sample column names
            normal_samples: List of normal sample column names

        Returns:
            True if both groups have samples, False otherwise
        """
        if len(cancer_samples) == 0 or len(normal_samples) == 0:
            logger.error(f"Insufficient samples: Cancer={len(cancer_samples)}, Normal={len(normal_samples)}")
            logger.error("Cannot perform statistical comparison without samples from both groups!")
            return False

        logger.info(f"Sample validation: Cancer={len(cancer_samples)}, Normal={len(normal_samples)}")
        return True

    @staticmethod
    def _hex_to_rgba(hex_color: str, alpha: float = LINK_ALPHA) -> str:
        """
        Convert hex color to RGBA string with transparency.

        Args:
            hex_color: Hex color string (e.g., "#FF0000")
            alpha: Transparency value (0.0 to 1.0)

        Returns:
            RGBA color string (e.g., "rgba(255, 0, 0, 0.3)")
        """
        r = int(hex_color[1:3], 16)
        g = int(hex_color[3:5], 16)
        b = int(hex_color[5:7], 16)
        return f"rgba({r}, {g}, {b}, {alpha})"

    @staticmethod
    def _calculate_regulation_status(
        df_with_stats: pd.DataFrame,
        log2fc_threshold: float = 1.0,
        fdr_threshold: float = 0.05
    ) -> pd.DataFrame:
        """
        Classify glycopeptides by regulation status (Cancer-centric).

        IMPORTANT: Requires BOTH fold change AND statistical significance.
        - Upregulated: Log2FC >= threshold AND FDR < threshold (Cancer > Normal)
        - Downregulated: Log2FC <= -threshold AND FDR < threshold (Cancer < Normal)
        - Unchanged: Either |Log2FC| < threshold OR FDR >= threshold OR FDR is NaN

        FDR=NaN Handling:
        - If FDR is NaN (missing), the glycopeptide is automatically classified as:
          * Regulation: "Unchanged" (cannot be Up/Down without valid FDR)
          * Significance: "Non-significant" (no statistical test result)
        - This ensures scientific validity: regulation requires statistical evidence

        Args:
            df_with_stats: DataFrame with Log2_Fold_Change and FDR columns
            log2fc_threshold: Fold change threshold (default: 1.0 = 2-fold)
            fdr_threshold: FDR significance threshold (default: 0.05)

        Returns:
            DataFrame with added 'Regulation' and 'Significance' columns
        """
        # Count FDR=NaN cases for logging
        fdr_nan_count = df_with_stats['FDR'].isna().sum()
        if fdr_nan_count > 0:
            logger.info(f"  FDR=NaN found: {fdr_nan_count} glycopeptides → treated as Non-significant/Unchanged")

        # Classify by regulation status (REQUIRES both fold change AND valid FDR < threshold)
        df_with_stats['Regulation'] = REGULATION_UNCHANGED

        # Upregulated: Log2FC >= threshold AND FDR < threshold AND FDR is not NaN
        upregulated_mask = (df_with_stats['Log2_Fold_Change'] >= log2fc_threshold) & \
                          (df_with_stats['FDR'] < fdr_threshold) & \
                          (df_with_stats['FDR'].notna())
        df_with_stats.loc[upregulated_mask, 'Regulation'] = REGULATION_UPREGULATED

        # Downregulated: Log2FC <= -threshold AND FDR < threshold AND FDR is not NaN
        downregulated_mask = (df_with_stats['Log2_Fold_Change'] <= -log2fc_threshold) & \
                            (df_with_stats['FDR'] < fdr_threshold) & \
                            (df_with_stats['FDR'].notna())
        df_with_stats.loc[downregulated_mask, 'Regulation'] = REGULATION_DOWNREGULATED

        # Classify by significance (FDR=NaN is treated as Non-significant)
        df_with_stats['Significance'] = SIGNIFICANCE_NON_SIGNIFICANT
        sig_mask = (df_with_stats['FDR'] < fdr_threshold) & (df_with_stats['FDR'].notna())
        df_with_stats.loc[sig_mask, 'Significance'] = SIGNIFICANCE_SIGNIFICANT

        return df_with_stats

    # ========================================================================
    # Main Visualization Methods (Public)
    # ========================================================================

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

        IMPORTANT - Regulation Definition (Cancer-centric):
        - "Upregulated" means Cancer > Normal (Log2FC ≥ +threshold AND FDR < threshold)
        - "Downregulated" means Cancer < Normal (Log2FC ≤ -threshold AND FDR < threshold)
        - "Unchanged" means either |Log2FC| < threshold OR FDR ≥ threshold

        Scientific Validity:
        - Requires BOTH fold change AND statistical significance (FDR < 0.05)
        - This follows standard proteomics practice (e.g., volcano plot methodology)
        - Prevents false positives from being classified as "regulated"

        Args:
            df: Annotated DataFrame with all samples
            vip_scores: VIP scores DataFrame from PLS-DA
            config: Data preparation configuration
            log2fc_threshold: Log2 fold change threshold for regulation (default: 1.0 = 2-fold)
            fdr_threshold: FDR threshold for significance (default: 0.05)
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

        # Get sample lists and validate
        cancer_samples, normal_samples = self._get_sample_lists(df_prepared)
        if not self._validate_samples(cancer_samples, normal_samples):
            return None

        # Calculate statistical significance
        logger.info("Calculating statistical significance for Sankey flows...")
        df_with_stats = calculate_statistical_significance(
            df_prep=df_prepared,
            cancer_samples=cancer_samples,
            normal_samples=normal_samples,
            method='mannwhitneyu',
            fdr_correction=True
        )

        # Classify by regulation status (uses helper method)
        df_with_stats = self._calculate_regulation_status(
            df_with_stats,
            log2fc_threshold=log2fc_threshold,
            fdr_threshold=fdr_threshold
        )

        logger.info(f"Total glycopeptides: {len(df_with_stats)}")

        # Define nodes (using constants)
        glycan_types = GLYCAN_TYPES
        regulation_states = [REGULATION_UPREGULATED, REGULATION_DOWNREGULATED, REGULATION_UNCHANGED]
        significance_states = [SIGNIFICANCE_SIGNIFICANT, SIGNIFICANCE_NON_SIGNIFICANT]

        # Create node labels and indices
        node_labels = glycan_types + regulation_states + significance_states
        node_dict = {label: idx for idx, label in enumerate(node_labels)}

        # Define node colors (using constants)
        node_colors = []
        # Glycan type colors
        for gt in glycan_types:
            node_colors.append(EXTENDED_CATEGORY_COLORS.get(gt, '#CCCCCC'))
        # Regulation colors
        node_colors.append(REGULATION_COLOR_UP)
        node_colors.append(REGULATION_COLOR_DOWN)
        node_colors.append(REGULATION_COLOR_UNCHANGED)
        # Significance colors
        node_colors.append(SIGNIFICANCE_COLOR_SIGNIFICANT)
        node_colors.append(SIGNIFICANCE_COLOR_NON_SIGNIFICANT)

        # Build flows
        source_indices = []
        target_indices = []
        values = []
        link_colors = []

        # Flow 1: Glycan Type → Regulation (2-way groupby to avoid double-counting)
        glycan_to_reg = df_with_stats.groupby(['GlycanTypeCategory', 'Regulation']).size().reset_index(name='Count')
        for _, row in glycan_to_reg.iterrows():
            glycan_type = row['GlycanTypeCategory']
            regulation = row['Regulation']
            count = row['Count']

            if glycan_type in node_dict and regulation in node_dict:
                # Glycan Type → Regulation
                source_indices.append(node_dict[glycan_type])
                target_indices.append(node_dict[regulation])
                values.append(count)

                # Color link by glycan type with transparency (using helper)
                base_color = EXTENDED_CATEGORY_COLORS.get(glycan_type, '#CCCCCC')
                link_colors.append(self._hex_to_rgba(base_color, LINK_ALPHA))

        logger.info(f"Flow 1 (Glycan → Regulation): {len(glycan_to_reg)} links")

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

                # Color link by regulation status with transparency (using constants and helper)
                if regulation == REGULATION_UPREGULATED:
                    link_colors.append(self._hex_to_rgba(REGULATION_COLOR_UP, LINK_ALPHA))
                elif regulation == REGULATION_DOWNREGULATED:
                    link_colors.append(self._hex_to_rgba(REGULATION_COLOR_DOWN, LINK_ALPHA))
                else:
                    link_colors.append(self._hex_to_rgba(REGULATION_COLOR_UNCHANGED, LINK_ALPHA))

        logger.info(f"Flow 2 (Regulation → Significance): {len(regulation_to_sig)} links")
        logger.info(f"Total flows: {len(source_indices)} links")

        # ====================================================================
        # Data Validation: Prevent Over-Counting
        # ====================================================================
        total_glycopeptides = len(df_with_stats)
        total_flow = sum(values)

        logger.info(f"Data validation:")
        logger.info(f"  Total glycopeptides: {total_glycopeptides}")
        logger.info(f"  Total flow (sum of link values): {total_flow}")
        logger.info(f"  Source node flow count: {len(values)}")

        # Check: Total glycopeptides should equal total flow from source nodes
        # (Each glycopeptide flows through exactly once: Glycan → Regulation → Significance)
        if total_glycopeptides != total_flow:
            logger.warning(f"⚠️ VALIDATION WARNING: Glycopeptides ({total_glycopeptides}) ≠ Total flow ({total_flow})")
            logger.warning(f"   This indicates possible over-counting or under-counting!")
        else:
            logger.info(f"✓ Validation passed: No over-counting detected")

        # ====================================================================
        # Boundary Value Check: Verify Classification at Thresholds
        # ====================================================================
        logger.info(f"Boundary value checks (thresholds: Log2FC={log2fc_threshold}, FDR={fdr_threshold}):")

        # Check glycopeptides at exact boundary values
        boundary_log2fc = df_with_stats[
            (df_with_stats['Log2_Fold_Change'].abs() == log2fc_threshold) &
            (df_with_stats['FDR'] < fdr_threshold)
        ]
        if len(boundary_log2fc) > 0:
            logger.info(f"  Boundary Log2FC (±{log2fc_threshold}): {len(boundary_log2fc)} glycopeptides")
            for _, row in boundary_log2fc.head(3).iterrows():
                logger.info(f"    {row['Peptide']}: Log2FC={row['Log2_Fold_Change']:.3f}, FDR={row['FDR']:.4f}, Reg={row['Regulation']}")

        boundary_fdr = df_with_stats[df_with_stats['FDR'] == fdr_threshold]
        if len(boundary_fdr) > 0:
            logger.info(f"  Boundary FDR ({fdr_threshold}): {len(boundary_fdr)} glycopeptides")
            for _, row in boundary_fdr.head(3).iterrows():
                logger.info(f"    {row['Peptide']}: Log2FC={row['Log2_Fold_Change']:.3f}, FDR={row['FDR']:.4f}, Reg={row['Regulation']}")

        # ====================================================================
        # Calculate Node Totals Directly for Accurate Hover Display
        # ====================================================================
        # Calculate total flow through each node by summing incoming/outgoing links
        node_totals = [0] * len(node_labels)

        for i in range(len(source_indices)):
            src = source_indices[i]
            tgt = target_indices[i]
            val = values[i]

            # Add to source node (outgoing flow)
            node_totals[src] += val
            # Note: For middle/target nodes, we only count incoming flow once
            # (to avoid double-counting at middle nodes which have both in and out)

        # For target nodes (Significance), add incoming flow
        for i in range(len(source_indices)):
            tgt = target_indices[i]
            val = values[i]
            # Only add if this is a target node (Significance nodes)
            if tgt >= len(glycan_types) + len(regulation_states):
                node_totals[tgt] = node_totals[tgt] if node_totals[tgt] > 0 else val

        # Better approach: Calculate node totals from the dataframe directly
        node_totals_corrected = []
        # Glycan types: count from df_with_stats
        for gt in glycan_types:
            count = len(df_with_stats[df_with_stats['GlycanTypeCategory'] == gt])
            node_totals_corrected.append(count)
        # Regulation states: count from df_with_stats
        for reg in regulation_states:
            count = len(df_with_stats[df_with_stats['Regulation'] == reg])
            node_totals_corrected.append(count)
        # Significance states: count from df_with_stats
        for sig in significance_states:
            count = len(df_with_stats[df_with_stats['Significance'] == sig])
            node_totals_corrected.append(count)

        logger.info(f"Node totals calculated: {len(node_totals_corrected)} nodes")

        # Create Sankey diagram with corrected node hover
        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=node_labels,
                color=node_colors,
                customdata=node_totals_corrected,  # Pass calculated totals
                hovertemplate='%{label}<br>Total: %{customdata}<extra></extra>'
            ),
            link=dict(
                source=source_indices,
                target=target_indices,
                value=values,
                color=link_colors,
                hovertemplate='%{source.label} → %{target.label}<br>Count: %{value}<extra></extra>'
            )
        )])

        # Update layout with improved styling and clear definitions
        fig.update_layout(
            title=dict(
                text=(
                    f"<b>Glycan Type Flow Analysis: Regulation & Statistical Significance</b><br>"
                    f"<sub><b>Regulation Definition (Cancer-centric):</b> Upregulated = Cancer > Normal, "
                    f"Downregulated = Cancer < Normal<br>"
                    f"<b>Criteria:</b> |Log2FC| ≥ {log2fc_threshold} (2-fold change) AND FDR < {fdr_threshold} (statistically significant)<br>"
                    f"<b>Note:</b> Unchanged = |Log2FC| < {log2fc_threshold} OR FDR ≥ {fdr_threshold} OR FDR=NaN</sub>"
                ),
                x=0.5,
                xanchor='center',
                font=dict(size=18, family='Arial, sans-serif', color='#2C3E50')
            ),
            font=dict(size=12, family='Arial, sans-serif', color='#2C3E50'),
            height=800,
            width=1400,
            plot_bgcolor='white',
            paper_bgcolor='#FAFAFA',  # Light gray background for better contrast
            margin=dict(l=50, r=50, t=150, b=50)  # Increased top margin for longer title
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
        logger.info("Sankey diagram generation complete:")
        logger.info(f"  Nodes: {len(node_labels)} ({len(glycan_types)} glycan types + {len(regulation_states)} regulations + {len(significance_states)} significance)")

        # Log breakdown by category
        logger.info("  Glycan type distribution:")
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
        split by regulation status

        Flow structure:
        Group (Cancer/Normal) → Glycan Type + Regulation (10 nodes: HM_Up, HM_Down, F_Up, F_Down, etc.)

        Features:
        - Vertical deployment of both groups on left side
        - 10 glycan type nodes on right side (5 types × 2 regulation states: Up/Down)
        - Links colored by glycan type (HM=green, F=red, S=pink, SF=orange, C/H=blue)
        - Legend showing glycan type color mapping

        IMPORTANT - Regulation Interpretation (Cancer-centric):
        - "Upregulated" means Cancer > Normal (Log2FC ≥ +threshold AND FDR < threshold)
        - "Downregulated" means Cancer < Normal (Log2FC ≤ -threshold AND FDR < threshold)
        - Example: "Normal → HM_Up" means "HM glycopeptide detected in Normal samples,
          but classified as upregulated in Cancer vs Normal comparison"

        Scientific Validity:
        - Requires BOTH fold change AND statistical significance (FDR < 0.05)
        - Only glycopeptides passing both thresholds are classified as regulated
        - This follows standard proteomics practice (similar to volcano plot)

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

        # Get sample lists and validate
        cancer_samples, normal_samples = self._get_sample_lists(df_prepared)
        if not self._validate_samples(cancer_samples, normal_samples):
            return None

        # Calculate statistical significance
        logger.info("Calculating statistical significance and regulation status...")
        df_with_stats = calculate_statistical_significance(
            df_prep=df_prepared,
            cancer_samples=cancer_samples,
            normal_samples=normal_samples,
            method='mannwhitneyu',
            fdr_correction=True
        )

        # Classify by regulation status (uses helper method)
        df_with_stats = self._calculate_regulation_status(
            df_with_stats,
            log2fc_threshold=log2fc_threshold,
            fdr_threshold=fdr_threshold
        )

        # Determine group presence for each glycopeptide (using constant)
        df_with_stats['Has_Cancer'] = df_with_stats['Cancer_Mean'].notna() & (df_with_stats['Cancer_Mean'] > MIN_INTENSITY_THRESHOLD)
        df_with_stats['Has_Normal'] = df_with_stats['Normal_Mean'].notna() & (df_with_stats['Normal_Mean'] > MIN_INTENSITY_THRESHOLD)

        # Define nodes - 10 glycan type nodes (5 types × 2 regulation states) (using constants)
        groups = [GROUP_CANCER, GROUP_NORMAL]
        glycan_types = GLYCAN_TYPES

        # Create glycan type nodes with regulation suffix (using constants)
        glycan_reg_nodes = []
        for gt in glycan_types:
            glycan_reg_nodes.append(f"{gt}_{REGULATION_SUFFIX_UP}")
            glycan_reg_nodes.append(f"{gt}_{REGULATION_SUFFIX_DOWN}")

        # Create node labels and indices
        node_labels = groups + glycan_reg_nodes
        node_dict = {label: idx for idx, label in enumerate(node_labels)}

        # Define node colors (using constants)
        node_colors = []
        # Group colors
        node_colors.append(GROUP_COLOR_CANCER)
        node_colors.append(GROUP_COLOR_NORMAL)

        # Glycan type colors (from plot_config) - pair for Up/Down
        for gt in glycan_types:
            base_color = EXTENDED_CATEGORY_COLORS.get(gt, '#CCCCCC')
            node_colors.append(base_color)  # Up - same color
            node_colors.append(base_color)  # Down - same color

        # Build flows from Cancer/Normal to Glycan Types (split by regulation)
        source_indices = []
        target_indices = []
        values = []
        link_colors = []
        link_labels = []

        # Process Cancer flows (using constants and helper)
        for gt in glycan_types:
            # Get base glycan type color and convert to RGBA (using helper)
            base_color = EXTENDED_CATEGORY_COLORS.get(gt, '#CCCCCC')
            rgba = self._hex_to_rgba(base_color, LINK_ALPHA_GROUP)

            # Upregulated flows
            gt_data_up = df_with_stats[
                (df_with_stats['GlycanTypeCategory'] == gt) &
                df_with_stats['Has_Cancer'] &
                (df_with_stats['Regulation'] == REGULATION_UPREGULATED)
            ]
            if len(gt_data_up) > 0:
                source_indices.append(node_dict[GROUP_CANCER])
                target_indices.append(node_dict[f'{gt}_{REGULATION_SUFFIX_UP}'])
                values.append(len(gt_data_up))
                link_colors.append(rgba)
                link_labels.append(f"{GROUP_CANCER} → {gt} {REGULATION_UPREGULATED}<br>Count: {len(gt_data_up)}")

            # Downregulated flows
            gt_data_down = df_with_stats[
                (df_with_stats['GlycanTypeCategory'] == gt) &
                df_with_stats['Has_Cancer'] &
                (df_with_stats['Regulation'] == REGULATION_DOWNREGULATED)
            ]
            if len(gt_data_down) > 0:
                source_indices.append(node_dict[GROUP_CANCER])
                target_indices.append(node_dict[f'{gt}_{REGULATION_SUFFIX_DOWN}'])
                values.append(len(gt_data_down))
                link_colors.append(rgba)
                link_labels.append(f"{GROUP_CANCER} → {gt} {REGULATION_DOWNREGULATED}<br>Count: {len(gt_data_down)}")

        # Process Normal flows (using constants and helper)
        for gt in glycan_types:
            # Get base glycan type color and convert to RGBA (using helper)
            base_color = EXTENDED_CATEGORY_COLORS.get(gt, '#CCCCCC')
            rgba = self._hex_to_rgba(base_color, LINK_ALPHA_GROUP)

            # Upregulated flows
            gt_data_up = df_with_stats[
                (df_with_stats['GlycanTypeCategory'] == gt) &
                df_with_stats['Has_Normal'] &
                (df_with_stats['Regulation'] == REGULATION_UPREGULATED)
            ]
            if len(gt_data_up) > 0:
                source_indices.append(node_dict[GROUP_NORMAL])
                target_indices.append(node_dict[f'{gt}_{REGULATION_SUFFIX_UP}'])
                values.append(len(gt_data_up))
                link_colors.append(rgba)
                link_labels.append(f"{GROUP_NORMAL} → {gt} {REGULATION_UPREGULATED}<br>Count: {len(gt_data_up)}")

            # Downregulated flows
            gt_data_down = df_with_stats[
                (df_with_stats['GlycanTypeCategory'] == gt) &
                df_with_stats['Has_Normal'] &
                (df_with_stats['Regulation'] == REGULATION_DOWNREGULATED)
            ]
            if len(gt_data_down) > 0:
                source_indices.append(node_dict[GROUP_NORMAL])
                target_indices.append(node_dict[f'{gt}_{REGULATION_SUFFIX_DOWN}'])
                values.append(len(gt_data_down))
                link_colors.append(rgba)
                link_labels.append(f"{GROUP_NORMAL} → {gt} {REGULATION_DOWNREGULATED}<br>Count: {len(gt_data_down)}")

        logger.info(f"Total flows created: {len(source_indices)} links")

        # ====================================================================
        # Data Validation: Check for Over-Counting (per group)
        # ====================================================================
        total_glycopeptides = len(df_with_stats)
        total_cancer_present = df_with_stats['Has_Cancer'].sum()
        total_normal_present = df_with_stats['Has_Normal'].sum()

        logger.info(f"Data validation (Group Sankey):")
        logger.info(f"  Total glycopeptides: {total_glycopeptides}")
        logger.info(f"  Present in Cancer: {total_cancer_present}")
        logger.info(f"  Present in Normal: {total_normal_present}")
        logger.info(f"  Total flows: {len(values)}")

        # Note: Total flow is NOT expected to equal total_glycopeptides here
        # because each glycopeptide can be present in Cancer, Normal, or both
        # So total flow = cancer flows + normal flows (may be > total_glycopeptides)

        # ====================================================================
        # Boundary Value Check
        # ====================================================================
        logger.info(f"Boundary value checks (thresholds: Log2FC={log2fc_threshold}, FDR={fdr_threshold}):")

        # Check boundary cases
        boundary_log2fc = df_with_stats[
            (df_with_stats['Log2_Fold_Change'].abs() == log2fc_threshold) &
            (df_with_stats['FDR'] < fdr_threshold)
        ]
        if len(boundary_log2fc) > 0:
            logger.info(f"  Boundary Log2FC (±{log2fc_threshold}): {len(boundary_log2fc)} glycopeptides")

        boundary_fdr = df_with_stats[df_with_stats['FDR'] == fdr_threshold]
        if len(boundary_fdr) > 0:
            logger.info(f"  Boundary FDR ({fdr_threshold}): {len(boundary_fdr)} glycopeptides")

        # ====================================================================
        # Calculate Node Totals Directly for Accurate Hover Display
        # ====================================================================
        node_totals_corrected = []

        # Group nodes: count present glycopeptides
        node_totals_corrected.append(int(total_cancer_present))  # Cancer
        node_totals_corrected.append(int(total_normal_present))  # Normal

        # Glycan type nodes (10 nodes: 5 types × 2 regulation states)
        for gt in glycan_types:
            # Up regulation: present in Cancer OR Normal with Upregulated status
            count_up = len(df_with_stats[
                (df_with_stats['GlycanTypeCategory'] == gt) &
                (df_with_stats['Regulation'] == REGULATION_UPREGULATED) &
                (df_with_stats['Has_Cancer'] | df_with_stats['Has_Normal'])
            ])
            node_totals_corrected.append(count_up)

            # Down regulation: present in Cancer OR Normal with Downregulated status
            count_down = len(df_with_stats[
                (df_with_stats['GlycanTypeCategory'] == gt) &
                (df_with_stats['Regulation'] == REGULATION_DOWNREGULATED) &
                (df_with_stats['Has_Cancer'] | df_with_stats['Has_Normal'])
            ])
            node_totals_corrected.append(count_down)

        logger.info(f"Node totals calculated: {len(node_totals_corrected)} nodes")

        # Create Sankey diagram with corrected node hover
        # Position: 2 groups on left + 10 glycan types on right
        x_coords = [0.01, 0.01] + [0.99] * 10  # Groups left, 10 glycan nodes right

        # Y-coordinates: distribute 10 nodes evenly on right side
        y_coords = [0.2, 0.8]  # Cancer and Normal
        y_coords.extend([0.05 + i * 0.09 for i in range(10)])  # 10 glycan nodes evenly spaced

        fig = go.Figure(data=[go.Sankey(
            arrangement='snap',  # Snap to grid for vertical alignment
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=1.0),
                label=node_labels,
                color=node_colors,
                x=x_coords,
                y=y_coords,
                customdata=node_totals_corrected,  # Pass calculated totals
                hovertemplate='%{label}<br>Total: %{customdata}<extra></extra>'
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

        # Update layout with improved styling and clear definitions
        fig.update_layout(
            title=dict(
                text=(
                    f"<b>Glycan Type Distribution: Cancer vs Normal Groups</b><br>"
                    f"<sub><b>Regulation Definition (Cancer-centric):</b> Up = Cancer > Normal, Down = Cancer < Normal<br>"
                    f"<b>Criteria:</b> |Log2FC| ≥ {log2fc_threshold} (2-fold) AND FDR < {fdr_threshold}<br>"
                    f"<b>Note:</b> Each glycan type split into Up/Down regulation states • FDR=NaN treated as Non-significant</sub>"
                ),
                x=0.5,
                xanchor='center',
                font=dict(size=18, family='Arial, sans-serif', color='#2C3E50')
            ),
            font=dict(size=13, family='Arial, sans-serif', color='#2C3E50'),
            height=1400,  # Increased to accommodate 10 glycan nodes (Up/Down split)
            width=1200,
            plot_bgcolor='white',
            paper_bgcolor='#FAFAFA',  # Light gray background for better contrast
            margin=dict(l=80, r=80, t=170, b=250),  # Increased top margin for longer title
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
                    text='<b>Glycan Types (Up/Down)</b>',
                    showarrow=False,
                    font=dict(size=14, family='Arial, sans-serif', color='#2C3E50'),
                    xanchor='center'
                ),
                # Add legend box underneath main diagram
                # y=-0.15 places it well below the Sankey flows (which end at y=0)
                dict(
                    x=0.5, y=-0.15,
                    xref='paper', yref='paper',
                    text='<b>Link Colors (Glycan Type):</b><br>'
                         '<span style="color:#00CC00">━━━</span> HM (High-Mannose) | '
                         '<span style="color:#FF0000">━━━</span> F (Fucosylated) | '
                         '<span style="color:#FF69B4">━━━</span> S (Sialylated) | '
                         '<span style="color:#FFA500">━━━</span> SF (Sialofucosylated) | '
                         '<span style="color:#0000FF">━━━</span> C/H (Complex/Hybrid)',
                    showarrow=False,
                    font=dict(size=13, family='Arial, sans-serif', color='#2C3E50'),
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
            fig.write_image(str(output_png), width=1200, height=1400, scale=2)
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
