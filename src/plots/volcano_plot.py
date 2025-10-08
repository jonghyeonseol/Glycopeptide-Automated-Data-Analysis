"""
Volcano Plot Module for pGlyco Auto Combine
Handles differential expression visualization

UPDATED: Now uses centralized data preparation for consistency
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
from adjustText import adjust_text
from ..utils import save_trace_data
from ..data_preparation import (
    DataPreparationConfig,
    prepare_visualization_data,
    calculate_statistical_significance
)
from .plot_config import (
    VOLCANO_POINT_SIZE, VOLCANO_POINT_ALPHA,
    VOLCANO_THRESHOLD_LINEWIDTH, VOLCANO_THRESHOLD_ALPHA,
    VOLCANO_POINT_EDGEWIDTH, VOLCANO_LABEL_FONTSIZE,
    VOLCANO_LABEL_WEIGHT, VOLCANO_LABEL_PADDING, VOLCANO_LABEL_LINEWIDTH,
    GLYCAN_COLORS,
    apply_standard_axis_style, apply_standard_legend,
    add_sample_size_annotation,  # Phase 2.2 enhancement
    get_regulation_style  # Phase 3 enhancement
)

logger = logging.getLogger(__name__)


class VolcanoPlotMixin:
    """Mixin class for Volcano plot visualization"""

    def _calculate_bootstrap_fold_change_ci(self, cancer_values: np.ndarray,
                                            normal_values: np.ndarray,
                                            n_bootstrap: int = 1000,
                                            confidence_level: float = 0.95) -> tuple:
        """
        Calculate bootstrap confidence interval for fold change

        Args:
            cancer_values: Cancer group intensities (filtered for non-zero)
            normal_values: Normal group intensities (filtered for non-zero)
            n_bootstrap: Number of bootstrap iterations
            confidence_level: Confidence level (default 0.95 for 95% CI)

        Returns:
            Tuple of (log2fc, log2fc_ci_lower, log2fc_ci_upper)
        """
        # Handle edge cases
        if len(cancer_values) < 2 or len(normal_values) < 2:
            return (np.nan, np.nan, np.nan)

        # Calculate observed fold change
        cancer_mean = np.mean(cancer_values)
        normal_mean = np.mean(normal_values)
        observed_log2fc = np.log2((cancer_mean + 1) / (normal_mean + 1))

        # Bootstrap resampling
        bootstrap_log2fcs = []
        np.random.seed(42)  # For reproducibility

        for _ in range(n_bootstrap):
            # Resample with replacement
            cancer_resample = np.random.choice(cancer_values, size=len(cancer_values), replace=True)
            normal_resample = np.random.choice(normal_values, size=len(normal_values), replace=True)

            # Calculate log2FC for this bootstrap sample
            cancer_mean_boot = np.mean(cancer_resample)
            normal_mean_boot = np.mean(normal_resample)
            log2fc_boot = np.log2((cancer_mean_boot + 1) / (normal_mean_boot + 1))

            bootstrap_log2fcs.append(log2fc_boot)

        bootstrap_log2fcs = np.array(bootstrap_log2fcs)

        # Calculate confidence interval (percentile method)
        alpha = 1 - confidence_level
        ci_lower = np.percentile(bootstrap_log2fcs, alpha / 2 * 100)
        ci_upper = np.percentile(bootstrap_log2fcs, (1 - alpha / 2) * 100)

        return (observed_log2fc, ci_lower, ci_upper)

    def plot_volcano(self, df: pd.DataFrame, vip_df: pd.DataFrame,
                     config: DataPreparationConfig = None,
                     fdr_threshold: float = 0.05, fc_threshold: float = 1.5,
                     figsize: tuple = (12, 10)):
        """
        Create volcano plot showing log2(fold change) vs -log10(FDR)

        UPDATED: Uses centralized data preparation for consistency

        Args:
            df: Annotated DataFrame with intensity data
            vip_df: DataFrame with VIP scores
            config: DataPreparationConfig (uses default if None)
            fdr_threshold: FDR significance threshold (default 0.05)
            fc_threshold: Fold change threshold (default 1.5)
            figsize: Figure size (width, height)
        """
        # Use default config if not provided
        if config is None:
            config = DataPreparationConfig(
                min_detection_pct=0.30,
                min_samples=5,
                missing_data_method='skipna'
            )

        logger.info("Creating volcano plot with STANDARDIZED data preparation...")

        # STANDARDIZED DATA PREPARATION (eliminates inline filtering)
        # Data is already filtered by DataPipeline in main.py - no need to refilter
        volcano_df = prepare_visualization_data(
            df=df,
            config=config,
            vip_scores=vip_df,
            merge_method='left',  # Keep all glycopeptides, add VIP where available
            apply_detection_filter=False,  # Data already filtered in main.py
            log_prefix="[Volcano Plot] "
        )

        # Get sample columns for statistical tests
        from ..utils import get_sample_columns
        cancer_samples, normal_samples = get_sample_columns(df)

        # STANDARDIZED STATISTICAL SIGNIFICANCE
        volcano_df = calculate_statistical_significance(
            df_prep=volcano_df,
            cancer_samples=cancer_samples,
            normal_samples=normal_samples,
            method='mannwhitneyu',
            fdr_correction=True
        )

        # Add glycopeptide identifier
        volcano_df['Glycopeptide'] = volcano_df['Peptide'] + '_' + volcano_df['GlycanComposition']

        # Rename for backward compatibility
        volcano_df['Log2FC'] = volcano_df['Log2_Fold_Change']
        volcano_df['-Log10FDR'] = -np.log10(volcano_df['FDR'])

        logger.info(f"Volcano plot data prepared: {len(volcano_df)} glycopeptides")

        # Phase 2.2: Calculate bootstrap confidence intervals for top significant features
        logger.info("Calculating bootstrap CIs for fold changes (top 20 significant features)...")

        # Select top 20 most significant glycopeptides for CI calculation
        top_significant = volcano_df[volcano_df['FDR'] < 0.05].nlargest(
            20, '-Log10FDR') if len(volcano_df[volcano_df['FDR'] < 0.05]) > 0 else pd.DataFrame()

        if len(top_significant) > 0:
            ci_results = []
            for idx in top_significant.index:
                # Get raw intensities for this glycopeptide
                cancer_vals = df.loc[idx, cancer_samples].replace('', 0).astype(float).values
                normal_vals = df.loc[idx, normal_samples].replace('', 0).astype(float).values

                # Filter non-zero values for CI calculation
                cancer_vals = cancer_vals[cancer_vals > 0]
                normal_vals = normal_vals[normal_vals > 0]

                # Calculate bootstrap CI
                log2fc, ci_lower, ci_upper = self._calculate_bootstrap_fold_change_ci(
                    cancer_vals, normal_vals, n_bootstrap=1000
                )

                ci_results.append({
                    'Index': idx,
                    'Log2FC_CI_Lower': ci_lower,
                    'Log2FC_CI_Upper': ci_upper
                })

            # Merge CI results back into volcano_df
            ci_df = pd.DataFrame(ci_results).set_index('Index')
            volcano_df = volcano_df.join(ci_df, how='left')

            logger.info(f"✓ Bootstrap CIs calculated for {len(ci_results)} glycopeptides")
        else:
            # Add empty CI columns if no significant features
            volcano_df['Log2FC_CI_Lower'] = np.nan
            volcano_df['Log2FC_CI_Upper'] = np.nan
            logger.info("No significant features found for CI calculation")

        # Classification
        volcano_df['Regulation'] = 'Non-significant'

        # Up-regulated: log2FC > threshold AND FDR < threshold
        up_mask = (volcano_df['Log2FC'] > np.log2(fc_threshold)) & (volcano_df['FDR'] < fdr_threshold)
        volcano_df.loc[up_mask, 'Regulation'] = 'Up in Cancer'

        # Down-regulated: log2FC < -threshold AND FDR < threshold
        down_mask = (volcano_df['Log2FC'] < -np.log2(fc_threshold)) & (volcano_df['FDR'] < fdr_threshold)
        volcano_df.loc[down_mask, 'Regulation'] = 'Down in Cancer'

        # Create plot
        fig, ax = plt.subplots(figsize=figsize)

        # Phase 3: Use color + shape encoding for colorblind accessibility
        # Plot points by regulation status
        for regulation in ['Up in Cancer', 'Down in Cancer', 'Non-significant']:
            mask = volcano_df['Regulation'] == regulation
            subset = volcano_df[mask]

            if len(subset) == 0:
                continue

            # Get color and marker (Phase 3: colorblind-safe)
            color, marker = get_regulation_style(regulation)

            # Size by VIP score (scale between 50 and 150 for better visibility)
            if len(subset) > 0 and subset['VIP_Score'].max() > 0:
                sizes = 50 + (subset['VIP_Score'] / subset['VIP_Score'].max() * 100)
            else:
                sizes = VOLCANO_POINT_SIZE

            # Shape symbols for legend clarity
            shape_symbol = {'Up in Cancer': '▲', 'Down in Cancer': '▼', 'Non-significant': '●'}[regulation]

            ax.scatter(subset['Log2FC'], subset['-Log10FDR'],
                       c=color, s=sizes, alpha=VOLCANO_POINT_ALPHA,
                       edgecolors='black', linewidths=VOLCANO_POINT_EDGEWIDTH,
                       marker=marker,  # Phase 3: Shape encoding
                       label=f"{regulation} ({shape_symbol}) n={len(subset)}", zorder=3)

        # Phase 2.2: Add confidence interval error bars for top significant features
        features_with_ci = volcano_df[~volcano_df['Log2FC_CI_Lower'].isna()]
        if len(features_with_ci) > 0:
            for idx, row in features_with_ci.iterrows():
                # Calculate error bar lengths (asymmetric)
                x_err_lower = row['Log2FC'] - row['Log2FC_CI_Lower']
                x_err_upper = row['Log2FC_CI_Upper'] - row['Log2FC']

                # Determine color based on regulation status
                if row['Regulation'] == 'Up in Cancer':
                    err_color = GROUP_PALETTE['Cancer']
                elif row['Regulation'] == 'Down in Cancer':
                    err_color = GROUP_PALETTE['Normal']
                else:
                    err_color = '#95A5A6'

                # Plot horizontal error bar (95% CI for Log2FC)
                ax.errorbar(
                    row['Log2FC'], row['-Log10FDR'],
                    xerr=[[x_err_lower], [x_err_upper]],
                    fmt='none',  # No marker
                    ecolor=err_color,
                    elinewidth=2,
                    capsize=4,
                    capthick=2,
                    alpha=0.7,
                    zorder=4
                )

            logger.info(f"Added 95% CI error bars for {len(features_with_ci)} significant features")

        # Add threshold lines using standardized styling
        ax.axhline(-np.log10(fdr_threshold), color='gray', linestyle='--',
                   linewidth=VOLCANO_THRESHOLD_LINEWIDTH, alpha=VOLCANO_THRESHOLD_ALPHA, zorder=1)
        ax.axvline(np.log2(fc_threshold), color='gray', linestyle='--',
                   linewidth=VOLCANO_THRESHOLD_LINEWIDTH, alpha=VOLCANO_THRESHOLD_ALPHA, zorder=1)
        ax.axvline(-np.log2(fc_threshold), color='gray', linestyle='--',
                   linewidth=VOLCANO_THRESHOLD_LINEWIDTH, alpha=VOLCANO_THRESHOLD_ALPHA, zorder=1)

        # Use standardized glycan type colors from plot_config
        glycan_type_colors = GLYCAN_COLORS

        # =============================================================================
        # INTELLIGENT LABEL SELECTION - Show only top VOLCANO_MAX_LABELS
        # =============================================================================
        # Strategy: Select the most significant AND most different glycopeptides
        # 1. Must pass FDR < 0.05
        # 2. Must have substantial effect size (|log2FC| > 0.5)
        # 3. Ranked by combined score: |log2FC| × -log10(FDR) × (1 + VIP)
        # 4. Spatially distributed (avoid clustering in one region)

        significant_candidates = volcano_df[
            (volcano_df['FDR'] < 0.05) & (abs(volcano_df['Log2FC']) > 0.5)
        ].copy()

        if len(significant_candidates) > 0:
            # Combined annotation priority score
            significant_candidates['AnnotationScore'] = (
                abs(significant_candidates['Log2FC']) *
                significant_candidates['-Log10FDR'] *
                (1 + significant_candidates['VIP_Score'])
            )

            # Select top VOLCANO_MAX_LABELS, ensuring spatial distribution
            # Pick 1 from left (down-regulated), 1 from top, 1 from right (up-regulated)
            down_reg = significant_candidates[significant_candidates['Log2FC'] < -1].nlargest(1, 'AnnotationScore')
            up_reg = significant_candidates[significant_candidates['Log2FC'] > 1].nlargest(1, 'AnnotationScore')
            top_sig = significant_candidates.nlargest(1, '-Log10FDR')

            # Combine and deduplicate
            to_annotate = pd.concat([down_reg, up_reg, top_sig]).drop_duplicates()
        else:
            to_annotate = pd.DataFrame()

        # =============================================================================
        # CREATE HIGH-QUALITY LABELS
        # =============================================================================
        texts = []
        for _, row in to_annotate.iterrows():
            peptide = row['Peptide']
            glycan = row['GlycanComposition']

            # Smart abbreviation: keep peptide readable
            if len(peptide) > 10:
                peptide_display = f"{peptide[:7]}..."
            else:
                peptide_display = peptide

            # Two-line label for clarity
            label = f"{peptide_display}\n{glycan}"

            # Color by glycan type for biological context
            glycan_type = row.get('GlycanTypeCategory', 'C/H')
            text_color = glycan_type_colors.get(glycan_type, '#000000')

            # Create label with MAXIMUM VISIBILITY settings
            texts.append(ax.text(
                row['Log2FC'], row['-Log10FDR'], label,
                fontsize=VOLCANO_LABEL_FONTSIZE,      # Extra large (14pt)
                fontweight=VOLCANO_LABEL_WEIGHT,       # Bold
                ha='center', va='center',
                color=text_color,
                bbox=dict(
                    boxstyle=f'round,pad={VOLCANO_LABEL_PADDING}',  # Extra padding
                    facecolor='white',
                    edgecolor=text_color,
                    alpha=0.98,                         # Nearly opaque
                    linewidth=VOLCANO_LABEL_LINEWIDTH   # Extra thick border (2.5)
                ),
                zorder=100  # Ensure labels are on top
            ))

        # =============================================================================
        # SMART LABEL POSITIONING - Avoid all overlaps
        # =============================================================================
        if len(texts) > 0:
            try:
                adjust_text(
                    texts,
                    arrowprops=dict(arrowstyle='->', color='black', lw=1.5, alpha=0.7),
                    expand_points=(2.0, 2.0),   # Large exclusion zone around points
                    expand_text=(2.0, 2.0),     # Large buffer between labels
                    force_points=0.7,            # Strong repulsion from data
                    force_text=0.7,              # Strong repulsion between labels
                    ax=ax,
                    lim=1000                     # More iterations for better placement
                )
            except Exception as e:
                logger.warning(f"adjustText optimization failed: {e}")
                # Fallback: simple offset positioning
                for i, text_obj in enumerate(texts):
                    text_obj.set_position((text_obj.get_position()[0], text_obj.get_position()[1] + (i * 0.3)))

        # Apply standardized styling
        apply_standard_axis_style(
            ax,
            xlabel='Log2 Fold Change (Cancer / Normal)',
            ylabel='-Log10 FDR',
            title='Volcano Plot: Differential Glycopeptide Expression\nCancer vs Normal',
            grid=True
        )

        # Apply standardized legend (positioned outside plot area)
        apply_standard_legend(ax)

        # Add statistics text
        n_up = len(volcano_df[volcano_df['Regulation'] == 'Up in Cancer'])
        n_down = len(volcano_df[volcano_df['Regulation'] == 'Down in Cancer'])
        n_ns = len(volcano_df[volcano_df['Regulation'] == 'Non-significant'])
        n_with_ci = len(features_with_ci) if len(features_with_ci) > 0 else 0

        stats_text = f"FC threshold: {fc_threshold}x | FDR < {fdr_threshold}\n"
        stats_text += f"Up: {n_up} | Down: {n_down} | NS: {n_ns}\n"
        if n_with_ci > 0:
            stats_text += f"95% CI shown for top {n_with_ci} features"

        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), zorder=10)

        # Add sample size annotation (Phase 2.2 enhancement)
        n_cancer = len(cancer_samples)
        n_normal = len(normal_samples)
        add_sample_size_annotation(ax, n_cancer=n_cancer, n_normal=n_normal,
                                   location='lower right', fontsize=10)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'volcano_plot.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved volcano plot to {output_file}")

        # Save trace data
        save_trace_data(volcano_df, self.output_dir, 'volcano_plot_data.csv')

        plt.close()

        return volcano_df
