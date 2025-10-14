"""
Volcano Plot Module for pGlyco Auto Combine
Handles differential expression visualization

UPDATED: Now uses centralized data preparation for consistency
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
from ..utils import save_trace_data
from ..data_preparation import (
    DataPreparationConfig,
    prepare_visualization_data,
    calculate_statistical_significance
)
from .plot_config import (
    VOLCANO_POINT_SIZE, VOLCANO_POINT_ALPHA, VOLCANO_DPI, VOLCANO_FIGSIZE,
    VOLCANO_THRESHOLD_LINEWIDTH, VOLCANO_THRESHOLD_ALPHA,
    VOLCANO_POINT_EDGEWIDTH, VOLCANO_LABEL_FONTSIZE,
    VOLCANO_LABEL_WEIGHT, VOLCANO_LABEL_PADDING, VOLCANO_LABEL_LINEWIDTH,
    PLOT_LINE_LINEWIDTH,  # Phase 10.3.3: Linewidth unification
    GLYCAN_COLORS,
    GROUP_PALETTE,
    apply_standard_axis_style, apply_standard_legend,
    add_sample_size_annotation,  # Phase 2.2 enhancement
    get_regulation_style,  # Phase 3 enhancement
    save_publication_figure,  # Phase 2.3: Optimized saving
    create_fancy_bbox, apply_publication_theme,  # ✨ Enhanced styling
    GRADIENT_ALPHA_START, GRADIENT_ALPHA_END,  # ✨ Gradient effects
    COLOR_CANCER, COLOR_NORMAL,  # Premium Material Design colors
    TITLE_SIZE, AXIS_LABEL_SIZE, ANNOTATION_SIZE,  # Enhanced typography
    DESIGN_SYSTEM_AVAILABLE,
    ALPHA_HIGH, ALPHA_VERY_HIGH,  # Alpha constants
    EDGE_COLOR_BLACK,  # Edge color standardization,
    # Linestyle constants (Phase 10.3.8)
    THRESHOLD_LINESTYLE,
    # Zorder constants (Phase 10.3.7)
    ZORDER_BACKGROUND, ZORDER_GRID, ZORDER_SEPARATOR,
    ZORDER_DATA_LOW, ZORDER_DATA_HIGH,
    ZORDER_THRESHOLD, ZORDER_ANNOTATION,
    ZORDER_OVERLAY, ZORDER_EFFECT,
    ZORDER_TOP, ZORDER_ABSOLUTE_TOP
)

# Import premium design system if available
if DESIGN_SYSTEM_AVAILABLE:
    from .design_system import VisualEffects, ColorSystem

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
                     fdr_threshold: float = 0.05, log2fc_threshold: float = 1.0,
                     figsize: tuple = None):
        """
        Create volcano plot showing log2(fold change) vs -log10(FDR)

        UPDATED: Uses centralized data preparation for consistency

        Args:
            df: Annotated DataFrame with intensity data
            vip_df: DataFrame with VIP scores
            config: DataPreparationConfig (uses default if None)
            fdr_threshold: FDR significance threshold (default 0.05)
            log2fc_threshold: Log2 fold change threshold (default 1.0 = 2-fold, common in literature)
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
                # Get Peptide and GlycanComposition for this row
                peptide = top_significant.loc[idx, 'Peptide']
                glycan = top_significant.loc[idx, 'GlycanComposition']

                # Find matching row in original df using Peptide+GlycanComposition
                mask = (df['Peptide'] == peptide) & (df['GlycanComposition'] == glycan)
                if mask.sum() == 0:
                    continue  # Skip if not found

                df_row = df[mask].iloc[0]

                # Get raw intensities for this glycopeptide
                cancer_vals = df_row[cancer_samples].replace('', 0).astype(float).values
                normal_vals = df_row[normal_samples].replace('', 0).astype(float).values

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

        # ========================================
        # PHASE 2.1 FIX: Explicit FDR=NaN handling with logging
        # ========================================
        # Count glycopeptides with NaN FDR (insufficient samples for test)
        fdr_nan_mask = volcano_df['FDR'].isna()
        n_fdr_nan = fdr_nan_mask.sum()

        if n_fdr_nan > 0:
            logger.warning(f"Found {n_fdr_nan} glycopeptides with FDR=NaN (insufficient samples <3 in one or both groups)")
            logger.info(f"  These will be classified as 'Non-significant' regardless of fold change")

        # Classification
        volcano_df['Regulation'] = 'Non-significant'

        # CRITICAL: Up-regulated requires VALID FDR (not NaN)
        # log2FC > threshold AND FDR < threshold AND FDR is not NaN
        up_mask = (
            (volcano_df['Log2FC'] > log2fc_threshold) &
            (volcano_df['FDR'] < fdr_threshold) &
            volcano_df['FDR'].notna()  # ← EXPLICIT: Exclude FDR=NaN
        )
        volcano_df.loc[up_mask, 'Regulation'] = 'Up in Cancer'

        # CRITICAL: Down-regulated requires VALID FDR (not NaN)
        # log2FC < -threshold AND FDR < threshold AND FDR is not NaN
        down_mask = (
            (volcano_df['Log2FC'] < -log2fc_threshold) &
            (volcano_df['FDR'] < fdr_threshold) &
            volcano_df['FDR'].notna()  # ← EXPLICIT: Exclude FDR=NaN
        )
        volcano_df.loc[down_mask, 'Regulation'] = 'Down in Cancer'

        # ========================================
        # PHASE 2.2 FIX: Data validation
        # ========================================
        # Validate: Total glycopeptides = sum of all regulation categories
        n_total = len(volcano_df)
        n_up = up_mask.sum()
        n_down = down_mask.sum()
        n_ns = (volcano_df['Regulation'] == 'Non-significant').sum()

        logger.info(f"Volcano plot regulation breakdown:")
        logger.info(f"  Total glycopeptides: {n_total}")
        logger.info(f"  Up-regulated (Cancer > Normal): {n_up}")
        logger.info(f"  Down-regulated (Cancer < Normal): {n_down}")
        logger.info(f"  Non-significant: {n_ns}")
        logger.info(f"  FDR=NaN (insufficient samples): {n_fdr_nan}")

        # VALIDATION: Check for data integrity
        if n_up + n_down + n_ns != n_total:
            logger.error(f"❌ DATA VALIDATION FAILED: Category sum ({n_up + n_down + n_ns}) ≠ Total ({n_total})")
            raise ValueError("Volcano plot data validation failed: categories don't sum to total")
        else:
            logger.info(f"✓ Data validation passed: {n_up} + {n_down} + {n_ns} = {n_total}")

        # VALIDATION: Check that all FDR=NaN are classified as Non-significant
        fdr_nan_but_significant = volcano_df[fdr_nan_mask & (volcano_df['Regulation'] != 'Non-significant')]
        if len(fdr_nan_but_significant) > 0:
            logger.error(f"❌ VALIDATION FAILED: {len(fdr_nan_but_significant)} glycopeptides with FDR=NaN incorrectly classified as significant!")
            raise ValueError("FDR=NaN glycopeptides incorrectly classified")
        else:
            logger.info(f"✓ FDR=NaN validation passed: All {n_fdr_nan} FDR=NaN glycopeptides are Non-significant")

        # Create plot with optimized figure size
        if figsize is None:
            figsize = VOLCANO_FIGSIZE
        fig, ax = plt.subplots(figsize=figsize)

        # ✨ PREMIUM: Add subtle gradient background for depth
        if DESIGN_SYSTEM_AVAILABLE:
            VisualEffects.add_gradient_background(ax, '#FFFFFF', '#F8F9FA', direction='vertical')

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
                       edgecolors=EDGE_COLOR_BLACK, linewidths=VOLCANO_POINT_EDGEWIDTH,
                       marker=marker,  # Phase 3: Shape encoding
                       label=f"{regulation} ({shape_symbol}) n={len(subset)}", zorder=ZORDER_DATA_LOW)

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
                    elinewidth=PLOT_LINE_LINEWIDTH,
                    capsize=4,
                    capthick=2,
                    alpha=ALPHA_HIGH,
                    zorder=ZORDER_DATA_HIGH
                )

            logger.info(f"Added 95% CI error bars for {len(features_with_ci)} significant features")

        # Add threshold lines using standardized styling
        ax.axhline(-np.log10(fdr_threshold), color='gray', linestyle=THRESHOLD_LINESTYLE,
                   linewidth=VOLCANO_THRESHOLD_LINEWIDTH, alpha=VOLCANO_THRESHOLD_ALPHA, zorder=ZORDER_SEPARATOR)
        ax.axvline(log2fc_threshold, color='gray', linestyle=THRESHOLD_LINESTYLE,
                   linewidth=VOLCANO_THRESHOLD_LINEWIDTH, alpha=VOLCANO_THRESHOLD_ALPHA, zorder=ZORDER_SEPARATOR)
        ax.axvline(-log2fc_threshold, color='gray', linestyle=THRESHOLD_LINESTYLE,
                   linewidth=VOLCANO_THRESHOLD_LINEWIDTH, alpha=VOLCANO_THRESHOLD_ALPHA, zorder=ZORDER_SEPARATOR)

        # Apply standardized styling
        apply_standard_axis_style(
            ax,
            xlabel='Log2 Fold Change (Cancer / Normal)',
            ylabel='-Log10 FDR',
            title='Volcano Plot: Differential Glycopeptide Expression\nCancer vs Normal',
            grid=True
        )

        # ✨ PREMIUM: Enhance title with better typography
        title = ax.get_title()
        ax.set_title(title, fontsize=TITLE_SIZE, fontweight='bold',
                    family='Inter', pad=20)

        # Apply standardized legend (positioned outside plot area)
        apply_standard_legend(ax)

        # Add statistics text
        n_up = len(volcano_df[volcano_df['Regulation'] == 'Up in Cancer'])
        n_down = len(volcano_df[volcano_df['Regulation'] == 'Down in Cancer'])
        n_ns = len(volcano_df[volcano_df['Regulation'] == 'Non-significant'])
        n_with_ci = len(features_with_ci) if len(features_with_ci) > 0 else 0

        stats_text = f"Log2 FC threshold: ±{log2fc_threshold} ({2**log2fc_threshold:.1f}-fold) | FDR < {fdr_threshold}\n"
        stats_text += f"Up: {n_up} | Down: {n_down} | NS: {n_ns}\n"
        if n_with_ci > 0:
            stats_text += f"95% CI shown for top {n_with_ci} features"

        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                fontsize=ANNOTATION_SIZE, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=ALPHA_VERY_HIGH), zorder=ZORDER_THRESHOLD)

        # Add sample size annotation (Phase 2.2 enhancement)
        n_cancer = len(cancer_samples)
        n_normal = len(normal_samples)
        add_sample_size_annotation(ax, n_cancer=n_cancer, n_normal=n_normal,
                                   location='lower right', fontsize=ANNOTATION_SIZE)

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        # Save plot with optimized settings (150 DPI for complex plot)
        output_file = self.output_dir / 'volcano_plot.png'
        save_publication_figure(fig, output_file, dpi=VOLCANO_DPI)
        logger.info(f"✨ Saved ENHANCED volcano plot to {output_file} (150 DPI, gradient glow effects)")

        # Save trace data
        save_trace_data(volcano_df, self.output_dir, 'volcano_plot_data.csv')

        plt.close()

        return volcano_df
