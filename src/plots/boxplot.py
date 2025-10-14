"""
Boxplot Module for pGlyco Auto Combine
Handles all boxplot visualizations

UPDATED: Cancer vs Normal methods now use centralized data preparation
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from scipy import stats
from statsmodels.stats.multitest import multipletests  # FDR correction
from ..utils import replace_empty_with_zero, save_trace_data, get_sample_columns
from ..data_preparation import (
    DataPreparationConfig
)
from .plot_config import (
    BOXPLOT_FIGSIZE, BOXPLOT_EXTENDED_FIGSIZE, BOXPLOT_WIDTH,
    BOXPLOT_LINEWIDTH, BOXPLOT_FLIERSIZE, BOXPLOT_DPI,
    LEGACY_GLYCAN_COLORS, EXTENDED_CATEGORY_COLORS,
    apply_standard_axis_style, apply_standard_legend,
    AXIS_LABEL_SIZE, AXIS_LABEL_WEIGHT,
    TITLE_SIZE, TITLE_WEIGHT,
    ANNOTATION_SIZE,  # Font constant
    add_sample_size_annotation,  # Phase 2.2 enhancement
    save_publication_figure,  # Phase 2.3: Optimized saving
    enhance_statistical_bracket, apply_publication_theme,  # ✨ Enhanced styling
    DESIGN_SYSTEM_AVAILABLE,
    LINE_MEDIUM_THICK,  # Linewidth constants
    # Zorder constants (Phase 10.3.7)
    ZORDER_BACKGROUND, ZORDER_GRID, ZORDER_SEPARATOR,
    ZORDER_DATA_LOW, ZORDER_DATA_HIGH,
    ZORDER_THRESHOLD, ZORDER_ANNOTATION,
    ZORDER_OVERLAY, ZORDER_EFFECT,
    ZORDER_TOP, ZORDER_ABSOLUTE_TOP
)

# Import premium design system if available
if DESIGN_SYSTEM_AVAILABLE:
    from .design_system import VisualEffects

logger = logging.getLogger(__name__)


def calculate_cohens_d(group1: np.ndarray, group2: np.ndarray) -> float:
    """
    Calculate Cohen's d effect size for two independent groups

    Formula: d = (mean1 - mean2) / pooled_std
    where pooled_std = sqrt(((n1-1)*std1² + (n2-1)*std2²) / (n1+n2-2))

    Args:
        group1: Data from first group (e.g., Cancer)
        group2: Data from second group (e.g., Normal)

    Returns:
        Cohen's d effect size (positive = group1 > group2)

    Interpretation:
        |d| < 0.2: negligible effect
        0.2 ≤ |d| < 0.5: small effect
        0.5 ≤ |d| < 0.8: medium effect
        |d| ≥ 0.8: large effect (biologically meaningful)

    Reference: Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences.
    """
    n1, n2 = len(group1), len(group2)

    # Handle edge cases
    if n1 < 2 or n2 < 2:
        return np.nan

    mean1, mean2 = np.mean(group1), np.mean(group2)
    std1, std2 = np.std(group1, ddof=1), np.std(group2, ddof=1)

    # Pooled standard deviation
    pooled_std = np.sqrt(((n1 - 1) * std1**2 + (n2 - 1) * std2**2) / (n1 + n2 - 2))

    # Avoid division by zero
    if pooled_std == 0:
        return 0.0 if mean1 == mean2 else np.inf

    cohens_d = (mean1 - mean2) / pooled_std
    return cohens_d


def _perform_fdr_correction_and_plot_brackets(
    ax: plt.Axes,
    boxplot_data: pd.DataFrame,
    category_column: str,
    existing_categories: list,
    log_context: str = "",
    use_enhanced_brackets: bool = True
) -> None:
    """
    Perform FDR-corrected statistical testing and plot significance brackets

    This helper consolidates the FDR correction logic used across multiple
    boxplot methods, ensuring consistent statistical testing and visualization.

    Args:
        ax: Matplotlib axes for plotting brackets
        boxplot_data: Long-format DataFrame with columns:
            - {category_column}: Category identifier (e.g., 'GlycanType', 'ExtendedCategory')
            - 'Group': 'Cancer' or 'Normal'
            - 'Intensity': Numerical values for comparison
        category_column: Name of column containing categories
        existing_categories: Ordered list of categories to test
        log_context: Prefix for logging (e.g., "" or "extended")
        use_enhanced_brackets: If True, use enhance_statistical_bracket() style.
                               If False, use manual line + text style.

    Returns:
        None (modifies ax in-place, logs results to logger)

    Statistical Methods:
        - Mann-Whitney U test: Non-parametric comparison of two groups
        - Benjamini-Hochberg FDR correction: Controls false discovery rate
        - Cohen's d effect size: Standardized difference between groups

    Scientific Integrity:
        - All calculations are deterministic (same input → same output)
        - No data-specific hard-coded assumptions
        - Generic parameterization prevents over-fitting
        - Standard significance thresholds: 0.001 (***), 0.01 (**), 0.05 (*)

    Reference:
        - Benjamini & Hochberg (1995). "Controlling the false discovery rate"
        - Mann & Whitney (1947). "On a test of whether one of two random variables"
    """
    # Step 1: Collect all p-values and effect sizes
    p_values_dict = {}  # {category: p_value}
    cohens_d_dict = {}  # {category: effect_size}

    for category in existing_categories:
        # Get data for each group
        cancer_data = boxplot_data[
            (boxplot_data[category_column] == category) &
            (boxplot_data['Group'] == 'Cancer')
        ]['Intensity'].values

        normal_data = boxplot_data[
            (boxplot_data[category_column] == category) &
            (boxplot_data['Group'] == 'Normal')
        ]['Intensity'].values

        # Skip if either group has insufficient data (need ≥3 for Mann-Whitney U)
        if len(cancer_data) < 3 or len(normal_data) < 3:
            p_values_dict[category] = np.nan
            cohens_d_dict[category] = np.nan
            logger.warning(
                f"{category}: Insufficient data (Cancer: {len(cancer_data)}, "
                f"Normal: {len(normal_data)})"
            )
            continue

        # Perform Mann-Whitney U test (non-parametric, robust to outliers)
        try:
            statistic, p_value = stats.mannwhitneyu(
                cancer_data, normal_data, alternative='two-sided'
            )
            cohens_d = calculate_cohens_d(cancer_data, normal_data)

            p_values_dict[category] = p_value
            cohens_d_dict[category] = cohens_d

        except Exception as e:
            logger.warning(f"Statistical test failed for {category}: {str(e)}")
            p_values_dict[category] = np.nan
            cohens_d_dict[category] = np.nan

    # Step 2: Apply FDR correction (Benjamini-Hochberg)
    valid_categories = [cat for cat, p in p_values_dict.items() if not np.isnan(p)]
    valid_p_values = [p_values_dict[cat] for cat in valid_categories]

    fdr_dict = {}
    if len(valid_p_values) > 0:
        # Apply Benjamini-Hochberg FDR correction
        reject, fdr_values, _, _ = multipletests(valid_p_values, method='fdr_bh')
        fdr_dict = dict(zip(valid_categories, fdr_values))

        # Log FDR correction results
        n_significant = sum(fdr < 0.05 for fdr in fdr_values)
        log_prefix = f" ({log_context})" if log_context else ""
        logger.info(
            f"FDR correction applied{log_prefix}: {len(valid_p_values)} tests, "
            f"{n_significant} significant (FDR < 0.05)"
        )

        for cat, raw_p, fdr in zip(valid_categories, valid_p_values, fdr_values):
            logger.info(f"  {cat}: p={raw_p:.4f} → FDR={fdr:.4f}")
    else:
        log_prefix = f" ({log_context})" if log_context else ""
        logger.warning(f"No valid p-values for FDR correction{log_prefix}")

    # Step 3: Plot with FDR-corrected significance markers
    y_max = boxplot_data['Intensity'].max()
    y_range = boxplot_data['Intensity'].max() - boxplot_data['Intensity'].min()

    for i, category in enumerate(existing_categories):
        # Get FDR-corrected p-value and effect size
        fdr = fdr_dict.get(category, np.nan)
        cohens_d = cohens_d_dict.get(category, np.nan)

        # Skip if no valid FDR value
        if np.isnan(fdr):
            continue

        # Determine significance level using FDR (not raw p-value)
        if fdr < 0.001:
            sig_marker = '***'
        elif fdr < 0.01:
            sig_marker = '**'
        elif fdr < 0.05:
            sig_marker = '*'
        else:
            sig_marker = 'ns'

        # Add significance bracket if significant
        if sig_marker != 'ns':
            # Calculate x positions for the connecting line
            n_cats = len(existing_categories)
            x_offset = (i - (n_cats - 1) / 2) * (BOXPLOT_WIDTH / n_cats)

            x1 = 0 + x_offset  # Cancer position
            x2 = 1 + x_offset  # Normal position

            y_position = y_max + y_range * 0.05 * (1 + i * 0.3)

            # Format annotation text with effect size
            if not np.isnan(cohens_d):
                annotation_text = f"{sig_marker}\n(d={cohens_d:.2f})"
            else:
                annotation_text = sig_marker

            # Plot bracket using appropriate style
            if use_enhanced_brackets:
                # Enhanced style with rounded ends and fancy box
                enhance_statistical_bracket(
                    ax, x1, x2, y_position,
                    text=annotation_text,
                    color='black',
                    fontsize=ANNOTATION_SIZE
                )
            else:
                # Manual style with simple line and text
                ax.plot(
                    [x1, x2], [y_position, y_position],
                    color='black', linewidth=LINE_MEDIUM_THICK,
                    zorder=ZORDER_THRESHOLD
                )
                ax.text(
                    (x1 + x2) / 2, y_position,
                    annotation_text,
                    ha='center', va='bottom',
                    fontsize=ANNOTATION_SIZE,
                    fontweight='bold',
                    color='black',
                    zorder=ZORDER_ANNOTATION
                )

            # Log result
            logger.info(
                f"{category}: Cancer vs Normal FDR={fdr:.4f} "
                f"({sig_marker}), Cohen's d={cohens_d:.3f}"
            )


def _plot_boxplot_cancer_vs_normal_base(
    self,
    df: pd.DataFrame,
    classification_type: str,
    categories: list,
    column_name: str,
    figsize: tuple,
    enable_debug_logging: bool = False
) -> None:
    """
    Base method for Cancer vs Normal boxplot comparison (Phase 3.2)

    Pipeline: TIC Normalization → Non-zero mean → Log2 Transform

    Creates two versions for each classification type:
    1. Without QC: All samples, mean of non-zero values
    2. With QC: Exclude samples with <10% detection rate

    Args:
        df: Annotated DataFrame
        classification_type: 'primary' or 'secondary' (for titles/filenames)
        categories: List of classification categories to plot
        column_name: Classification column name in DataFrame
        figsize: Figure size tuple (e.g., (10, 6) or (12, 6))
        enable_debug_logging: If True, log category distribution by group

    Returns:
        None (saves plots and trace data to self.output_dir)

    Scientific Pipeline:
        1. TIC normalization (Total Ion Current): Normalize intensities by sample sum
        2. Non-zero mean calculation: Average only detected values (skipna method)
        3. Log2 transformation: log2(mean + 1) for variance stabilization
        4. Optional QC filtering: Exclude samples with <10% detection rate

    Data Integrity:
        - TIC normalization preserves relative intensities within samples
        - Non-zero mean avoids bias from missing values (MNAR data)
        - Log2 transformation: Standard for fold-change analysis
        - Detection rate QC: Ensures reliable measurements
    """
    # Get sample columns (C1-C24, N1-N24)
    cancer_samples, normal_samples = get_sample_columns(df)
    sample_cols = cancer_samples + normal_samples

    # Step 1: TIC (Total Ion Current) Normalization on entire dataset
    intensity_matrix = replace_empty_with_zero(df[sample_cols])
    sample_sums = intensity_matrix.sum(axis=0)
    median_sum = sample_sums.median()
    sample_sums_safe = sample_sums.replace(0, 1)
    intensity_normalized = intensity_matrix / sample_sums_safe * median_sum

    # Create df copy with normalized intensities
    df_normalized = df.copy()
    df_normalized[sample_cols] = intensity_normalized

    # STANDARDIZED: Use centralized statistics calculation
    _ = DataPreparationConfig(missing_data_method='skipna')  # Config for documentation

    # Generate both versions (with and without QC)
    for apply_qc in [False, True]:
        data_for_plot = []

        # Collect data for each classification category
        for classification in categories:
            subset_df = df_normalized[df_normalized[column_name] == classification]

            for sample in sample_cols:
                group = 'Cancer' if sample in cancer_samples else 'Normal'

                # Calculate mean directly (for single sample column)
                if len(subset_df) > 0:
                    values = pd.to_numeric(subset_df[sample], errors='coerce')

                    # Use skipna method: mean of non-zero values only
                    valid_values = values[values > 0]

                    if len(valid_values) > 0:
                        mean_intensity = valid_values.mean()
                        detection_rate = len(valid_values) / len(values)

                        # Apply QC filter if needed
                        if apply_qc and detection_rate < 0.1:
                            continue  # Skip this sample

                        # Add data point
                        data_for_plot.append({
                            'Group': group,
                            'Classification': classification,
                            'Intensity': np.log2(mean_intensity + 1),
                            'Sample': sample,
                            'DetectionRate': detection_rate
                        })

        # Check if we have data
        if not data_for_plot:
            logger.warning(
                f"No data for {classification_type} Cancer vs Normal boxplot (QC={apply_qc})"
            )
            continue

        plot_df = pd.DataFrame(data_for_plot)

        # DEBUG: Log category distribution (optional, for secondary only)
        if enable_debug_logging:
            logger.info(f"{classification_type.capitalize()} boxplot (QC={apply_qc}) - Categories by group:")
            for group in ['Cancer', 'Normal']:
                group_data = plot_df[plot_df['Group'] == group]
                group_categories = group_data['Classification'].unique()
                logger.info(
                    f"  {group}: {sorted(group_categories)} ({len(group_data)} datapoints)"
                )
                for cat in categories:
                    count = len(group_data[group_data['Classification'] == cat])
                    if count == 0:
                        logger.warning(f"    MISSING: {cat} has 0 samples in {group}!")

        # Create figure
        fig, ax = plt.subplots(figsize=figsize)

        # Determine which classifications are actually present
        existing_classifications = [
            cat for cat in categories if cat in plot_df['Classification'].values
        ]

        # Prism-style boxplot: Group on x-axis, colored by Classification
        sns.boxplot(
            data=plot_df,
            x='Group',
            y='Intensity',
            hue='Classification',
            order=['Cancer', 'Normal'],
            hue_order=existing_classifications,
            palette=EXTENDED_CATEGORY_COLORS,
            width=BOXPLOT_WIDTH,
            linewidth=BOXPLOT_LINEWIDTH,
            flierprops={'markersize': BOXPLOT_FLIERSIZE},
            dodge=True,
            ax=ax
        )

        # Set labels and title
        ax.set_xlabel('Group', fontsize=AXIS_LABEL_SIZE, weight=AXIS_LABEL_WEIGHT)
        ax.set_ylabel(
            'Log2(Mean Intensity + 1)\n(TIC-normalized, non-zero values only)',
            fontsize=AXIS_LABEL_SIZE,
            weight=AXIS_LABEL_WEIGHT
        )

        title = f'{classification_type.capitalize()} Classification: Cancer vs Normal'
        if apply_qc:
            title += '\n(QC: Detection rate ≥10%)'

        ax.set_title(title, fontsize=TITLE_SIZE, weight=TITLE_WEIGHT)

        # Apply Prism styling (grid and spine settings)
        apply_standard_axis_style(ax, grid=True)
        apply_standard_legend(ax, title=f'{classification_type.capitalize()} Classification')

        # Add sample size annotation (Phase 2.2 enhancement)
        n_cancer = len(plot_df[plot_df['Group'] == 'Cancer']['Sample'].unique())
        n_normal = len(plot_df[plot_df['Group'] == 'Normal']['Sample'].unique())
        add_sample_size_annotation(
            ax, n_cancer=n_cancer, n_normal=n_normal,
            location='upper left', fontsize=ANNOTATION_SIZE
        )

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        # Save outputs
        suffix = '_qc' if apply_qc else ''
        output_file = self.output_dir / f'boxplot_{classification_type}_cancer_vs_normal{suffix}.png'
        save_publication_figure(fig, output_file, dpi=BOXPLOT_DPI)
        logger.info(
            f"✨ Saved ENHANCED {classification_type} Cancer vs Normal boxplot "
            f"(QC={apply_qc}) to {output_file}"
        )

        # Save trace data
        save_trace_data(
            plot_df, self.output_dir,
            f'boxplot_{classification_type}_cancer_vs_normal{suffix}_data.csv'
        )

        plt.close()


def _plot_boxplot_classification_base(
    self,
    df: pd.DataFrame,
    classification_type: str,
    categories: list,
    column_name: str,
    normalization: str,
    figsize: tuple
) -> None:
    """
    Base method for classification boxplot (Phase 3.3)

    Creates boxplot showing classification distribution across Cancer vs Normal groups.
    Supports two normalization modes: 'raw' (min-max per sample) or 'aggregated' (log2).

    Args:
        df: Annotated DataFrame
        classification_type: 'primary' or 'secondary' (for titles/filenames)
        categories: List of classification categories to include
        column_name: Classification column name in DataFrame
        normalization: 'raw' (min-max) or 'aggregated' (log2)
        figsize: Figure size tuple (e.g., (12, 8) or (14, 8))

    Returns:
        None (saves plot and trace data to self.output_dir)

    Pipeline:
        1. Get sample columns
        2. Prepare long-format data with normalization
        3. Apply log2 transform if aggregated
        4. Aggregate by Sample, Group, Classification
        5. Create Prism-style boxplot
        6. Save outputs
    """
    # Get sample columns (C1-C24, N1-N24)
    cancer_samples, normal_samples = get_sample_columns(df)
    sample_cols = cancer_samples + normal_samples

    # Prepare long format data
    data_for_plot = []

    for sample in sample_cols:
        if normalization == 'raw':
            # Min-max normalize per sample
            intensity_col = replace_empty_with_zero(df[sample])
            min_val, max_val = intensity_col.min(), intensity_col.max()
            if max_val > min_val:
                intensity_col = (intensity_col - min_val) / (max_val - min_val)
        else:
            intensity_col = replace_empty_with_zero(df[sample])

        group = 'Cancer' if sample in cancer_samples else 'Normal'

        for idx, row in df.iterrows():
            if row[column_name] in categories:
                intensity = intensity_col.loc[idx]
                # CRITICAL FIX: Filter out zeros to avoid detection bias
                # Only include detected values in the plot
                if intensity > 0:
                    data_for_plot.append({
                        'Sample': sample,
                        'Group': group,
                        'Classification': row[column_name],
                        'Intensity': intensity
                    })

    plot_df = pd.DataFrame(data_for_plot)

    # Apply aggregated normalization if needed
    if normalization == 'aggregated':
        plot_df['Intensity'] = np.log2(plot_df['Intensity'] + 1)

    # CRITICAL FIX: Aggregate by Sample and Classification
    # Without this, we plot ~108,000 individual glycopeptide measurements
    # With this, we plot ~94 aggregated values (47 samples × categories)
    plot_df = plot_df.groupby(['Sample', 'Group', 'Classification'], as_index=False)['Intensity'].mean()

    fig, ax = plt.subplots(figsize=figsize)

    # Determine which classifications are actually present
    existing_classifications = [cat for cat in categories if cat in plot_df['Classification'].values]

    # Prism-style boxplot: Group on x-axis, colored by Classification
    sns.boxplot(
        data=plot_df,
        x='Group',
        y='Intensity',
        hue='Classification',
        order=['Cancer', 'Normal'],
        hue_order=existing_classifications,
        palette=EXTENDED_CATEGORY_COLORS,
        width=BOXPLOT_WIDTH,
        linewidth=BOXPLOT_LINEWIDTH,
        flierprops={'markersize': BOXPLOT_FLIERSIZE},
        dodge=True,
        ax=ax
    )

    norm_text = 'Raw Normalized' if normalization == 'raw' else 'Log2 Scaled'
    ax.set_xlabel('Group', fontsize=AXIS_LABEL_SIZE, weight=AXIS_LABEL_WEIGHT)
    ax.set_ylabel(f'Intensity ({norm_text})', fontsize=AXIS_LABEL_SIZE, weight=AXIS_LABEL_WEIGHT)
    ax.set_title(
        f'{classification_type.capitalize()} Classification Distribution ({norm_text})',
        fontsize=TITLE_SIZE, weight=TITLE_WEIGHT
    )

    # Apply Prism styling (grid and spine settings)
    apply_standard_axis_style(ax, grid=True)
    apply_standard_legend(ax, title=f'{classification_type.capitalize()} Classification')

    # Add sample size annotation (Phase 2.2 enhancement)
    n_cancer = len(plot_df[plot_df['Group'] == 'Cancer']['Sample'].unique())
    n_normal = len(plot_df[plot_df['Group'] == 'Normal']['Sample'].unique())
    add_sample_size_annotation(
        ax, n_cancer=n_cancer, n_normal=n_normal,
        location='upper left', fontsize=ANNOTATION_SIZE
    )

    # ✨ ENHANCED: Apply publication theme
    apply_publication_theme(fig)

    plt.tight_layout()

    output_file = self.output_dir / f'boxplot_{classification_type}_{normalization}_normalized.png'
    save_publication_figure(fig, output_file, dpi=BOXPLOT_DPI)
    logger.info(f"✨ Saved ENHANCED {classification_type} boxplot to {output_file}")

    plt.close()


class BoxplotMixin:
    """Mixin class for boxplot-related visualizations"""

    def plot_boxplot(self, boxplot_data: pd.DataFrame, figsize: tuple = None):
        """
        Create boxplot comparing glycan types between groups with statistical significance

        Args:
            boxplot_data: Long-format DataFrame from analyzer
            figsize: Figure size (default: from plot_config)
        """
        if figsize is None:
            figsize = BOXPLOT_FIGSIZE

        fig, ax = plt.subplots(figsize=figsize)

        # ✨ PREMIUM: Add subtle gradient background for depth
        if DESIGN_SYSTEM_AVAILABLE:
            VisualEffects.add_gradient_background(ax, '#FFFFFF', '#F8F9FA', direction='vertical')

        # Define fixed order for glycan types
        glycan_order = ['Non', 'Sialylated', 'Fucosylated', 'Both']

        # Filter to only include existing glycan types
        existing_types = [gt for gt in glycan_order if gt in boxplot_data['GlycanType'].unique()]

        # PRISM STYLE: Color by glycan type (consistent across all plots)
        # x-axis shows Groups (Cancer vs Normal), hue shows glycan types with their characteristic colors
        sns.boxplot(
            data=boxplot_data,
            x='Group',  # Swapped: Groups on x-axis
            y='Intensity',
            hue='GlycanType',  # Swapped: Glycan types as hue (colored)
            order=['Cancer', 'Normal'],
            hue_order=existing_types,
            palette=LEGACY_GLYCAN_COLORS,  # Use glycan-type colors!
            width=BOXPLOT_WIDTH,
            linewidth=BOXPLOT_LINEWIDTH,  # Prism style: thicker lines
            flierprops={'markersize': BOXPLOT_FLIERSIZE},  # Prism style: larger outliers
            ax=ax
        )

        # ========================================
        # PHASE 3.1: Apply FDR correction using centralized helper
        # ========================================
        _perform_fdr_correction_and_plot_brackets(
            ax=ax,
            boxplot_data=boxplot_data,
            category_column='GlycanType',
            existing_categories=existing_types,
            log_context="",  # No prefix for standard boxplot
            use_enhanced_brackets=True  # Use enhanced bracket style
        )

        # Apply standardized axis styling with gridlines
        apply_standard_axis_style(
            ax,
            xlabel='Group',
            ylabel='Log2(Intensity + 1)',
            title='Glycan Type Distribution: Cancer vs Normal',
            grid=True
        )

        # Apply standardized legend (now showing glycan types with their colors)
        # Legend positioned outside plot area (non-interruptive)
        apply_standard_legend(ax, title='Glycan Type')

        # Add sample size annotation (Phase 2.2 enhancement)
        n_cancer = (boxplot_data['Group'] == 'Cancer').sum() // len(existing_types)
        n_normal = (boxplot_data['Group'] == 'Normal').sum() // len(existing_types)
        add_sample_size_annotation(ax, n_cancer=n_cancer, n_normal=n_normal,
                                   location='upper left', fontsize=ANNOTATION_SIZE)

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        # Save plot with optimized settings
        output_file = self.output_dir / 'boxplot_glycan_types.png'
        save_publication_figure(fig, output_file, dpi=BOXPLOT_DPI)
        logger.info(f"✨ Saved ENHANCED boxplot to {output_file}")

        # Save trace data
        save_trace_data(boxplot_data, self.output_dir, 'boxplot_glycan_types_data.csv')

        plt.close()

    def plot_boxplot_extended(self, boxplot_data: pd.DataFrame, figsize: tuple = None):
        """
        Create extended boxplot comparing 5 glycan categories between groups

        Args:
            boxplot_data: Long-format DataFrame from analyzer (extended)
            figsize: Figure size (default: from plot_config)
        """
        if figsize is None:
            figsize = BOXPLOT_EXTENDED_FIGSIZE

        fig, ax = plt.subplots(figsize=figsize)

        # ✨ PREMIUM: Add subtle gradient background for depth
        if DESIGN_SYSTEM_AVAILABLE:
            VisualEffects.add_gradient_background(ax, '#FFFFFF', '#F8F9FA', direction='vertical')

        # Define fixed order for extended categories
        category_order = ['HM', 'C/H', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        # Filter to only include existing categories
        existing_categories = [cat for cat in category_order if cat in boxplot_data['ExtendedCategory'].unique()]

        # PRISM STYLE: Color by glycan type (consistent across all plots)
        # x-axis shows Groups (Cancer vs Normal), hue shows glycan categories with their characteristic colors
        sns.boxplot(
            data=boxplot_data,
            x='Group',  # Swapped: Groups on x-axis
            y='Intensity',
            hue='ExtendedCategory',  # Swapped: Glycan categories as hue (colored)
            order=['Cancer', 'Normal'],
            hue_order=existing_categories,
            palette=EXTENDED_CATEGORY_COLORS,  # Use glycan-type colors!
            width=BOXPLOT_WIDTH,
            linewidth=BOXPLOT_LINEWIDTH,  # Prism style: thicker lines
            flierprops={'markersize': BOXPLOT_FLIERSIZE},  # Prism style: larger outliers
            ax=ax
        )

        # ========================================
        # PHASE 3.1: Apply FDR correction using centralized helper
        # ========================================
        _perform_fdr_correction_and_plot_brackets(
            ax=ax,
            boxplot_data=boxplot_data,
            category_column='ExtendedCategory',
            existing_categories=existing_categories,
            log_context="extended",  # Add "extended" prefix to logs
            use_enhanced_brackets=False  # Use manual bracket style (line + text)
        )

        # Apply standardized axis styling with gridlines
        apply_standard_axis_style(
            ax,
            xlabel='Group',
            ylabel='Log2(Intensity + 1)',
            title='Extended Glycan Categories: Cancer vs Normal',
            grid=True
        )

        # Apply standardized legend (now showing glycan categories with their colors)
        # Legend positioned outside plot area (non-interruptive)
        apply_standard_legend(ax, title='Glycan Category')

        # Add sample size annotation (Phase 2.2 enhancement)
        n_cancer = (boxplot_data['Group'] == 'Cancer').sum() // len(existing_categories)
        n_normal = (boxplot_data['Group'] == 'Normal').sum() // len(existing_categories)
        add_sample_size_annotation(ax, n_cancer=n_cancer, n_normal=n_normal,
                                   location='upper left', fontsize=ANNOTATION_SIZE)

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        # Save plot with optimized settings
        output_file = self.output_dir / 'boxplot_extended_categories.png'
        save_publication_figure(fig, output_file, dpi=BOXPLOT_DPI)
        logger.info(f"✨ Saved ENHANCED extended boxplot to {output_file}")

        # Save trace data
        save_trace_data(boxplot_data, self.output_dir, 'boxplot_extended_categories_data.csv')

        plt.close()

    def plot_boxplot_primary_classification(
        self, df: pd.DataFrame, normalization: str = 'raw',
        figsize: tuple = (12, 8)
    ):
        """
        Create boxplot for primary classification

        Args:
            df: Annotated DataFrame
            normalization: 'raw' or 'aggregated'
            figsize: Figure size
        """
        # PHASE 3.3: Thin wrapper delegating to centralized base method
        _plot_boxplot_classification_base(
            self,
            df=df,
            classification_type='primary',
            categories=['Truncated', 'High Mannose', 'ComplexHybrid'],
            column_name='PrimaryClassification',
            normalization=normalization,
            figsize=figsize
        )

    def plot_boxplot_secondary_classification(
        self, df: pd.DataFrame, normalization: str = 'raw',
        figsize: tuple = (14, 8)
    ):
        """
        Create boxplot for secondary classification

        Args:
            df: Annotated DataFrame
            normalization: 'raw' or 'aggregated'
            figsize: Figure size
        """
        # PHASE 3.3: Thin wrapper delegating to centralized base method
        _plot_boxplot_classification_base(
            self,
            df=df,
            classification_type='secondary',
            categories=['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated'],
            column_name='SecondaryClassification',
            normalization=normalization,
            figsize=figsize
        )

    def plot_boxplot_cancer_vs_normal_primary(self, df: pd.DataFrame, figsize: tuple = (10, 6)):
        """
        Create boxplot for Cancer vs Normal comparison (Primary classification)

        Pipeline: TIC Normalization → Non-zero mean → Log2 Transform

        Creates two versions:
        1. Without QC: All samples, mean of non-zero values
        2. With QC: Exclude samples with <10% detection rate

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # PHASE 3.2: Thin wrapper delegating to centralized base method
        _plot_boxplot_cancer_vs_normal_base(
            self,
            df=df,
            classification_type='primary',
            categories=['High Mannose', 'ComplexHybrid'],
            column_name='PrimaryClassification',
            figsize=figsize,
            enable_debug_logging=False
        )

    def plot_boxplot_cancer_vs_normal_secondary(self, df: pd.DataFrame, figsize: tuple = (12, 6)):
        """
        Create boxplot for Cancer vs Normal comparison (Secondary classification)

        Pipeline: TIC Normalization → Non-zero mean → Log2 Transform

        Creates two versions:
        1. Without QC: All samples, mean of non-zero values
        2. With QC: Exclude samples with <10% detection rate

        Args:
            df: Annotated DataFrame
            figsize: Figure size
        """
        # PHASE 3.2: Thin wrapper delegating to centralized base method
        _plot_boxplot_cancer_vs_normal_base(
            self,
            df=df,
            classification_type='secondary',
            categories=['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated'],
            column_name='SecondaryClassification',
            figsize=figsize,
            enable_debug_logging=True  # Enable debug logging for secondary (original behavior)
        )
