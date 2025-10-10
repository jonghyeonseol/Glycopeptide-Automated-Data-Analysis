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
    add_sample_size_annotation,  # Phase 2.2 enhancement
    save_publication_figure,  # Phase 2.3: Optimized saving
    enhance_statistical_bracket, apply_publication_theme  # ✨ Enhanced styling
)

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

        # Perform statistical tests and add significance markers
        # Now comparing Cancer vs Normal for each glycan type
        y_max = boxplot_data['Intensity'].max()
        y_range = boxplot_data['Intensity'].max() - boxplot_data['Intensity'].min()

        # Add significance markers between Cancer and Normal for each glycan type
        for i, glycan_type in enumerate(existing_types):
            # Get data for each group
            cancer_data = boxplot_data[
                (boxplot_data['GlycanType'] == glycan_type) &
                (boxplot_data['Group'] == 'Cancer')
            ]['Intensity'].values

            normal_data = boxplot_data[
                (boxplot_data['GlycanType'] == glycan_type) &
                (boxplot_data['Group'] == 'Normal')
            ]['Intensity'].values

            # Skip if either group has insufficient data
            if len(cancer_data) < 3 or len(normal_data) < 3:
                continue

            # Perform Mann-Whitney U test (non-parametric)
            try:
                statistic, p_value = stats.mannwhitneyu(cancer_data, normal_data, alternative='two-sided')

                # Calculate Cohen's d effect size (Phase 1.1 enhancement)
                cohens_d = calculate_cohens_d(cancer_data, normal_data)

                # Determine significance level
                if p_value < 0.001:
                    sig_marker = '***'
                elif p_value < 0.01:
                    sig_marker = '**'
                elif p_value < 0.05:
                    sig_marker = '*'
                else:
                    sig_marker = 'ns'

                # ✨ ENHANCED: Add significance bracket if significant
                if sig_marker != 'ns':
                    # Calculate x positions for the connecting line
                    n_types = len(existing_types)
                    x_offset = (i - (n_types - 1) / 2) * (BOXPLOT_WIDTH / n_types)

                    x1 = 0 + x_offset  # Cancer position
                    x2 = 1 + x_offset  # Normal position

                    y_position = y_max + y_range * 0.05 * (1 + i * 0.3)

                    # Format annotation text with effect size
                    if not np.isnan(cohens_d):
                        annotation_text = f"{sig_marker}\n(d={cohens_d:.2f})"
                    else:
                        annotation_text = sig_marker

                    # ✨ Use enhanced statistical bracket (rounded ends, fancy box)
                    enhance_statistical_bracket(
                        ax, x1, x2, y_position,
                        text=annotation_text,
                        color='black',
                        fontsize=11
                    )

                    logger.info(
                        f"{glycan_type}: Cancer vs Normal p={p_value:.4f} "
                        f"({sig_marker}), Cohen's d={cohens_d:.3f}"
                    )

            except Exception as e:
                logger.warning(f"Statistical test failed for {glycan_type}: {str(e)}")

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
                                   location='upper left', fontsize=10)

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

        # Perform statistical tests and add significance markers
        # Now comparing Cancer vs Normal for each glycan category
        y_max = boxplot_data['Intensity'].max()
        y_range = boxplot_data['Intensity'].max() - boxplot_data['Intensity'].min()

        # Add significance markers between Cancer and Normal for each category
        for i, category in enumerate(existing_categories):
            # Get data for each group
            cancer_data = boxplot_data[
                (boxplot_data['ExtendedCategory'] == category) &
                (boxplot_data['Group'] == 'Cancer')
            ]['Intensity'].values

            normal_data = boxplot_data[
                (boxplot_data['ExtendedCategory'] == category) &
                (boxplot_data['Group'] == 'Normal')
            ]['Intensity'].values

            # Skip if either group has insufficient data
            if len(cancer_data) < 3 or len(normal_data) < 3:
                continue

            # Perform Mann-Whitney U test (non-parametric)
            try:
                statistic, p_value = stats.mannwhitneyu(cancer_data, normal_data, alternative='two-sided')

                # Calculate Cohen's d effect size (Phase 1.1 enhancement)
                cohens_d = calculate_cohens_d(cancer_data, normal_data)

                # Determine significance level
                if p_value < 0.001:
                    sig_marker = '***'
                elif p_value < 0.01:
                    sig_marker = '**'
                elif p_value < 0.05:
                    sig_marker = '*'
                else:
                    sig_marker = 'ns'

                # Add significance marker if significant (connecting Cancer and Normal groups)
                if sig_marker != 'ns':
                    # Calculate x positions for the connecting line
                    n_cats = len(existing_categories)
                    x_offset = (i - (n_cats - 1) / 2) * (BOXPLOT_WIDTH / n_cats)

                    x1 = 0 + x_offset  # Cancer position
                    x2 = 1 + x_offset  # Normal position

                    y_position = y_max + y_range * 0.05 * (1 + i * 0.3)

                    # Draw line connecting the two groups
                    ax.plot([x1, x2], [y_position, y_position],
                            color='black', linewidth=1.5, zorder=10)

                    # Add significance marker WITH effect size (Phase 1.1)
                    if not np.isnan(cohens_d):
                        annotation_text = f"{sig_marker}\n(d={cohens_d:.2f})"
                    else:
                        annotation_text = sig_marker

                    ax.text(
                        (x1 + x2) / 2,
                        y_position,
                        annotation_text,
                        ha='center',
                        va='bottom',
                        fontsize=11,  # Slightly smaller for two-line text
                        fontweight='bold',
                        color='black',
                        zorder=11
                    )

                    logger.info(
                        f"{category}: Cancer vs Normal p={p_value:.4f} "
                        f"({sig_marker}), Cohen's d={cohens_d:.3f}"
                    )

            except Exception as e:
                logger.warning(f"Statistical test failed for {category}: {str(e)}")

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
                                   location='upper left', fontsize=10)

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
        # Identify sample columns
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Primary classification categories
        primary_categories = ['Truncated', 'High Mannose', 'ComplexHybrid']

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

            group = 'Cancer' if sample.startswith('C') else 'Normal'

            for idx, row in df.iterrows():
                if row['PrimaryClassification'] in primary_categories:
                    intensity = intensity_col.loc[idx]
                    # CRITICAL FIX: Filter out zeros to avoid detection bias
                    # Only include detected values in the plot
                    if intensity > 0:
                        data_for_plot.append({
                            'Sample': sample,
                            'Group': group,
                            'Classification': row['PrimaryClassification'],
                            'Intensity': intensity
                        })

        plot_df = pd.DataFrame(data_for_plot)

        # Apply aggregated normalization if needed
        if normalization == 'aggregated':
            plot_df['Intensity'] = np.log2(plot_df['Intensity'] + 1)

        # CRITICAL FIX: Aggregate by Sample and Classification
        # Without this, we plot ~108,000 individual glycopeptide measurements
        # With this, we plot ~94 aggregated values (47 samples × 2 categories)
        plot_df = plot_df.groupby(['Sample', 'Group', 'Classification'], as_index=False)['Intensity'].mean()

        fig, ax = plt.subplots(figsize=figsize)

        # Determine which classifications are actually present
        existing_classifications = [cat for cat in primary_categories if cat in plot_df['Classification'].values]

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
        ax.set_title(f'Primary Classification Distribution ({norm_text})', fontsize=TITLE_SIZE, weight=TITLE_WEIGHT)

        # Apply Prism styling (grid and spine settings)
        apply_standard_axis_style(ax, grid=True)
        apply_standard_legend(ax, title='Primary Classification')

        # Add sample size annotation (Phase 2.2 enhancement)
        n_cancer = len(plot_df[plot_df['Group'] == 'Cancer']['Sample'].unique())
        n_normal = len(plot_df[plot_df['Group'] == 'Normal']['Sample'].unique())
        add_sample_size_annotation(ax, n_cancer=n_cancer, n_normal=n_normal,
                                   location='upper left', fontsize=10)

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        output_file = self.output_dir / f'boxplot_primary_{normalization}_normalized.png'
        save_publication_figure(fig, output_file, dpi=BOXPLOT_DPI)
        logger.info(f"✨ Saved ENHANCED primary boxplot to {output_file}")

        plt.close()

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
        # Identify sample columns
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Secondary classification categories
        secondary_categories = ['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

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

            group = 'Cancer' if sample.startswith('C') else 'Normal'

            for idx, row in df.iterrows():
                if row['SecondaryClassification'] in secondary_categories:
                    intensity = intensity_col.loc[idx]
                    # CRITICAL FIX: Filter out zeros to avoid detection bias
                    # Only include detected values in the plot
                    if intensity > 0:
                        data_for_plot.append({
                            'Sample': sample,
                            'Group': group,
                            'Classification': row['SecondaryClassification'],
                            'Intensity': intensity
                        })

        plot_df = pd.DataFrame(data_for_plot)

        # Apply aggregated normalization if needed
        if normalization == 'aggregated':
            plot_df['Intensity'] = np.log2(plot_df['Intensity'] + 1)

        # CRITICAL FIX: Aggregate by Sample and Classification
        # Without this, we plot ~108,000 individual glycopeptide measurements
        # With this, we plot ~235 aggregated values (47 samples × 5 categories)
        plot_df = plot_df.groupby(['Sample', 'Group', 'Classification'], as_index=False)['Intensity'].mean()

        fig, ax = plt.subplots(figsize=figsize)

        # Determine which classifications are actually present
        existing_classifications = [cat for cat in secondary_categories if cat in plot_df['Classification'].values]

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
        ax.set_title(f'Secondary Classification Distribution ({norm_text})', fontsize=TITLE_SIZE, weight=TITLE_WEIGHT)

        # Apply Prism styling (grid and spine settings)
        apply_standard_axis_style(ax, grid=True)
        apply_standard_legend(ax, title='Secondary Classification')

        # Add sample size annotation (Phase 2.2 enhancement)
        n_cancer = len(plot_df[plot_df['Group'] == 'Cancer']['Sample'].unique())
        n_normal = len(plot_df[plot_df['Group'] == 'Normal']['Sample'].unique())
        add_sample_size_annotation(ax, n_cancer=n_cancer, n_normal=n_normal,
                                   location='upper left', fontsize=10)

        # ✨ ENHANCED: Apply publication theme
        apply_publication_theme(fig)

        plt.tight_layout()

        output_file = self.output_dir / f'boxplot_secondary_{normalization}_normalized.png'
        save_publication_figure(fig, output_file, dpi=BOXPLOT_DPI)
        logger.info(f"✨ Saved ENHANCED secondary boxplot to {output_file}")

        plt.close()

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

        primary_categories = ['High Mannose', 'ComplexHybrid']

        # STANDARDIZED: Use centralized statistics calculation
        _ = DataPreparationConfig(missing_data_method='skipna')  # Config for documentation

        # Generate both versions
        for apply_qc in [False, True]:
            data_for_plot = []

            for classification in primary_categories:
                subset_df = df_normalized[df_normalized['PrimaryClassification'] == classification]

                for sample in sample_cols:
                    group = 'Cancer' if sample.startswith('C') else 'Normal'

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

            if not data_for_plot:
                logger.warning(f"No data for primary Cancer vs Normal boxplot (QC={apply_qc})")
                continue

            plot_df = pd.DataFrame(data_for_plot)

            fig, ax = plt.subplots(figsize=figsize)

            # Determine which classifications are actually present
            existing_classifications = [cat for cat in primary_categories if cat in plot_df['Classification'].values]

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

            ax.set_xlabel('Group', fontsize=AXIS_LABEL_SIZE, weight=AXIS_LABEL_WEIGHT)
            ax.set_ylabel('Log2(Mean Intensity + 1)\n(TIC-normalized, non-zero values only)',
                          fontsize=AXIS_LABEL_SIZE, weight=AXIS_LABEL_WEIGHT)

            title = 'Primary Classification: Cancer vs Normal'
            if apply_qc:
                title += '\n(QC: Detection rate ≥10%)'

            ax.set_title(title, fontsize=TITLE_SIZE, weight=TITLE_WEIGHT)

            # Apply Prism styling (grid and spine settings)
            apply_standard_axis_style(ax, grid=True)
            apply_standard_legend(ax, title='Primary Classification')

            # Add sample size annotation (Phase 2.2 enhancement)
            n_cancer = len(plot_df[plot_df['Group'] == 'Cancer']['Sample'].unique())
            n_normal = len(plot_df[plot_df['Group'] == 'Normal']['Sample'].unique())
            add_sample_size_annotation(ax, n_cancer=n_cancer, n_normal=n_normal,
                                       location='upper left', fontsize=10)

            # ✨ ENHANCED: Apply publication theme
            apply_publication_theme(fig)

            plt.tight_layout()

            suffix = '_qc' if apply_qc else ''
            output_file = self.output_dir / f'boxplot_primary_cancer_vs_normal{suffix}.png'
            save_publication_figure(fig, output_file, dpi=BOXPLOT_DPI)
            logger.info(f"✨ Saved ENHANCED primary Cancer vs Normal boxplot (QC={apply_qc}) to {output_file}")

            # Save trace data
            save_trace_data(plot_df, self.output_dir, f'boxplot_primary_cancer_vs_normal{suffix}_data.csv')

            plt.close()

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

        secondary_categories = ['High Mannose', 'Complex/Hybrid', 'Fucosylated', 'Sialylated', 'Sialofucosylated']

        # STANDARDIZED: Use centralized statistics calculation
        _ = DataPreparationConfig(missing_data_method='skipna')  # Config for documentation

        # Generate both versions
        for apply_qc in [False, True]:
            data_for_plot = []

            for classification in secondary_categories:
                subset_df = df_normalized[df_normalized['SecondaryClassification'] == classification]

                for sample in sample_cols:
                    group = 'Cancer' if sample.startswith('C') else 'Normal'

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

            if not data_for_plot:
                logger.warning(f"No data for secondary Cancer vs Normal boxplot (QC={apply_qc})")
                continue

            plot_df = pd.DataFrame(data_for_plot)

            # DEBUG: Log categories by group
            logger.info(f"Secondary boxplot (QC={apply_qc}) - Categories by group:")
            for group in ['Cancer', 'Normal']:
                group_data = plot_df[plot_df['Group'] == group]
                categories = group_data['Classification'].unique()
                logger.info(f"  {group}: {sorted(categories)} ({len(group_data)} datapoints)")
                for cat in secondary_categories:
                    count = len(group_data[group_data['Classification'] == cat])
                    if count == 0:
                        logger.warning(f"    MISSING: {cat} has 0 samples in {group}!")

            fig, ax = plt.subplots(figsize=figsize)

            # Determine which classifications are actually present
            existing_classifications = [cat for cat in secondary_categories if cat in plot_df['Classification'].values]

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

            ax.set_xlabel('Group', fontsize=AXIS_LABEL_SIZE, weight=AXIS_LABEL_WEIGHT)
            ax.set_ylabel('Log2(Mean Intensity + 1)\n(TIC-normalized, non-zero values only)',
                          fontsize=AXIS_LABEL_SIZE, weight=AXIS_LABEL_WEIGHT)

            title = 'Secondary Classification: Cancer vs Normal'
            if apply_qc:
                title += '\n(QC: Detection rate ≥10%)'

            ax.set_title(title, fontsize=TITLE_SIZE, weight=TITLE_WEIGHT)

            # Apply Prism styling (grid and spine settings)
            apply_standard_axis_style(ax, grid=True)
            apply_standard_legend(ax, title='Secondary Classification')

            # Add sample size annotation (Phase 2.2 enhancement)
            n_cancer = len(plot_df[plot_df['Group'] == 'Cancer']['Sample'].unique())
            n_normal = len(plot_df[plot_df['Group'] == 'Normal']['Sample'].unique())
            add_sample_size_annotation(ax, n_cancer=n_cancer, n_normal=n_normal,
                                       location='upper left', fontsize=10)

            # ✨ ENHANCED: Apply publication theme
            apply_publication_theme(fig)

            plt.tight_layout()

            suffix = '_qc' if apply_qc else ''
            output_file = self.output_dir / f'boxplot_secondary_cancer_vs_normal{suffix}.png'
            save_publication_figure(fig, output_file, dpi=BOXPLOT_DPI)
            logger.info(f"✨ Saved ENHANCED secondary Cancer vs Normal boxplot (QC={apply_qc}) to {output_file}")

            # Save trace data
            save_trace_data(plot_df, self.output_dir, f'boxplot_secondary_cancer_vs_normal{suffix}_data.csv')

            plt.close()
