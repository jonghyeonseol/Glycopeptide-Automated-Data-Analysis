"""
Data Pipeline Module for pGlyco Auto Combine
Single source of truth for data filtering and preparation

CRITICAL: This module ensures ALL visualizations use the same filtered dataset,
guaranteeing consistent glycan-type ratios across all outputs.

Architecture:
    Raw Data → DataPipeline.filter_dataset() → Filtered Data → ALL Visualizations
"""

import pandas as pd
from pathlib import Path
import logging

from .data_preparation import (
    DataPreparationConfig,
    filter_by_detection_frequency,
    get_sample_columns
)
from .metadata_collector import get_metadata_collector

logger = logging.getLogger(__name__)


class DataPipeline:
    """
    Centralized data filtering pipeline

    Ensures consistency across all visualizations by:
    - Applying detection filter ONCE
    - Tracking before/after statistics
    - Providing detailed filtering reports
    - Saving both raw and filtered datasets
    """

    def __init__(self, config: DataPreparationConfig):
        """
        Initialize DataPipeline

        Args:
            config: DataPreparationConfig with filtering parameters
        """
        self.config = config
        self.filtering_stats = {}
        self.glycan_type_stats = {}

    def filter_dataset(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Apply detection frequency filter to dataset

        This is the SINGLE POINT where filtering occurs in the entire pipeline.
        All downstream analyses and visualizations use the filtered dataset.

        Args:
            df: Annotated DataFrame (raw, unfiltered)

        Returns:
            Filtered DataFrame
        """
        logger.info("=" * 80)
        logger.info("DATA PIPELINE: Applying Detection Filter (Single Source of Truth)")
        logger.info("=" * 80)

        # Store raw dataset statistics
        self._calculate_raw_statistics(df)

        # Apply detection filter
        df_filtered = filter_by_detection_frequency(df, self.config)

        # Store filtered dataset statistics
        self._calculate_filtered_statistics(df, df_filtered)

        # Log glycan-type ratio changes
        self._log_glycan_ratio_changes()

        logger.info("=" * 80)
        logger.info("GUARANTEE: All visualizations will use this filtered dataset")
        logger.info("RESULT: Consistent glycan-type ratios across all outputs")
        logger.info("=" * 80)

        return df_filtered

    def _calculate_raw_statistics(self, df: pd.DataFrame):
        """Calculate statistics for raw (unfiltered) dataset"""
        self.filtering_stats['total_before'] = len(df)

        # Glycan type distribution (raw)
        if 'GlycanType' in df.columns:
            glycan_counts_raw = df['GlycanType'].value_counts().to_dict()
            glycan_pct_raw = {k: v / len(df) * 100 for k, v in glycan_counts_raw.items()}

            self.glycan_type_stats['raw'] = {
                'counts': glycan_counts_raw,
                'percentages': glycan_pct_raw
            }

            logger.info(f"\nGlycan Type Distribution (RAW, n={len(df)}):")
            for glycan_type in sorted(glycan_counts_raw.keys()):
                count = glycan_counts_raw[glycan_type]
                pct = glycan_pct_raw[glycan_type]
                logger.info(f"  {glycan_type:15s}: {count:5d} ({pct:5.1f}%)")

    def _calculate_filtered_statistics(self, df_raw: pd.DataFrame, df_filtered: pd.DataFrame):
        """Calculate statistics for filtered dataset"""
        total_before = len(df_raw)
        total_after = len(df_filtered)
        total_removed = total_before - total_after
        pct_removed = total_removed / total_before * 100

        self.filtering_stats['total_after'] = total_after
        self.filtering_stats['total_removed'] = total_removed
        self.filtering_stats['pct_removed'] = pct_removed

        logger.info("\nDetection Filtering Results:")
        logger.info(f"  Threshold: ≥{self.config.min_detection_pct * 100:.0f}% detection in at least one group")
        logger.info(f"  Before:    {total_before:5d} glycopeptides")
        logger.info(f"  After:     {total_after:5d} glycopeptides")
        logger.info(f"  Removed:   {total_removed:5d} ({pct_removed:5.1f}%)")

        # Glycan type distribution (filtered)
        if 'GlycanType' in df_filtered.columns:
            glycan_counts_filtered = df_filtered['GlycanType'].value_counts().to_dict()
            glycan_pct_filtered = {k: v / len(df_filtered) * 100 for k, v in glycan_counts_filtered.items()}

            self.glycan_type_stats['filtered'] = {
                'counts': glycan_counts_filtered,
                'percentages': glycan_pct_filtered
            }

            logger.info(f"\nGlycan Type Distribution (FILTERED, n={len(df_filtered)}):")
            for glycan_type in sorted(glycan_counts_filtered.keys()):
                count = glycan_counts_filtered[glycan_type]
                pct = glycan_pct_filtered[glycan_type]
                logger.info(f"  {glycan_type:15s}: {count:5d} ({pct:5.1f}%)")

    def _log_glycan_ratio_changes(self):
        """Log how glycan-type ratios changed due to filtering"""
        if 'raw' not in self.glycan_type_stats or 'filtered' not in self.glycan_type_stats:
            return

        raw_pct = self.glycan_type_stats['raw']['percentages']
        filtered_pct = self.glycan_type_stats['filtered']['percentages']

        # Find all unique glycan types
        all_types = set(raw_pct.keys()) | set(filtered_pct.keys())

        logger.info("\nGlycan-Type Ratio Changes (Raw → Filtered):")
        logger.info(f"  {'Type':15s}  {'Raw %':>8s}  {'Filtered %':>12s}  {'Change':>10s}")
        logger.info(f"  {'-' * 15}  {'-' * 8}  {'-' * 12}  {'-' * 10}")

        for glycan_type in sorted(all_types):
            raw_p = raw_pct.get(glycan_type, 0)
            filt_p = filtered_pct.get(glycan_type, 0)
            change = filt_p - raw_p
            change_str = f"{change:+.1f}%"

            logger.info(f"  {glycan_type:15s}  {raw_p:7.1f}%  {filt_p:11.1f}%  {change_str:>10s}")

    def get_filtering_report(self) -> str:
        """
        Generate detailed filtering report for analysis summary

        Returns:
            Markdown-formatted report string
        """
        report = []
        report.append("\n" + "=" * 80)
        report.append("DATA FILTERING REPORT")
        report.append("=" * 80)
        report.append("")

        # Filter settings
        report.append("Filter Settings:")
        report.append(f"  - Detection threshold: ≥{self.config.min_detection_pct * 100:.0f}% in at least one group")
        report.append(f"  - Minimum samples: {self.config.min_samples}")
        report.append(f"  - Missing data method: {self.config.missing_data_method}")
        report.append("")

        # Filtering results
        if 'total_before' in self.filtering_stats:
            report.append("Filtering Results:")
            report.append(f"  - Before filtering: {self.filtering_stats['total_before']:,} glycopeptides")
            report.append(f"  - After filtering:  {self.filtering_stats['total_after']:,} glycopeptides")
            removed_count = self.filtering_stats['total_removed']
            removed_pct = self.filtering_stats['pct_removed']
            report.append(f"  - Removed:          {removed_count:,} ({removed_pct:.1f}%)")
            report.append(f"  - Reason: <{self.config.min_detection_pct * 100:.0f}% detection in both groups")
            report.append("")

        # Glycan type distributions
        if 'filtered' in self.glycan_type_stats:
            report.append("Glycan Type Distribution (FILTERED DATA - used in all visualizations):")

            filtered_counts = self.glycan_type_stats['filtered']['counts']
            filtered_pct = self.glycan_type_stats['filtered']['percentages']
            total = sum(filtered_counts.values())

            for glycan_type in sorted(filtered_counts.keys()):
                count = filtered_counts[glycan_type]
                pct = filtered_pct[glycan_type]
                report.append(f"  - {glycan_type:15s}: {count:5d} ({pct:5.1f}%)")
            report.append(f"  - {'TOTAL':15s}: {total:5d} (100.0%)")
            report.append("")

        # Data consistency guarantee
        report.append("Data Consistency Guarantee:")
        report.append("  ✓ All statistical tests use filtered data")
        report.append("  ✓ All visualizations use filtered data")
        report.append("  ✓ Glycan-type ratios are IDENTICAL across all outputs")
        report.append("  ✓ No double-filtering or inconsistent subsets")
        report.append("")
        report.append("=" * 80)

        return "\n".join(report)

    def save_datasets(self,
                      df_raw: pd.DataFrame,
                      df_filtered: pd.DataFrame,
                      output_dir: Path,
                      raw_filename: str = 'integrated.csv',
                      filtered_filename: str = 'integrated_filtered.csv'):
        """
        Save both raw and filtered datasets with ALCOA++ metadata headers

        Args:
            df_raw: Raw unfiltered dataset
            df_filtered: Filtered dataset
            output_dir: Output directory
            raw_filename: Filename for raw dataset
            filtered_filename: Filename for filtered dataset
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Get metadata collector
        metadata_collector = get_metadata_collector()

        # Save raw dataset with metadata header
        raw_path = output_dir / raw_filename
        with open(raw_path, 'w') as f:
            # Write metadata header
            f.write(metadata_collector.get_metadata_header_lines())
            f.write("# Dataset: RAW (Unfiltered)\n")
            f.write(f"# Glycopeptides: {len(df_raw)}\n")
            f.write("# Description: Complete integrated dataset before filtering\n")
            f.write("#\n")
            # Write DataFrame
            df_raw.to_csv(f, index=False)
        logger.info(f"Saved raw dataset: {raw_path} ({len(df_raw)} glycopeptides)")

        # Save filtered dataset with metadata header
        filtered_path = output_dir / filtered_filename
        with open(filtered_path, 'w') as f:
            # Write metadata header
            f.write(metadata_collector.get_metadata_header_lines())
            f.write("# Dataset: FILTERED (Used in all analyses)\n")
            f.write(f"# Glycopeptides: {len(df_filtered)}\n")
            f.write(f"# Filter: ≥{self.config.min_detection_pct * 100:.0f}% detection in at least one group\n")
            f.write("# Description: Quality-controlled dataset used in ALL downstream analyses\n")
            f.write("#\n")
            # Write DataFrame
            df_filtered.to_csv(f, index=False)
        logger.info(f"Saved filtered dataset: {filtered_path} ({len(df_filtered)} glycopeptides)")

        # Save filtering report with metadata
        report_path = output_dir / 'filtering_report.txt'
        with open(report_path, 'w') as f:
            # Write metadata header
            f.write(metadata_collector.get_metadata_header_lines())
            f.write("\n")
            # Write filtering report
            f.write(self.get_filtering_report())
        logger.info(f"Saved filtering report: {report_path}")

    def validate_filtering(self, df_raw: pd.DataFrame, df_filtered: pd.DataFrame) -> bool:
        """
        Validate that filtering was applied correctly

        Args:
            df_raw: Raw unfiltered dataset
            df_filtered: Filtered dataset

        Returns:
            True if validation passes

        Raises:
            ValueError if validation fails
        """
        # Check 1: Filtered should be subset of raw
        if len(df_filtered) > len(df_raw):
            raise ValueError(f"Filtered dataset ({len(df_filtered)}) is larger than raw ({len(df_raw)})")

        # Check 2: All peptide+glycan combinations in filtered should exist in raw
        raw_combos = set(zip(df_raw['Peptide'], df_raw['GlycanComposition']))
        filtered_combos = set(zip(df_filtered['Peptide'], df_filtered['GlycanComposition']))

        if not filtered_combos.issubset(raw_combos):
            raise ValueError("Filtered dataset contains combinations not in raw dataset")

        # Check 3: Each glycopeptide must meet filter criteria (detection % OR sample count)
        cancer_samples, normal_samples = get_sample_columns(df_filtered)

        cancer_detection = (df_filtered[cancer_samples] != '').sum(axis=1) / len(cancer_samples)
        normal_detection = (df_filtered[normal_samples] != '').sum(axis=1) / len(normal_samples)
        cancer_count = (df_filtered[cancer_samples] != '').sum(axis=1)
        normal_count = (df_filtered[normal_samples] != '').sum(axis=1)

        for idx in df_filtered.index:
            c_det = cancer_detection.loc[idx]
            n_det = normal_detection.loc[idx]
            c_cnt = cancer_count.loc[idx]
            n_cnt = normal_count.loc[idx]

            # Must meet: (≥threshold detection OR ≥min_samples) in at least one group
            cancer_ok = (c_det >= self.config.min_detection_pct - 0.001) or (c_cnt >= self.config.min_samples)
            normal_ok = (n_det >= self.config.min_detection_pct - 0.001) or (n_cnt >= self.config.min_samples)

            if not (cancer_ok or normal_ok):
                raise ValueError(
                    f"Glycopeptide at index {idx} doesn't meet filter criteria:\n"
                    f"  Cancer: {c_det:.1%} detection ({c_cnt} samples), "
                    f"Normal: {n_det:.1%} detection ({n_cnt} samples)"
                )

        logger.info("✓ Filtering validation passed")
        return True
