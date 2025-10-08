"""
Data Validation Module for pGlyco Auto Combine
Ensures data consistency across all visualizations

CRITICAL: Run validation checks after preparing data for multiple visualizations
to catch any inconsistencies before they affect scientific conclusions.
"""

import pandas as pd
import numpy as np
from typing import List, Dict
from .logger_config import get_logger
from .exceptions import ValidationError

logger = get_logger(__name__)


class DataConsistencyValidator:
    """Validator for checking data consistency across visualizations"""

    def __init__(self, tolerance: float = 1e-6):
        """
        Initialize validator

        Args:
            tolerance: Numerical tolerance for floating-point comparisons
        """
        self.tolerance = tolerance
        self.validation_results = []

    def validate_glycopeptide_overlap(self,
                                      df1: pd.DataFrame,
                                      df2: pd.DataFrame,
                                      name1: str = "Dataset1",
                                      name2: str = "Dataset2",
                                      min_overlap_pct: float = 0.80) -> bool:
        """
        Validate that two datasets have sufficient glycopeptide overlap

        Args:
            df1: First DataFrame
            df2: Second DataFrame
            name1: Name of first dataset (for reporting)
            name2: Name of second dataset (for reporting)
            min_overlap_pct: Minimum required overlap percentage

        Returns:
            True if validation passes

        Raises:
            ValidationError if overlap is below threshold
        """
        # Create glycopeptide identifiers
        gp1 = set(df1['Peptide'] + '|' + df1['GlycanComposition'])
        gp2 = set(df2['Peptide'] + '|' + df2['GlycanComposition'])

        overlap = gp1.intersection(gp2)
        only_in_1 = gp1 - gp2
        only_in_2 = gp2 - gp1

        overlap_pct_1 = len(overlap) / len(gp1) if len(gp1) > 0 else 0
        overlap_pct_2 = len(overlap) / len(gp2) if len(gp2) > 0 else 0

        result = {
            'check': 'glycopeptide_overlap',
            'dataset1': name1,
            'dataset2': name2,
            'total_in_1': len(gp1),
            'total_in_2': len(gp2),
            'overlap': len(overlap),
            'only_in_1': len(only_in_1),
            'only_in_2': len(only_in_2),
            'overlap_pct_1': overlap_pct_1,
            'overlap_pct_2': overlap_pct_2,
            'passed': overlap_pct_1 >= min_overlap_pct and overlap_pct_2 >= min_overlap_pct
        }

        self.validation_results.append(result)

        logger.info(f"Glycopeptide overlap check: {name1} vs {name2}")
        logger.info(f"  {name1}: {len(gp1)} glycopeptides")
        logger.info(f"  {name2}: {len(gp2)} glycopeptides")
        logger.info(f"  Overlap: {len(overlap)} ({overlap_pct_1 * 100:.1f}% of {name1}, "
                    f"{overlap_pct_2 * 100:.1f}% of {name2})")
        logger.info(f"  Only in {name1}: {len(only_in_1)}")
        logger.info(f"  Only in {name2}: {len(only_in_2)}")

        if not result['passed']:
            error_msg = (f"Glycopeptide overlap too low: {overlap_pct_1 * 100:.1f}% ({name1}) and "
                         f"{overlap_pct_2 * 100:.1f}% ({name2}), required: {min_overlap_pct * 100:.1f}%")
            logger.error(error_msg)
            raise ValidationError(error_msg)

        logger.info(f"  ✓ PASS: Overlap sufficient (≥{min_overlap_pct * 100:.0f}%)")
        return True

    def validate_intensity_consistency(self,
                                       df1: pd.DataFrame,
                                       df2: pd.DataFrame,
                                       intensity_col: str,
                                       name1: str = "Dataset1",
                                       name2: str = "Dataset2",
                                       max_diff_pct: float = 0.01) -> bool:
        """
        Validate that intensity values are consistent between datasets

        Args:
            df1: First DataFrame
            df2: Second DataFrame
            intensity_col: Column name to compare (e.g., 'Cancer_Mean')
            name1: Name of first dataset
            name2: Name of second dataset
            max_diff_pct: Maximum allowed difference percentage

        Returns:
            True if validation passes

        Raises:
            ValidationError if differences exceed threshold
        """
        # Merge on Peptide + GlycanComposition
        merged = df1[['Peptide', 'GlycanComposition', intensity_col]].merge(
            df2[['Peptide', 'GlycanComposition', intensity_col]],
            on=['Peptide', 'GlycanComposition'],
            suffixes=('_1', '_2'),
            how='inner'
        )

        if len(merged) == 0:
            logger.warning(f"No common glycopeptides between {name1} and {name2}")
            return True

        val1 = merged[f'{intensity_col}_1']
        val2 = merged[f'{intensity_col}_2']

        # Calculate relative differences (handling zeros/NaN)
        # Use absolute difference for values near zero
        abs_diff = np.abs(val1 - val2)
        avg_val = (np.abs(val1) + np.abs(val2)) / 2
        rel_diff = np.where(avg_val > self.tolerance,
                            abs_diff / avg_val,
                            abs_diff)

        max_rel_diff = np.nanmax(rel_diff)
        mean_rel_diff = np.nanmean(rel_diff)
        n_exceed = (rel_diff > max_diff_pct).sum()

        result = {
            'check': 'intensity_consistency',
            'column': intensity_col,
            'dataset1': name1,
            'dataset2': name2,
            'n_compared': len(merged),
            'max_relative_diff': max_rel_diff,
            'mean_relative_diff': mean_rel_diff,
            'n_exceed_threshold': n_exceed,
            'passed': max_rel_diff <= max_diff_pct
        }

        self.validation_results.append(result)

        logger.info(f"Intensity consistency check: {intensity_col}")
        logger.info(f"  Compared: {len(merged)} glycopeptides")
        logger.info(f"  Max relative difference: {max_rel_diff * 100:.3f}%")
        logger.info(f"  Mean relative difference: {mean_rel_diff * 100:.3f}%")
        logger.info(f"  Exceeding threshold: {n_exceed}")

        if not result['passed']:
            # Find glycopeptides with largest differences
            worst_idx = np.argmax(rel_diff)
            worst_row = merged.iloc[worst_idx]

            error_msg = (f"Intensity values inconsistent for {intensity_col}: "
                         f"max difference {max_rel_diff * 100:.3f}%, threshold {max_diff_pct * 100:.3f}%. "
                         f"Worst case: {worst_row['Peptide']}|{worst_row['GlycanComposition']} "
                         f"({val1.iloc[worst_idx]:.2e} vs {val2.iloc[worst_idx]:.2e})")
            logger.error(error_msg)
            raise ValidationError(error_msg)

        logger.info("  ✓ PASS: Intensities consistent")
        return True

    def validate_detection_statistics(self,
                                      df1: pd.DataFrame,
                                      df2: pd.DataFrame,
                                      detection_col: str,
                                      name1: str = "Dataset1",
                                      name2: str = "Dataset2") -> bool:
        """
        Validate that detection statistics are identical

        Args:
            df1: First DataFrame
            df2: Second DataFrame
            detection_col: Column name (e.g., 'Cancer_Detection_Pct')
            name1: Name of first dataset
            name2: Name of second dataset

        Returns:
            True if validation passes

        Raises:
            ValidationError if detection stats differ
        """
        # Merge on Peptide + GlycanComposition
        merged = df1[['Peptide', 'GlycanComposition', detection_col]].merge(
            df2[['Peptide', 'GlycanComposition', detection_col]],
            on=['Peptide', 'GlycanComposition'],
            suffixes=('_1', '_2'),
            how='inner'
        )

        if len(merged) == 0:
            logger.warning(f"No common glycopeptides between {name1} and {name2}")
            return True

        val1 = merged[f'{detection_col}_1']
        val2 = merged[f'{detection_col}_2']

        # Detection percentages should be EXACTLY equal (integer counts / fixed sample size)
        abs_diff = np.abs(val1 - val2)
        max_diff = np.max(abs_diff)
        n_diff = (abs_diff > self.tolerance).sum()

        result = {
            'check': 'detection_statistics',
            'column': detection_col,
            'dataset1': name1,
            'dataset2': name2,
            'n_compared': len(merged),
            'max_difference': max_diff,
            'n_different': n_diff,
            'passed': max_diff <= self.tolerance
        }

        self.validation_results.append(result)

        logger.info(f"Detection statistics check: {detection_col}")
        logger.info(f"  Compared: {len(merged)} glycopeptides")
        logger.info(f"  Max difference: {max_diff:.6f}")
        logger.info(f"  Different values: {n_diff}")

        if not result['passed']:
            error_msg = (f"Detection statistics differ for {detection_col}: "
                         f"max difference {max_diff:.6f}, expected 0")
            logger.error(error_msg)
            raise ValidationError(error_msg)

        logger.info("  ✓ PASS: Detection statistics identical")
        return True

    def validate_all_visualizations(self,
                                    datasets: Dict[str, pd.DataFrame],
                                    required_columns: List[str] = None) -> bool:
        """
        Comprehensive validation across all visualization datasets

        Args:
            datasets: Dictionary mapping visualization names to DataFrames
            required_columns: List of columns that must be present and consistent

        Returns:
            True if all validations pass

        Raises:
            ValidationError if any validation fails
        """
        if required_columns is None:
            required_columns = ['Cancer_Mean', 'Normal_Mean',
                                'Cancer_Detection_Pct', 'Normal_Detection_Pct']

        logger.info("=" * 80)
        logger.info("COMPREHENSIVE DATA CONSISTENCY VALIDATION")
        logger.info("=" * 80)

        # Get list of dataset names
        names = list(datasets.keys())

        if len(names) < 2:
            logger.warning("Need at least 2 datasets for comparison")
            return True

        # Pairwise validation
        all_passed = True

        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                name1 = names[i]
                name2 = names[j]
                df1 = datasets[name1]
                df2 = datasets[name2]

                logger.info(f"\n--- Comparing: {name1} vs {name2} ---")

                try:
                    # Check 1: Glycopeptide overlap
                    self.validate_glycopeptide_overlap(df1, df2, name1, name2,
                                                       min_overlap_pct=0.80)

                    # Check 2: Intensity consistency
                    for col in required_columns:
                        if col in df1.columns and col in df2.columns:
                            if 'Mean' in col or 'Sum' in col:
                                # Allow 1% tolerance for means/sums (numerical precision)
                                self.validate_intensity_consistency(
                                    df1, df2, col, name1, name2, max_diff_pct=0.01
                                )
                            elif 'Detection_Pct' in col:
                                # Detection percentages should be exact
                                self.validate_detection_statistics(
                                    df1, df2, col, name1, name2
                                )

                except ValidationError as e:
                    logger.error(f"VALIDATION FAILED: {name1} vs {name2}")
                    logger.error(f"  {str(e)}")
                    all_passed = False

        logger.info("=" * 80)
        if all_passed:
            logger.info("✓ ALL VALIDATIONS PASSED")
        else:
            logger.error("✗ VALIDATION FAILURES DETECTED")
        logger.info("=" * 80)

        return all_passed

    def get_validation_summary(self) -> pd.DataFrame:
        """
        Get summary of all validation checks performed

        Returns:
            DataFrame with validation results
        """
        if not self.validation_results:
            return pd.DataFrame()

        return pd.DataFrame(self.validation_results)

    def save_validation_report(self, output_path: str):
        """
        Save validation report to file

        Args:
            output_path: Path to save validation report
        """
        summary = self.get_validation_summary()

        if summary.empty:
            logger.warning("No validation results to save")
            return

        summary.to_csv(output_path, index=False)
        logger.info(f"Validation report saved to {output_path}")

        # Also create human-readable text report
        txt_path = output_path.replace('.csv', '.txt')
        with open(txt_path, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("DATA CONSISTENCY VALIDATION REPORT\n")
            f.write("=" * 80 + "\n\n")

            for result in self.validation_results:
                f.write(f"Check: {result['check']}\n")
                for key, value in result.items():
                    if key != 'check':
                        f.write(f"  {key}: {value}\n")
                f.write("\n")

            # Summary
            total_checks = len(self.validation_results)
            passed_checks = sum(1 for r in self.validation_results if r.get('passed', False))

            f.write("=" * 80 + "\n")
            f.write(f"SUMMARY: {passed_checks}/{total_checks} checks passed\n")
            f.write("=" * 80 + "\n")

        logger.info(f"Text report saved to {txt_path}")


def quick_consistency_check(df1: pd.DataFrame,
                            df2: pd.DataFrame,
                            name1: str = "Dataset1",
                            name2: str = "Dataset2") -> bool:
    """
    Quick consistency check between two datasets

    Args:
        df1: First DataFrame
        df2: Second DataFrame
        name1: Name of first dataset
        name2: Name of second dataset

    Returns:
        True if datasets are consistent
    """
    validator = DataConsistencyValidator()

    try:
        # Basic checks
        validator.validate_glycopeptide_overlap(df1, df2, name1, name2, min_overlap_pct=0.95)

        # Intensity checks
        for col in ['Cancer_Mean', 'Normal_Mean']:
            if col in df1.columns and col in df2.columns:
                validator.validate_intensity_consistency(df1, df2, col, name1, name2)

        return True

    except ValidationError as e:
        logger.error(f"Consistency check failed: {str(e)}")
        return False
