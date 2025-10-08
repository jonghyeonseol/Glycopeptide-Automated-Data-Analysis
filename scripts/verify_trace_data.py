#!/usr/bin/env python3
"""
Verification Script for Glycopeptide Comparison Heatmap Trace Data
This script performs comprehensive checks to ensure data integrity and correctness.
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path


def verify_heatmap_data():
    """Comprehensive verification of trace data"""

    print("="*80)
    print("Verification of Glycopeptide Comparison Heatmap Trace Data")
    print("="*80)

    # Check if files exist
    summary_path = Path('Results/Trace/glycopeptide_comparison_heatmap_summary.csv')
    full_data_path = Path('Results/Trace/glycopeptide_comparison_heatmap_data.csv')

    if not summary_path.exists():
        print(f"\n❌ ERROR: Summary file not found at {summary_path}")
        print("   Run the visualization first: python3 test_comparison_heatmap.py")
        return False

    if not full_data_path.exists():
        print(f"\n❌ ERROR: Full data file not found at {full_data_path}")
        print("   Run the visualization first: python3 test_comparison_heatmap.py")
        return False

    # Load data
    print("\nLoading data files...")
    try:
        summary = pd.read_csv(summary_path)
        full_data = pd.read_csv(full_data_path)
    except Exception as e:
        print(f"\n❌ ERROR: Failed to load data files: {e}")
        return False

    print("✓ Data loaded successfully")
    print(f"  - Summary rows: {len(summary)}")
    print(f"  - Full data rows: {len(full_data)}")
    print(f"  - Summary columns: {len(summary.columns)}")
    print(f"  - Full data columns: {len(full_data.columns)}")

    all_checks_passed = True

    # ==================== Check 1: Verify Cancer Mean Calculations ====================
    print(f"\n{'='*80}")
    print("CHECK 1: Verifying Cancer_Mean calculations")
    print(f"{'='*80}")

    cancer_cols = [f'C{i}' for i in range(1, 25)]
    available_cancer_cols = [col for col in cancer_cols if col in full_data.columns]

    print(f"Found {len(available_cancer_cols)} cancer sample columns")

    mismatches = 0
    max_error = 0

    for idx, row in full_data.iterrows():
        cancer_values = row[available_cancer_cols].fillna(0).astype(float)
        manual_mean = cancer_values.mean()
        saved_mean = row['Cancer_Mean']

        error = abs(manual_mean - saved_mean)
        max_error = max(max_error, error)

        if error > 0.01:
            mismatches += 1
            if mismatches <= 3:  # Show first 3 mismatches
                print(f"  Row {idx}: Manual={manual_mean:.2f}, Saved={saved_mean:.2f}, Error={error:.2f}")

    print("\nResults:")
    print(f"  - Total rows checked: {len(full_data)}")
    print(f"  - Mismatches (error > 0.01): {mismatches}")
    print(f"  - Maximum error: {max_error:.6f}")

    if mismatches == 0:
        print("  ✓ PASS: All Cancer_Mean values verified")
    else:
        print(f"  ✗ FAIL: {mismatches} mismatches found")
        all_checks_passed = False

    # ==================== Check 2: Verify Normal Mean Calculations ====================
    print(f"\n{'='*80}")
    print("CHECK 2: Verifying Normal_Mean calculations")
    print(f"{'='*80}")

    normal_cols = [f'N{i}' for i in list(range(1, 19)) + list(range(20, 25))]  # N19 missing
    available_normal_cols = [col for col in normal_cols if col in full_data.columns]

    print(f"Found {len(available_normal_cols)} normal sample columns")

    mismatches = 0
    max_error = 0

    for idx, row in full_data.iterrows():
        normal_values = row[available_normal_cols].fillna(0).astype(float)
        manual_mean = normal_values.mean()
        saved_mean = row['Normal_Mean']

        error = abs(manual_mean - saved_mean)
        max_error = max(max_error, error)

        if error > 0.01:
            mismatches += 1
            if mismatches <= 3:
                print(f"  Row {idx}: Manual={manual_mean:.2f}, Saved={saved_mean:.2f}, Error={error:.2f}")

    print("\nResults:")
    print(f"  - Total rows checked: {len(full_data)}")
    print(f"  - Mismatches (error > 0.01): {mismatches}")
    print(f"  - Maximum error: {max_error:.6f}")

    if mismatches == 0:
        print("  ✓ PASS: All Normal_Mean values verified")
    else:
        print(f"  ✗ FAIL: {mismatches} mismatches found")
        all_checks_passed = False

    # ==================== Check 3: Verify VIP Score Sorting ====================
    print(f"\n{'='*80}")
    print("CHECK 3: Verifying VIP score sorting (Y-axis)")
    print(f"{'='*80}")

    vip_by_pos = summary.groupby('Plot_Y_Position')['PeptideVIP'].mean().sort_index()

    print("VIP scores by Y position (should be descending):")
    for pos, vip in vip_by_pos.head(10).items():
        print(f"  Y={pos}: VIP={vip:.4f}")

    is_sorted = all(vip_by_pos.iloc[i] >= vip_by_pos.iloc[i+1]
                   for i in range(len(vip_by_pos)-1))

    print("\nResults:")
    print(f"  - Total Y positions: {len(vip_by_pos)}")
    print(f"  - Properly sorted (descending): {is_sorted}")

    if is_sorted:
        print("  ✓ PASS: VIP scores properly sorted")
    else:
        print("  ✗ FAIL: VIP scores not properly sorted")
        all_checks_passed = False

    # ==================== Check 4: Verify Glycan Type Grouping ====================
    print(f"\n{'='*80}")
    print("CHECK 4: Verifying glycan type grouping (X-axis)")
    print(f"{'='*80}")

    glycan_grouping = summary.sort_values('Plot_X_Position')[['Plot_X_Position', 'GlycanTypeCategory']].drop_duplicates()

    print("Glycan type distribution:")
    for gtype in ['HM', 'F', 'S', 'SF', 'C/H']:
        positions = summary[summary['GlycanTypeCategory'] == gtype]['Plot_X_Position'].unique()
        positions = sorted(positions)

        if len(positions) > 0:
            is_contiguous = len(positions) == 1 or all(positions[i+1] == positions[i]+1 for i in range(len(positions)-1))
            status = "✓ Contiguous" if is_contiguous else "✗ Non-contiguous"
            print(f"  {gtype}: {len(positions)} glycans at positions {positions[0]}-{positions[-1]} {status}")

            if not is_contiguous:
                all_checks_passed = False

    # ==================== Check 5: Verify Fold Change Calculations ====================
    print(f"\n{'='*80}")
    print("CHECK 5: Verifying fold change calculations")
    print(f"{'='*80}")

    fc_mismatches = 0
    log2fc_mismatches = 0

    for idx, row in summary.iterrows():
        if row['Normal_Mean'] > 0:
            # Verify Fold_Change
            manual_fc = row['Cancer_Mean'] / row['Normal_Mean']
            if not np.isinf(row['Fold_Change']) and abs(manual_fc - row['Fold_Change']) > 0.01:
                fc_mismatches += 1

            # Verify Log2_Fold_Change
            if row['Cancer_Mean'] > 0:
                manual_log2fc = np.log2(manual_fc)
                if not np.isnan(row['Log2_Fold_Change']) and abs(manual_log2fc - row['Log2_Fold_Change']) > 0.01:
                    log2fc_mismatches += 1

    print("\nResults:")
    print(f"  - Fold_Change mismatches: {fc_mismatches}")
    print(f"  - Log2_Fold_Change mismatches: {log2fc_mismatches}")

    if fc_mismatches == 0 and log2fc_mismatches == 0:
        print("  ✓ PASS: All fold change calculations verified")
    else:
        print("  ✗ FAIL: Found calculation errors")
        all_checks_passed = False

    # ==================== Check 6: Verify Plot Flags ====================
    print(f"\n{'='*80}")
    print("CHECK 6: Verifying plot flags (Cancer_Dot_Plotted, Normal_Dot_Plotted)")
    print(f"{'='*80}")

    cancer_flag_errors = sum((summary['Cancer_Mean'] > 0) != summary['Cancer_Dot_Plotted'])
    normal_flag_errors = sum((summary['Normal_Mean'] > 0) != summary['Normal_Dot_Plotted'])

    print("\nResults:")
    print(f"  - Cancer flag errors: {cancer_flag_errors}")
    print(f"  - Normal flag errors: {normal_flag_errors}")

    if cancer_flag_errors == 0 and normal_flag_errors == 0:
        print("  ✓ PASS: All plot flags verified")
    else:
        print("  ✗ FAIL: Found flag errors")
        all_checks_passed = False

    # ==================== Check 7: Verify Alpha Values ====================
    print(f"\n{'='*80}")
    print("CHECK 7: Verifying alpha (transparency) values")
    print(f"{'='*80}")

    cancer_alpha_errors = sum((summary['Cancer_Alpha'] < 0) | (summary['Cancer_Alpha'] > 1))
    normal_alpha_errors = sum((summary['Normal_Alpha'] < 0) | (summary['Normal_Alpha'] > 1))

    print("\nResults:")
    print(f"  - Cancer alpha values out of range [0, 1]: {cancer_alpha_errors}")
    print(f"  - Normal alpha values out of range [0, 1]: {normal_alpha_errors}")

    if cancer_alpha_errors == 0 and normal_alpha_errors == 0:
        print("  ✓ PASS: All alpha values in valid range")
    else:
        print("  ✗ FAIL: Found invalid alpha values")
        all_checks_passed = False

    # ==================== Check 8: Verify Sample Counts ====================
    print(f"\n{'='*80}")
    print("CHECK 8: Verifying sample counts")
    print(f"{'='*80}")

    max_cancer_samples = len(available_cancer_cols)
    max_normal_samples = len(available_normal_cols)

    cancer_count_errors = sum(summary['Cancer_SampleCount'] > max_cancer_samples)
    normal_count_errors = sum(summary['Normal_SampleCount'] > max_normal_samples)

    print("\nResults:")
    print(f"  - Max cancer samples: {max_cancer_samples}")
    print(f"  - Max normal samples: {max_normal_samples}")
    print(f"  - Cancer sample count errors: {cancer_count_errors}")
    print(f"  - Normal sample count errors: {normal_count_errors}")

    if cancer_count_errors == 0 and normal_count_errors == 0:
        print("  ✓ PASS: All sample counts valid")
    else:
        print("  ✗ FAIL: Found invalid sample counts")
        all_checks_passed = False

    # ==================== Final Summary ====================
    print(f"\n{'='*80}")
    print("VERIFICATION SUMMARY")
    print(f"{'='*80}")

    if all_checks_passed:
        print("\n✓✓✓ ALL CHECKS PASSED ✓✓✓")
        print("\nThe trace data is RELIABLE and can be used for manual verification.")
        print("Every value in the visualization can be traced back to the original data.")
    else:
        print("\n✗✗✗ SOME CHECKS FAILED ✗✗✗")
        print("\nPlease review the detailed output above for specific issues.")
        print("Contact support if problems persist.")

    print(f"\n{'='*80}\n")

    return all_checks_passed


if __name__ == "__main__":
    success = verify_heatmap_data()
    sys.exit(0 if success else 1)
