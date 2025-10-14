"""
Helper Functions for Glycopeptide Comparison Heatmaps
Phase 2.1 Refactoring - Extract Common Functions

This module contains pure functions used by multiple heatmap methods
to eliminate code duplication.
"""

import re
import numpy as np
import pandas as pd
from typing import Tuple, List, Dict


def glycan_sort_key(glycan_comp: str) -> tuple:
    """
    Extract numbers from glycan composition for natural sorting

    Sorts glycan compositions by monosaccharide type and count.
    Order: H (Hexose), N (HexNAc), A (NeuAc), F (Fucose), G (GlcNAc)

    Args:
        glycan_comp: Glycan composition string (e.g., "H(5)N(4)A(2)F(1)")

    Returns:
        Tuple of (monosaccharide_order, count) pairs for sorting

    Example:
        >>> glycan_sort_key("H(5)N(4)A(2)")
        ((0, 5), (1, 4), (2, 2))
        >>> glycan_sort_key("H(6)N(5)")
        ((0, 6), (1, 5))
    """
    # Extract all monosaccharide types and their counts
    # Pattern: Letter(number)
    parts = re.findall(r'([A-Z]+)\((\d+)\)', glycan_comp)

    # Create a tuple of (monosaccharide, count as int) for sorting
    # Sort by: H, N, A, F, G order (typical glycan structure)
    monosaccharide_order = {'H': 0, 'N': 1, 'A': 2, 'F': 3, 'G': 4}

    sort_tuple = []
    for mono, count in parts:
        order = monosaccharide_order.get(mono, 99)
        sort_tuple.append((order, int(count)))

    return tuple(sort_tuple)


def get_symbol_info(
    row: pd.Series,
    linewidth_bold: float,
    linewidth_normal: float
) -> Tuple[str, str, float, str]:
    """
    Determine symbol markers and linewidth based on group presence

    Analyzes Cancer_Mean and Normal_Mean to determine:
    - Which symbols to plot (× for Cancer, + for Normal)
    - Linewidth (bold for qualitative differences, normal otherwise)
    - Group presence status

    Args:
        row: DataFrame row with 'Cancer_Mean' and 'Normal_Mean' columns
        linewidth_bold: Linewidth for qualitative differences (e.g., 5.0, 3.0)
        linewidth_normal: Linewidth when both groups present (e.g., 2.5)

    Returns:
        Tuple of (cancer_symbol, normal_symbol, linewidth, group_presence_label)
        - cancer_symbol: 'x' if Cancer present, None otherwise
        - normal_symbol: '+' if Normal present, None otherwise
        - linewidth: Bold or normal based on qualitative difference
        - group_presence_label: 'Both groups', 'Cancer only', 'Normal only', or 'Neither'

    Example:
        >>> row = pd.Series({'Cancer_Mean': 100, 'Normal_Mean': np.nan})
        >>> get_symbol_info(row, 5.0, 2.5)
        ('x', None, 5.0, 'Cancer only')
    """
    has_cancer = not np.isnan(row['Cancer_Mean']) and row['Cancer_Mean'] > 0
    has_normal = not np.isnan(row['Normal_Mean']) and row['Normal_Mean'] > 0

    # Qualitative difference: only one group present
    is_qualitatively_different = (has_cancer and not has_normal) or (has_normal and not has_cancer)

    # Set linewidth based on qualitative difference
    if is_qualitatively_different:
        linewidth = linewidth_bold  # Bold for qualitative difference
    else:
        linewidth = linewidth_normal  # Normal for both groups present

    # Determine symbols and group presence label
    if has_cancer and has_normal:
        return 'x', '+', linewidth, 'Both groups'
    elif has_cancer:
        return 'x', None, linewidth, 'Cancer only'
    elif has_normal:
        return None, '+', linewidth, 'Normal only'
    else:
        return None, None, 0, 'Neither'


def create_peptide_order(df: pd.DataFrame, vip_col: str = 'VIP_Score') -> List[str]:
    """
    Create peptide order sorted by VIP score (descending)

    Groups glycopeptides by peptide, takes max VIP score per peptide,
    and returns peptides sorted by VIP score in descending order.

    Args:
        df: DataFrame with 'Peptide' and VIP score columns
        vip_col: Name of VIP score column (default: 'VIP_Score')

    Returns:
        List of peptide names sorted by max VIP score (highest first)

    Example:
        >>> df = pd.DataFrame({
        ...     'Peptide': ['PEPTIDEA', 'PEPTIDEA', 'PEPTIDEB'],
        ...     'VIP_Score': [2.5, 3.0, 1.5]
        ... })
        >>> create_peptide_order(df)
        ['PEPTIDEA', 'PEPTIDEB']  # PEPTIDEA has max VIP 3.0 > PEPTIDEB 1.5
    """
    peptide_max_vip = df.groupby('Peptide')[vip_col].max().sort_values(ascending=False)
    return peptide_max_vip.index.tolist()


def create_glycan_order(
    df: pd.DataFrame,
    glycan_types: List[str],
    include_positions: bool = True
) -> Tuple[List[str], Dict[str, Dict[str, int]]]:
    """
    Create glycan order grouped by type with natural sorting

    Groups glycans by type (HM, F, S, SF, C/H), sorts each group
    naturally using glycan_sort_key(), and concatenates them.

    Args:
        df: DataFrame with 'GlycanComposition' and 'GlycanTypeCategory' columns
        glycan_types: Ordered list of glycan types (e.g., ['HM', 'F', 'S', 'SF', 'C/H'])
        include_positions: If True, track position of each glycan type (default: True)

    Returns:
        Tuple of (glycan_order, glycan_type_positions):
        - glycan_order: List of glycan compositions sorted by type and composition
        - glycan_type_positions: Dict mapping glycan_type -> {'start': int, 'end': int}
          Empty dict if include_positions=False

    Example:
        >>> df = pd.DataFrame({
        ...     'GlycanComposition': ['H(5)N(4)', 'H(6)N(5)', 'H(4)N(4)A(1)'],
        ...     'GlycanTypeCategory': ['HM', 'HM', 'S']
        ... })
        >>> glycan_order, positions = create_glycan_order(df, ['HM', 'S'])
        >>> glycan_order
        ['H(5)N(4)', 'H(6)N(5)', 'H(4)N(4)A(1)']
        >>> positions
        {'HM': {'start': 0, 'end': 2}, 'S': {'start': 2, 'end': 3}}
    """
    glycan_order = []
    glycan_type_positions = {}
    current_pos = 0

    for glycan_type in glycan_types:
        type_glycans = df[df['GlycanTypeCategory'] == glycan_type]['GlycanComposition'].unique()

        if len(type_glycans) > 0:
            # Sort glycan compositions by numeric values (natural sorting)
            type_glycans_sorted = sorted(type_glycans, key=glycan_sort_key)

            if include_positions:
                glycan_type_positions[glycan_type] = {
                    'start': current_pos,
                    'end': current_pos + len(type_glycans_sorted)
                }

            glycan_order.extend(type_glycans_sorted)
            current_pos += len(type_glycans_sorted)

    return glycan_order, glycan_type_positions


def create_index_mappings(
    peptide_order: List[str],
    glycan_order: List[str]
) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    Create index mapping dictionaries for peptides and glycans

    Maps peptide/glycan names to their integer positions in the heatmap.

    Args:
        peptide_order: Ordered list of peptide names (y-axis)
        glycan_order: Ordered list of glycan compositions (x-axis)

    Returns:
        Tuple of (peptide_to_idx, glycan_to_idx) dictionaries

    Example:
        >>> peptides = ['PEPTIDEA', 'PEPTIDEB']
        >>> glycans = ['H(5)N(4)', 'H(6)N(5)']
        >>> pep_idx, gly_idx = create_index_mappings(peptides, glycans)
        >>> pep_idx
        {'PEPTIDEA': 0, 'PEPTIDEB': 1}
        >>> gly_idx
        {'H(5)N(4)': 0, 'H(6)N(5)': 1}
    """
    peptide_to_idx = {p: i for i, p in enumerate(peptide_order)}
    glycan_to_idx = {g: i for i, g in enumerate(glycan_order)}

    return peptide_to_idx, glycan_to_idx
