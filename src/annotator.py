"""
Annotator Module for pGlyco Auto Combine
Handles Fucosylation and Sialylation annotation
"""

import pandas as pd
import re
from functools import lru_cache
from typing import Optional

from .constants import (
    DEFAULT_SIALYLATION_MARKER,
    DEFAULT_FUCOSYLATION_MARKER,
    MONOSACCHARIDE_H,
    MONOSACCHARIDE_N,
    MONOSACCHARIDE_A,
    MONOSACCHARIDE_F,
    MONOSACCHARIDE_G,
    GLYCAN_TYPE_HM,
    GLYCAN_TYPE_F,
    GLYCAN_TYPE_S,
    GLYCAN_TYPE_SF,
    GLYCAN_TYPE_CH,
    GLYCAN_TYPE_UNKNOWN,
    PRIMARY_TRUNCATED,
    PRIMARY_HIGH_MANNOSE,
    PRIMARY_COMPLEX_HYBRID,
    PRIMARY_OUTLIER,
    PRIMARY_UNKNOWN,
    SECONDARY_TRUNCATED,
    SECONDARY_HIGH_MANNOSE,
    SECONDARY_COMPLEX_HYBRID,
    SECONDARY_FUCOSYLATED,
    SECONDARY_SIALYLATED,
    SECONDARY_SIALOFUCOSYLATED,
    SECONDARY_OUTLIER,
    SECONDARY_UNKNOWN,
    LEGACY_NON,
    LEGACY_SIALYLATED,
    LEGACY_FUCOSYLATED,
    LEGACY_BOTH,
    HIGH_MANNOSE_MIN_H,
    HIGH_MANNOSE_EXACT_N,
    COMPLEX_HYBRID_MIN_N,
    GLYCAN_COMPOSITION_PATTERN
)
from .exceptions import InvalidGlycanCompositionError, AnnotationError
from .logger_config import get_logger

logger = get_logger(__name__)


class GlycanAnnotator:
    """Annotate glycopeptides based on glycan composition"""

    def __init__(self,
                 sialylation_marker: str = DEFAULT_SIALYLATION_MARKER,
                 fucosylation_marker: str = DEFAULT_FUCOSYLATION_MARKER):
        """
        Initialize GlycanAnnotator

        Args:
            sialylation_marker: Marker for sialylation (default: 'A' for NeuAc)
            fucosylation_marker: Marker for fucosylation (default: 'F')
        """
        self.sialylation_marker = sialylation_marker
        self.fucosylation_marker = fucosylation_marker

        # Create cached extraction method with current markers
        self._extract_cached = lru_cache(maxsize=1024)(self._extract_monosaccharide_impl)

    def _extract_monosaccharide_impl(self, glycan_composition: str, monosaccharide: str) -> int:
        """
        Implementation of monosaccharide extraction (cached)

        Args:
            glycan_composition: Glycan composition string
            monosaccharide: Monosaccharide symbol

        Returns:
            Count of the monosaccharide
        """
        # Pattern to match monosaccharide with count: e.g., A(2), F(1)
        pattern = rf'{monosaccharide}\((\d+)\)'
        match = re.search(pattern, glycan_composition)

        if match:
            return int(match.group(1))
        return 0

    def extract_monosaccharide_count(self, glycan_composition: str, monosaccharide: str) -> int:
        """
        Extract the count of a specific monosaccharide from glycan composition
        Results are cached for performance.

        Args:
            glycan_composition: Glycan composition string (e.g., "H(5)N(4)A(2)F(1)")
            monosaccharide: Monosaccharide symbol (e.g., "A", "F")

        Returns:
            Count of the monosaccharide (0 if not present)
        """
        if pd.isna(glycan_composition) or glycan_composition == "":
            return 0

        return self._extract_cached(str(glycan_composition), monosaccharide)

    def is_sialylated(self, glycan_composition: str) -> bool:
        """
        Check if glycan is sialylated (contains A for NeuAc)

        Args:
            glycan_composition: Glycan composition string

        Returns:
            True if sialylated, False otherwise
        """
        count = self.extract_monosaccharide_count(glycan_composition, self.sialylation_marker)
        return count > 0

    def is_fucosylated(self, glycan_composition: str) -> bool:
        """
        Check if glycan is fucosylated (contains F)

        Args:
            glycan_composition: Glycan composition string

        Returns:
            True if fucosylated, False otherwise
        """
        count = self.extract_monosaccharide_count(glycan_composition, self.fucosylation_marker)
        return count > 0

    def get_sialylation_count(self, glycan_composition: str) -> int:
        """Get the number of sialic acid residues"""
        return self.extract_monosaccharide_count(glycan_composition, self.sialylation_marker)

    def get_fucosylation_count(self, glycan_composition: str) -> int:
        """Get the number of fucose residues"""
        return self.extract_monosaccharide_count(glycan_composition, self.fucosylation_marker)

    def is_high_mannose(self, glycan_composition: str) -> bool:
        """
        Check if glycan is high mannose type
        Criteria: H≥5, N=2, no other modifications (A, F, G)

        Args:
            glycan_composition: Glycan composition string

        Returns:
            True if high mannose, False otherwise
        """
        if pd.isna(glycan_composition) or glycan_composition == "":
            return False

        # Must not have F, A, or G
        has_f = self.is_fucosylated(glycan_composition)
        has_a = self.is_sialylated(glycan_composition)
        has_g = self.extract_monosaccharide_count(glycan_composition, MONOSACCHARIDE_G) > 0

        if has_f or has_a or has_g:
            return False

        # Must have H≥5 and N=2 (using constants)
        h_count = self.extract_monosaccharide_count(glycan_composition, MONOSACCHARIDE_H)
        n_count = self.extract_monosaccharide_count(glycan_composition, MONOSACCHARIDE_N)

        return h_count >= HIGH_MANNOSE_MIN_H and n_count == HIGH_MANNOSE_EXACT_N

    def is_complex_hybrid(self, glycan_composition: str) -> bool:
        """
        Check if glycan is complex or hybrid type
        Criteria: No F, No A, and N >= 3

        Args:
            glycan_composition: Glycan composition string

        Returns:
            True if complex/hybrid, False otherwise
        """
        if pd.isna(glycan_composition) or glycan_composition == "":
            return False

        # Must not have F or A
        has_f = self.is_fucosylated(glycan_composition)
        has_a = self.is_sialylated(glycan_composition)

        if has_f or has_a:
            return False

        # Must have N >= 3 (using constant)
        n_count = self.extract_monosaccharide_count(glycan_composition, MONOSACCHARIDE_N)
        return n_count >= COMPLEX_HYBRID_MIN_N

    def get_glycan_type_category(self, glycan_composition: str) -> str:
        """
        Get glycan type category for visualization
        Categories: HM (High-mannose), F (Fucosylated), S (Sialylated),
                    SF (Sialofucosylated), C/H (Complex/Hybrid)

        Args:
            glycan_composition: Glycan composition string

        Returns:
            Glycan type category string
        """
        if pd.isna(glycan_composition) or glycan_composition == "":
            return GLYCAN_TYPE_UNKNOWN

        # Check for high mannose first (H≥5, N=2, no A/F/G)
        if self.is_high_mannose(glycan_composition):
            return GLYCAN_TYPE_HM

        # Check for sialylation and fucosylation
        has_a = self.is_sialylated(glycan_composition)
        has_f = self.is_fucosylated(glycan_composition)

        # Determine category based on modifications
        if has_a and has_f:
            return GLYCAN_TYPE_SF
        elif has_a:
            return GLYCAN_TYPE_S
        elif has_f:
            return GLYCAN_TYPE_F
        else:
            # Everything else is Complex/Hybrid
            return GLYCAN_TYPE_CH

    def get_primary_classification(self, glycan_composition: str) -> str:
        """
        Determine primary classification based on N (HexNAc) count

        Args:
            glycan_composition: Glycan composition string

        Returns:
            Primary classification: 'Truncated', 'High Mannose', 'ComplexHybrid', or 'Outlier'
        """
        if pd.isna(glycan_composition) or glycan_composition == "":
            return PRIMARY_UNKNOWN

        n_count = self.extract_monosaccharide_count(glycan_composition, MONOSACCHARIDE_N)
        has_f = self.is_fucosylated(glycan_composition)
        has_a = self.is_sialylated(glycan_composition)

        if n_count < 2:
            return PRIMARY_TRUNCATED
        elif n_count == 2:
            if has_f or has_a:
                return PRIMARY_OUTLIER
            else:
                return PRIMARY_HIGH_MANNOSE
        else:  # n_count >= 3
            return PRIMARY_COMPLEX_HYBRID

    def get_secondary_classification(self, glycan_composition: str) -> str:
        """
        Determine secondary classification based on primary type and F/A content

        Args:
            glycan_composition: Glycan composition string

        Returns:
            Secondary classification string
        """
        if pd.isna(glycan_composition) or glycan_composition == "":
            return SECONDARY_UNKNOWN

        primary = self.get_primary_classification(glycan_composition)
        has_f = self.is_fucosylated(glycan_composition)
        has_a = self.is_sialylated(glycan_composition)

        if primary == PRIMARY_TRUNCATED:
            return SECONDARY_TRUNCATED
        elif primary == PRIMARY_OUTLIER:
            return SECONDARY_OUTLIER
        elif primary == PRIMARY_HIGH_MANNOSE:
            return SECONDARY_HIGH_MANNOSE
        elif primary == PRIMARY_COMPLEX_HYBRID:
            if not has_f and not has_a:
                return SECONDARY_COMPLEX_HYBRID
            elif has_f and not has_a:
                return SECONDARY_FUCOSYLATED
            elif not has_f and has_a:
                return SECONDARY_SIALYLATED
            else:  # has_f and has_a
                return SECONDARY_SIALOFUCOSYLATED
        else:
            return SECONDARY_UNKNOWN

    def annotate_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Annotate DataFrame with primary and secondary glycan classifications

        Args:
            df: DataFrame with 'GlycanComposition' column

        Returns:
            DataFrame with added annotation columns:
            - PrimaryClassification: Truncated/High Mannose/ComplexHybrid/Outlier
            - SecondaryClassification: Detailed classification based on F/A content
            - (Legacy columns maintained for compatibility)
        """
        if 'GlycanComposition' not in df.columns:
            raise ValueError("DataFrame must contain 'GlycanComposition' column")

        df_annotated = df.copy()

        # Extract N count for classification
        df_annotated['N_count'] = df_annotated['GlycanComposition'].apply(
            lambda x: self.extract_monosaccharide_count(x, 'N')
        )

        # Check sialylation
        df_annotated['IsSialylated'] = df_annotated['GlycanComposition'].apply(self.is_sialylated)
        df_annotated['Sialylation'] = df_annotated['IsSialylated'].apply(
            lambda x: 'Sialylated' if x else 'Non-sialylated'
        )
        df_annotated['SialylationCount'] = df_annotated['GlycanComposition'].apply(
            self.get_sialylation_count
        )

        # Check fucosylation
        df_annotated['IsFucosylated'] = df_annotated['GlycanComposition'].apply(self.is_fucosylated)
        df_annotated['Fucosylation'] = df_annotated['IsFucosylated'].apply(
            lambda x: 'Fucosylated' if x else 'Non-fucosylated'
        )
        df_annotated['FucosylationCount'] = df_annotated['GlycanComposition'].apply(
            self.get_fucosylation_count
        )

        # Primary classification
        df_annotated['PrimaryClassification'] = df_annotated['GlycanComposition'].apply(
            self.get_primary_classification
        )

        # Secondary classification (final column)
        df_annotated['SecondaryClassification'] = df_annotated['GlycanComposition'].apply(
            self.get_secondary_classification
        )

        # Add GlycanTypeCategory for new visualization
        df_annotated['GlycanTypeCategory'] = df_annotated['GlycanComposition'].apply(
            self.get_glycan_type_category
        )

        # Legacy columns for backward compatibility (using constants)
        def determine_glycan_type(row):
            sia = row['IsSialylated']
            fuc = row['IsFucosylated']

            if sia and fuc:
                return LEGACY_BOTH
            elif sia:
                return LEGACY_SIALYLATED
            elif fuc:
                return LEGACY_FUCOSYLATED
            else:
                return LEGACY_NON

        df_annotated['GlycanType'] = df_annotated.apply(determine_glycan_type, axis=1)

        # Check High Mannose (legacy)
        df_annotated['IsHighMannose'] = df_annotated['GlycanComposition'].apply(self.is_high_mannose)
        df_annotated['HighMannose'] = df_annotated['IsHighMannose'].apply(
            lambda x: 'High mannose' if x else 'Not high mannose'
        )

        # Check Complex/Hybrid (legacy)
        df_annotated['IsComplexHybrid'] = df_annotated['GlycanComposition'].apply(self.is_complex_hybrid)
        df_annotated['ComplexHybrid'] = df_annotated['IsComplexHybrid'].apply(
            lambda x: 'C/H' if x else 'Not C/H'
        )

        # Log statistics
        logger.info(f"Annotation complete:")
        logger.info(f"\nPrimary Classification:")
        logger.info(df_annotated['PrimaryClassification'].value_counts())
        logger.info(f"\nSecondary Classification:")
        logger.info(df_annotated['SecondaryClassification'].value_counts())
        logger.info(f"\nLegacy Statistics:")
        logger.info(f"  - Sialylated: {df_annotated['IsSialylated'].sum()} ({df_annotated['IsSialylated'].sum()/len(df_annotated)*100:.1f}%)")
        logger.info(f"  - Fucosylated: {df_annotated['IsFucosylated'].sum()} ({df_annotated['IsFucosylated'].sum()/len(df_annotated)*100:.1f}%)")

        return df_annotated


if __name__ == "__main__":
    # Test the annotator
    test_data = pd.DataFrame({
        'Peptide': ['PEPTIDE1', 'PEPTIDE2', 'PEPTIDE3', 'PEPTIDE4'],
        'GlycanComposition': ['H(5)N(4)A(2)', 'H(5)N(4)A(2)F(1)', 'H(5)N(4)', 'H(5)N(4)F(1)']
    })

    annotator = GlycanAnnotator()
    annotated_data = annotator.annotate_dataframe(test_data)
    print(annotated_data[['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation', 'GlycanType']])
