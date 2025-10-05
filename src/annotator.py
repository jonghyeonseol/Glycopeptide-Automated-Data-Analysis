"""
Annotator Module for pGlyco Auto Combine
Handles Fucosylation and Sialylation annotation
"""

import pandas as pd
import re
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class GlycanAnnotator:
    """Annotate glycopeptides based on glycan composition"""

    def __init__(self, sialylation_marker: str = "A", fucosylation_marker: str = "F"):
        """
        Initialize GlycanAnnotator

        Args:
            sialylation_marker: Marker for sialylation (default: 'A' for NeuAc)
            fucosylation_marker: Marker for fucosylation (default: 'F')
        """
        self.sialylation_marker = sialylation_marker
        self.fucosylation_marker = fucosylation_marker

    def extract_monosaccharide_count(self, glycan_composition: str, monosaccharide: str) -> int:
        """
        Extract the count of a specific monosaccharide from glycan composition

        Args:
            glycan_composition: Glycan composition string (e.g., "H(5)N(4)A(2)F(1)")
            monosaccharide: Monosaccharide symbol (e.g., "A", "F")

        Returns:
            Count of the monosaccharide (0 if not present)
        """
        if pd.isna(glycan_composition) or glycan_composition == "":
            return 0

        # Pattern to match monosaccharide with count: e.g., A(2), F(1)
        pattern = rf'{monosaccharide}\((\d+)\)'
        match = re.search(pattern, str(glycan_composition))

        if match:
            return int(match.group(1))
        return 0

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
        has_g = self.extract_monosaccharide_count(glycan_composition, 'G') > 0

        if has_f or has_a or has_g:
            return False

        # Must have H≥5 and N=2
        h_count = self.extract_monosaccharide_count(glycan_composition, 'H')
        n_count = self.extract_monosaccharide_count(glycan_composition, 'N')

        return h_count >= 5 and n_count == 2

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

        # Must have N >= 3
        n_count = self.extract_monosaccharide_count(glycan_composition, 'N')
        return n_count >= 3

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
            return 'Unknown'

        # Check for high mannose first (H≥5, N=2, no A/F/G)
        if self.is_high_mannose(glycan_composition):
            return 'HM'

        # Check for sialylation and fucosylation
        has_a = self.is_sialylated(glycan_composition)
        has_f = self.is_fucosylated(glycan_composition)

        # Determine category based on modifications
        if has_a and has_f:
            return 'SF'
        elif has_a:
            return 'S'
        elif has_f:
            return 'F'
        else:
            # Everything else is Complex/Hybrid
            return 'C/H'

    def get_primary_classification(self, glycan_composition: str) -> str:
        """
        Determine primary classification based on N (HexNAc) count

        Args:
            glycan_composition: Glycan composition string

        Returns:
            Primary classification: 'Truncated', 'High Mannose', 'ComplexHybrid', or 'Outlier'
        """
        if pd.isna(glycan_composition) or glycan_composition == "":
            return 'Unknown'

        n_count = self.extract_monosaccharide_count(glycan_composition, 'N')
        has_f = self.is_fucosylated(glycan_composition)
        has_a = self.is_sialylated(glycan_composition)

        if n_count < 2:
            return 'Truncated'
        elif n_count == 2:
            if has_f or has_a:
                return 'Outlier'
            else:
                return 'High Mannose'
        else:  # n_count >= 3
            return 'ComplexHybrid'

    def get_secondary_classification(self, glycan_composition: str) -> str:
        """
        Determine secondary classification based on primary type and F/A content

        Args:
            glycan_composition: Glycan composition string

        Returns:
            Secondary classification string
        """
        if pd.isna(glycan_composition) or glycan_composition == "":
            return 'Unknown'

        primary = self.get_primary_classification(glycan_composition)
        has_f = self.is_fucosylated(glycan_composition)
        has_a = self.is_sialylated(glycan_composition)

        if primary == 'Truncated':
            return 'Truncated'
        elif primary == 'Outlier':
            return 'Outlier'
        elif primary == 'High Mannose':
            return 'High Mannose'
        elif primary == 'ComplexHybrid':
            if not has_f and not has_a:
                return 'Complex/Hybrid'
            elif has_f and not has_a:
                return 'Fucosylated'
            elif not has_f and has_a:
                return 'Sialylated'
            else:  # has_f and has_a
                return 'Sialofucosylated'
        else:
            return 'Unknown'

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

        # Legacy columns for backward compatibility
        def determine_glycan_type(row):
            sia = row['IsSialylated']
            fuc = row['IsFucosylated']

            if sia and fuc:
                return 'Both'
            elif sia:
                return 'Sialylated'
            elif fuc:
                return 'Fucosylated'
            else:
                return 'Non'

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
