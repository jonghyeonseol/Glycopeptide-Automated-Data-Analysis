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

    def annotate_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Annotate DataFrame with Sialylation and Fucosylation columns

        Args:
            df: DataFrame with 'GlycanComposition' column

        Returns:
            DataFrame with added annotation columns:
            - Sialylation: 'Sialylated' or 'Non-sialylated'
            - Fucosylation: 'Fucosylated' or 'Non-fucosylated'
            - SialylationCount: Number of sialic acid residues
            - FucosylationCount: Number of fucose residues
            - GlycanType: Combined annotation (Non/Sialylated/Fucosylated/Both)
        """
        if 'GlycanComposition' not in df.columns:
            raise ValueError("DataFrame must contain 'GlycanComposition' column")

        df_annotated = df.copy()

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

        # Combined glycan type
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

        # Log statistics
        logger.info(f"Annotation complete:")
        logger.info(f"  - Sialylated: {df_annotated['IsSialylated'].sum()} ({df_annotated['IsSialylated'].sum()/len(df_annotated)*100:.1f}%)")
        logger.info(f"  - Fucosylated: {df_annotated['IsFucosylated'].sum()} ({df_annotated['IsFucosylated'].sum()/len(df_annotated)*100:.1f}%)")
        logger.info(f"\nGlycan Type Distribution:")
        logger.info(df_annotated['GlycanType'].value_counts())

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
