"""
Data Loader Module for pGlyco Auto Combine
Handles CSV integration from multiple input files
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re
from typing import List, Dict
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DataLoader:
    """Load and integrate glycoproteomics data from multiple CSV files"""

    def __init__(self, dataset_dir: str, required_columns: List[str]):
        """
        Initialize DataLoader

        Args:
            dataset_dir: Path to directory containing CSV files
            required_columns: List of required column names
        """
        self.dataset_dir = Path(dataset_dir)
        self.required_columns = required_columns

    def get_csv_files(self) -> List[Path]:
        """Get all CSV files in dataset directory"""
        csv_files = sorted(self.dataset_dir.glob("*.csv"))
        logger.info(f"Found {len(csv_files)} CSV files in {self.dataset_dir}")
        return csv_files

    def extract_sample_id(self, filename: str) -> str:
        """
        Extract sample ID from filename (e.g., C_01.csv -> C1)

        Args:
            filename: Name of the CSV file

        Returns:
            Sample ID (e.g., 'C1', 'N10')
        """
        match = re.search(r'([CN])_(\d+)', filename)
        if match:
            prefix = match.group(1)
            number = int(match.group(2))
            return f"{prefix}{number}"
        return filename

    def load_single_file(self, file_path: Path) -> pd.DataFrame:
        """
        Load a single CSV file and extract required columns

        Args:
            file_path: Path to CSV file

        Returns:
            DataFrame with required columns
        """
        try:
            df = pd.read_csv(file_path)

            # Check if required columns exist
            missing_cols = [col for col in self.required_columns if col not in df.columns]
            if missing_cols:
                logger.warning(f"Missing columns in {file_path.name}: {missing_cols}")
                return None

            # Extract required columns
            df_subset = df[self.required_columns].copy()

            # Add sample ID
            sample_id = self.extract_sample_id(file_path.name)
            df_subset['SampleID'] = sample_id

            logger.info(f"Loaded {len(df_subset)} rows from {file_path.name} (Sample: {sample_id})")
            return df_subset

        except Exception as e:
            logger.error(f"Error loading {file_path.name}: {str(e)}")
            return None

    def integrate_data(self, qc_filters: Dict = None) -> pd.DataFrame:
        """
        Integrate data from all CSV files into a single wide-format DataFrame

        Args:
            qc_filters: Dictionary of quality control filters

        Returns:
            Integrated DataFrame with structure:
            Peptide | GlycanComposition | C1 | C2 | ... | N1 | N2 | ...
        """
        csv_files = self.get_csv_files()

        if not csv_files:
            raise ValueError(f"No CSV files found in {self.dataset_dir}")

        # Load all files
        all_data = []
        for csv_file in csv_files:
            df = self.load_single_file(csv_file)
            if df is not None:
                all_data.append(df)

        if not all_data:
            raise ValueError("No valid data loaded from CSV files")

        # Concatenate all data
        combined_df = pd.concat(all_data, ignore_index=True)
        logger.info(f"Combined data shape: {combined_df.shape}")

        # Apply QC filters
        if qc_filters:
            min_area = qc_filters.get('min_isotope_area', 0)
            before_count = len(combined_df)
            combined_df = combined_df[combined_df['IsotopeArea'] >= min_area]
            logger.info(f"QC filter: Removed {before_count - len(combined_df)} rows with IsotopeArea < {min_area}")

        # Pivot to wide format
        # Group by Peptide, GlycanComposition, and Proteins, with SampleID as columns
        integrated_df = combined_df.pivot_table(
            index=['Peptide', 'GlycanComposition', 'Proteins'],
            columns='SampleID',
            values='IsotopeArea',
            aggfunc='sum'  # Sum if there are duplicates
        ).reset_index()

        # Sort columns: C1, C2, ..., C24, N1, N2, ..., N24
        sample_cols = [col for col in integrated_df.columns if col not in ['Peptide', 'GlycanComposition', 'Proteins']]

        # Separate C and N samples
        c_samples = sorted([col for col in sample_cols if col.startswith('C')],
                          key=lambda x: int(x[1:]))
        n_samples = sorted([col for col in sample_cols if col.startswith('N')],
                          key=lambda x: int(x[1:]))

        # Reorder columns: Peptide, GlycanComposition, sample columns, Proteins (at the end)
        ordered_cols = ['Peptide', 'GlycanComposition'] + c_samples + n_samples + ['Proteins']
        integrated_df = integrated_df[ordered_cols]

        # Fill NaN with empty string for missing values
        integrated_df = integrated_df.fillna('')

        logger.info(f"Integrated data shape: {integrated_df.shape}")
        logger.info(f"Sample columns: {c_samples + n_samples}")

        return integrated_df

    def save_integrated_data(self, df: pd.DataFrame, output_path: str):
        """
        Save integrated data to CSV

        Args:
            df: Integrated DataFrame
            output_path: Path to output CSV file
        """
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        df.to_csv(output_file, index=False)
        logger.info(f"Saved integrated data to {output_file}")


if __name__ == "__main__":
    # Test the data loader
    loader = DataLoader(
        dataset_dir="Dataset",
        required_columns=["Peptide", "GlycanComposition", "IsotopeArea"]
    )

    integrated_data = loader.integrate_data(qc_filters={'min_isotope_area': 0})
    print(integrated_data.head())
    print(f"\nShape: {integrated_data.shape}")
