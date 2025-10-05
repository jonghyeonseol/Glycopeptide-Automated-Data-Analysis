"""
Data Loader Module for pGlyco Auto Combine
Handles CSV integration from multiple input files
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional

from .constants import (
    REQUIRED_INPUT_COLUMNS,
    CSV_FILE_PATTERN,
    CANCER_PREFIX,
    NORMAL_PREFIX
)
from .exceptions import (
    NoDataFilesError,
    MissingColumnError,
    EmptyDataError,
    DataLoadError
)
from .utils import extract_sample_id, ensure_directory
from .logger_config import get_logger

logger = get_logger(__name__)


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
        """
        Get all CSV files in dataset directory

        Returns:
            List of CSV file paths

        Raises:
            NoDataFilesError: If no CSV files are found
        """
        csv_files = sorted(self.dataset_dir.glob(CSV_FILE_PATTERN))
        logger.info(f"Found {len(csv_files)} CSV files in {self.dataset_dir}")

        if not csv_files:
            raise NoDataFilesError(str(self.dataset_dir))

        return csv_files

    def load_single_file(self, file_path: Path) -> Optional[pd.DataFrame]:
        """
        Load a single CSV file and extract required columns

        Args:
            file_path: Path to CSV file

        Returns:
            DataFrame with required columns, or None if loading fails

        Raises:
            MissingColumnError: If required columns are missing
            DataLoadError: If file cannot be loaded
        """
        try:
            df = pd.read_csv(file_path)

            # Check if required columns exist
            missing_cols = [col for col in self.required_columns if col not in df.columns]
            if missing_cols:
                raise MissingColumnError(file_path.name, missing_cols)

            # Extract required columns
            df_subset = df[self.required_columns].copy()

            # Add sample ID
            sample_id = extract_sample_id(file_path.name)
            df_subset['SampleID'] = sample_id

            logger.info(f"Loaded {len(df_subset)} rows from {file_path.name} (Sample: {sample_id})")
            return df_subset

        except MissingColumnError:
            # Re-raise MissingColumnError
            raise
        except Exception as e:
            raise DataLoadError(f"Failed to load {file_path.name}: {str(e)}")

    def integrate_data(self, qc_filters: Optional[Dict] = None) -> pd.DataFrame:
        """
        Integrate data from all CSV files into a single wide-format DataFrame

        Args:
            qc_filters: Dictionary of quality control filters

        Returns:
            Integrated DataFrame with structure:
            Peptide | GlycanComposition | C1 | C2 | ... | N1 | N2 | ...

        Raises:
            NoDataFilesError: If no CSV files are found
            EmptyDataError: If no valid data is loaded
        """
        csv_files = self.get_csv_files()

        # Load all files
        all_data = []
        for csv_file in csv_files:
            try:
                df = self.load_single_file(csv_file)
                if df is not None:
                    all_data.append(df)
            except (MissingColumnError, DataLoadError) as e:
                logger.warning(f"Skipping {csv_file.name}: {str(e)}")
                continue

        if not all_data:
            raise EmptyDataError("No valid data loaded from CSV files")

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
        c_samples = sorted([col for col in sample_cols if col.startswith(CANCER_PREFIX)],
                          key=lambda x: int(x[1:]))
        n_samples = sorted([col for col in sample_cols if col.startswith(NORMAL_PREFIX)],
                          key=lambda x: int(x[1:]))

        # Reorder columns: Peptide, GlycanComposition, sample columns, Proteins (at the end)
        ordered_cols = ['Peptide', 'GlycanComposition'] + c_samples + n_samples + ['Proteins']
        integrated_df = integrated_df[ordered_cols]

        # Fill NaN with empty string for missing values
        integrated_df = integrated_df.fillna('')

        logger.info(f"Integrated data shape: {integrated_df.shape}")
        logger.info(f"Sample columns: {c_samples + n_samples}")

        return integrated_df

    def save_integrated_data(self, df: pd.DataFrame, output_path: str) -> Path:
        """
        Save integrated data to CSV

        Args:
            df: Integrated DataFrame
            output_path: Path to output CSV file

        Returns:
            Path to the saved file

        Raises:
            FileOperationError: If file cannot be saved
        """
        output_file = Path(output_path)
        ensure_directory(output_file.parent)

        df.to_csv(output_file, index=False)
        logger.info(f"Saved integrated data to {output_file}")

        return output_file


if __name__ == "__main__":
    # Test the data loader
    loader = DataLoader(
        dataset_dir="Dataset",
        required_columns=["Peptide", "GlycanComposition", "IsotopeArea"]
    )

    integrated_data = loader.integrate_data(qc_filters={'min_isotope_area': 0})
    print(integrated_data.head())
    print(f"\nShape: {integrated_data.shape}")
