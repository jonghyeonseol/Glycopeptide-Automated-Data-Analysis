"""
Analyzer Module for pGlyco Auto Combine
Handles PCA and statistical analysis
"""

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy import stats
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class GlycanAnalyzer:
    """Perform statistical analysis on glycoproteomics data"""

    def __init__(self, n_components: int = 2, log_transform: bool = True):
        """
        Initialize GlycanAnalyzer

        Args:
            n_components: Number of PCA components
            log_transform: Whether to log-transform intensity values
        """
        self.n_components = n_components
        self.log_transform = log_transform
        self.pca = None
        self.scaler = None

    def prepare_intensity_matrix(self, df: pd.DataFrame) -> tuple:
        """
        Prepare intensity matrix for analysis

        Args:
            df: Annotated DataFrame with sample columns

        Returns:
            Tuple of (intensity_matrix, sample_names, feature_names)
        """
        # Identify sample columns (C1, C2, ..., N1, N2, ...)
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid']

        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Extract intensity matrix
        intensity_matrix = df[sample_cols].copy()

        # Replace empty strings with 0
        intensity_matrix = intensity_matrix.replace('', 0)

        # Convert to numeric
        intensity_matrix = intensity_matrix.apply(pd.to_numeric, errors='coerce').fillna(0)

        # Transpose: rows = samples, columns = features
        intensity_matrix_t = intensity_matrix.T

        # Log transform if specified
        if self.log_transform:
            # Add pseudocount to avoid log(0)
            intensity_matrix_t = np.log2(intensity_matrix_t + 1)

        logger.info(f"Intensity matrix shape: {intensity_matrix_t.shape} (samples x features)")

        return intensity_matrix_t, sample_cols, df[['Peptide', 'GlycanComposition']].values

    def perform_pca(self, df: pd.DataFrame) -> dict:
        """
        Perform PCA analysis

        Args:
            df: Annotated DataFrame

        Returns:
            Dictionary containing PCA results
        """
        intensity_matrix, sample_names, feature_names = self.prepare_intensity_matrix(df)

        # Standardize the data
        self.scaler = StandardScaler()
        intensity_scaled = self.scaler.fit_transform(intensity_matrix)

        # Perform PCA
        self.pca = PCA(n_components=self.n_components)
        pca_coords = self.pca.fit_transform(intensity_scaled)

        # Create results DataFrame
        pca_df = pd.DataFrame(
            pca_coords,
            columns=[f'PC{i+1}' for i in range(self.n_components)],
            index=sample_names
        )

        # Add sample group information (C or N)
        pca_df['Group'] = pca_df.index.map(lambda x: 'Cancer' if x.startswith('C') else 'Normal')

        # Calculate explained variance
        explained_variance = self.pca.explained_variance_ratio_

        logger.info(f"PCA completed:")
        logger.info(f"  - PC1 explained variance: {explained_variance[0]*100:.2f}%")
        logger.info(f"  - PC2 explained variance: {explained_variance[1]*100:.2f}%")
        logger.info(f"  - Total explained variance: {explained_variance.sum()*100:.2f}%")

        return {
            'pca_df': pca_df,
            'explained_variance': explained_variance,
            'loadings': self.pca.components_,
            'feature_names': feature_names
        }

    def calculate_statistics_by_glycan_type(self, df: pd.DataFrame, group_col: str = 'GlycanType') -> pd.DataFrame:
        """
        Calculate statistics for each glycan type across samples

        Args:
            df: Annotated DataFrame
            group_col: Column to group by (default: 'GlycanType')

        Returns:
            DataFrame with statistics for each glycan type
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid']

        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Separate C and N samples
        c_samples = [col for col in sample_cols if col.startswith('C')]
        n_samples = [col for col in sample_cols if col.startswith('N')]

        # Calculate statistics for each glycan type
        stats_list = []

        for glycan_type in df[group_col].unique():
            subset = df[df[group_col] == glycan_type]

            # Convert to numeric and replace empty with 0
            c_values = subset[c_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)
            n_values = subset[n_samples].replace('', 0).apply(pd.to_numeric, errors='coerce').fillna(0)

            # Calculate mean intensity for each sample
            c_mean = c_values.sum(axis=0).mean()
            n_mean = n_values.sum(axis=0).mean()

            # Calculate total intensity
            c_total = c_values.sum().sum()
            n_total = n_values.sum().sum()

            # Count glycopeptides
            count = len(subset)

            stats_list.append({
                'GlycanType': glycan_type,
                'Count': count,
                'Cancer_Mean': c_mean,
                'Normal_Mean': n_mean,
                'Cancer_Total': c_total,
                'Normal_Total': n_total,
                'Fold_Change': c_mean / n_mean if n_mean > 0 else np.inf,
                'Log2_Fold_Change': np.log2(c_mean / n_mean) if (c_mean > 0 and n_mean > 0) else np.nan
            })

        stats_df = pd.DataFrame(stats_list)

        logger.info(f"\nStatistics by Glycan Type:")
        logger.info(stats_df.to_string(index=False))

        return stats_df

    def prepare_boxplot_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Prepare data for boxplot visualization

        Args:
            df: Annotated DataFrame

        Returns:
            Long-format DataFrame for plotting
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid']

        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Melt to long format
        df_long = df.melt(
            id_vars=['GlycanType'],
            value_vars=sample_cols,
            var_name='Sample',
            value_name='Intensity'
        )

        # Convert intensity to numeric
        df_long['Intensity'] = pd.to_numeric(df_long['Intensity'].replace('', 0), errors='coerce').fillna(0)

        # Add group information
        df_long['Group'] = df_long['Sample'].apply(lambda x: 'Cancer' if x.startswith('C') else 'Normal')

        # Log transform
        if self.log_transform:
            df_long['Intensity'] = np.log2(df_long['Intensity'] + 1)

        # Remove zero intensities for better visualization
        df_long = df_long[df_long['Intensity'] > 0]

        logger.info(f"Boxplot data prepared: {len(df_long)} observations")

        return df_long

    def prepare_boxplot_data_extended(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Prepare data for extended boxplot visualization with 5 categories

        Args:
            df: Annotated DataFrame

        Returns:
            Long-format DataFrame for plotting with extended categories
        """
        # Identify sample columns
        metadata_cols = ['Peptide', 'GlycanComposition', 'Sialylation', 'Fucosylation',
                        'IsSialylated', 'IsFucosylated', 'SialylationCount',
                        'FucosylationCount', 'GlycanType', 'HighMannose', 'ComplexHybrid',
                        'IsHighMannose', 'IsComplexHybrid']

        sample_cols = [col for col in df.columns if col not in metadata_cols]

        # Create extended category based on all annotations
        def determine_extended_category(row):
            if row['IsHighMannose']:
                return 'HM'
            elif row['IsComplexHybrid']:
                return 'C/H'
            elif row['IsSialylated'] and row['IsFucosylated']:
                return 'Sialofucosylated'
            elif row['IsSialylated']:
                return 'Sialylated'
            elif row['IsFucosylated']:
                return 'Fucosylated'
            else:
                return 'Other'

        df_copy = df.copy()
        df_copy['ExtendedCategory'] = df_copy.apply(determine_extended_category, axis=1)

        # Melt to long format
        df_long = df_copy.melt(
            id_vars=['ExtendedCategory'],
            value_vars=sample_cols,
            var_name='Sample',
            value_name='Intensity'
        )

        # Convert intensity to numeric
        df_long['Intensity'] = pd.to_numeric(df_long['Intensity'].replace('', 0), errors='coerce').fillna(0)

        # Add group information
        df_long['Group'] = df_long['Sample'].apply(lambda x: 'Cancer' if x.startswith('C') else 'Normal')

        # Log transform
        if self.log_transform:
            df_long['Intensity'] = np.log2(df_long['Intensity'] + 1)

        # Remove zero intensities for better visualization
        df_long = df_long[df_long['Intensity'] > 0]

        logger.info(f"Extended boxplot data prepared: {len(df_long)} observations")

        return df_long


if __name__ == "__main__":
    # Test with sample data
    pass
