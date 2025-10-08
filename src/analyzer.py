"""
Analyzer Module for pGlyco Auto Combine
Handles PCA and statistical analysis
"""

import pandas as pd
import numpy as np
from typing import Dict, Tuple, Optional
from sklearn.decomposition import PCA
from sklearn.preprocessing import RobustScaler
from sklearn.cross_decomposition import PLSRegression

from .constants import (
    DEFAULT_PCA_COMPONENTS,
    DEFAULT_LOG_TRANSFORM,
    DEFAULT_PLSDA_COMPONENTS,
    LOG_TRANSFORM_PSEUDOCOUNT,
    GROUP_CANCER
)
from .exceptions import (
    InsufficientDataError
)
from .utils import (
    replace_empty_with_zero,
    to_numeric_safe,
    get_sample_columns,
    get_all_sample_columns,
    get_sample_group,
    log_transform as utils_log_transform,
    calculate_fold_change
)
from .logger_config import get_logger

logger = get_logger(__name__)


class GlycanAnalyzer:
    """Perform statistical analysis on glycoproteomics data"""

    def __init__(self,
                 n_components: int = DEFAULT_PCA_COMPONENTS,
                 log_transform: bool = DEFAULT_LOG_TRANSFORM):
        """
        Initialize GlycanAnalyzer

        Args:
            n_components: Number of PCA components (default: 2)
            log_transform: Whether to log-transform intensity values (default: True)
        """
        self.n_components = n_components
        self.log_transform = log_transform
        self.pca: Optional[PCA] = None
        self.scaler: Optional[RobustScaler] = None
        self.plsda: Optional[PLSRegression] = None
        self.vip_scores: Optional[pd.DataFrame] = None

    def prepare_intensity_matrix(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, list, np.ndarray]:
        """
        Prepare intensity matrix for analysis with TIC normalization

        Args:
            df: Annotated DataFrame with sample columns

        Returns:
            Tuple of (intensity_matrix, sample_names, feature_names)

        Raises:
            InsufficientDataError: If insufficient samples are present
        """
        # Get sample columns using utility function (eliminates duplication)
        sample_cols = get_all_sample_columns(df)

        if len(sample_cols) < 2:
            raise InsufficientDataError("PCA", 2, len(sample_cols))

        # Extract intensity matrix
        intensity_matrix = df[sample_cols].copy()

        # Replace empty strings with 0 and convert to numeric
        intensity_matrix = replace_empty_with_zero(intensity_matrix)

        # Transpose: rows = samples, columns = features
        intensity_matrix_t = intensity_matrix.T

        # Step 1: TIC (Total Ion Current) Normalization - normalize each sample
        logger.info("Applying TIC normalization...")
        sample_sums = intensity_matrix_t.sum(axis=1)
        median_sum = sample_sums.median()

        # Avoid division by zero
        sample_sums_safe = sample_sums.replace(0, 1)
        intensity_matrix_t = intensity_matrix_t.div(sample_sums_safe, axis=0) * median_sum

        logger.info(f"  - Sample sum range before: {sample_sums.min():.2e} - {sample_sums.max():.2e}")
        logger.info(f"  - Median sum (target): {median_sum:.2e}")

        # Step 2: Log2 transform if specified
        if self.log_transform:
            logger.info("Applying Log2 transformation...")
            # Use utility function with constant pseudocount
            intensity_matrix_t = utils_log_transform(intensity_matrix_t, LOG_TRANSFORM_PSEUDOCOUNT)

        logger.info(f"Intensity matrix shape: {intensity_matrix_t.shape} (samples x features)")

        return intensity_matrix_t, sample_cols, df[['Peptide', 'GlycanComposition']].values

    def perform_pca(self, df: pd.DataFrame) -> dict:
        """
        Perform PCA analysis with RobustScaler
        Pipeline: TIC Normalization → Log2 Transform → RobustScaler → PCA

        Args:
            df: Annotated DataFrame

        Returns:
            Dictionary containing PCA results
        """
        intensity_matrix, sample_names, feature_names = self.prepare_intensity_matrix(df)

        # Step 3: RobustScaler - scale features using median and IQR (robust to outliers)
        logger.info("Applying RobustScaler (median and IQR-based scaling)...")
        self.scaler = RobustScaler()
        intensity_scaled = self.scaler.fit_transform(intensity_matrix)

        # Step 4: Perform PCA
        logger.info("Performing PCA...")
        self.pca = PCA(n_components=self.n_components)
        pca_coords = self.pca.fit_transform(intensity_scaled)

        # Create results DataFrame
        pca_df = pd.DataFrame(
            pca_coords,
            columns =[f'PC{i + 1}' for i in range(self.n_components)],
            index=sample_names
        )

        # Add sample group information (C or N) using utility function
        pca_df['Group'] = pca_df.index.map(get_sample_group)

        # Calculate explained variance
        explained_variance = self.pca.explained_variance_ratio_

        logger.info("PCA completed:")
        logger.info(f"  - PC1 explained variance: {explained_variance[0] * 100:.2f}%")
        logger.info(f"  - PC2 explained variance: {explained_variance[1] * 100:.2f}%")
        logger.info(f"  - Total explained variance: {explained_variance.sum() * 100:.2f}%")

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
        # Get sample columns using utility function (eliminates duplication #2)
        c_samples, n_samples = get_sample_columns(df)

        # Calculate statistics for each glycan type
        stats_list = []

        for glycan_type in df[group_col].unique():
            subset = df[df[group_col] == glycan_type]

            # Convert to numeric and replace empty with 0
            c_values = replace_empty_with_zero(subset[c_samples])
            n_values = replace_empty_with_zero(subset[n_samples])

            # Calculate mean intensity for each sample
            c_mean = c_values.sum(axis=0).mean()
            n_mean = n_values.sum(axis=0).mean()

            # Calculate total intensity
            c_total = c_values.sum().sum()
            n_total = n_values.sum().sum()

            # Count glycopeptides
            count = len(subset)

            # Use utility function for fold change calculation
            stats_list.append({
                'GlycanType': glycan_type,
                'Count': count,
                'Cancer_Mean': c_mean,
                'Normal_Mean': n_mean,
                'Cancer_Total': c_total,
                'Normal_Total': n_total,
                'Fold_Change': calculate_fold_change(c_mean, n_mean, log_scale=False),
                'Log2_Fold_Change': calculate_fold_change(c_mean, n_mean, log_scale=True)
            })

        stats_df = pd.DataFrame(stats_list)

        logger.info("\nStatistics by Glycan Type:")
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
        # Get sample columns using utility function (eliminates duplication #3)
        sample_cols = get_all_sample_columns(df)

        # Melt to long format
        df_long = df.melt(
            id_vars=['GlycanType'],
            value_vars=sample_cols,
            var_name='Sample',
            value_name='Intensity'
        )

        # Convert intensity to numeric
        df_long['Intensity'] = to_numeric_safe(df_long['Intensity'].replace('', np.nan))

        # Add group information using utility function
        df_long['Group'] = df_long['Sample'].apply(get_sample_group)

        # Log transform
        if self.log_transform:
            df_long['Intensity'] = np.log2(df_long['Intensity'] + 1)

        # Remove zero intensities for better visualization
        df_long = df_long[df_long['Intensity'] > 0]

        # Aggregate by Sample and GlycanType (mean of all glycopeptides within each category)
        df_long = df_long.groupby(['Sample', 'GlycanType', 'Group'], as_index=False)['Intensity'].mean()

        logger.info(f"Boxplot data prepared: {len(df_long)} aggregated observations")

        return df_long

    def prepare_boxplot_data_extended(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Prepare data for extended boxplot visualization with 5 categories

        Args:
            df: Annotated DataFrame

        Returns:
            Long-format DataFrame for plotting with extended categories
        """
        # Get sample columns using utility function (eliminates duplication #4)
        sample_cols = get_all_sample_columns(df)

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
        df_long['Intensity'] = to_numeric_safe(df_long['Intensity'].replace('', np.nan))

        # Add group information using utility function
        df_long['Group'] = df_long['Sample'].apply(get_sample_group)

        # Log transform
        if self.log_transform:
            df_long['Intensity'] = np.log2(df_long['Intensity'] + 1)

        # Remove zero intensities for better visualization
        df_long = df_long[df_long['Intensity'] > 0]

        # Aggregate by Sample and ExtendedCategory (mean of all glycopeptides within each category)
        df_long = df_long.groupby(['Sample', 'ExtendedCategory', 'Group'], as_index=False)['Intensity'].mean()

        logger.info(f"Extended boxplot data prepared: {len(df_long)} aggregated observations")

        return df_long

    def perform_plsda(self, df: pd.DataFrame, n_components: int = DEFAULT_PLSDA_COMPONENTS) -> Dict:
        """
        Perform PLS-DA (Partial Least Squares Discriminant Analysis)
        for Cancer vs Normal classification

        UPDATED: Now accepts PRE-FILTERED data from DataPipeline.
        Filtering is applied ONCE in main.py via DataPipeline to ensure consistency.

        Args:
            df: Annotated DataFrame (ALREADY FILTERED by DataPipeline)
            n_components: Number of PLS components (default: 2)

        Returns:
            Dictionary containing PLS-DA results and VIP scores
        """
        # NO FILTERING HERE - Data is pre-filtered by DataPipeline
        logger.info(f"Performing PLS-DA analysis on {len(df)} pre-filtered glycopeptides...")

        intensity_matrix, sample_names, feature_info = self.prepare_intensity_matrix(df)

        # Create binary labels (0 = Normal, 1 = Cancer) using utility function
        y = np.array([1 if get_sample_group(name) == GROUP_CANCER else 0 for name in sample_names])

        # Scale data using RobustScaler
        logger.info("Performing PLS-DA analysis...")
        self.scaler = RobustScaler()
        X_scaled = self.scaler.fit_transform(intensity_matrix)

        # Fit PLS-DA model
        self.plsda = PLSRegression(n_components=n_components)
        self.plsda.fit(X_scaled, y)

        # Get PLS scores
        X_scores = self.plsda.x_scores_

        # Calculate VIP scores
        vip_scores = self._calculate_vip_scores(X_scaled, y, self.plsda)

        # Create VIP results DataFrame
        vip_df = pd.DataFrame({
            'Peptide': feature_info[:, 0],
            'GlycanComposition': feature_info[:, 1],
            'VIP_Score': vip_scores
        })

        # Sort by VIP score
        vip_df = vip_df.sort_values('VIP_Score', ascending=False).reset_index(drop=True)

        logger.info(f"PLS-DA completed with {n_components} components")
        logger.info(f"Top 10 VIP scores:\n{vip_df.head(10)}")

        return {
            'vip_scores': vip_df,
            'plsda_model': self.plsda,
            'X_scores': X_scores,
            'y_labels': y,
            'sample_names': sample_names
        }

    def _calculate_vip_scores(self, X, y, pls_model):
        """
        Calculate VIP (Variable Importance in Projection) scores

        Args:
            X: Scaled feature matrix
            y: Target labels
            pls_model: Fitted PLSRegression model

        Returns:
            Array of VIP scores for each feature
        """
        t = pls_model.x_scores_  # PLS scores
        w = pls_model.x_weights_  # PLS weights
        q = pls_model.y_loadings_  # Y loadings

        p, h = w.shape  # p = number of features, h = number of components

        # Calculate sum of squares explained by each component
        s = np.diag(t.T @ t @ q.T @ q).reshape(h, -1)
        total_s = np.sum(s)

        # Calculate VIP scores
        vip = np.zeros((p,))
        for i in range(p):
            weight = np.array([(w[i, j] ** 2) * s[j] for j in range(h)])
            vip[i] = np.sqrt(p * np.sum(weight) / total_s)

        return vip

    def validate_vip_with_bootstrap(self, X: np.ndarray, y: np.ndarray,
                                    n_iterations: int = 100, n_components: int = 2) -> Dict:
        """
        Validate VIP scores using bootstrap resampling

        SCIENTIFIC VALIDITY: Bootstrap validation assesses stability of VIP scores.
        Features with stable VIP scores across resampled datasets are more reliable biomarkers.

        Args:
            X: Scaled feature matrix (samples × features)
            y: Binary labels (Cancer=1, Normal=0)
            n_iterations: Number of bootstrap iterations (default: 100)
            n_components: Number of PLS components

        Returns:
            Dictionary with:
            - vip_mean: Mean VIP score across iterations
            - vip_std: Standard deviation of VIP scores
            - vip_ci_lower: Lower 95% confidence interval
            - vip_ci_upper: Upper 95% confidence interval
            - stability: Boolean array (True if CI doesn't include 1.0)
        """
        n_samples, n_features = X.shape
        vip_matrix = np.zeros((n_iterations, n_features))

        logger.info(f"Running VIP bootstrap validation ({n_iterations} iterations)...")

        for i in range(n_iterations):
            # Resample with replacement
            indices = np.random.choice(n_samples, n_samples, replace=True)
            X_boot = X[indices]
            y_boot = y[indices]

            # Fit PLS-DA model
            pls_boot = PLSRegression(n_components=n_components)
            pls_boot.fit(X_boot, y_boot)

            # Calculate VIP scores for this iteration
            vip_boot = self._calculate_vip_scores(X_boot, y_boot, pls_boot)
            vip_matrix[i, :] = vip_boot

        # Calculate statistics across iterations
        vip_mean = np.mean(vip_matrix, axis=0)
        vip_std = np.std(vip_matrix, axis=0)
        vip_ci_lower = np.percentile(vip_matrix, 2.5, axis=0)
        vip_ci_upper = np.percentile(vip_matrix, 97.5, axis=0)

        # Stability: VIP > 1 if lower CI > 1.0 (consistently important)
        stability = vip_ci_lower > 1.0

        logger.info("Bootstrap validation complete:")
        logger.info(f"  - {np.sum(stability)} features have stable VIP > 1.0 (CI lower > 1.0)")
        logger.info(f"  - Mean VIP score: {np.mean(vip_mean):.3f} ± {np.mean(vip_std):.3f}")

        return {
            'vip_mean': vip_mean,
            'vip_std': vip_std,
            'vip_ci_lower': vip_ci_lower,
            'vip_ci_upper': vip_ci_upper,
            'stability': stability
        }

    def get_top_vip_by_glycopeptide(self, df: pd.DataFrame, plsda_results: dict, top_n: int = 10) -> pd.DataFrame:
        """
        Get top N glycopeptides by VIP score

        Args:
            df: Annotated DataFrame
            plsda_results: Results from perform_plsda()
            top_n: Number of top features to return

        Returns:
            DataFrame with top N glycopeptides
        """
        return plsda_results['vip_scores'].head(top_n)

    def get_top_vip_by_glycan_type(self, df: pd.DataFrame, plsda_results: dict, top_n: int = 10,
                                   classification_col: str = 'SecondaryClassification') -> pd.DataFrame:
        """
        Get top N glycan types by aggregated VIP score

        Args:
            df: Annotated DataFrame
            plsda_results: Results from perform_plsda()
            top_n: Number of top glycan types to return
            classification_col: Column to group by (PrimaryClassification or SecondaryClassification)

        Returns:
            DataFrame with top N glycan types
        """
        vip_df = plsda_results['vip_scores'].copy()

        # Merge with glycan type information
        df_info = df[['Peptide', 'GlycanComposition', classification_col]].copy()
        vip_with_type = vip_df.merge(df_info, on=['Peptide', 'GlycanComposition'], how='left')

        # Aggregate VIP scores by glycan type (mean)
        glycan_type_vip = vip_with_type.groupby(classification_col)['VIP_Score'].agg(['mean', 'count']).reset_index()
        glycan_type_vip.columns = ['GlycanType', 'Mean_VIP_Score', 'Count']
        glycan_type_vip = glycan_type_vip.sort_values('Mean_VIP_Score', ascending=False).reset_index(drop=True)

        logger.info(f"Top {top_n} Glycan Types by VIP Score:\n{glycan_type_vip.head(top_n)}")

        return glycan_type_vip.head(top_n)

    def get_top_vip_by_peptide(self, df: pd.DataFrame, plsda_results: dict, top_n: int = 10) -> pd.DataFrame:
        """
        Get top N peptides by aggregated VIP score

        Args:
            df: Annotated DataFrame
            plsda_results: Results from perform_plsda()
            top_n: Number of top peptides to return

        Returns:
            DataFrame with top N peptides
        """
        vip_df = plsda_results['vip_scores'].copy()

        # Aggregate VIP scores by peptide (mean across all glycan modifications)
        peptide_vip = vip_df.groupby('Peptide')['VIP_Score'].agg(['mean', 'count']).reset_index()
        peptide_vip.columns = ['Peptide', 'Mean_VIP_Score', 'Count']
        peptide_vip = peptide_vip.sort_values('Mean_VIP_Score', ascending=False).reset_index(drop=True)

        logger.info(f"Top {top_n} Peptides by VIP Score:\n{peptide_vip.head(top_n)}")

        return peptide_vip.head(top_n)


if __name__ == "__main__":
    # Test with sample data
    pass
