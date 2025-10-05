"""
PLS-DA Analyzer Module
Handles PLS-DA analysis and VIP score calculation
"""

import pandas as pd
import numpy as np
from typing import Dict
from sklearn.cross_decomposition import PLSRegression

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.constants import DEFAULT_PLSDA_COMPONENTS, GROUP_CANCER
from src.utils import get_sample_group
from src.logger_config import get_logger
from .base_analyzer import BaseAnalyzer

logger = get_logger(__name__)


class PLSDAAnalyzer(BaseAnalyzer):
    """
    PLS-DA (Partial Least Squares Discriminant Analysis) Analyzer

    Performs supervised analysis for Cancer vs Normal classification
    and calculates VIP (Variable Importance in Projection) scores
    """

    def __init__(
        self,
        n_components: int = DEFAULT_PLSDA_COMPONENTS,
        log_transform: bool = True
    ):
        """
        Initialize PLSDAAnalyzer

        Args:
            n_components: Number of PLS components
            log_transform: Whether to apply log2 transformation
        """
        super().__init__(log_transform)
        self.n_components = n_components
        self.plsda: PLSRegression = None

    def perform_plsda(self, df: pd.DataFrame, n_components: int = None) -> Dict:
        """
        Perform PLS-DA analysis for Cancer vs Normal classification

        Args:
            df: Annotated DataFrame (pre-filtered)
            n_components: Number of components (uses default if None)

        Returns:
            Dictionary with PLS-DA results and VIP scores:
            - vip_scores: DataFrame with VIP scores
            - plsda_model: Fitted PLS model
            - X_scores: PLS scores
            - y_labels: Binary labels
            - sample_names: Sample identifiers
        """
        if n_components is None:
            n_components = self.n_components

        logger.info(f"Performing PLS-DA on {len(df)} pre-filtered glycopeptides...")

        # Prepare data
        intensity_matrix, sample_names, feature_info = self.prepare_intensity_matrix(df)

        # Create binary labels (0 = Normal, 1 = Cancer)
        y = np.array([
            1 if get_sample_group(name) == GROUP_CANCER else 0
            for name in sample_names
        ])

        # Scale data
        X_scaled = self.scale_features(intensity_matrix)

        # Fit PLS-DA model
        logger.info("Fitting PLS-DA model...")
        self.plsda = PLSRegression(n_components=n_components)
        self.plsda.fit(X_scaled, y)

        # Get PLS scores
        X_scores = self.plsda.x_scores_

        # Calculate VIP scores
        vip_scores = self._calculate_vip_scores(X_scaled, y)

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

    def _calculate_vip_scores(self, X: np.ndarray, y: np.ndarray) -> np.ndarray:
        """
        Calculate VIP (Variable Importance in Projection) scores

        Args:
            X: Scaled feature matrix
            y: Target labels

        Returns:
            Array of VIP scores for each feature
        """
        if self.plsda is None:
            raise ValueError("PLS-DA model not fitted. Call perform_plsda() first.")

        t = self.plsda.x_scores_  # PLS scores
        w = self.plsda.x_weights_  # PLS weights
        q = self.plsda.y_loadings_  # Y loadings

        p, h = w.shape  # p = features, h = components

        # Calculate sum of squares explained by each component
        s = np.diag(t.T @ t @ q.T @ q).reshape(h, -1)
        total_s = np.sum(s)

        # Calculate VIP scores
        vip = np.zeros((p,))
        for i in range(p):
            weight = np.array([(w[i, j] ** 2) * s[j] for j in range(h)])
            vip[i] = np.sqrt(p * np.sum(weight) / total_s)

        return vip

    def validate_vip_with_bootstrap(
        self,
        X: np.ndarray,
        y: np.ndarray,
        n_iterations: int = 100
    ) -> Dict:
        """
        Validate VIP scores using bootstrap resampling

        Assesses stability of VIP scores across resampled datasets.
        Features with stable VIP scores are more reliable biomarkers.

        Args:
            X: Scaled feature matrix
            y: Binary labels
            n_iterations: Number of bootstrap iterations

        Returns:
            Dictionary with bootstrap validation results:
            - mean_vip: Mean VIP scores
            - std_vip: Standard deviation of VIP scores
            - cv_vip: Coefficient of variation
        """
        n_samples, n_features = X.shape
        vip_bootstrap = np.zeros((n_iterations, n_features))

        logger.info(f"Bootstrap VIP validation: {n_iterations} iterations...")

        for i in range(n_iterations):
            # Resample with replacement
            indices = np.random.choice(n_samples, n_samples, replace=True)
            X_boot = X[indices]
            y_boot = y[indices]

            # Fit PLS model
            pls_boot = PLSRegression(n_components=self.n_components)
            pls_boot.fit(X_boot, y_boot)

            # Calculate VIP scores
            t = pls_boot.x_scores_
            w = pls_boot.x_weights_
            q = pls_boot.y_loadings_

            p, h = w.shape
            s = np.diag(t.T @ t @ q.T @ q).reshape(h, -1)
            total_s = np.sum(s)

            for j in range(p):
                weight = np.array([(w[j, k] ** 2) * s[k] for k in range(h)])
                vip_bootstrap[i, j] = np.sqrt(p * np.sum(weight) / total_s)

        # Calculate statistics
        mean_vip = np.mean(vip_bootstrap, axis=0)
        std_vip = np.std(vip_bootstrap, axis=0)
        cv_vip = std_vip / (mean_vip + 1e-10)  # Avoid division by zero

        logger.info("Bootstrap validation completed")

        return {
            'mean_vip': mean_vip,
            'std_vip': std_vip,
            'cv_vip': cv_vip,
            'bootstrap_vips': vip_bootstrap
        }

    def get_top_vip_by_category(
        self,
        df: pd.DataFrame,
        vip_df: pd.DataFrame,
        category_col: str = 'GlycanType',
        top_n: int = 10
    ) -> pd.DataFrame:
        """
        Get top VIP scores grouped by category

        Args:
            df: Annotated DataFrame
            vip_df: VIP scores DataFrame
            category_col: Column to group by
            top_n: Number of top features per category

        Returns:
            DataFrame with top VIP scores per category
        """
        # Merge VIP scores with annotations
        merged = vip_df.merge(
            df[['Peptide', 'GlycanComposition', category_col]],
            on=['Peptide', 'GlycanComposition'],
            how='left'
        )

        # Get top N per category
        top_per_category = (
            merged.groupby(category_col)
            .apply(lambda x: x.nlargest(top_n, 'VIP_Score'))
            .reset_index(drop=True)
        )

        return top_per_category
