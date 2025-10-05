"""
PCA Analyzer Module
Handles Principal Component Analysis
"""

import pandas as pd
import numpy as np
from typing import Dict
from sklearn.decomposition import PCA

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.constants import DEFAULT_PCA_COMPONENTS
from src.utils import get_sample_group
from src.logger_config import get_logger
from .base_analyzer import BaseAnalyzer

logger = get_logger(__name__)


class PCAAnalyzer(BaseAnalyzer):
    """
    Principal Component Analysis

    Performs PCA on glycoproteomics data with proper normalization pipeline:
    TIC Normalization → Log2 Transform → RobustScaler → PCA
    """

    def __init__(
        self,
        n_components: int = DEFAULT_PCA_COMPONENTS,
        log_transform: bool = True
    ):
        """
        Initialize PCAAnalyzer

        Args:
            n_components: Number of PCA components
            log_transform: Whether to apply log2 transformation
        """
        super().__init__(log_transform)
        self.n_components = n_components
        self.pca: PCA = None

    def perform_pca(self, df: pd.DataFrame) -> Dict:
        """
        Perform PCA analysis

        Pipeline: TIC Normalization → Log2 Transform → RobustScaler → PCA

        Args:
            df: Annotated DataFrame

        Returns:
            Dictionary with PCA results:
            - pca_df: DataFrame with PC coordinates and groups
            - explained_variance: Explained variance ratios
            - loadings: Feature loadings on PCs
            - pca_model: Fitted PCA object
        """
        # Prepare data
        intensity_matrix, sample_names, feature_names = self.prepare_intensity_matrix(df)

        # Scale features
        intensity_scaled = self.scale_features(intensity_matrix)

        # Perform PCA
        logger.info("Performing PCA...")
        self.pca = PCA(n_components=self.n_components)
        pca_coords = self.pca.fit_transform(intensity_scaled)

        # Create results DataFrame
        pca_df = pd.DataFrame(
            pca_coords,
            columns=[f'PC{i+1}' for i in range(self.n_components)],
            index=sample_names
        )

        # Add sample group information
        pca_df['Group'] = pca_df.index.map(get_sample_group)

        # Calculate explained variance
        explained_variance = self.pca.explained_variance_ratio_

        logger.info(f"PCA completed:")
        for i, var in enumerate(explained_variance):
            logger.info(f"  PC{i+1}: {var*100:.2f}% variance explained")

        # Get loadings (feature contributions)
        loadings = pd.DataFrame(
            self.pca.components_.T,
            columns=[f'PC{i+1}' for i in range(self.n_components)],
            index=range(len(feature_names))
        )

        return {
            'pca_df': pca_df,
            'explained_variance': explained_variance,
            'loadings': loadings,
            'pca_model': self.pca
        }
