"""
Base Analyzer Module
Provides shared functionality for all analyzers
"""

import pandas as pd
import numpy as np
from typing import Tuple
from sklearn.preprocessing import RobustScaler

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.constants import LOG_TRANSFORM_PSEUDOCOUNT
from src.exceptions import InsufficientDataError
from src.utils import (
    replace_empty_with_zero,
    get_all_sample_columns,
    log_transform as utils_log_transform,
    detect_log_transform  # Phase 1.3: Auto-detect log transformation
)
from src.logger_config import get_logger

logger = get_logger(__name__)


class BaseAnalyzer:
    """
    Base class for analyzers

    Provides common functionality:
    - Intensity matrix preparation
    - TIC normalization
    - Log transformation
    - Scaling
    """

    def __init__(self, log_transform: bool = True):
        """
        Initialize BaseAnalyzer

        Args:
            log_transform: Whether to apply log2 transformation
        """
        self.log_transform = log_transform
        self.scaler: RobustScaler = None

    def prepare_intensity_matrix(
        self,
        df: pd.DataFrame,
        min_samples: int = 2
    ) -> Tuple[pd.DataFrame, list, np.ndarray]:
        """
        Prepare intensity matrix for analysis with TIC normalization

        Pipeline: Extract → TIC Normalize → Log Transform

        Args:
            df: Annotated DataFrame with sample columns
            min_samples: Minimum number of samples required

        Returns:
            Tuple of (intensity_matrix_transposed, sample_names, feature_names)

        Raises:
            InsufficientDataError: If insufficient samples
        """
        # Get sample columns
        sample_cols = get_all_sample_columns(df)

        if len(sample_cols) < min_samples:
            raise InsufficientDataError("Analysis", min_samples, len(sample_cols))

        # Extract intensity matrix
        intensity_matrix = df[sample_cols].copy()

        # Replace empty strings with 0 and convert to numeric
        intensity_matrix = replace_empty_with_zero(intensity_matrix)

        # Transpose: rows = samples, columns = features
        intensity_matrix_t = intensity_matrix.T

        # Step 1: TIC (Total Ion Current) Normalization
        logger.info("Applying TIC normalization...")
        sample_sums = intensity_matrix_t.sum(axis=1)
        median_sum = sample_sums.median()

        # Avoid division by zero
        sample_sums_safe = sample_sums.replace(0, 1)
        intensity_matrix_t = intensity_matrix_t.div(sample_sums_safe, axis=0) * median_sum

        logger.info(f"  - Sample sum range: {sample_sums.min():.2e} - {sample_sums.max():.2e}")
        logger.info(f"  - Median sum: {median_sum:.2e}")

        # Step 2: Log2 transform if specified
        if self.log_transform:
            # Phase 1.3: Check if data is already log-transformed to prevent double transformation
            detection_result = detect_log_transform(intensity_matrix_t)

            if detection_result['is_log_transformed']:
                logger.warning(
                    f"Data appears to be already log-transformed "
                    f"(confidence: {detection_result['confidence']:.2f}, "
                    f"method: {detection_result['method']}). "
                    f"Skipping log transformation to prevent log(log(x)) error."
                )
                self.log_transform = False  # Disable to prevent future attempts
            else:
                logger.info(
                    f"Applying Log2 transformation... "
                    f"(data confirmed as raw, confidence: {1 - detection_result['confidence']:.2f})"
                )
                intensity_matrix_t = utils_log_transform(
                    intensity_matrix_t,
                    LOG_TRANSFORM_PSEUDOCOUNT
                )

        logger.info(f"Intensity matrix shape: {intensity_matrix_t.shape} (samples x features)")

        # Feature names
        feature_names = df[['Peptide', 'GlycanComposition']].values

        return intensity_matrix_t, sample_cols, feature_names

    def scale_features(self, data: np.ndarray) -> np.ndarray:
        """
        Scale features using RobustScaler

        Args:
            data: Data to scale (samples x features)

        Returns:
            Scaled data
        """
        logger.info("Applying RobustScaler (median and IQR-based scaling)...")
        self.scaler = RobustScaler()
        return self.scaler.fit_transform(data)
