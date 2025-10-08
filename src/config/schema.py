"""
Configuration Schema
Defines and validates configuration structure
"""

from typing import Dict, Any
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.exceptions import ConfigurationError


class ConfigSchema:
    """
    Configuration schema validator

    Validates structure and types of configuration
    """

    def __init__(self):
        self.required_keys = [
            'paths',
            'processing',
            'annotation',
            'analysis',
            'visualization'
        ]

        self.required_paths = [
            'dataset_dir',
            'results_dir',
            'output_file'
        ]

        self.required_processing = [
            'required_columns',
            'qc_filters'
        ]

    def validate(self, config: Dict[str, Any]) -> None:
        """
        Validate configuration

        Args:
            config: Configuration dictionary

        Raises:
            ConfigurationError: If validation fails
        """
        # Check top-level keys
        for key in self.required_keys:
            if key not in config:
                raise ConfigurationError(f"Missing required config section: {key}")

        # Validate paths
        self._validate_paths(config.get('paths', {}))

        # Validate processing
        self._validate_processing(config.get('processing', {}))

        # Validate annotation
        self._validate_annotation(config.get('annotation', {}))

        # Validate analysis
        self._validate_analysis(config.get('analysis', {}))

    def _validate_paths(self, paths: Dict) -> None:
        """Validate paths section"""
        for key in self.required_paths:
            if key not in paths:
                raise ConfigurationError(f"Missing required path: {key}")

    def _validate_processing(self, processing: Dict) -> None:
        """Validate processing section"""
        for key in self.required_processing:
            if key not in processing:
                raise ConfigurationError(f"Missing required processing config: {key}")

        # Validate required_columns
        required_cols = processing.get('required_columns', [])
        if not isinstance(required_cols, list) or len(required_cols) == 0:
            raise ConfigurationError("required_columns must be a non-empty list")

    def _validate_annotation(self, annotation: Dict) -> None:
        """Validate annotation section"""
        required = ['sialylation_marker', 'fucosylation_marker']
        for key in required:
            if key not in annotation:
                raise ConfigurationError(f"Missing required annotation config: {key}")

    def _validate_analysis(self, analysis: Dict) -> None:
        """Validate analysis section"""
        if 'pca' not in analysis:
            raise ConfigurationError("Missing required analysis config: pca")

        pca = analysis['pca']
        if 'n_components' not in pca:
            raise ConfigurationError("Missing PCA config: n_components")
