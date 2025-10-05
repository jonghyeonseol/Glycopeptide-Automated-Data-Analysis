"""
Configuration Validation for pGlyco Auto Combine
Validates config.yaml structure and values
"""

from pathlib import Path
from typing import Any, Dict, List
import yaml

from .exceptions import (
    ConfigurationError,
    MissingConfigKeyError,
    InvalidConfigValueError
)


class ConfigValidator:
    """Validates configuration dictionary"""

    # Required top-level keys
    REQUIRED_KEYS = ['paths', 'processing', 'annotation', 'analysis', 'visualization']

    # Required nested keys
    REQUIRED_NESTED_KEYS = {
        'paths': ['dataset_dir', 'results_dir', 'output_file'],
        'processing': ['required_columns', 'qc_filters'],
        'annotation': ['sialylation_marker', 'fucosylation_marker'],
        'analysis': ['pca', 'statistical_tests'],
        'visualization': ['figsize', 'dpi', 'colors']
    }

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize validator with config dictionary

        Args:
            config: Configuration dictionary loaded from YAML
        """
        self.config = config
        self.errors: List[str] = []

    def validate(self) -> bool:
        """
        Validate configuration

        Returns:
            True if valid, False otherwise

        Raises:
            ConfigurationError: If validation fails with details
        """
        self._validate_required_keys()
        self._validate_paths()
        self._validate_processing()
        self._validate_annotation()
        self._validate_analysis()
        self._validate_visualization()

        if self.errors:
            error_message = "Configuration validation failed:\n" + "\n".join(
                f"  - {error}" for error in self.errors
            )
            raise ConfigurationError(error_message)

        return True

    def _validate_required_keys(self) -> None:
        """Validate that all required top-level keys exist"""
        for key in self.REQUIRED_KEYS:
            if key not in self.config:
                self.errors.append(f"Missing required key: '{key}'")

        # Validate nested keys
        for parent_key, required_nested in self.REQUIRED_NESTED_KEYS.items():
            if parent_key not in self.config:
                continue  # Already reported above

            for nested_key in required_nested:
                if nested_key not in self.config[parent_key]:
                    self.errors.append(
                        f"Missing required key: '{parent_key}.{nested_key}'"
                    )

    def _validate_paths(self) -> None:
        """Validate paths configuration"""
        if 'paths' not in self.config:
            return

        paths = self.config['paths']

        # Validate dataset_dir
        if 'dataset_dir' in paths:
            if not isinstance(paths['dataset_dir'], str):
                self.errors.append("paths.dataset_dir must be a string")

        # Validate results_dir
        if 'results_dir' in paths:
            if not isinstance(paths['results_dir'], str):
                self.errors.append("paths.results_dir must be a string")

        # Validate output_file
        if 'output_file' in paths:
            if not isinstance(paths['output_file'], str):
                self.errors.append("paths.output_file must be a string")
            elif not paths['output_file'].endswith('.csv'):
                self.errors.append("paths.output_file must end with '.csv'")

    def _validate_processing(self) -> None:
        """Validate processing configuration"""
        if 'processing' not in self.config:
            return

        processing = self.config['processing']

        # Validate required_columns
        if 'required_columns' in processing:
            if not isinstance(processing['required_columns'], list):
                self.errors.append("processing.required_columns must be a list")
            elif len(processing['required_columns']) < 3:
                self.errors.append(
                    "processing.required_columns must contain at least 3 columns "
                    "(Peptide, GlycanComposition, IsotopeArea)"
                )

        # Validate qc_filters
        if 'qc_filters' in processing:
            if not isinstance(processing['qc_filters'], dict):
                self.errors.append("processing.qc_filters must be a dictionary")
            elif 'min_isotope_area' in processing['qc_filters']:
                min_area = processing['qc_filters']['min_isotope_area']
                if not isinstance(min_area, (int, float)) or min_area < 0:
                    self.errors.append(
                        "processing.qc_filters.min_isotope_area must be >= 0"
                    )

    def _validate_annotation(self) -> None:
        """Validate annotation configuration"""
        if 'annotation' not in self.config:
            return

        annotation = self.config['annotation']

        # Validate markers
        for marker_type in ['sialylation_marker', 'fucosylation_marker']:
            if marker_type in annotation:
                marker = annotation[marker_type]
                if not isinstance(marker, str) or len(marker) != 1:
                    self.errors.append(
                        f"annotation.{marker_type} must be a single character string"
                    )

    def _validate_analysis(self) -> None:
        """Validate analysis configuration"""
        if 'analysis' not in self.config:
            return

        analysis = self.config['analysis']

        # Validate PCA settings
        if 'pca' in analysis:
            pca = analysis['pca']

            if 'n_components' in pca:
                n_comp = pca['n_components']
                if not isinstance(n_comp, int) or n_comp < 1:
                    self.errors.append(
                        "analysis.pca.n_components must be a positive integer"
                    )

            if 'log_transform' in pca:
                if not isinstance(pca['log_transform'], bool):
                    self.errors.append(
                        "analysis.pca.log_transform must be a boolean"
                    )

        # Validate statistical_tests settings
        if 'statistical_tests' in analysis:
            stat_tests = analysis['statistical_tests']

            if 'alpha' in stat_tests:
                alpha = stat_tests['alpha']
                if not isinstance(alpha, (int, float)) or not (0 < alpha < 1):
                    self.errors.append(
                        "analysis.statistical_tests.alpha must be between 0 and 1"
                    )

    def _validate_visualization(self) -> None:
        """Validate visualization configuration"""
        if 'visualization' not in self.config:
            return

        viz = self.config['visualization']

        # Validate DPI
        if 'dpi' in viz:
            dpi = viz['dpi']
            if not isinstance(dpi, int) or dpi < 72:
                self.errors.append(
                    "visualization.dpi must be an integer >= 72"
                )

        # Validate figsize
        if 'figsize' in viz:
            if not isinstance(viz['figsize'], dict):
                self.errors.append("visualization.figsize must be a dictionary")
            else:
                for plot_type, size in viz['figsize'].items():
                    if not isinstance(size, list) or len(size) != 2:
                        self.errors.append(
                            f"visualization.figsize.{plot_type} must be a list of 2 numbers [width, height]"
                        )
                    elif not all(isinstance(x, (int, float)) and x > 0 for x in size):
                        self.errors.append(
                            f"visualization.figsize.{plot_type} values must be positive numbers"
                        )

        # Validate colors
        if 'colors' in viz:
            if not isinstance(viz['colors'], dict):
                self.errors.append("visualization.colors must be a dictionary")

        # Validate glycopeptide_comparison settings
        if 'glycopeptide_comparison' in viz:
            comp = viz['glycopeptide_comparison']

            if 'enabled' in comp:
                if not isinstance(comp['enabled'], bool):
                    self.errors.append(
                        "visualization.glycopeptide_comparison.enabled must be a boolean"
                    )

            if 'max_peptides' in comp:
                max_pep = comp['max_peptides']
                if not isinstance(max_pep, int) or max_pep < 1:
                    self.errors.append(
                        "visualization.glycopeptide_comparison.max_peptides must be a positive integer"
                    )

            if 'max_glycans_per_type' in comp:
                max_gly = comp['max_glycans_per_type']
                if not isinstance(max_gly, int) or max_gly < 1:
                    self.errors.append(
                        "visualization.glycopeptide_comparison.max_glycans_per_type must be a positive integer"
                    )


def load_and_validate_config(config_path: str) -> Dict[str, Any]:
    """
    Load configuration from YAML file and validate it

    Args:
        config_path: Path to config.yaml file

    Returns:
        Validated configuration dictionary

    Raises:
        ConfigurationError: If config file is invalid or missing
    """
    config_file = Path(config_path)

    if not config_file.exists():
        raise ConfigurationError(f"Configuration file not found: {config_path}")

    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ConfigurationError(f"Invalid YAML syntax in {config_path}: {str(e)}")
    except Exception as e:
        raise ConfigurationError(f"Failed to load {config_path}: {str(e)}")

    if config is None:
        raise ConfigurationError(f"Configuration file is empty: {config_path}")

    # Validate the loaded config
    validator = ConfigValidator(config)
    validator.validate()

    return config
