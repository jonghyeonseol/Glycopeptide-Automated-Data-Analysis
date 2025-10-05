"""
Configuration Manager
Unified configuration loading and validation
"""

from typing import Dict, Any, Optional
from pathlib import Path
import yaml

import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from .schema import ConfigSchema
from .defaults import DEFAULT_CONFIG
from src.exceptions import ConfigurationError
from src.logger_config import get_logger

logger = get_logger(__name__)


class ConfigManager:
    """
    Unified configuration manager

    Features:
    - YAML loading
    - Schema validation
    - Default values
    - Config merging
    - Type checking
    """

    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize config manager

        Args:
            config_path: Path to config.yaml (optional)
        """
        self.config_path = config_path
        self.config: Dict[str, Any] = {}
        self.schema = ConfigSchema()

        if config_path:
            self.load(config_path)

    def load(self, config_path: str) -> Dict[str, Any]:
        """
        Load configuration from YAML file

        Args:
            config_path: Path to config file

        Returns:
            Configuration dictionary

        Raises:
            ConfigurationError: If config is invalid
        """
        config_file = Path(config_path)

        if not config_file.exists():
            raise ConfigurationError(f"Config file not found: {config_path}")

        try:
            with open(config_file) as f:
                loaded_config = yaml.safe_load(f)

            # Merge with defaults
            self.config = self._merge_with_defaults(loaded_config)

            # Validate
            self.schema.validate(self.config)

            logger.info(f"Loaded configuration from {config_path}")
            return self.config

        except yaml.YAMLError as e:
            raise ConfigurationError(f"Invalid YAML: {str(e)}")
        except Exception as e:
            raise ConfigurationError(f"Failed to load config: {str(e)}")

    def _merge_with_defaults(self, loaded: Dict) -> Dict:
        """
        Merge loaded config with defaults

        Args:
            loaded: Loaded configuration

        Returns:
            Merged configuration
        """
        merged = DEFAULT_CONFIG.copy()

        # Deep merge
        for key, value in loaded.items():
            if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
                merged[key] = {**merged[key], **value}
            else:
                merged[key] = value

        return merged

    def get(self, key: str, default: Any = None) -> Any:
        """
        Get configuration value

        Args:
            key: Configuration key (supports dot notation: 'paths.dataset_dir')
            default: Default value if key not found

        Returns:
            Configuration value
        """
        keys = key.split('.')
        value = self.config

        for k in keys:
            if isinstance(value, dict) and k in value:
                value = value[k]
            else:
                return default

        return value

    def set(self, key: str, value: Any) -> None:
        """
        Set configuration value

        Args:
            key: Configuration key (supports dot notation)
            value: Value to set
        """
        keys = key.split('.')
        config = self.config

        for k in keys[:-1]:
            if k not in config:
                config[k] = {}
            config = config[k]

        config[keys[-1]] = value

    def to_dict(self) -> Dict[str, Any]:
        """
        Get configuration as dictionary

        Returns:
            Configuration dictionary
        """
        return self.config.copy()

    def save(self, output_path: str) -> None:
        """
        Save configuration to YAML file

        Args:
            output_path: Path to save config
        """
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, 'w') as f:
            yaml.dump(self.config, f, default_flow_style=False, sort_keys=False)

        logger.info(f"Saved configuration to {output_path}")

    @staticmethod
    def from_dict(config_dict: Dict[str, Any]) -> 'ConfigManager':
        """
        Create ConfigManager from dictionary

        Args:
            config_dict: Configuration dictionary

        Returns:
            ConfigManager instance
        """
        manager = ConfigManager()
        manager.config = config_dict
        manager.schema.validate(config_dict)
        return manager
