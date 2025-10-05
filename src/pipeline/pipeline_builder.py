"""
Pipeline Builder
Provides fluent API for building pipelines
"""

from typing import Dict, Any, List, Optional
from pathlib import Path

from .glyco_pipeline import GlycoPipeline
from ..config_validator import load_and_validate_config
from ..logger_config import setup_logging, get_logger

logger = get_logger(__name__)


class GlycoPipelineBuilder:
    """
    Builder for GlycoPipeline

    Provides fluent interface for pipeline construction:
        pipeline = (GlycoPipelineBuilder()
            .with_config('config.yaml')
            .build())
    """

    def __init__(self):
        self._config: Optional[Dict[str, Any]] = None
        self._config_path: Optional[str] = None
        self._log_configured = False

    def with_config(self, config_path: str) -> 'GlycoPipelineBuilder':
        """
        Load configuration from file

        Args:
            config_path: Path to config.yaml

        Returns:
            Self for chaining
        """
        self._config_path = config_path
        return self

    def with_config_dict(self, config: Dict[str, Any]) -> 'GlycoPipelineBuilder':
        """
        Use provided configuration dictionary

        Args:
            config: Configuration dictionary

        Returns:
            Self for chaining
        """
        self._config = config
        return self

    def with_logging(self) -> 'GlycoPipelineBuilder':
        """
        Setup logging

        Returns:
            Self for chaining
        """
        if not self._log_configured:
            setup_logging()
            self._log_configured = True
        return self

    def build(self) -> GlycoPipeline:
        """
        Build the pipeline

        Returns:
            Configured GlycoPipeline instance

        Raises:
            ValueError: If no configuration provided
        """
        # Setup logging if not already done
        if not self._log_configured:
            self.with_logging()

        # Load config if path provided
        if self._config_path and not self._config:
            logger.info(f"Loading configuration from {self._config_path}")
            self._config = load_and_validate_config(self._config_path)

        if not self._config:
            raise ValueError("No configuration provided. Use with_config() or with_config_dict()")

        logger.info("Building GlycoPipeline...")
        return GlycoPipeline(self._config)

    def build_and_run(self) -> Any:
        """
        Build and immediately run the pipeline

        Returns:
            Pipeline state with results
        """
        pipeline = self.build()
        return pipeline.run()
