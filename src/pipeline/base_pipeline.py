"""
Base Pipeline Abstract Class
Defines the interface for all pipeline implementations
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, List, Optional
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class PipelineState:
    """
    Container for pipeline state and intermediate results
    Tracks data flow through pipeline stages
    """

    def __init__(self):
        self.raw_data: Optional[Any] = None
        self.integrated_data: Optional[Any] = None
        self.annotated_data_raw: Optional[Any] = None
        self.filtered_data: Optional[Any] = None
        self.analysis_results: Dict[str, Any] = {}
        self.visualization_paths: List[Path] = []
        self.reports: Dict[str, Path] = {}

    def get(self, key: str, default: Any = None) -> Any:
        """Get a value from state"""
        return getattr(self, key, default)

    def set(self, key: str, value: Any) -> None:
        """Set a value in state"""
        setattr(self, key, value)

    def update(self, **kwargs) -> None:
        """Update multiple state values"""
        for key, value in kwargs.items():
            setattr(self, key, value)


class BasePipeline(ABC):
    """
    Abstract base class for pipeline implementations

    Defines the interface and common functionality for all pipelines.
    Subclasses must implement run() and can override hooks for customization.
    """

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize pipeline

        Args:
            config: Configuration dictionary
        """
        self.config = config
        self.state = PipelineState()
        self.workflows: List = []

    @abstractmethod
    def run(self) -> PipelineState:
        """
        Execute the pipeline

        Returns:
            PipelineState with results
        """

    def add_workflow(self, workflow) -> 'BasePipeline':
        """
        Add a workflow to the pipeline

        Args:
            workflow: Workflow instance to add

        Returns:
            Self for chaining
        """
        self.workflows.append(workflow)
        return self

    def on_start(self) -> None:
        """Hook called before pipeline starts"""
        logger.info("=" * 80)
        logger.info("Pipeline Starting")
        logger.info("=" * 80)

    def on_complete(self) -> None:
        """Hook called after pipeline completes"""
        logger.info("=" * 80)
        logger.info("Pipeline Completed Successfully")
        logger.info("=" * 80)

    def on_error(self, error: Exception) -> None:
        """Hook called when pipeline encounters an error"""
        logger.error(f"Pipeline Failed: {str(error)}", exc_info=True)

    def validate_config(self) -> None:
        """Validate pipeline configuration"""
        required_keys = ['paths', 'processing', 'annotation', 'analysis', 'visualization']
        for key in required_keys:
            if key not in self.config:
                raise ValueError(f"Missing required config key: {key}")
