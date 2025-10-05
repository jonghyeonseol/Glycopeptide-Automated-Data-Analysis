"""
Base Plot Module
Provides abstract base class for all visualizations
"""

from abc import ABC, abstractmethod
from enum import Enum
from typing import Any, Dict, Optional
from pathlib import Path
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.logger_config import get_logger

logger = get_logger(__name__)


class PlotCategory(Enum):
    """Plot categories for organization"""
    CORE = "core"           # Essential plots (always generated)
    ADVANCED = "advanced"   # Optional/configurable plots
    DIAGNOSTIC = "diagnostic"  # QC and diagnostic plots


class BasePlot(ABC):
    """
    Abstract base class for all visualizations

    Provides:
    - Consistent interface
    - Auto-registration with registry
    - Standardized saving
    - Metadata tracking
    """

    def __init__(
        self,
        name: str,
        category: PlotCategory,
        description: str = "",
        enabled: bool = True
    ):
        """
        Initialize BasePlot

        Args:
            name: Plot identifier
            category: Plot category
            description: Human-readable description
            enabled: Whether plot is enabled by default
        """
        self.name = name
        self.category = category
        self.description = description
        self.enabled = enabled
        self.output_path: Optional[Path] = None

    @abstractmethod
    def generate(
        self,
        data: pd.DataFrame,
        config: Dict[str, Any],
        output_dir: Path,
        **kwargs
    ) -> Path:
        """
        Generate the visualization

        Args:
            data: Input DataFrame
            config: Configuration dictionary
            output_dir: Output directory
            **kwargs: Additional plot-specific arguments

        Returns:
            Path to saved plot file

        Raises:
            PlotGenerationError: If plot generation fails
        """
        pass

    def is_enabled(self, config: Dict[str, Any]) -> bool:
        """
        Check if plot is enabled based on config

        Args:
            config: Configuration dictionary

        Returns:
            True if plot should be generated
        """
        # Can be overridden for plot-specific logic
        return self.enabled

    def get_dependencies(self) -> list:
        """
        Get list of required data/analysis dependencies

        Returns:
            List of dependency names (e.g., ['pca_results', 'vip_scores'])
        """
        return []

    def save_metadata(self, output_dir: Path) -> None:
        """
        Save plot metadata for reproducibility

        Args:
            output_dir: Directory to save metadata
        """
        metadata = {
            'name': self.name,
            'category': self.category.value,
            'description': self.description,
            'output_path': str(self.output_path) if self.output_path else None
        }

        metadata_file = output_dir / f"{self.name}_metadata.json"
        import json
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}: {self.name} ({self.category.value})>"
