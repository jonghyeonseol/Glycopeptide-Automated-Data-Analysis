"""
Visualization Coordinator
Orchestrates plot generation using the registry pattern
"""

from typing import Dict, Any, List
from pathlib import Path
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from .base_plot import BasePlot, PlotCategory
from .registry import VisualizationRegistry
from src.logger_config import get_logger

logger = get_logger(__name__)


class VisualizationCoordinator:
    """
    Coordinates visualization generation

    Features:
    - Uses registry to discover plots
    - Category-based generation (core, advanced, etc.)
    - Dependency resolution
    - Error handling per plot
    - Progress tracking
    """

    def __init__(
        self,
        output_dir: Path,
        config: Dict[str, Any],
        registry: VisualizationRegistry = None
    ):
        """
        Initialize coordinator

        Args:
            output_dir: Directory for output plots
            config: Configuration dictionary
            registry: VisualizationRegistry instance (uses singleton if None)
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.config = config
        self.registry = registry or VisualizationRegistry()
        self.generated_plots: List[Path] = []
        self.failed_plots: List[str] = []

    def generate_all(
        self,
        data: pd.DataFrame,
        **kwargs
    ) -> List[Path]:
        """
        Generate all enabled plots

        Args:
            data: Input DataFrame
            **kwargs: Additional data/results to pass to plots

        Returns:
            List of paths to generated plots
        """
        logger.info("Generating all visualizations...")

        enabled_plots = self.registry.get_enabled(self.config)
        logger.info(f"Found {len(enabled_plots)} enabled plots")

        for plot in enabled_plots:
            self._generate_plot(plot, data, **kwargs)

        logger.info(
            f"\nVisualization Summary:\n"
            f"  Generated: {len(self.generated_plots)}\n"
            f"  Failed: {len(self.failed_plots)}"
        )

        return self.generated_plots

    def generate_category(
        self,
        category: PlotCategory,
        data: pd.DataFrame,
        **kwargs
    ) -> List[Path]:
        """
        Generate plots in a specific category

        Args:
            category: Plot category to generate
            data: Input DataFrame
            **kwargs: Additional data/results

        Returns:
            List of paths to generated plots
        """
        logger.info(f"Generating {category.value} visualizations...")

        plots = self.registry.get_category(category)
        enabled_plots = [p for p in plots if p.is_enabled(self.config)]

        for plot in enabled_plots:
            self._generate_plot(plot, data, **kwargs)

        return self.generated_plots

    def _generate_plot(
        self,
        plot: BasePlot,
        data: pd.DataFrame,
        **kwargs
    ) -> None:
        """
        Generate a single plot with error handling

        Args:
            plot: Plot instance to generate
            data: Input DataFrame
            **kwargs: Additional data/results
        """
        try:
            logger.info(f"  Generating: {plot.name}...")

            # Check dependencies
            missing_deps = self._check_dependencies(plot, kwargs)
            if missing_deps:
                logger.warning(
                    f"  Skipping {plot.name}: Missing dependencies {missing_deps}"
                )
                self.failed_plots.append(plot.name)
                return

            # Generate plot
            output_path = plot.generate(
                data=data,
                config=self.config,
                output_dir=self.output_dir,
                **kwargs
            )

            self.generated_plots.append(output_path)
            logger.info(f"    ✓ Saved to: {output_path.name}")

        except Exception as e:
            logger.error(f"    ✗ Failed to generate {plot.name}: {str(e)}")
            self.failed_plots.append(plot.name)

    def _check_dependencies(
        self,
        plot: BasePlot,
        available_data: Dict
    ) -> List[str]:
        """
        Check if plot dependencies are satisfied

        Args:
            plot: Plot instance
            available_data: Available data/results

        Returns:
            List of missing dependencies
        """
        required_deps = plot.get_dependencies()
        missing = [dep for dep in required_deps if dep not in available_data]
        return missing

    def generate_by_names(
        self,
        plot_names: List[str],
        data: pd.DataFrame,
        **kwargs
    ) -> List[Path]:
        """
        Generate specific plots by name

        Args:
            plot_names: List of plot names to generate
            data: Input DataFrame
            **kwargs: Additional data/results

        Returns:
            List of paths to generated plots
        """
        for name in plot_names:
            plot = self.registry.get(name)
            if plot:
                self._generate_plot(plot, data, **kwargs)
            else:
                logger.warning(f"Plot '{name}' not found in registry")

        return self.generated_plots

    def list_available_plots(self) -> None:
        """List all available plots"""
        self.registry.list_plots()
