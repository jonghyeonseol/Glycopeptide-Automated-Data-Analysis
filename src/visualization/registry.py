"""
Visualization Registry
Manages plot registration and retrieval
"""

from typing import Dict, List, Optional, Type
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from .base_plot import BasePlot, PlotCategory
from src.logger_config import get_logger

logger = get_logger(__name__)


class VisualizationRegistry:
    """
    Registry for visualization plots

    Provides:
    - Auto-registration of plots
    - Category-based retrieval
    - Plot discovery
    - Dependency resolution
    """

    _instance = None
    _plots: Dict[str, BasePlot] = {}

    def __new__(cls):
        """Singleton pattern"""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    @classmethod
    def register(
        cls,
        name: str,
        category: PlotCategory = PlotCategory.CORE,
        description: str = "",
        enabled: bool = True
    ):
        """
        Decorator to register a plot class

        Usage:
            @VisualizationRegistry.register('pca_plot', PlotCategory.CORE)
            class PCAPlotWrapper(BasePlot):
                ...

        Args:
            name: Plot identifier
            category: Plot category
            description: Human-readable description
            enabled: Default enabled state

        Returns:
            Decorator function
        """
        def decorator(plot_class: Type[BasePlot]):
            # Instantiate the plot
            plot_instance = plot_class(
                name=name,
                category=category,
                description=description,
                enabled=enabled
            )
            cls._plots[name] = plot_instance
            logger.debug(f"Registered plot: {name} ({category.value})")
            return plot_class

        return decorator

    @classmethod
    def get(cls, name: str) -> Optional[BasePlot]:
        """
        Get a plot by name

        Args:
            name: Plot identifier

        Returns:
            Plot instance or None
        """
        return cls._plots.get(name)

    @classmethod
    def get_category(cls, category: PlotCategory) -> List[BasePlot]:
        """
        Get all plots in a category

        Args:
            category: Plot category

        Returns:
            List of plot instances
        """
        return [
            plot for plot in cls._plots.values()
            if plot.category == category
        ]

    @classmethod
    def get_all(cls) -> List[BasePlot]:
        """
        Get all registered plots

        Returns:
            List of all plot instances
        """
        return list(cls._plots.values())

    @classmethod
    def get_enabled(cls, config: Dict) -> List[BasePlot]:
        """
        Get all enabled plots based on config

        Args:
            config: Configuration dictionary

        Returns:
            List of enabled plot instances
        """
        return [
            plot for plot in cls._plots.values()
            if plot.is_enabled(config)
        ]

    @classmethod
    def list_plots(cls) -> None:
        """Print all registered plots"""
        logger.info("Registered Plots:")
        for category in PlotCategory:
            plots = cls.get_category(category)
            if plots:
                logger.info(f"\n  {category.value.upper()}:")
                for plot in plots:
                    status = "✓" if plot.enabled else "✗"
                    logger.info(f"    [{status}] {plot.name} - {plot.description}")

    @classmethod
    def clear(cls) -> None:
        """Clear all registered plots (for testing)"""
        cls._plots.clear()
