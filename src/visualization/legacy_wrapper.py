"""
Legacy Plot Wrapper
Wraps existing GlycanVisualizer plots for use with new registry

This provides a bridge between old plot modules and new registry pattern
without rewriting all 18+ plot modules.
"""

from pathlib import Path
from typing import Dict, Any
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from .base_plot import BasePlot, PlotCategory
from .registry import VisualizationRegistry
from src.visualizer import GlycanVisualizer
from src.logger_config import get_logger

logger = get_logger(__name__)


class LegacyPlotWrapper(BasePlot):
    """
    Wrapper for existing GlycanVisualizer methods

    Allows gradual migration to registry pattern while
    maintaining compatibility with existing plot code.
    """

    def __init__(
        self,
        name: str,
        category: PlotCategory,
        plot_method: str,
        description: str = "",
        enabled: bool = True,
        dependencies: list = None
    ):
        """
        Initialize wrapper

        Args:
            name: Plot identifier
            category: Plot category
            plot_method: Method name in GlycanVisualizer
            description: Human-readable description
            enabled: Default enabled state
            dependencies: List of required dependencies
        """
        super().__init__(name, category, description, enabled)
        self.plot_method = plot_method
        self.dependencies = dependencies or []

    def generate(
        self,
        data: pd.DataFrame,
        config: Dict[str, Any],
        output_dir: Path,
        **kwargs
    ) -> Path:
        """
        Generate plot using legacy visualizer

        Args:
            data: Input DataFrame
            config: Configuration
            output_dir: Output directory
            **kwargs: Additional arguments (pca_results, vip_scores, etc.)

        Returns:
            Path to generated plot
        """
        # Create visualizer instance
        visualizer = GlycanVisualizer(
            output_dir=str(output_dir),
            dpi=config.get('visualization', {}).get('dpi', 300),
            colors=config.get('visualization', {}).get('colors', {})
        )

        # Get the plot method
        method = getattr(visualizer, self.plot_method, None)
        if not method:
            raise AttributeError(f"GlycanVisualizer has no method '{self.plot_method}'")

        # Call the method with appropriate arguments
        method_args = self._prepare_arguments(data, kwargs)
        method(**method_args)

        # Return the expected output path
        output_path = output_dir / f"{self.name}.png"
        self.output_path = output_path
        return output_path

    def _prepare_arguments(
        self,
        data: pd.DataFrame,
        kwargs: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Prepare arguments for legacy plot method

        Args:
            data: Input DataFrame
            kwargs: Additional arguments

        Returns:
            Dictionary of arguments for plot method
        """
        # Map common argument patterns
        args = {'df': data}

        # Add common optional arguments if available
        if 'pca_results' in kwargs:
            args['pca_results'] = kwargs['pca_results']
        if 'vip_scores' in kwargs:
            args['vip_scores'] = kwargs['vip_scores']
        if 'boxplot_data' in kwargs:
            args['boxplot_data'] = kwargs['boxplot_data']
        if 'data_prep_config' in kwargs:
            args['config'] = kwargs['data_prep_config']

        return args

    def get_dependencies(self) -> list:
        """Return list of dependencies"""
        return self.dependencies


# ==================== Register Core Plots ====================

@VisualizationRegistry.register(
    'pca_plot',
    PlotCategory.CORE,
    'PCA plot showing sample separation',
    enabled=True
)
class PCAPlotWrapper(LegacyPlotWrapper):
    def __init__(self, **kwargs):
        super().__init__(
            name='pca_plot',
            category=PlotCategory.CORE,
            plot_method='plot_pca',
            description='PCA plot showing sample separation',
            dependencies=['pca_results']
        )


@VisualizationRegistry.register(
    'boxplot',
    PlotCategory.CORE,
    'Boxplot of glycan type intensities',
    enabled=True
)
class BoxplotWrapper(LegacyPlotWrapper):
    def __init__(self, **kwargs):
        super().__init__(
            name='boxplot_glycan_types',
            category=PlotCategory.CORE,
            plot_method='plot_boxplot',
            description='Boxplot of glycan type intensities',
            dependencies=['boxplot_data']
        )


@VisualizationRegistry.register(
    'heatmap',
    PlotCategory.CORE,
    'Heatmap of top glycopeptides',
    enabled=True
)
class HeatmapWrapper(LegacyPlotWrapper):
    def __init__(self, **kwargs):
        super().__init__(
            name='heatmap_top_glycopeptides',
            category=PlotCategory.CORE,
            plot_method='plot_heatmap',
            description='Heatmap of top glycopeptides'
        )


# ==================== Register Advanced Plots ====================

@VisualizationRegistry.register(
    'volcano_plot',
    PlotCategory.ADVANCED,
    'Volcano plot for differential expression',
    enabled=True
)
class VolcanoPlotWrapper(LegacyPlotWrapper):
    def __init__(self, **kwargs):
        super().__init__(
            name='volcano_plot',
            category=PlotCategory.ADVANCED,
            plot_method='plot_volcano',
            description='Volcano plot for differential expression',
            dependencies=['vip_scores', 'data_prep_config']
        )


@VisualizationRegistry.register(
    'vip_scores_r',
    PlotCategory.ADVANCED,
    'VIP scores visualization (R/ggplot2)',
    enabled=True
)
class VIPScoresRWrapper(LegacyPlotWrapper):
    def __init__(self, **kwargs):
        super().__init__(
            name='vip_score_glycopeptide_r',
            category=PlotCategory.ADVANCED,
            plot_method='plot_vip_scores_glycopeptide_r',
            description='VIP scores by glycopeptide (R)',
            dependencies=['vip_scores']
        )
