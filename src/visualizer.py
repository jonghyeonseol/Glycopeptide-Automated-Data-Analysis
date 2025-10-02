"""
Visualizer Module for pGlyco Auto Combine
Handles data visualization (PCA, boxplot, heatmap, histogram, VIP scores, distribution)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

# Import all mixin classes
from src.plots.pca_plot import PCAPlotMixin
from src.plots.boxplot import BoxplotMixin
from src.plots.heatmap import HeatmapMixin
from src.plots.histogram import HistogramMixin
from src.plots.vip_score_plot import VIPScorePlotMixin
from src.plots.vip_score_plot_r import VIPScorePlotRMixin
from src.plots.distribution_plot import DistributionPlotMixin
from src.plots.volcano_plot import VolcanoPlotMixin
from src.plots.site_specific_heatmap import SiteSpecificHeatmapMixin
from src.plots.cv_distribution_plot import CVDistributionPlotMixin
from src.plots.correlation_matrix_plot import CorrelationMatrixPlotMixin
from src.plots.venn_diagram_plot import VennDiagramPlotMixin
from src.plots.radar_chart_plot import RadarChartPlotMixin

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class GlycanVisualizer(
    PCAPlotMixin,
    BoxplotMixin,
    HeatmapMixin,
    HistogramMixin,
    VIPScorePlotMixin,
    VIPScorePlotRMixin,
    DistributionPlotMixin,
    VolcanoPlotMixin,
    SiteSpecificHeatmapMixin,
    CVDistributionPlotMixin,
    CorrelationMatrixPlotMixin,
    VennDiagramPlotMixin,
    RadarChartPlotMixin
):
    """Create visualizations for glycoproteomics data"""

    def __init__(self, output_dir: str, dpi: int = 300, colors: dict = None):
        """
        Initialize GlycanVisualizer

        Args:
            output_dir: Directory to save plots
            dpi: Resolution for saved figures
            colors: Color scheme for glycan types
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.dpi = dpi

        # Default color scheme
        self.colors = colors or {
            'Non': '#CCCCCC',
            'Sialylated': "#FF00EA",
            'Fucosylated': "#DB3434",
            'Both': "#FFB700"
        }

        # Set style
        sns.set_style("whitegrid")
        plt.rcParams['font.size'] = 10
        plt.rcParams['axes.labelsize'] = 12
        plt.rcParams['axes.titlesize'] = 14

    def plot_all(self, df: pd.DataFrame, pca_results: dict, boxplot_data: pd.DataFrame, boxplot_data_extended: pd.DataFrame = None):
        """
        Generate all plots

        Args:
            df: Annotated DataFrame
            pca_results: PCA results from analyzer
            boxplot_data: Boxplot data from analyzer
            boxplot_data_extended: Extended boxplot data from analyzer (optional)
        """
        logger.info("Generating all visualizations...")

        self.plot_pca(pca_results)
        self.plot_pca_by_glycan_type(df, pca_results)
        self.plot_boxplot(boxplot_data)

        # Plot extended boxplot if data is provided
        if boxplot_data_extended is not None:
            self.plot_boxplot_extended(boxplot_data_extended)

        self.plot_glycan_type_distribution(df)
        self.plot_heatmap(df)
        self.plot_heatmap_full_profile(df)
        self.plot_histogram_normalized(df)

        logger.info(f"All visualizations saved to {self.output_dir}")


if __name__ == "__main__":
    # Test
    pass
