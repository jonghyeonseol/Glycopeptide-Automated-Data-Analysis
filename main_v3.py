#!/usr/bin/env python3
"""
pGlyco Auto Combine - Main Pipeline (v3.0)
Redesigned architecture with workflow-based pipeline

This is a SIMPLIFIED version using the new pipeline architecture.
Compare with main.py (380 lines) - this is only ~100 lines!
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from src.pipeline import GlycoPipelineBuilder
from src.exceptions import PGlycoAutoError
from src.logger_config import get_logger
from src.visualizer import GlycanVisualizer
from src.data_preparation import get_standard_config_from_dict

logger = get_logger(__name__)


def run_visualizations(state, config):
    """
    Run all visualizations

    Args:
        state: Pipeline state with analysis results
        config: Configuration dictionary
    """
    logger.info("\n[5/6] Creating visualizations...")

    results_dir = config['paths']['results_dir']
    visualizer = GlycanVisualizer(
        output_dir=results_dir,
        dpi=config['visualization']['dpi'],
        colors=config['visualization']['colors']
    )

    # Core visualizations
    visualizer.plot_all(
        df=state.filtered_data,
        pca_results=state.analysis_results['pca'],
        boxplot_data=state.analysis_results['boxplot_data'],
        boxplot_data_extended=state.analysis_results['boxplot_data_extended']
    )

    # Primary/Secondary classification plots
    logger.info("Creating classification visualizations...")
    visualizer.plot_histogram_primary_classification(state.filtered_data, normalization='raw')
    visualizer.plot_histogram_primary_classification(state.filtered_data, normalization='aggregated')
    visualizer.plot_histogram_secondary_classification(state.filtered_data, normalization='raw')
    visualizer.plot_histogram_secondary_classification(state.filtered_data, normalization='aggregated')
    visualizer.plot_histogram_cancer_vs_normal_primary(state.filtered_data)
    visualizer.plot_histogram_cancer_vs_normal_secondary(state.filtered_data)

    # VIP score visualizations (R/ggplot2)
    logger.info("Creating VIP score plots with R/ggplot2...")
    vip_scores = state.analysis_results['plsda']['vip_scores']
    visualizer.plot_vip_scores_glycopeptide_r(state.filtered_data, vip_scores)
    visualizer.plot_vip_scores_glycan_composition_r(state.filtered_data, vip_scores)
    visualizer.plot_vip_scores_peptide_r(state.filtered_data, vip_scores)
    visualizer.plot_vip_scores_peptide_grouped_r(state.filtered_data, vip_scores)

    # Box plots
    logger.info("Creating box plots...")
    visualizer.plot_boxplot_primary_classification(state.filtered_data, normalization='raw')
    visualizer.plot_boxplot_primary_classification(state.filtered_data, normalization='aggregated')
    visualizer.plot_boxplot_secondary_classification(state.filtered_data, normalization='raw')
    visualizer.plot_boxplot_cancer_vs_normal_primary(state.filtered_data)
    visualizer.plot_boxplot_cancer_vs_normal_secondary(state.filtered_data)

    # Advanced visualizations
    logger.info("Creating advanced visualizations...")
    data_prep_config = get_standard_config_from_dict(config)

    visualizer.plot_volcano(state.filtered_data, vip_scores, config=data_prep_config)
    visualizer.plot_site_specific_heatmap(state.filtered_data, vip_scores)
    visualizer.plot_cv_distribution(state.filtered_data)
    visualizer.plot_correlation_matrix(state.filtered_data)
    visualizer.plot_correlation_clustermap(state.filtered_data)
    visualizer.plot_glycan_venn_diagram(state.filtered_data)
    visualizer.plot_radar_chart(state.filtered_data)
    visualizer.plot_pie_chart_glycan_types(state.filtered_data)
    visualizer.plot_pie_chart_primary_classification(state.filtered_data)
    visualizer.plot_pie_chart_secondary_classification(state.filtered_data)

    if config['visualization']['glycopeptide_comparison']['enabled']:
        visualizer.plot_glycopeptide_comparison_heatmap(
            df=state.filtered_data,
            vip_scores=vip_scores,
            config=data_prep_config,
            figsize=tuple(config['visualization']['glycopeptide_comparison']['figsize']),
            max_peptides=config['visualization']['glycopeptide_comparison']['max_peptides'],
            max_glycans_per_type=config['visualization']['glycopeptide_comparison']['max_glycans_per_type']
        )


def generate_summary_report(state, config):
    """
    Generate summary report

    Args:
        state: Pipeline state
        config: Configuration
    """
    logger.info("\n[6/6] Generating summary report...")

    # This will be moved to reporting/ module in future phase
    # For now, use existing summary generation logic
    from main import generate_summary  # Import from old main
    # TODO: Replace with ReportBuilder pattern in Phase 4


def main():
    """Main entry point"""
    logger.info("="*80)
    logger.info("pGlyco Auto Combine v3.0 - Workflow-Based Pipeline")
    logger.info("="*80)

    try:
        # Build and run pipeline (simplified!)
        state = (GlycoPipelineBuilder()
                .with_config('config.yaml')
                .with_logging()
                .build_and_run())

        # Run visualizations (will be moved to workflow in next phase)
        run_visualizations(state, state.config)

        # Generate reports (will be moved to ReportBuilder in Phase 4)
        # generate_summary_report(state, config)

        logger.info("\n" + "="*80)
        logger.info("Pipeline completed successfully!")
        logger.info("="*80)

        print("\n" + "="*80)
        print("Analysis Complete!")
        print("="*80)
        print(f"\nResults saved to: Results/")
        print("\n" + "="*80 + "\n")

    except PGlycoAutoError as e:
        logger.error(f"Pipeline failed: {str(e)}")
        sys.exit(1)
    except KeyboardInterrupt:
        logger.warning("\nPipeline interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
