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

    visualizer.plot_volcano(state.filtered_data, vip_scores, config=data_prep_config, log2fc_threshold=2.0)
    visualizer.plot_site_specific_heatmap(state.filtered_data, vip_scores)
    visualizer.plot_cv_distribution(state.filtered_data)

    # Correlation visualizations (separate + combined)
    logger.info("Creating correlation visualizations...")
    visualizer.plot_correlation_matrix(state.filtered_data)  # Separate Cancer/Normal
    visualizer.plot_correlation_clustermap(state.filtered_data)  # Separate clustermaps
    visualizer.plot_correlation_matrix_combined(state.filtered_data)  # Combined 47x47 matrix
    visualizer.plot_correlation_cross_group(state.filtered_data)  # Cancer vs Normal only
    visualizer.plot_correlation_clustermap_combined(state.filtered_data)  # Combined clustermap

    visualizer.plot_glycan_venn_diagram(state.filtered_data)
    visualizer.plot_radar_chart(state.filtered_data)
    visualizer.plot_pie_chart_glycan_types(state.filtered_data)
    visualizer.plot_pie_chart_primary_classification(state.filtered_data)
    visualizer.plot_pie_chart_secondary_classification(state.filtered_data)

    # Significance-filtered pie chart (|Log2 FC| >= 2)
    logger.info("Creating significance-filtered pie chart (highly differential features)...")
    visualizer.plot_pie_chart_significant_glycan_types(
        df=state.filtered_data,
        vip_df=vip_scores,
        config=data_prep_config,
        log2fc_threshold=2.0,  # 4-fold change
        fdr_threshold=0.05
    )

    if config['visualization']['glycopeptide_comparison']['enabled']:
        visualizer.plot_glycopeptide_comparison_heatmap(
            df=state.filtered_data,
            vip_scores=vip_scores,
            config=data_prep_config,
            figsize=tuple(config['visualization']['glycopeptide_comparison']['figsize']),
            max_peptides=config['visualization']['glycopeptide_comparison']['max_peptides'],
            max_glycans_per_type=config['visualization']['glycopeptide_comparison']['max_glycans_per_type']
        )

        # Also generate full-scale landscape view (ALL glycopeptides)
        logger.info("Creating full-scale glycopeptide comparison heatmap (complete landscape)...")
        visualizer.plot_glycopeptide_comparison_heatmap_full(
            df=state.filtered_data,
            vip_scores=vip_scores,
            config=data_prep_config
        )

        # Generate glycan-type-specific heatmaps (5 separate heatmaps)
        logger.info("Creating glycan-type-specific comparison heatmaps...")
        for glycan_type in ['HM', 'F', 'S', 'SF', 'C/H']:
            visualizer.plot_glycopeptide_comparison_heatmap_by_type(
                df=state.filtered_data,
                vip_scores=vip_scores,
                glycan_type=glycan_type,
                config=data_prep_config
            )

    # Generate Sankey diagrams
    logger.info("Creating Sankey diagram for glycan type flows...")
    visualizer.plot_glycan_type_sankey(
        df=state.filtered_data,
        vip_scores=vip_scores,
        config=data_prep_config,
        log2fc_threshold=1.0,  # 2-fold change
        fdr_threshold=0.05
    )

    logger.info("Creating Sankey diagram for Group → Glycan Type distribution...")
    visualizer.plot_group_to_glycan_sankey(
        df=state.filtered_data,
        vip_scores=vip_scores,
        config=data_prep_config,
        log2fc_threshold=1.0,  # 2-fold change
        fdr_threshold=0.05
    )


def generate_summary_report(state, config):
    """
    Generate summary report

    Args:
        state: Pipeline state
        config: Configuration
    """
    logger.info("\n[6/6] Generating summary report...")

    from src.reporting import SummaryReport
    from pathlib import Path

    # Read filtering report
    filtering_report_path = Path(config['paths']['results_dir']) / 'filtering_report.txt'
    filtering_report = ''
    if filtering_report_path.exists():
        with open(filtering_report_path, 'r') as f:
            filtering_report = f.read()

    # Generate summary report
    report = SummaryReport()
    summary_content = report.generate({
        'pipeline_state': state,
        'config': config,
        'filtering_report': filtering_report
    })

    # Save to file
    output_path = Path(config['paths']['results_dir']) / 'analysis_summary.txt'
    with open(output_path, 'w') as f:
        f.write(summary_content)

    logger.info(f"Saved summary report to {output_path}")


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

        # Generate summary report
        generate_summary_report(state, state.config)

        logger.info("\n" + "="*80)
        logger.info("Pipeline completed successfully!")
        logger.info("="*80)

        print("\n" + "="*80)
        print("Analysis Complete!")
        print("="*80)
        print("\nResults saved to: Results/")
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
