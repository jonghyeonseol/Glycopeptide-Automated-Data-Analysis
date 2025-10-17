#!/usr/bin/env python3
"""
pGlyco Auto Combine - Main Pipeline (v3.0)
Redesigned architecture with workflow-based pipeline

This is a SIMPLIFIED version using the new pipeline architecture.
Compare with main.py (380 lines) - this is only ~100 lines!
"""

import sys
import time
from pathlib import Path
from contextlib import contextmanager

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from src.pipeline import GlycoPipelineBuilder
from src.exceptions import PGlycoAutoError
from src.logger_config import get_logger
from src.visualizer import GlycanVisualizer
from src.data_preparation import get_standard_config_from_dict
from src.constants import (
    LOG2FC_THRESHOLD_STRICT,
    LOG2FC_THRESHOLD_MODERATE,
    FDR_THRESHOLD_DEFAULT,
    LOG_TRANSFORM_PSEUDOCOUNT,
    DETECTION_FILTER_DEFAULT_PCT,
    DETECTION_FILTER_DEFAULT_SAMPLES
)

logger = get_logger(__name__)


@contextmanager
def timed_operation(operation_name: str, config: dict = None):
    """
    Context manager for timing pipeline operations

    Args:
        operation_name: Name of the operation being timed
        config: Optional config dict to check if timing is enabled

    Yields:
        None

    Usage:
        with timed_operation("Data Integration"):
            # ... perform operation ...
    """
    # Check if timing is enabled (default to True if not specified)
    timing_enabled = True
    if config:
        timing_enabled = config.get('execution', {}).get('timing', {}).get('enabled', True)

    if timing_enabled:
        logger.info(f"▶ Starting: {operation_name}")
        start_time = time.time()

    try:
        yield
    finally:
        if timing_enabled:
            duration = time.time() - start_time
            logger.info(f"✓ Completed: {operation_name} (Duration: {duration:.2f}s)")



def run_visualizations(state, config):
    """
    Run all visualizations (Phase 11.1: with progress reporting)

    Args:
        state: Pipeline state with analysis results
        config: Configuration dictionary
    """
    logger.info("\n[5/6] Creating visualizations...")

    # Phase 11.1: Check if progress reporting is enabled
    progress_enabled = config.get('execution', {}).get('progress_reporting', {}).get('visualizations', True)

    # Count total visualizations for progress tracking
    total_plots = 39  # Approximate count of all visualizations
    current_plot = 0

    def log_progress(plot_name: str):
        """Log progress for individual plots"""
        nonlocal current_plot
        current_plot += 1
        if progress_enabled:
            logger.info(f"  [{current_plot}/{total_plots}] Creating {plot_name}")

    results_dir = config['paths']['results_dir']
    visualizer = GlycanVisualizer(
        output_dir=results_dir,
        dpi=config['visualization']['dpi'],
        colors=config['visualization']['colors']
    )

    # Core visualizations
    log_progress("core visualizations (PCA, heatmaps, histograms)")
    visualizer.plot_all(
        df=state.filtered_data,
        pca_results=state.analysis_results['pca'],
        boxplot_data=state.analysis_results['boxplot_data'],
        boxplot_data_extended=state.analysis_results['boxplot_data_extended']
    )

    # Primary/Secondary classification plots
    log_progress("histogram - primary classification (raw)")
    visualizer.plot_histogram_primary_classification(state.filtered_data, normalization='raw')
    log_progress("histogram - primary classification (aggregated)")
    visualizer.plot_histogram_primary_classification(state.filtered_data, normalization='aggregated')
    log_progress("histogram - secondary classification (raw)")
    visualizer.plot_histogram_secondary_classification(state.filtered_data, normalization='raw')
    log_progress("histogram - secondary classification (aggregated)")
    visualizer.plot_histogram_secondary_classification(state.filtered_data, normalization='aggregated')
    log_progress("histogram - cancer vs normal (primary)")
    visualizer.plot_histogram_cancer_vs_normal_primary(state.filtered_data)
    log_progress("histogram - cancer vs normal (secondary)")
    visualizer.plot_histogram_cancer_vs_normal_secondary(state.filtered_data)

    # VIP score visualizations (R/ggplot2)
    vip_scores = state.analysis_results['plsda']['vip_scores']
    log_progress("VIP scores - glycopeptide (R/ggplot2)")
    visualizer.plot_vip_scores_glycopeptide_r(state.filtered_data, vip_scores)
    log_progress("VIP scores - glycan composition (R/ggplot2)")
    visualizer.plot_vip_scores_glycan_composition_r(state.filtered_data, vip_scores)
    log_progress("VIP scores - peptide (R/ggplot2)")
    visualizer.plot_vip_scores_peptide_r(state.filtered_data, vip_scores)
    log_progress("VIP scores - peptide grouped (R/ggplot2)")
    visualizer.plot_vip_scores_peptide_grouped_r(state.filtered_data, vip_scores)

    # Box plots
    log_progress("boxplot - primary classification (raw)")
    visualizer.plot_boxplot_primary_classification(state.filtered_data, normalization='raw')
    log_progress("boxplot - primary classification (aggregated)")
    visualizer.plot_boxplot_primary_classification(state.filtered_data, normalization='aggregated')
    log_progress("boxplot - secondary classification (raw)")
    visualizer.plot_boxplot_secondary_classification(state.filtered_data, normalization='raw')
    log_progress("boxplot - cancer vs normal (primary)")
    visualizer.plot_boxplot_cancer_vs_normal_primary(state.filtered_data)
    log_progress("boxplot - cancer vs normal (secondary)")
    visualizer.plot_boxplot_cancer_vs_normal_secondary(state.filtered_data)

    # Advanced visualizations
    data_prep_config = get_standard_config_from_dict(config)

    # Get thresholds from config (Phase 11.1: use centralized thresholds)
    log2fc_strict = config.get('analysis', {}).get('differential_expression', {}).get('log2fc_strict', LOG2FC_THRESHOLD_STRICT)

    log_progress("volcano plot")
    visualizer.plot_volcano(state.filtered_data, vip_scores, config=data_prep_config, log2fc_threshold=log2fc_strict)
    log_progress("site-specific heatmap")
    visualizer.plot_site_specific_heatmap(state.filtered_data, vip_scores)
    log_progress("CV distribution")
    visualizer.plot_cv_distribution(state.filtered_data)

    # Correlation visualizations (separate + combined)
    log_progress("correlation matrix (separate)")
    visualizer.plot_correlation_matrix(state.filtered_data)  # Separate Cancer/Normal
    log_progress("correlation clustermap (separate)")
    visualizer.plot_correlation_clustermap(state.filtered_data)  # Separate clustermaps
    log_progress("correlation matrix (combined 47x47)")
    visualizer.plot_correlation_matrix_combined(state.filtered_data)  # Combined 47x47 matrix
    log_progress("correlation cross-group")
    visualizer.plot_correlation_cross_group(state.filtered_data)  # Cancer vs Normal only
    log_progress("correlation clustermap (combined)")
    visualizer.plot_correlation_clustermap_combined(state.filtered_data)  # Combined clustermap

    log_progress("Venn diagram")
    visualizer.plot_glycan_venn_diagram(state.filtered_data)
    log_progress("radar chart")
    visualizer.plot_radar_chart(state.filtered_data)
    log_progress("pie chart - glycan types")
    visualizer.plot_pie_chart_glycan_types(state.filtered_data)
    log_progress("pie chart - primary classification")
    visualizer.plot_pie_chart_primary_classification(state.filtered_data)
    log_progress("pie chart - secondary classification")
    visualizer.plot_pie_chart_secondary_classification(state.filtered_data)

    # Get FDR threshold from config (Phase 11.1: use centralized thresholds)
    fdr_threshold = config.get('analysis', {}).get('differential_expression', {}).get('fdr_threshold', FDR_THRESHOLD_DEFAULT)

    # Significance-filtered pie chart (|Log2 FC| >= 2)
    log_progress("pie chart - significant features (highly differential)")
    visualizer.plot_pie_chart_significant_glycan_types(
        df=state.filtered_data,
        vip_df=vip_scores,
        config=data_prep_config,
        log2fc_threshold=log2fc_strict,  # Use strict threshold from config
        fdr_threshold=fdr_threshold
    )

    if config['visualization']['glycopeptide_comparison']['enabled']:
        log_progress("glycopeptide comparison heatmap")
        visualizer.plot_glycopeptide_comparison_heatmap(
            df=state.filtered_data,
            vip_scores=vip_scores,
            config=data_prep_config,
            figsize=tuple(config['visualization']['glycopeptide_comparison']['figsize']),
            max_peptides=config['visualization']['glycopeptide_comparison']['max_peptides'],
            max_glycans_per_type=config['visualization']['glycopeptide_comparison']['max_glycans_per_type']
        )

        # Also generate full-scale landscape view (ALL glycopeptides)
        log_progress("glycopeptide comparison heatmap (full landscape)")
        visualizer.plot_glycopeptide_comparison_heatmap_full(
            df=state.filtered_data,
            vip_scores=vip_scores,
            config=data_prep_config
        )

        # Generate glycan-type-specific heatmaps (5 separate heatmaps)
        for glycan_type in ['HM', 'F', 'S', 'SF', 'C/H']:
            log_progress(f"glycopeptide comparison heatmap ({glycan_type})")
            visualizer.plot_glycopeptide_comparison_heatmap_by_type(
                df=state.filtered_data,
                vip_scores=vip_scores,
                glycan_type=glycan_type,
                config=data_prep_config
            )

    # Get moderate threshold for Sankey diagrams (Phase 11.1: use centralized thresholds)
    log2fc_moderate = config.get('analysis', {}).get('differential_expression', {}).get('log2fc_moderate', LOG2FC_THRESHOLD_MODERATE)

    # Generate Sankey diagrams
    log_progress("Sankey diagram - glycan type flows")
    visualizer.plot_glycan_type_sankey(
        df=state.filtered_data,
        vip_scores=vip_scores,
        config=data_prep_config,
        log2fc_threshold=log2fc_moderate,  # Use moderate threshold from config
        fdr_threshold=fdr_threshold
    )

    log_progress("Sankey diagram - Group → Glycan Type distribution")
    visualizer.plot_group_to_glycan_sankey(
        df=state.filtered_data,
        vip_scores=vip_scores,
        config=data_prep_config,
        log2fc_threshold=log2fc_moderate,  # Use moderate threshold from config
        fdr_threshold=fdr_threshold
    )

    logger.info(f"✓ Visualization generation complete: {current_plot} plots created")


def generate_summary_report(state, config):
    """
    Generate summary report

    Args:
        state: Pipeline state
        config: Configuration
    """
    logger.info("\n[6/6] Generating summary report...")

    from src.reporting import SummaryReport
    from src.preprocessing_tracker import PreprocessingTracker
    from pathlib import Path
    import json

    # Read filtering report
    filtering_report_path = Path(config['paths']['results_dir']) / 'filtering_report.txt'
    filtering_report = ''
    if filtering_report_path.exists():
        with open(filtering_report_path, 'r') as f:
            filtering_report = f.read()

    # Phase 1.1: Generate and save preprocessing state
    logger.info("Generating preprocessing state for reproducibility...")
    tracker = PreprocessingTracker()

    # Extract preprocessing information from config and pipeline state
    tracker.mark_tic_normalized("tic")
    # Phase 11.1: Use constant instead of hardcoded value
    tracker.mark_log_transformed(
        pseudocount=LOG_TRANSFORM_PSEUDOCOUNT,
        base=2
    )
    tracker.mark_scaled("RobustScaler")

    # Set data statistics
    cancer_samples = [col for col in state.filtered_data.columns if col.startswith('C') and col[1:].isdigit()]
    normal_samples = [col for col in state.filtered_data.columns if col.startswith('N') and col[1:].isdigit()]
    tracker.set_sample_counts(len(cancer_samples), len(normal_samples))

    # Get glycopeptide counts from filtering report if available
    if hasattr(state, 'raw_data') and state.raw_data is not None:
        tracker.set_glycopeptide_counts(len(state.raw_data), len(state.filtered_data))
    else:
        # Only filtered count available (raw_data not stored in state)
        tracker.set_glycopeptide_counts(len(state.filtered_data))

    # Record filtering configuration (Phase 11.1: use constants as defaults)
    detection_config = config.get('analysis', {}).get('detection_filter', {})
    tracker.mark_filtered(
        min_detection_pct=detection_config.get('min_detection_pct', DETECTION_FILTER_DEFAULT_PCT),
        min_samples=detection_config.get('min_samples', DETECTION_FILTER_DEFAULT_SAMPLES)
    )

    # Save preprocessing state as JSON
    preprocessing_path = Path(config['paths']['results_dir']) / 'preprocessing_state.json'
    tracker.save(str(preprocessing_path))
    logger.info(f"Saved preprocessing state to {preprocessing_path}")

    # Generate summary report
    report = SummaryReport()
    summary_content = report.generate({
        'pipeline_state': state,
        'config': config,
        'filtering_report': filtering_report,
        'preprocessing_tracker': tracker  # Pass tracker for summary inclusion
    })

    # Save to file
    output_path = Path(config['paths']['results_dir']) / 'analysis_summary.txt'
    with open(output_path, 'w') as f:
        f.write(summary_content)

    logger.info(f"Saved summary report to {output_path}")


def main():
    """Main entry point"""
    logger.info("="*80)
    logger.info("pGlyco Auto Combine v3.0 - Workflow-Based Pipeline (Phase 11.1)")
    logger.info("="*80)

    overall_start_time = time.time()

    try:
        # Build and run pipeline (Phase 11.1: with execution timing)
        with timed_operation("Core Pipeline (Data + Analysis)"):
            state = (GlycoPipelineBuilder()
                    .with_config('config.yaml')
                    .with_logging()
                    .build_and_run())

        # Run visualizations with timing (Phase 11.1)
        with timed_operation("Visualization Generation", state.config):
            run_visualizations(state, state.config)

        # Generate summary report with timing (Phase 11.1)
        with timed_operation("Summary Report Generation", state.config):
            generate_summary_report(state, state.config)

        overall_duration = time.time() - overall_start_time

        logger.info("\n" + "="*80)
        logger.info("Pipeline completed successfully!")
        logger.info(f"Total execution time: {overall_duration:.2f}s ({overall_duration/60:.2f} minutes)")
        logger.info("="*80)

        print("\n" + "="*80)
        print("Analysis Complete!")
        print("="*80)
        print("\nResults saved to: Results/")
        print(f"\nTotal execution time: {overall_duration:.2f}s ({overall_duration/60:.2f} minutes)")
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
