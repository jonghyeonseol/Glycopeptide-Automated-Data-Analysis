#!/usr/bin/env python3
"""
pGlyco Auto Combine - Main Pipeline
Automated glycoproteomics data integration and analysis

ALCOA++ Compliance: This pipeline implements regulatory-grade data integrity
with full audit trails, metadata collection, and data verification.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

# Import refactored modules
from src.config_validator import load_and_validate_config
from src.logger_config import setup_logging, get_logger
from src.data_loader import DataLoader
from src.annotator import GlycanAnnotator
from src.analyzer import GlycanAnalyzer
from src.visualizer import GlycanVisualizer
from src.exceptions import PGlycoAutoError
from src.data_preparation import get_standard_config_from_dict
from src.data_pipeline import DataPipeline
from src.constants import (
    SEPARATOR_LINE, BOOTSTRAP_ITERATIONS, BOOTSTRAP_CONFIDENCE_LEVEL,
    STABILITY_THRESHOLD, COHENS_D_LARGE_EFFECT, VIP_THRESHOLD,
    CV_N_FOLDS, PERMUTATION_N_ITERATIONS,
    FILENAME_GLYCAN_STATS, FILENAME_VIP_SCORES,
    FILENAME_BOOTSTRAP_VALIDATION, FILENAME_STABLE_BIOMARKERS,
    FILENAME_CV_RESULTS, FILENAME_COHENS_D, FILENAME_PCA_PERMUTATION,
    FILENAME_ANALYSIS_SUMMARY, FILENAME_INPUT_MANIFEST,
    FILENAME_OUTPUT_MANIFEST, FILENAME_EXECUTION_METADATA
)

# ALCOA++ Compliance modules
from src.metadata_collector import get_metadata_collector
from src.audit_logger import get_audit_logger, EventType
from src.data_integrity import DataIntegrityManager

# Setup logging (single configuration for entire application)
setup_logging()
logger = get_logger(__name__)


# ==============================================================================
# HELPER FUNCTIONS FOR PIPELINE STAGES
# ==============================================================================

def initialize_compliance_systems(config):
    """
    Initialize ALCOA++ compliance systems and audit logging.

    Args:
        config: Configuration dictionary

    Returns:
        tuple: (metadata_collector, audit_logger, data_integrity, data_prep_config)
    """
    import os

    metadata_collector = get_metadata_collector()
    audit_logger = get_audit_logger()
    data_integrity = DataIntegrityManager()

    # Set researcher ID if available
    researcher_id = os.environ.get('PGLYCO_RESEARCHER_ID')
    if researcher_id:
        metadata_collector.set_researcher_id(researcher_id)
        logger.info(f"Researcher ID: {researcher_id}")

    # Initialize audit log
    results_dir_path = Path(config['paths']['results_dir'])
    results_dir_path.mkdir(parents=True, exist_ok=True)
    audit_logger.initialize(results_dir_path)

    audit_logger.log_event(
        EventType.CONFIG_LOADED,
        "Configuration loaded and validated",
        data={'config_file': 'config.yaml'}
    )

    # Create data preparation config
    data_prep_config = get_standard_config_from_dict(config)
    logger.info(
        f"Data preparation config loaded: "
        f"detection={data_prep_config.min_detection_pct * 100:.0f}%, "
        f"min_samples={data_prep_config.min_samples}, "
        f"method={data_prep_config.missing_data_method}"
    )

    return metadata_collector, audit_logger, data_integrity, data_prep_config


def load_and_integrate_dataset(config, audit_logger, data_integrity):
    """
    Load and integrate CSV files from dataset directory.

    Args:
        config: Configuration dictionary
        audit_logger: Audit logger instance
        data_integrity: Data integrity manager

    Returns:
        pd.DataFrame: Integrated dataset
    """
    logger.info("\n[2/7] Loading and integrating CSV files...")
    dataset_dir = config['paths']['dataset_dir']
    results_dir = config['paths']['results_dir']

    audit_logger.log_event(
        EventType.DATA_LOAD_START,
        f"Starting data loading from {dataset_dir}"
    )

    loader = DataLoader(
        dataset_dir=dataset_dir,
        required_columns=config['processing']['required_columns']
    )

    integrated_data = loader.integrate_data(
        qc_filters=config['processing']['qc_filters']
    )

    logger.info(f"Integrated data shape: {integrated_data.shape}")

    # Log data integration
    n_samples = len([col for col in integrated_data.columns if col.startswith(('C', 'N'))])
    audit_logger.log_data_integration(
        n_files=len(list(Path(dataset_dir).glob('*.csv'))),
        n_glycopeptides=len(integrated_data),
        n_samples=n_samples
    )

    # Create input data integrity manifest
    logger.info("\nCreating input data integrity manifest...")
    input_manifest_path = Path(results_dir) / FILENAME_INPUT_MANIFEST
    data_integrity.create_input_manifest(
        dataset_dir=Path(dataset_dir),
        output_path=input_manifest_path
    )

    return integrated_data


def annotate_dataset(integrated_data, config):
    """
    Annotate glycan compositions with modification information.

    Args:
        integrated_data: Integrated DataFrame
        config: Configuration dictionary

    Returns:
        pd.DataFrame: Annotated dataset (unfiltered)
    """
    logger.info("\n[3/7] Annotating glycan compositions...")

    annotator = GlycanAnnotator(
        sialylation_marker=config['annotation']['sialylation_marker'],
        fucosylation_marker=config['annotation']['fucosylation_marker']
    )

    annotated_data_raw = annotator.annotate_dataframe(integrated_data)
    logger.info(f"Annotated data shape: {annotated_data_raw.shape}")

    return annotated_data_raw


def apply_data_filter(annotated_data_raw, data_prep_config, config):
    """
    Apply detection frequency filter via DataPipeline.

    Args:
        annotated_data_raw: Unfiltered annotated DataFrame
        data_prep_config: Data preparation configuration
        config: Configuration dictionary

    Returns:
        tuple: (filtered_data, pipeline, filtering_report)
    """
    logger.info("\n[4/7] Applying detection filter via DataPipeline...")

    pipeline = DataPipeline(data_prep_config)
    annotated_data = pipeline.filter_dataset(annotated_data_raw)

    logger.info(f"Filtered data shape: {annotated_data.shape}")
    logger.info(f"Glycopeptides removed: {len(annotated_data_raw) - len(annotated_data)}")

    # Get filtering report
    filtering_report = pipeline.get_filtering_report()

    # Save datasets
    results_dir = Path(config['paths']['results_dir'])
    pipeline.save_datasets(
        raw_data=annotated_data_raw,
        filtered_data=annotated_data,
        results_dir=results_dir
    )

    # Validate filtering
    if not pipeline.validate_filtering(annotated_data):
        raise ValueError("Filtering validation failed!")

    return annotated_data, pipeline, filtering_report


def perform_analysis(annotated_data, config):
    """
    Perform statistical analysis (PCA, PLS-DA, statistics).

    Args:
        annotated_data: Filtered annotated DataFrame
        config: Configuration dictionary

    Returns:
        tuple: (analyzer, pca_results, plsda_results, statistics_df)
    """
    logger.info("\n[5/7] Performing statistical analysis...")

    analyzer = GlycanAnalyzer(
        log_transform=config['analysis']['pca']['log_transform'],
        n_components=config['analysis']['pca']['n_components']
    )

    # PCA
    logger.info("Performing PCA...")
    pca_results = analyzer.perform_pca(annotated_data)

    # PLS-DA
    logger.info("Performing PLS-DA...")
    plsda_results = analyzer.perform_plsda(annotated_data)

    # Statistics by glycan type
    logger.info("Calculating statistics by glycan type...")
    statistics_df = analyzer.calculate_glycan_type_statistics(annotated_data)

    # Save results
    results_dir = Path(config['paths']['results_dir'])
    stats_file = results_dir / FILENAME_GLYCAN_STATS
    statistics_df.to_csv(stats_file, index=False)
    logger.info(f"Saved statistics to {stats_file}")

    vip_file = results_dir / FILENAME_VIP_SCORES
    plsda_results['vip_scores'].to_csv(vip_file, index=False)
    logger.info(f"Saved VIP scores to {vip_file}")

    return analyzer, pca_results, plsda_results, statistics_df


def main():
    """Main pipeline execution with ALCOA++ compliance"""
    logger.info(SEPARATOR_LINE)
    logger.info("pGlyco Auto Combine - Glycoproteomics Data Analysis Pipeline")
    logger.info("ALCOA++ Compliant: Full audit trail and data integrity verification")
    logger.info(SEPARATOR_LINE)

    # [1/7] Load and validate configuration
    logger.info("\n[1/7] Loading and validating configuration...")
    try:
        config = load_and_validate_config('config.yaml')
        logger.info("Configuration validated successfully")
    except PGlycoAutoError as e:
        logger.error(f"Configuration error: {str(e)}")
        sys.exit(1)

    # Initialize compliance systems
    metadata_collector, audit_logger, data_integrity, data_prep_config = \
        initialize_compliance_systems(config)

    # Load and integrate data
    integrated_data = load_and_integrate_dataset(config, audit_logger, data_integrity)

    # Annotate data
    annotated_data_raw = annotate_dataset(integrated_data, config)

    # Log annotation statistics
    annotation_stats = {
        'sialylated': int(annotated_data_raw['IsSialylated'].sum()),
        'fucosylated': int(annotated_data_raw['IsFucosylated'].sum()),
        'glycan_types': annotated_data_raw['GlycanType'].value_counts().to_dict()
    }
    audit_logger.log_annotation(
        n_glycopeptides=len(annotated_data_raw),
        annotation_stats=annotation_stats
    )

    # Apply detection filter
    annotated_data, pipeline, filtering_report = apply_data_filter(
        annotated_data_raw, data_prep_config, config
    )

    # Prepare clean datasets for report
    sample_columns = [col for col in annotated_data.columns
                      if col.startswith('C') or col.startswith('N')]
    output_columns = (
        ['Peptide', 'GlycanComposition'] + sample_columns +
        ['Sialylation', 'Fucosylation', 'HighMannose', 'ComplexHybrid',
         'PrimaryClassification', 'SecondaryClassification',
         'GlycanTypeCategory', 'Proteins']
    )
    clean_raw = annotated_data_raw[output_columns].copy()
    clean_filtered = annotated_data[output_columns].copy()

    # Log filtering to audit trail
    n_removed = len(annotated_data_raw) - len(annotated_data)
    pct_removed = n_removed / len(annotated_data_raw) * 100 if len(annotated_data_raw) > 0 else 0
    audit_logger.log_filtering(
        n_before=len(annotated_data_raw),
        n_after=len(annotated_data),
        n_removed=n_removed,
        pct_removed=pct_removed,
        filter_criteria={
            'min_detection_pct': data_prep_config.min_detection_pct,
            'min_samples': data_prep_config.min_samples,
            'method': 'detection_frequency_30pct_OR_5samples'
        }
    )

    # Validate sample sizes
    from src.utils import validate_statistical_power
    cancer_samples = [col for col in annotated_data.columns if col.startswith('C') and col[1:].isdigit()]
    normal_samples = [col for col in annotated_data.columns if col.startswith('N') and col[1:].isdigit()]
    validate_statistical_power(cancer_samples, normal_samples, min_n=5)

    # Perform statistical analysis
    analyzer, pca_results, plsda_results, stats_df = perform_analysis(annotated_data, config)

    # Define file paths for later reference
    results_dir = config['paths']['results_dir']
    stats_file = Path(results_dir) / FILENAME_GLYCAN_STATS
    vip_file = Path(results_dir) / FILENAME_VIP_SCORES
    input_manifest_path = Path(results_dir) / FILENAME_INPUT_MANIFEST

    # Log to audit trail
    audit_logger.log_pca(
        n_components=config['analysis']['pca']['n_components'],
        variance_explained=pca_results['explained_variance'].tolist()
    )
    audit_logger.log_plsda(n_components=2, n_features=len(annotated_data))
    audit_logger.log_vip_scores(
        n_glycopeptides=len(plsda_results['vip_scores']),
        vip_threshold=1.0
    )

    # Prepare boxplot data
    logger.info("Preparing boxplot data...")
    boxplot_data = analyzer.prepare_boxplot_data(annotated_data)
    boxplot_data_extended = analyzer.prepare_boxplot_data_extended(annotated_data)

    # Get top VIP scores
    vip_glycopeptide = analyzer.get_top_vip_by_glycopeptide(annotated_data, plsda_results, top_n=10)

    # ==========================================================================
    # PHASE 2: PUBLICATION-QUALITY STATISTICAL VALIDATION
    # ==========================================================================
    logger.info("\n[PHASE 2] Running publication-quality statistical validation...")
    logger.info(SEPARATOR_LINE)

    from src.statistical_validation import StatisticalValidator
    validator = StatisticalValidator(random_state=42)

    # Get intensity matrix and labels for validation
    intensity_matrix, sample_names, feature_info = analyzer.prepare_intensity_matrix(annotated_data)
    from src.utils import get_sample_group
    from src.constants import GROUP_CANCER
    y_labels = np.array([1 if get_sample_group(name) == GROUP_CANCER else 0 for name in sample_names])
    feature_names = [f"{feat[0]}_{feat[1]}" for feat in feature_info]

    # 1. Bootstrap VIP Validation
    logger.info(f"\n1. Bootstrap VIP Validation ({BOOTSTRAP_ITERATIONS} iterations)...")
    bootstrap_results = validator.bootstrap_vip_validation(
        X=intensity_matrix.values,
        y=y_labels,
        feature_names=feature_names,
        n_iterations=BOOTSTRAP_ITERATIONS,
        n_components=2,
        confidence_level=BOOTSTRAP_CONFIDENCE_LEVEL
    )

    # Save bootstrap results
    bootstrap_file = Path(results_dir) / FILENAME_BOOTSTRAP_VALIDATION
    bootstrap_df = pd.DataFrame({
        'Peptide': [name.split('_')[0] for name in feature_names],
        'GlycanComposition': [name.split('_')[1] for name in feature_names],
        'VIP_Mean': bootstrap_results.vip_mean,
        'VIP_Std': bootstrap_results.vip_std,
        'VIP_CI_Lower': bootstrap_results.vip_ci_lower,
        'VIP_CI_Upper': bootstrap_results.vip_ci_upper,
        'Stability_Score': bootstrap_results.stability_score
    })
    bootstrap_df.to_csv(bootstrap_file, index=False)
    logger.info(f"✓ Bootstrap validation saved to {bootstrap_file}")

    # Get stable biomarkers
    stable_biomarkers = bootstrap_results.get_stable_biomarkers(
        stability_threshold=STABILITY_THRESHOLD,
        vip_threshold=VIP_THRESHOLD
    )
    stable_file = Path(results_dir) / FILENAME_STABLE_BIOMARKERS
    stable_biomarkers.to_csv(stable_file, index=False)
    logger.info(f"✓ Found {len(stable_biomarkers)} stable biomarkers (saved to {stable_file})")

    # 2. PLS-DA Cross-Validation
    logger.info(f"\n2. PLS-DA Cross-Validation ({CV_N_FOLDS}-fold)...")
    cv_results = validator.cross_validate_plsda(
        X=intensity_matrix.values,
        y=y_labels,
        n_components=2,
        n_folds=CV_N_FOLDS
    )

    # Save CV results
    cv_file = Path(results_dir) / FILENAME_CV_RESULTS
    with open(cv_file, 'w') as f:
        f.write(metadata_collector.get_metadata_header_lines())
        f.write("\n")
        f.write(cv_results.summary())
    logger.info(f"✓ Cross-validation saved to {cv_file}")

    # 3. Cohen's d Effect Size Calculations
    logger.info("\n3. Cohen's d Effect Size Calculations...")
    cancer_samples = [col for col in annotated_data.columns if col.startswith('C') and col[1:].isdigit()]
    normal_samples = [col for col in annotated_data.columns if col.startswith('N') and col[1:].isdigit()]

    cancer_matrix = intensity_matrix.loc[cancer_samples].values
    normal_matrix = intensity_matrix.loc[normal_samples].values

    effect_size_results = validator.calculate_cohens_d(
        group1=cancer_matrix,
        group2=normal_matrix,
        feature_names=feature_names
    )

    # Save effect sizes
    effect_size_file = Path(results_dir) / FILENAME_COHENS_D
    effect_size_df = pd.DataFrame({
        'Peptide': [name.split('_')[0] for name in feature_names],
        'GlycanComposition': [name.split('_')[1] for name in feature_names],
        'Cohens_d': effect_size_results.cohens_d,
        'Effect_Magnitude': effect_size_results.effect_magnitude,
        'Cancer_Mean': effect_size_results.group1_mean,
        'Normal_Mean': effect_size_results.group2_mean,
        'Pooled_Std': effect_size_results.pooled_std
    })
    effect_size_df.to_csv(effect_size_file, index=False)
    logger.info(f"✓ Effect sizes saved to {effect_size_file}")

    # Get large effect sizes (|d| >= 0.8)
    large_effects = effect_size_results.get_large_effects(threshold=COHENS_D_LARGE_EFFECT)
    large_effects_file = Path(results_dir) / 'large_effect_sizes.csv'
    large_effects.to_csv(large_effects_file, index=False)
    logger.info(f"✓ Found {len(large_effects)} features with large effect sizes")

    # 4. PCA Permutation Test
    logger.info(f"\n4. PCA Permutation Test ({PERMUTATION_N_ITERATIONS} permutations)...")
    pca_perm_results = validator.permutation_test_pca(
        X=intensity_matrix.values,
        y=y_labels,
        n_permutations=PERMUTATION_N_ITERATIONS
    )

    # Save permutation test results
    perm_file = Path(results_dir) / FILENAME_PCA_PERMUTATION
    with open(perm_file, 'w') as f:
        f.write(metadata_collector.get_metadata_header_lines())
        f.write("\n")
        f.write("PCA Permutation Test Results\n")
        f.write(SEPARATOR_LINE + "\n\n")
        f.write(f"Observed Separation: {pca_perm_results.observed_statistic:.4f}\n")
        f.write(f"P-value: {pca_perm_results.p_value:.4f}\n")
        f.write(f"Number of Permutations: {pca_perm_results.n_permutations}\n\n")
        if pca_perm_results.is_significant(alpha=0.05):
            f.write("✓ Result: SIGNIFICANT (p < 0.05)\n")
            f.write("The observed group separation in PCA space is significantly\n")
            f.write("greater than expected by chance.\n")
        else:
            f.write("✗ Result: NOT SIGNIFICANT (p >= 0.05)\n")
            f.write("The observed group separation could be due to chance.\n")
    logger.info(f"✓ Permutation test saved to {perm_file}")

    logger.info("\n" + SEPARATOR_LINE)
    logger.info("PHASE 2 Statistical Validation Complete!")
    logger.info(SEPARATOR_LINE)

    # Step 5: Create visualizations
    logger.info("\n[6/7] Creating visualizations...")
    audit_logger.log_event(
        EventType.VISUALIZATION_START,
        "Starting visualization generation"
    )

    visualizer = GlycanVisualizer(
        output_dir=results_dir,
        dpi=config['visualization']['dpi'],
        colors=config['visualization']['colors']
    )

    visualizer.plot_all(
        df=annotated_data,
        pca_results=pca_results,
        boxplot_data=boxplot_data,
        boxplot_data_extended=boxplot_data_extended
    )

    # Additional visualizations for primary/secondary classification
    logger.info("Creating primary/secondary classification histograms...")
    visualizer.plot_histogram_primary_classification(annotated_data, normalization='raw')
    visualizer.plot_histogram_primary_classification(annotated_data, normalization='aggregated')
    visualizer.plot_histogram_secondary_classification(annotated_data, normalization='raw')
    visualizer.plot_histogram_secondary_classification(annotated_data, normalization='aggregated')
    visualizer.plot_histogram_cancer_vs_normal_primary(annotated_data)
    visualizer.plot_histogram_cancer_vs_normal_secondary(annotated_data)

    # VIP score visualizations (with R/ggplot2)
    logger.info("Creating VIP score plots with R/ggplot2...")
    visualizer.plot_vip_scores_glycopeptide_r(annotated_data, plsda_results['vip_scores'])
    visualizer.plot_vip_scores_glycan_composition_r(annotated_data, plsda_results['vip_scores'])
    visualizer.plot_vip_scores_peptide_r(annotated_data, plsda_results['vip_scores'])
    visualizer.plot_vip_scores_peptide_grouped_r(annotated_data, plsda_results['vip_scores'])

    # Box plots corresponding to histograms
    logger.info("Creating box plots...")
    visualizer.plot_boxplot_primary_classification(annotated_data, normalization='raw')
    visualizer.plot_boxplot_primary_classification(annotated_data, normalization='aggregated')
    visualizer.plot_boxplot_secondary_classification(annotated_data, normalization='raw')
    # Removed: boxplot_secondary_aggregated_normalized (user request - not useful)
    # visualizer.plot_boxplot_secondary_classification(annotated_data, normalization='aggregated')
    visualizer.plot_boxplot_cancer_vs_normal_primary(annotated_data)
    visualizer.plot_boxplot_cancer_vs_normal_secondary(annotated_data)

    # Advanced visualizations (evidence-based from literature)
    logger.info("Creating advanced visualizations...")
    logger.info("  - Volcano plot (differential expression)...")
    volcano_data = visualizer.plot_volcano(annotated_data, plsda_results['vip_scores'], config=data_prep_config)

    logger.info("  - Site-specific glycosylation heatmap...")
    visualizer.plot_site_specific_heatmap(annotated_data, plsda_results['vip_scores'])

    logger.info("  - CV distribution plots...")
    visualizer.plot_cv_distribution(annotated_data)

    logger.info("  - Sample correlation matrices...")
    visualizer.plot_correlation_matrix(annotated_data)
    visualizer.plot_correlation_clustermap(annotated_data)

    logger.info("  - Venn diagram (glycan type overlap)...")
    visualizer.plot_glycan_venn_diagram(annotated_data)

    logger.info("  - Radar chart (glycan profile comparison)...")
    visualizer.plot_radar_chart(annotated_data)

    logger.info("  - Pie charts (glycan distribution)...")
    visualizer.plot_pie_chart_glycan_types(annotated_data)
    visualizer.plot_pie_chart_primary_classification(annotated_data)
    visualizer.plot_pie_chart_secondary_classification(annotated_data)

    logger.info("  - Glycopeptide comparison heatmap (Cancer vs Normal)...")
    if config['visualization']['glycopeptide_comparison']['enabled']:
        visualizer.plot_glycopeptide_comparison_heatmap(
            df=annotated_data,
            vip_scores=plsda_results['vip_scores'],
            config=data_prep_config,  # Pass standardized config
            figsize=tuple(config['visualization']['glycopeptide_comparison']['figsize']),
            max_peptides=config['visualization']['glycopeptide_comparison']['max_peptides'],
            max_glycans_per_type=config['visualization']['glycopeptide_comparison']['max_glycans_per_type']
        )

    # Phase 4.1: Missing data visualization (data integrity)
    logger.info("  - Missing data matrix (data integrity validation)...")
    visualizer.plot_missing_data_matrix(annotated_data_raw)  # Use RAW data to show full picture

    # Phase 4.2: PLS-DA diagnostic plots (model validation)
    logger.info("  - PLS-DA diagnostic plots (model validation)...")
    visualizer.plot_plsda_diagnostics(plsda_results, annotated_data)

    # Phase 4.3: Per-sample QC dashboard
    logger.info("  - Per-sample QC dashboard...")
    visualizer.plot_sample_qc_dashboard(annotated_data)

    # Mark visualizations complete
    audit_logger.log_event(
        EventType.VISUALIZATION_COMPLETE,
        "All visualizations generated successfully"
    )

    # ==========================================================================
    # PHASE 2.1: INTERACTIVE PLOTLY DASHBOARD
    # ==========================================================================
    logger.info("\n[PHASE 2.1] Generating interactive Plotly visualizations...")
    logger.info(SEPARATOR_LINE)

    from src.interactive_dashboard import InteractiveDashboard
    try:
        dashboard = InteractiveDashboard(
            output_dir=results_dir,
            colors=config['visualization']['colors']
        )

        # Generate all interactive plots with embedded static PNGs (Phase 2.3)
        dashboard.generate_all_interactive_plots(
            pca_results=pca_results,
            df=annotated_data,
            vip_df=plsda_results['vip_scores'],
            volcano_df=volcano_data,
            results_dir=Path(results_dir)  # For base64 embedding
        )

        logger.info("✓ Interactive dashboard generated successfully")
        logger.info(f"✓ Open {Path(results_dir) / 'interactive' / 'index.html'} to view dashboard")
    except Exception as e:
        logger.warning(f"Interactive dashboard generation encountered an issue: {e}")
        logger.warning("Continuing with pipeline completion...")

    logger.info(SEPARATOR_LINE)

    # Step 6: Summary report
    logger.info("\n[7/7] Generating summary report...")

    summary_file = Path(results_dir) / FILENAME_ANALYSIS_SUMMARY
    with open(summary_file, 'w') as f:
        # Add metadata header
        f.write(metadata_collector.get_metadata_header_lines())
        f.write("\n")
        f.write(SEPARATOR_LINE + "\n")
        f.write("pGlyco Auto Combine - Analysis Summary\n")
        f.write(SEPARATOR_LINE + "\n\n")

        # Add filtering report
        f.write(pipeline.get_filtering_report())
        f.write("\n")

        f.write("Data Integration:\n")
        f.write(f"  - Total glycopeptides: {len(annotated_data)}\n")
        f.write(f"  - Total samples: {len([col for col in annotated_data.columns if col.startswith(('C', 'N'))])}\n")
        f.write(f"  - Cancer samples: {len([col for col in annotated_data.columns if col.startswith('C')])}\n")
        f.write(f"  - Normal samples: {len([col for col in annotated_data.columns if col.startswith('N')])}\n\n")

        f.write("Glycan Annotation:\n")
        sialylated_count = annotated_data['IsSialylated'].sum()
        sialylated_pct = sialylated_count / len(annotated_data) * 100
        f.write(f"  - Sialylated: {sialylated_count} ({sialylated_pct:.1f}%)\n")
        fucosylated_count = annotated_data['IsFucosylated'].sum()
        fucosylated_pct = fucosylated_count / len(annotated_data) * 100
        f.write(f"  - Fucosylated: {fucosylated_count} ({fucosylated_pct:.1f}%)\n\n")

        f.write("Glycan Type Category Distribution (for Comparison Heatmap):\n")
        for glycan_cat, count in annotated_data['GlycanTypeCategory'].value_counts().items():
            f.write(f"  - {glycan_cat}: {count} ({count / len(annotated_data) * 100:.1f}%)\n")
        f.write("\n")

        f.write("Primary Classification Distribution:\n")
        for primary_class, count in annotated_data['PrimaryClassification'].value_counts().items():
            f.write(f"  - {primary_class}: {count} ({count / len(annotated_data) * 100:.1f}%)\n")

        f.write("\nSecondary Classification Distribution:\n")
        for secondary_class, count in annotated_data['SecondaryClassification'].value_counts().items():
            f.write(f"  - {secondary_class}: {count} ({count / len(annotated_data) * 100:.1f}%)\n")

        f.write("\nGlycan Type Distribution (Legacy):\n")
        for glycan_type, count in annotated_data['GlycanType'].value_counts().items():
            f.write(f"  - {glycan_type}: {count} ({count / len(annotated_data) * 100:.1f}%)\n")

        f.write("\nPCA Results:\n")
        f.write(f"  - PC1 explained variance: {pca_results['explained_variance'][0] * 100:.2f}%\n")
        f.write(f"  - PC2 explained variance: {pca_results['explained_variance'][1] * 100:.2f}%\n")
        f.write(f"  - Total explained variance: {pca_results['explained_variance'].sum() * 100:.2f}%\n")

        f.write("\nStatistics by Glycan Type:\n")
        f.write(stats_df.to_string(index=False))

        f.write("\n\nTop 10 VIP Scores (Glycopeptide):\n")
        f.write(vip_glycopeptide[['Peptide', 'GlycanComposition', 'VIP_Score']].to_string(index=False))
        f.write("\n\nNote: VIP scores by Glycan Type and Peptide are visualized in the corresponding PNG files.")

        f.write("\n\n" + SEPARATOR_LINE + "\n")
        f.write("Output Files:\n")
        f.write(f"  - Integrated data (RAW): integrated.csv ({len(clean_raw)} glycopeptides)\n")
        f.write(
            f"  - Integrated data (FILTERED): integrated_filtered.csv "
            f"({len(clean_filtered)} glycopeptides - USED IN ALL ANALYSES)\n"
        )
        f.write("  - Filtering report: filtering_report.txt\n")
        f.write(f"  - Statistics: {stats_file}\n")
        f.write(f"  - VIP Scores (all): {vip_file}\n")
        f.write("\nVisualization Files:\n")
        f.write("  - PCA plot: pca_plot.png\n")
        f.write("  - PCA samples: pca_samples.png\n")
        f.write("  - Boxplot: boxplot_glycan_types.png\n")
        f.write("  - Boxplot (extended): boxplot_extended_categories.png\n")
        f.write("  - Heatmap (Top 50): heatmap_top_glycopeptides.png\n")
        f.write("  - Heatmap (Full Profile): heatmap_full_glycan_profile.png\n")
        f.write("  - Distribution: glycan_type_distribution.png\n")
        f.write("  - Histogram (normalized): histogram_glycan_types_by_sample_normalized.png\n")
        f.write("  - Histogram (primary, raw norm): histogram_primary_raw_normalized.png\n")
        f.write("  - Histogram (primary, agg norm): histogram_primary_aggregated_normalized.png\n")
        f.write("  - Histogram (secondary, raw norm): histogram_secondary_raw_normalized.png\n")
        f.write("  - Histogram (secondary, agg norm): histogram_secondary_aggregated_normalized.png\n")
        f.write("  - Histogram (Cancer vs Normal, primary): histogram_primary_cancer_vs_normal.png\n")
        f.write("  - Histogram (Cancer vs Normal, secondary): histogram_secondary_cancer_vs_normal.png\n")
        f.write("  - Boxplot (primary, raw norm): boxplot_primary_raw_normalized.png\n")
        f.write("  - Boxplot (primary, agg norm): boxplot_primary_aggregated_normalized.png\n")
        f.write("  - Boxplot (secondary, raw norm): boxplot_secondary_raw_normalized.png\n")
        f.write("  - Boxplot (secondary, agg norm): boxplot_secondary_aggregated_normalized.png\n")
        f.write("  - Boxplot (Cancer vs Normal, primary): boxplot_primary_cancer_vs_normal.png\n")
        f.write("  - Boxplot (Cancer vs Normal, secondary): boxplot_secondary_cancer_vs_normal.png\n")
        f.write("  - VIP scores (glycopeptide, R/ggplot2): vip_score_glycopeptide_r.png\n")
        f.write("  - VIP scores (glycan composition, R/ggplot2): vip_score_glycan_composition_r.png\n")
        f.write("  - VIP scores (peptide, R/ggplot2): vip_score_peptide_r.png\n")
        f.write("  - VIP scores (peptide grouped, R/ggplot2): vip_score_peptide_grouped_r.png\n")
        f.write("\nAdvanced Visualizations (Evidence-based):\n")
        f.write("  - Volcano plot: volcano_plot.png\n")
        f.write("  - Site-specific glycosylation: site_specific_glycosylation_heatmap.png\n")
        f.write("  - CV distribution: cv_distribution.png\n")
        f.write("  - Correlation matrix (Cancer): correlation_matrix_cancer.png\n")
        f.write("  - Correlation matrix (Normal): correlation_matrix_normal.png\n")
        f.write("  - Correlation clustermap (Cancer): correlation_clustermap_cancer.png\n")
        f.write("  - Correlation clustermap (Normal): correlation_clustermap_normal.png\n")
        f.write("  - Venn diagram: venn_diagram_glycan_types.png\n")
        f.write("  - Radar chart: radar_chart_glycan_profile.png\n")
        f.write("  - Pie chart (glycan types): pie_chart_glycan_types.png\n")
        f.write("  - Pie chart (primary classification): pie_chart_primary_classification.png\n")
        f.write("  - Pie chart (secondary classification): pie_chart_secondary_classification.png\n")
        f.write("  - Glycopeptide comparison heatmap (Cancer vs Normal): glycopeptide_comparison_heatmap.png\n")
        f.write("\nPhase 4: Data Integrity & Model Validation:\n")
        f.write("  - Missing data matrix: missing_data_matrix.png\n")
        f.write("  - PLS-DA diagnostics: plsda_diagnostics.png\n")
        f.write("  - Sample QC dashboard: sample_qc_dashboard.png\n")
        f.write("\n")
        f.write("Data Traceability:\n")
        f.write("  All visualization source data is exported to Results/Trace/ folder as CSV files.\n")
        f.write("  Each PNG visualization has a corresponding *_data.csv file containing the exact\n")
        f.write("  data used to generate that plot, ensuring full reproducibility and transparency.\n")
        f.write(SEPARATOR_LINE + "\n")

    logger.info(f"Saved summary report to {summary_file}")

    # ==========================================================================
    # PHASE 2.3: PUBLICATION-READY REPORTING MATERIALS
    # ==========================================================================
    logger.info("\n[PHASE 2.3] Generating publication-ready reporting materials...")
    logger.info(SEPARATOR_LINE)

    from src.publication_report import generate_publication_report
    try:
        generate_publication_report(Path(results_dir))
    except Exception as e:
        logger.warning(f"Publication report generation encountered an issue: {e}")
        logger.warning("Continuing with pipeline completion...")

    # Create output data integrity manifest
    logger.info("\nCreating output data integrity manifest...")
    output_manifest_path = Path(results_dir) / FILENAME_OUTPUT_MANIFEST
    data_integrity.create_output_manifest(
        results_dir=Path(results_dir),
        output_path=output_manifest_path
    )

    # Save execution metadata
    logger.info("Saving execution metadata...")
    metadata_file = Path(results_dir) / FILENAME_EXECUTION_METADATA
    metadata_collector.save_metadata_json(metadata_file)

    # Finalize audit log
    logger.info("Finalizing audit trail...")
    audit_logger.finalize()

    # Print audit summary
    audit_summary = audit_logger.get_audit_summary()
    logger.info("\nAudit Trail Summary:")
    logger.info(f"  Total events recorded: {audit_summary['total_events']}")
    logger.info(f"  Audit log: {audit_summary['audit_file']}")

    logger.info("\n" + SEPARATOR_LINE)
    logger.info("Pipeline completed successfully with ALCOA++ compliance!")
    logger.info(SEPARATOR_LINE)
    logger.info("\nALCOA++ Compliance Outputs:")
    logger.info(f"  - Audit log: {audit_summary['audit_file']}")
    logger.info(f"  - Execution metadata: {metadata_file}")
    logger.info(f"  - Input manifest: {input_manifest_path}")
    logger.info(f"  - Output manifest: {output_manifest_path}")
    logger.info(SEPARATOR_LINE + "\n")

    print("\n" + SEPARATOR_LINE)
    print("Analysis Complete!")
    print(SEPARATOR_LINE)
    print(f"\nResults saved to: {results_dir}/")
    print("\nKey outputs:")
    print("  1. integrated_filtered.csv - Integrated and annotated data (FILTERED)")
    print("  2. analysis_summary.txt - Comprehensive summary")
    print("  3. Multiple visualization plots (PNG files)")
    print("\n" + SEPARATOR_LINE + "\n")


if __name__ == "__main__":
    try:
        main()
    except PGlycoAutoError as e:
        # Custom exceptions have user-friendly messages
        logger.error(f"Pipeline failed: {str(e)}")
        sys.exit(1)
    except KeyboardInterrupt:
        logger.warning("\nPipeline interrupted by user")
        sys.exit(130)
    except Exception as e:
        # Unexpected errors - show full traceback
        logger.error(f"Unexpected error: {str(e)}", exc_info=True)
        sys.exit(1)
