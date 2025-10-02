#!/usr/bin/env python3
"""
pGlyco Auto Combine - Main Pipeline
Automated glycoproteomics data integration and analysis
"""

import yaml
import sys
import logging
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from data_loader import DataLoader
from annotator import GlycanAnnotator
from analyzer import GlycanAnalyzer
from visualizer import GlycanVisualizer

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_config(config_path: str = 'config.yaml') -> dict:
    """Load configuration from YAML file"""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def main():
    """Main pipeline execution"""
    logger.info("="*80)
    logger.info("pGlyco Auto Combine - Glycoproteomics Data Analysis Pipeline")
    logger.info("="*80)

    # Load configuration
    logger.info("\n[1/6] Loading configuration...")
    config = load_config('config.yaml')

    # Paths
    dataset_dir = config['paths']['dataset_dir']
    results_dir = config['paths']['results_dir']
    output_file = Path(results_dir) / config['paths']['output_file']

    # Step 1: Load and integrate data
    logger.info("\n[2/6] Loading and integrating CSV files...")
    loader = DataLoader(
        dataset_dir=dataset_dir,
        required_columns=config['processing']['required_columns']
    )

    integrated_data = loader.integrate_data(
        qc_filters=config['processing']['qc_filters']
    )

    logger.info(f"Integrated data shape: {integrated_data.shape}")
    logger.info(f"Columns: {list(integrated_data.columns)}")

    # Step 2: Annotate data
    logger.info("\n[3/6] Annotating glycan compositions...")
    annotator = GlycanAnnotator(
        sialylation_marker=config['annotation']['sialylation_marker'],
        fucosylation_marker=config['annotation']['fucosylation_marker']
    )

    annotated_data = annotator.annotate_dataframe(integrated_data)

    # Prepare clean integrated data for output (only essential columns)
    # Keep: Peptide, GlycanComposition, all sample columns (C1-C24, N1-N24), Sialylation, Fucosylation, HighMannose, ComplexHybrid, PrimaryClassification, SecondaryClassification
    sample_columns = [col for col in annotated_data.columns
                     if col.startswith('C') or col.startswith('N')]

    output_columns = ['Peptide', 'GlycanComposition'] + sample_columns + ['Sialylation', 'Fucosylation', 'HighMannose', 'ComplexHybrid', 'PrimaryClassification', 'SecondaryClassification']
    clean_integrated_data = annotated_data[output_columns].copy()

    # Save clean integrated data
    logger.info(f"\nSaving integrated data to {output_file}...")
    loader.save_integrated_data(clean_integrated_data, output_file)

    # Step 3: Perform statistical analysis
    logger.info("\n[4/6] Performing statistical analysis...")
    analyzer = GlycanAnalyzer(
        n_components=config['analysis']['pca']['n_components'],
        log_transform=config['analysis']['pca']['log_transform']
    )

    # PCA analysis
    logger.info("Running PCA...")
    pca_results = analyzer.perform_pca(annotated_data)

    # Calculate statistics by glycan type
    logger.info("Calculating statistics by glycan type...")
    stats_df = analyzer.calculate_statistics_by_glycan_type(annotated_data)

    # Save statistics
    stats_file = Path(results_dir) / 'glycan_type_statistics.csv'
    stats_df.to_csv(stats_file, index=False)
    logger.info(f"Saved statistics to {stats_file}")

    # Prepare boxplot data
    logger.info("Preparing boxplot data...")
    boxplot_data = analyzer.prepare_boxplot_data(annotated_data)

    logger.info("Preparing extended boxplot data...")
    boxplot_data_extended = analyzer.prepare_boxplot_data_extended(annotated_data)

    # PLS-DA analysis and VIP score calculation
    logger.info("Performing PLS-DA analysis...")
    plsda_results = analyzer.perform_plsda(annotated_data, n_components=2)

    # Get VIP scores for summary report (TOP 10 only)
    logger.info("Calculating VIP scores...")
    vip_glycopeptide = analyzer.get_top_vip_by_glycopeptide(annotated_data, plsda_results, top_n=10)

    # Save all VIP scores
    vip_file = Path(results_dir) / 'vip_scores_all.csv'
    plsda_results['vip_scores'].to_csv(vip_file, index=False)
    logger.info(f"Saved all VIP scores to {vip_file}")

    # Step 4: Create visualizations
    logger.info("\n[5/6] Creating visualizations...")
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

    # VIP score visualizations (dot plots with Cancer/Normal coloring)
    logger.info("Creating VIP score plots...")
    visualizer.plot_vip_scores_glycopeptide(annotated_data, plsda_results['vip_scores'])
    visualizer.plot_vip_scores_glycan_type(annotated_data, plsda_results['vip_scores'])
    visualizer.plot_vip_scores_peptide_top10(annotated_data, plsda_results['vip_scores'])
    visualizer.plot_vip_scores_peptide_all(annotated_data, plsda_results['vip_scores'])

    # Box plots corresponding to histograms
    logger.info("Creating box plots...")
    visualizer.plot_boxplot_primary_classification(annotated_data, normalization='raw')
    visualizer.plot_boxplot_primary_classification(annotated_data, normalization='aggregated')
    visualizer.plot_boxplot_secondary_classification(annotated_data, normalization='raw')
    visualizer.plot_boxplot_secondary_classification(annotated_data, normalization='aggregated')
    visualizer.plot_boxplot_cancer_vs_normal_primary(annotated_data)
    visualizer.plot_boxplot_cancer_vs_normal_secondary(annotated_data)

    # Step 5: Summary report
    logger.info("\n[6/6] Generating summary report...")

    summary_file = Path(results_dir) / 'analysis_summary.txt'
    with open(summary_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("pGlyco Auto Combine - Analysis Summary\n")
        f.write("="*80 + "\n\n")

        f.write("Data Integration:\n")
        f.write(f"  - Total glycopeptides: {len(annotated_data)}\n")
        f.write(f"  - Total samples: {len([col for col in annotated_data.columns if col.startswith(('C', 'N'))])}\n")
        f.write(f"  - Cancer samples: {len([col for col in annotated_data.columns if col.startswith('C')])}\n")
        f.write(f"  - Normal samples: {len([col for col in annotated_data.columns if col.startswith('N')])}\n\n")

        f.write("Glycan Annotation:\n")
        f.write(f"  - Sialylated: {annotated_data['IsSialylated'].sum()} ({annotated_data['IsSialylated'].sum()/len(annotated_data)*100:.1f}%)\n")
        f.write(f"  - Fucosylated: {annotated_data['IsFucosylated'].sum()} ({annotated_data['IsFucosylated'].sum()/len(annotated_data)*100:.1f}%)\n\n")

        f.write("Primary Classification Distribution:\n")
        for primary_class, count in annotated_data['PrimaryClassification'].value_counts().items():
            f.write(f"  - {primary_class}: {count} ({count/len(annotated_data)*100:.1f}%)\n")

        f.write("\nSecondary Classification Distribution:\n")
        for secondary_class, count in annotated_data['SecondaryClassification'].value_counts().items():
            f.write(f"  - {secondary_class}: {count} ({count/len(annotated_data)*100:.1f}%)\n")

        f.write("\nGlycan Type Distribution (Legacy):\n")
        for glycan_type, count in annotated_data['GlycanType'].value_counts().items():
            f.write(f"  - {glycan_type}: {count} ({count/len(annotated_data)*100:.1f}%)\n")

        f.write("\nPCA Results:\n")
        f.write(f"  - PC1 explained variance: {pca_results['explained_variance'][0]*100:.2f}%\n")
        f.write(f"  - PC2 explained variance: {pca_results['explained_variance'][1]*100:.2f}%\n")
        f.write(f"  - Total explained variance: {pca_results['explained_variance'].sum()*100:.2f}%\n")

        f.write("\nStatistics by Glycan Type:\n")
        f.write(stats_df.to_string(index=False))

        f.write("\n\nTop 10 VIP Scores (Glycopeptide):\n")
        f.write(vip_glycopeptide[['Peptide', 'GlycanComposition', 'VIP_Score']].to_string(index=False))
        f.write("\n\nNote: VIP scores by Glycan Type and Peptide are visualized in the corresponding PNG files.")

        f.write("\n\n" + "="*80 + "\n")
        f.write("Output Files:\n")
        f.write(f"  - Integrated data: {output_file}\n")
        f.write(f"  - Statistics: {stats_file}\n")
        f.write(f"  - VIP Scores (all): {vip_file}\n")
        f.write(f"\nVisualization Files:\n")
        f.write(f"  - PCA plot: pca_plot.png\n")
        f.write(f"  - PCA samples: pca_samples.png\n")
        f.write(f"  - Boxplot: boxplot_glycan_types.png\n")
        f.write(f"  - Boxplot (extended): boxplot_extended_categories.png\n")
        f.write(f"  - Heatmap (Top 50): heatmap_top_glycopeptides.png\n")
        f.write(f"  - Heatmap (Full Profile): heatmap_full_glycan_profile.png\n")
        f.write(f"  - Distribution: glycan_type_distribution.png\n")
        f.write(f"  - Histogram (original): histogram_glycan_types_by_sample.png\n")
        f.write(f"  - Histogram (normalized): histogram_glycan_types_by_sample_normalized.png\n")
        f.write(f"  - Histogram (primary, raw norm): histogram_primary_raw_normalized.png\n")
        f.write(f"  - Histogram (primary, agg norm): histogram_primary_aggregated_normalized.png\n")
        f.write(f"  - Histogram (secondary, raw norm): histogram_secondary_raw_normalized.png\n")
        f.write(f"  - Histogram (secondary, agg norm): histogram_secondary_aggregated_normalized.png\n")
        f.write(f"  - Histogram (Cancer vs Normal, primary): histogram_primary_cancer_vs_normal.png\n")
        f.write(f"  - Histogram (Cancer vs Normal, secondary): histogram_secondary_cancer_vs_normal.png\n")
        f.write(f"  - Boxplot (primary, raw norm): boxplot_primary_raw_normalized.png\n")
        f.write(f"  - Boxplot (primary, agg norm): boxplot_primary_aggregated_normalized.png\n")
        f.write(f"  - Boxplot (secondary, raw norm): boxplot_secondary_raw_normalized.png\n")
        f.write(f"  - Boxplot (secondary, agg norm): boxplot_secondary_aggregated_normalized.png\n")
        f.write(f"  - Boxplot (Cancer vs Normal, primary): boxplot_primary_cancer_vs_normal.png\n")
        f.write(f"  - Boxplot (Cancer vs Normal, secondary): boxplot_secondary_cancer_vs_normal.png\n")
        f.write(f"  - VIP scores (glycopeptide): vip_score_glycopeptide.png\n")
        f.write(f"  - VIP scores (glycan type): vip_score_glycan_type.png\n")
        f.write(f"  - VIP scores (peptide, top 10): vip_score_peptide_top10.png\n")
        f.write(f"  - VIP scores (peptide, all): vip_score_peptide_all.png\n")
        f.write("="*80 + "\n")

    logger.info(f"Saved summary report to {summary_file}")

    logger.info("\n" + "="*80)
    logger.info("Pipeline completed successfully!")
    logger.info("="*80 + "\n")

    print("\n" + "="*80)
    print("Analysis Complete!")
    print("="*80)
    print(f"\nResults saved to: {results_dir}/")
    print(f"\nKey outputs:")
    print(f"  1. {output_file.name} - Integrated and annotated data")
    print(f"  2. analysis_summary.txt - Comprehensive summary")
    print(f"  3. Multiple visualization plots (PNG files)")
    print("\n" + "="*80 + "\n")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.error(f"Pipeline failed with error: {str(e)}", exc_info=True)
        sys.exit(1)
