"""
Glycoproteomics Pipeline Implementation
Main pipeline for pGlyco Auto Combine
"""

from pathlib import Path
from typing import Dict, Any

from .base_pipeline import BasePipeline, PipelineState
from .workflow import Workflow, WorkflowStep

# Import existing modules
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.data_loader import DataLoader
from src.annotator import GlycanAnnotator
from src.data_preparation import get_standard_config_from_dict
from src.data_pipeline import DataPipeline
from src.utils import validate_statistical_power
from src.logger_config import get_logger

# NEW: Use focused analyzers instead of monolithic GlycanAnalyzer
from src.analysis import PCAAnalyzer, PLSDAAnalyzer, StatisticsAnalyzer

logger = get_logger(__name__)


# ==================== Data Ingestion Workflow ====================

class LoadDataStep(WorkflowStep):
    """Load and integrate CSV files"""

    def __init__(self):
        super().__init__("Load Data", "Load and integrate CSV files from Dataset/")

    def execute(self, state: PipelineState, config: Dict[str, Any]) -> None:
        dataset_dir = config['paths']['dataset_dir']
        required_columns = config['processing']['required_columns']

        loader = DataLoader(
            dataset_dir=dataset_dir,
            required_columns=required_columns
        )

        integrated_data = loader.integrate_data(
            qc_filters=config['processing']['qc_filters']
        )

        state.integrated_data = integrated_data
        logger.info(f"Integrated data shape: {integrated_data.shape}")


# ==================== Annotation Workflow ====================

class AnnotateDataStep(WorkflowStep):
    """Annotate glycan compositions"""

    def __init__(self):
        super().__init__("Annotate Data", "Classify glycans by composition")

    def execute(self, state: PipelineState, config: Dict[str, Any]) -> None:
        annotator = GlycanAnnotator(
            sialylation_marker=config['annotation']['sialylation_marker'],
            fucosylation_marker=config['annotation']['fucosylation_marker']
        )

        state.annotated_data_raw = annotator.annotate_dataframe(state.integrated_data)
        logger.info(f"Annotated {len(state.annotated_data_raw)} glycopeptides")


# ==================== Filtering Workflow ====================

class FilterDataStep(WorkflowStep):
    """Apply detection frequency filter"""

    def __init__(self):
        super().__init__("Filter Data", "Apply detection filter (single source of truth)")

    def execute(self, state: PipelineState, config: Dict[str, Any]) -> None:
        data_prep_config = get_standard_config_from_dict(config)

        pipeline = DataPipeline(data_prep_config)
        state.filtered_data = pipeline.filter_dataset(state.annotated_data_raw)

        # Validate filtering
        pipeline.validate_filtering(state.annotated_data_raw, state.filtered_data)

        # Save datasets
        results_dir = config['paths']['results_dir']
        sample_columns = [col for col in state.filtered_data.columns
                          if col.startswith('C') or col.startswith('N')]

        output_columns = ['Peptide', 'GlycanComposition'] + sample_columns + [
            'Sialylation', 'Fucosylation', 'HighMannose', 'ComplexHybrid',
            'PrimaryClassification', 'SecondaryClassification',
            'GlycanTypeCategory', 'Proteins'
        ]

        clean_raw = state.annotated_data_raw[output_columns].copy()
        clean_filtered = state.filtered_data[output_columns].copy()

        pipeline.save_datasets(
            clean_raw,
            clean_filtered,
            results_dir,
            raw_filename='integrated.csv',
            filtered_filename='integrated_filtered.csv'
        )

        state.data_prep_config = data_prep_config


class ValidateStatisticalPowerStep(WorkflowStep):
    """Validate sample sizes for statistical tests"""

    def __init__(self):
        super().__init__("Validate Power", "Validate statistical power")

    def execute(self, state: PipelineState, config: Dict[str, Any]) -> None:
        cancer_samples = [col for col in state.filtered_data.columns
                          if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in state.filtered_data.columns
                          if col.startswith('N') and col[1:].isdigit()]

        validate_statistical_power(cancer_samples, normal_samples, min_n=5)


# ==================== Analysis Workflow ====================

class PCAAnalysisStep(WorkflowStep):
    """Perform PCA analysis"""

    def __init__(self):
        super().__init__("PCA Analysis", "Principal Component Analysis")

    def execute(self, state: PipelineState, config: Dict[str, Any]) -> None:
        # Use focused PCA analyzer
        pca_analyzer = PCAAnalyzer(
            n_components=config['analysis']['pca']['n_components'],
            log_transform=config['analysis']['pca']['log_transform']
        )

        state.analysis_results['pca'] = pca_analyzer.perform_pca(state.filtered_data)
        state.pca_analyzer = pca_analyzer


class StatisticsStep(WorkflowStep):
    """Calculate statistics by glycan type"""

    def __init__(self):
        super().__init__("Statistics", "Calculate statistics by glycan type")

    def execute(self, state: PipelineState, config: Dict[str, Any]) -> None:
        # Use focused statistics analyzer
        stats_analyzer = StatisticsAnalyzer(
            log_transform=config['analysis']['pca']['log_transform']
        )

        stats_df = stats_analyzer.calculate_statistics_by_glycan_type(state.filtered_data)
        state.analysis_results['statistics'] = stats_df
        state.stats_analyzer = stats_analyzer

        # Save statistics
        results_dir = config['paths']['results_dir']
        stats_file = Path(results_dir) / 'glycan_type_statistics.csv'
        stats_df.to_csv(stats_file, index=False)
        logger.info(f"Saved statistics to {stats_file}")


class BoxplotDataStep(WorkflowStep):
    """Prepare boxplot data"""

    def __init__(self):
        super().__init__("Boxplot Data", "Prepare boxplot visualization data")

    def execute(self, state: PipelineState, config: Dict[str, Any]) -> None:
        # Use statistics analyzer for boxplot data
        stats_analyzer = state.stats_analyzer

        state.analysis_results['boxplot_data'] = stats_analyzer.prepare_boxplot_data(
            state.filtered_data
        )
        state.analysis_results['boxplot_data_extended'] = stats_analyzer.prepare_boxplot_data_extended(
            state.filtered_data
        )


class PLSDAStep(WorkflowStep):
    """Perform PLS-DA and calculate VIP scores"""

    def __init__(self):
        super().__init__("PLS-DA", "PLS-DA analysis and VIP scores")

    def execute(self, state: PipelineState, config: Dict[str, Any]) -> None:
        # Use focused PLS-DA analyzer
        plsda_analyzer = PLSDAAnalyzer(
            n_components=2,
            log_transform=config['analysis']['pca']['log_transform']
        )

        plsda_results = plsda_analyzer.perform_plsda(state.filtered_data)
        state.analysis_results['plsda'] = plsda_results
        state.plsda_analyzer = plsda_analyzer

        # Save VIP scores
        results_dir = config['paths']['results_dir']
        vip_file = Path(results_dir) / 'vip_scores_all.csv'
        plsda_results['vip_scores'].to_csv(vip_file, index=False)
        logger.info(f"Saved VIP scores to {vip_file}")


# ==================== Main Pipeline ====================

class GlycoPipeline(BasePipeline):
    """
    Main glycoproteomics analysis pipeline

    Orchestrates the complete workflow from data loading to report generation
    """

    def __init__(self, config: Dict[str, Any]):
        super().__init__(config)
        self._build_workflows()

    def _build_workflows(self) -> None:
        """Build all workflow steps"""

        # Data Ingestion Workflow
        data_ingestion = Workflow("Data Ingestion")
        data_ingestion.add_step(LoadDataStep())

        # Annotation Workflow
        annotation = Workflow("Annotation")
        annotation.add_step(AnnotateDataStep())

        # Filtering Workflow
        filtering = Workflow("Filtering")
        filtering.add_step(FilterDataStep())
        filtering.add_step(ValidateStatisticalPowerStep())

        # Analysis Workflow
        analysis = Workflow("Analysis")
        analysis.add_step(PCAAnalysisStep())
        analysis.add_step(StatisticsStep())
        analysis.add_step(BoxplotDataStep())
        analysis.add_step(PLSDAStep())

        # Add all workflows
        self.workflows = [
            data_ingestion,
            annotation,
            filtering,
            analysis
        ]

    def run(self) -> PipelineState:
        """
        Execute the complete pipeline

        Returns:
            PipelineState with all results
        """
        try:
            self.on_start()
            self.validate_config()

            # Store config in state for access by downstream code
            self.state.config = self.config

            # Execute all workflows
            for workflow in self.workflows:
                workflow.execute(self.state, self.config)

            self.on_complete()
            return self.state

        except Exception as e:
            self.on_error(e)
            raise
