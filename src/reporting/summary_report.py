"""
Summary Report Generator
Generates comprehensive analysis summary
"""

from typing import Dict, Any
from pathlib import Path
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from .base_report import BaseReport
from src.logger_config import get_logger

logger = get_logger(__name__)


class SummaryReport(BaseReport):
    """
    Comprehensive analysis summary report

    Generates sections for:
    - Data filtering
    - Sample information
    - Glycan annotation
    - PCA results
    - Statistics
    - VIP scores
    - Output files
    """

    def __init__(self):
        super().__init__(title="pGlyco Auto Combine - Analysis Summary")

    def generate(self, data: Dict[str, Any]) -> str:
        """
        Generate summary report from pipeline state

        Args:
            data: Dictionary containing:
                - pipeline_state: PipelineState object
                - config: Configuration dict
                - filtering_report: Filtering report text
                - preprocessing_tracker: PreprocessingTracker object (optional)

        Returns:
            Complete summary report
        """
        state = data.get('pipeline_state')
        config = data.get('config')
        filtering_report = data.get('filtering_report', '')
        preprocessing_tracker = data.get('preprocessing_tracker')

        # Clear any existing sections
        self.clear()

        # Add sections in order
        self.add_section(
            "DATA FILTERING REPORT",
            filtering_report,
            order=0
        )

        # Phase 1.1: Add preprocessing section for reproducibility
        if preprocessing_tracker:
            self.add_section(
                "PREPROCESSING STATE (Reproducibility)",
                self._generate_preprocessing_info(preprocessing_tracker),
                order=1
            )

        self.add_section(
            "Data Integration",
            self._generate_data_info(state),
            order=2
        )

        self.add_section(
            "Glycan Annotation",
            self._generate_annotation_info(state.filtered_data),
            order=3
        )

        self.add_section(
            "PCA Results",
            self._generate_pca_info(state.analysis_results.get('pca')),
            order=4
        )

        self.add_section(
            "Statistics by Glycan Type",
            self._generate_statistics_info(state.analysis_results.get('statistics')),
            order=5
        )

        self.add_section(
            "Top 10 VIP Scores (Glycopeptide)",
            self._generate_vip_info(state.analysis_results.get('plsda')),
            order=6
        )

        self.add_section(
            "Output Files",
            self._generate_output_files_info(config, state),
            order=7
        )

        return self.render()

    def _generate_preprocessing_info(self, tracker) -> str:
        """
        Generate preprocessing state section for reproducibility

        Phase 1.1: Added to ensure complete documentation of all data transformations
        """
        lines = []
        lines.append("Transformations Applied:")

        state_dict = tracker.get_dict()
        transformations = state_dict.get('transformations', {})
        parameters = state_dict.get('parameters', {})
        statistics = state_dict.get('statistics', {})

        # List transformations with checkmarks
        trans_list = [
            ("TIC Normalization", transformations.get('tic_normalized', False), parameters.get('normalization_method', '')),
            ("Log2 Transformation", transformations.get('log2_transformed', False), f"pseudocount={parameters.get('log_transform_pseudocount', 1.0)}"),
            ("Scaling", transformations.get('scaled', False), parameters.get('scaler_type', '')),
            ("Detection Filtering", transformations.get('filtered', False), f"≥{parameters.get('min_detection_pct', 0.30)*100:.0f}%"),
        ]

        for name, applied, detail in trans_list:
            status = "✓" if applied else "✗"
            detail_str = f" ({detail})" if detail else ""
            lines.append(f"  [{status}] {name}{detail_str}")

        lines.append("")
        lines.append("Data Statistics:")
        lines.append(f"  - Original glycopeptides: {statistics.get('n_glycopeptides_original', 'N/A')}")
        lines.append(f"  - Filtered glycopeptides: {statistics.get('n_glycopeptides_filtered', 'N/A')}")
        lines.append(f"  - Cancer samples: {statistics.get('n_samples_cancer', 'N/A')}")
        lines.append(f"  - Normal samples: {statistics.get('n_samples_normal', 'N/A')}")

        lines.append("")
        lines.append("Missing Data Handling:")
        lines.append(f"  - Method: {parameters.get('missing_data_method', 'skipna')}")
        lines.append("  - Scientifically appropriate for MNAR (Missing Not At Random) data")

        lines.append("")
        lines.append("Reproducibility:")
        lines.append("  - Complete state saved to: Results/preprocessing_state.json")
        lines.append("  - Compatible with alphapeptstats preprocessing_info format")

        return "\n".join(lines)

    def _generate_data_info(self, state) -> str:
        """Generate data integration section"""
        df = state.filtered_data
        sample_cols = [col for col in df.columns if col.startswith(('C', 'N'))]
        cancer_cols = [col for col in sample_cols if col.startswith('C')]
        normal_cols = [col for col in sample_cols if col.startswith('N')]

        return (
            f"  - Total glycopeptides: {len(df)}\n"
            f"  - Total samples: {len(sample_cols)}\n"
            f"  - Cancer samples: {len(cancer_cols)}\n"
            f"  - Normal samples: {len(normal_cols)}"
        )

    def _generate_annotation_info(self, df: pd.DataFrame) -> str:
        """Generate glycan annotation section"""
        lines = []

        # Sialylation/Fucosylation
        sialylated = df['IsSialylated'].sum()
        fucosylated = df['IsFucosylated'].sum()
        lines.append(f"  - Sialylated: {sialylated} ({sialylated / len(df) * 100:.1f}%)")
        lines.append(f"  - Fucosylated: {fucosylated} ({fucosylated / len(df) * 100:.1f}%)")

        # Glycan Type Distribution
        lines.append("\nGlycan Type Distribution (Legacy):")
        for glycan_type, count in df['GlycanType'].value_counts().items():
            lines.append(f"  - {glycan_type}: {count} ({count / len(df) * 100:.1f}%)")

        return "\n".join(lines)

    def _generate_pca_info(self, pca_results: Dict) -> str:
        """Generate PCA results section"""
        if not pca_results:
            return "  No PCA results available"

        explained_var = pca_results['explained_variance']

        return (
            f"  - PC1 explained variance: {explained_var[0] * 100:.2f}%\n"
            f"  - PC2 explained variance: {explained_var[1] * 100:.2f}%\n"
            f"  - Total explained variance: {explained_var.sum() * 100:.2f}%"
        )

    def _generate_statistics_info(self, stats_df: pd.DataFrame) -> str:
        """Generate statistics section"""
        if stats_df is None or stats_df.empty:
            return "  No statistics available"

        return stats_df.to_string(index=False)

    def _generate_vip_info(self, plsda_results: Dict) -> str:
        """Generate VIP scores section"""
        if not plsda_results:
            return "  No VIP scores available"

        vip_df = plsda_results['vip_scores'].head(10)

        return (
            vip_df[['Peptide', 'GlycanComposition', 'VIP_Score']].to_string(index=False) +
            "\n\nNote: VIP scores by Glycan Type and Peptide are visualized in the corresponding PNG files."
        )

    def _generate_output_files_info(self, config: Dict, state) -> str:
        """Generate output files section"""
        results_dir = config['paths']['results_dir']

        lines = [
            "  - Integrated data (RAW): integrated.csv",
            "  - Integrated data (FILTERED): integrated_filtered.csv (USED IN ALL ANALYSES)",
            "  - Filtering report: filtering_report.txt",
            f"  - Statistics: {results_dir}/glycan_type_statistics.csv",
            f"  - VIP Scores (all): {results_dir}/vip_scores_all.csv",
            "",
            "Visualization Files:",
            "  - All PNG files saved to Results/ directory",
            "  - 300 DPI publication quality",
            "",
            "Data Traceability:",
            "  All visualization source data is exported to Results/Trace/ folder as CSV files.",
            "  Each PNG visualization has a corresponding *_data.csv file containing the exact",
            "  data used to generate that plot, ensuring full reproducibility and transparency."
        ]

        return "\n".join(lines)
