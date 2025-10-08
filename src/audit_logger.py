"""
Audit Trail Logging Module for pGlyco Auto Combine
Implements ALCOA++ contemporaneous logging requirements

Provides:
- Timestamped event logging (ISO 8601 format)
- JSON Lines format (append-only, tamper-evident)
- Critical decision point tracking
- Data processing audit trail
- Regulatory compliance documentation

All events are logged with:
- Timestamp (UTC, ISO 8601)
- Event type
- Event description
- Associated data (counts, parameters, etc.)
- User/system context
"""

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Any, Optional
from enum import Enum

from .logger_config import get_logger
from .metadata_collector import get_metadata_collector

logger = get_logger(__name__)


class EventType(Enum):
    """Enumeration of audit event types"""

    # Pipeline execution events
    PIPELINE_START = "pipeline_start"
    PIPELINE_COMPLETE = "pipeline_complete"
    PIPELINE_ERROR = "pipeline_error"

    # Data loading events
    DATA_LOAD_START = "data_load_start"
    DATA_LOAD_COMPLETE = "data_load_complete"
    DATA_FILE_LOADED = "data_file_loaded"

    # Data processing events
    DATA_INTEGRATION = "data_integration"
    DATA_ANNOTATION = "data_annotation"
    DATA_FILTERING = "data_filtering"
    DATA_VALIDATION = "data_validation"

    # Statistical analysis events
    ANALYSIS_PCA = "analysis_pca"
    ANALYSIS_PLSDA = "analysis_plsda"
    ANALYSIS_VIP = "analysis_vip_scores"
    ANALYSIS_STATISTICS = "analysis_statistics"

    # Visualization events
    VISUALIZATION_START = "visualization_start"
    VISUALIZATION_COMPLETE = "visualization_complete"
    PLOT_GENERATED = "plot_generated"

    # Configuration events
    CONFIG_LOADED = "config_loaded"
    CONFIG_VALIDATED = "config_validated"

    # Output events
    OUTPUT_SAVED = "output_saved"
    REPORT_GENERATED = "report_generated"

    # Quality control events
    QC_CHECK = "qc_check"
    QC_WARNING = "qc_warning"
    QC_FAILURE = "qc_failure"


class AuditLogger:
    """
    Audit trail logger for regulatory compliance

    Implements append-only logging in JSON Lines format for tamper-evidence.
    """

    _instance: Optional['AuditLogger'] = None
    _initialized: bool = False

    def __new__(cls):
        """Singleton pattern - only one instance allowed"""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        """Initialize audit logger (only once)"""
        if self._initialized:
            return

        self.audit_file: Optional[Path] = None
        self.event_count = 0
        self.metadata_collector = get_metadata_collector()

        self._initialized = True

    def initialize(self, output_dir: Path):
        """
        Initialize audit log file

        Args:
            output_dir: Directory to save audit log
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Create audit log with timestamp in filename
        execution_id = self.metadata_collector.execution_id
        self.audit_file = output_dir / f"audit_log_{execution_id}.jsonl"

        logger.info(f"Audit log initialized: {self.audit_file}")

        # Log initialization event
        self.log_event(
            EventType.PIPELINE_START,
            "Pipeline execution started",
            metadata=self.metadata_collector.get_metadata_dict()
        )

    def log_event(
        self,
        event_type: EventType,
        description: str,
        data: Optional[Dict[str, Any]] = None,
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        Log an audit event

        Args:
            event_type: Type of event (from EventType enum)
            description: Human-readable event description
            data: Event-specific data (counts, parameters, etc.)
            metadata: Additional metadata (optional)
        """
        # Create event record
        event = {
            'event_id': self.event_count,
            'timestamp': datetime.now(timezone.utc).isoformat(),
            'event_type': event_type.value,
            'description': description,
            'execution_id': self.metadata_collector.execution_id,
        }

        # Add data if provided
        if data:
            event['data'] = data

        # Add metadata if provided
        if metadata:
            event['metadata'] = metadata

        # Write to audit log (JSON Lines format - one JSON object per line)
        if self.audit_file:
            try:
                with open(self.audit_file, 'a') as f:
                    f.write(json.dumps(event) + '\n')
            except OSError as e:
                logger.error(f"Failed to write to audit log: {e}")

        # Increment event counter
        self.event_count += 1

        # Also log to standard logger for visibility
        logger.info(f"[AUDIT] {event_type.value}: {description}")

    def log_data_load(self, file_path: Path, row_count: int, column_count: int):
        """Log data file loading event"""
        self.log_event(
            EventType.DATA_FILE_LOADED,
            f"Loaded data file: {file_path.name}",
            data={
                'file_path': str(file_path),
                'file_name': file_path.name,
                'rows': row_count,
                'columns': column_count,
            }
        )

    def log_data_integration(self, n_files: int, n_glycopeptides: int, n_samples: int):
        """Log data integration event"""
        self.log_event(
            EventType.DATA_INTEGRATION,
            f"Integrated {n_files} files into {n_glycopeptides} glycopeptides × {n_samples} samples",
            data={
                'n_files': n_files,
                'n_glycopeptides': n_glycopeptides,
                'n_samples': n_samples,
            }
        )

    def log_annotation(self, n_glycopeptides: int, annotation_stats: Dict[str, int]):
        """Log annotation event"""
        self.log_event(
            EventType.DATA_ANNOTATION,
            f"Annotated {n_glycopeptides} glycopeptides",
            data={
                'n_glycopeptides': n_glycopeptides,
                'annotation_counts': annotation_stats,
            }
        )

    def log_filtering(
        self,
        n_before: int,
        n_after: int,
        n_removed: int,
        pct_removed: float,
        filter_criteria: Dict[str, Any]
    ):
        """Log data filtering event (CRITICAL for ALCOA++)"""
        self.log_event(
            EventType.DATA_FILTERING,
            f"Applied detection filter: {n_before} → {n_after} glycopeptides ({pct_removed:.1f}% removed)",
            data={
                'n_before': n_before,
                'n_after': n_after,
                'n_removed': n_removed,
                'pct_removed': pct_removed,
                'filter_criteria': filter_criteria,
            }
        )

    def log_pca(self, n_components: int, variance_explained: list):
        """Log PCA analysis event"""
        self.log_event(
            EventType.ANALYSIS_PCA,
            f"PCA analysis completed with {n_components} components",
            data={
                'n_components': n_components,
                'variance_explained': variance_explained,
                'total_variance': sum(variance_explained),
            }
        )

    def log_plsda(self, n_components: int, n_features: int):
        """Log PLS-DA analysis event"""
        self.log_event(
            EventType.ANALYSIS_PLSDA,
            f"PLS-DA analysis completed with {n_components} components on {n_features} features",
            data={
                'n_components': n_components,
                'n_features': n_features,
            }
        )

    def log_vip_scores(self, n_glycopeptides: int, vip_threshold: float = 1.0):
        """Log VIP score calculation event"""
        self.log_event(
            EventType.ANALYSIS_VIP,
            f"VIP scores calculated for {n_glycopeptides} glycopeptides",
            data={
                'n_glycopeptides': n_glycopeptides,
                'vip_threshold': vip_threshold,
            }
        )

    def log_statistical_test(
        self,
        test_name: str,
        n_tests: int,
        n_significant: int,
        alpha: float,
        correction_method: Optional[str] = None
    ):
        """Log statistical testing event"""
        self.log_event(
            EventType.ANALYSIS_STATISTICS,
            f"{test_name}: {n_significant}/{n_tests} significant (α={alpha})",
            data={
                'test_name': test_name,
                'n_tests': n_tests,
                'n_significant': n_significant,
                'alpha': alpha,
                'correction_method': correction_method or 'none',
            }
        )

    def log_plot_generated(self, plot_name: str, output_path: Path, dpi: int = 300):
        """Log plot generation event"""
        self.log_event(
            EventType.PLOT_GENERATED,
            f"Generated plot: {plot_name}",
            data={
                'plot_name': plot_name,
                'output_path': str(output_path),
                'dpi': dpi,
            }
        )

    def log_file_saved(self, file_type: str, file_path: Path, file_size_bytes: int):
        """Log file output event"""
        self.log_event(
            EventType.OUTPUT_SAVED,
            f"Saved {file_type}: {file_path.name}",
            data={
                'file_type': file_type,
                'file_path': str(file_path),
                'file_name': file_path.name,
                'file_size_bytes': file_size_bytes,
                'file_size_mb': file_size_bytes / (1024 * 1024),
            }
        )

    def log_qc_check(self, check_name: str, result: str, details: Optional[Dict[str, Any]] = None):
        """Log quality control check"""
        self.log_event(
            EventType.QC_CHECK,
            f"QC Check - {check_name}: {result}",
            data={
                'check_name': check_name,
                'result': result,
                'details': details or {},
            }
        )

    def log_qc_warning(self, warning_message: str, details: Optional[Dict[str, Any]] = None):
        """Log quality control warning"""
        self.log_event(
            EventType.QC_WARNING,
            warning_message,
            data=details
        )

    def finalize(self):
        """Finalize audit log at pipeline completion"""
        # Get final metadata
        final_metadata = self.metadata_collector.finalize()

        # Log completion
        self.log_event(
            EventType.PIPELINE_COMPLETE,
            "Pipeline execution completed successfully",
            data={
                'total_events': self.event_count,
                'duration_seconds': final_metadata['execution']['duration_seconds'],
                'duration_formatted': final_metadata['execution']['duration_formatted'],
            },
            metadata=final_metadata
        )

        logger.info(f"Audit log finalized: {self.event_count} events recorded")

    def get_audit_summary(self) -> Dict[str, Any]:
        """
        Get summary of audit log

        Returns:
            Dictionary with audit summary statistics
        """
        # Count events by type
        event_counts = {}

        if self.audit_file and self.audit_file.exists():
            with open(self.audit_file) as f:
                for line in f:
                    event = json.loads(line)
                    event_type = event['event_type']
                    event_counts[event_type] = event_counts.get(event_type, 0) + 1

        return {
            'total_events': self.event_count,
            'audit_file': str(self.audit_file) if self.audit_file else None,
            'event_counts_by_type': event_counts,
        }


# Global singleton instance accessor
def get_audit_logger() -> AuditLogger:
    """
    Get the global AuditLogger instance

    Returns:
        AuditLogger singleton instance
    """
    return AuditLogger()
