"""
Preprocessing State Tracker for pGlyco Auto Combine

Tracks full preprocessing history for reproducibility.
Inspired by alphapeptstats' preprocessing_info dictionary pattern.

Author: pGlyco Auto Combine Pipeline
Created: 2025-10-14
"""

from dataclasses import dataclass, field, asdict
from typing import Dict, List, Any, Optional
from datetime import datetime
import json
from pathlib import Path
import numpy as np


@dataclass
class TransformationStep:
    """Record of a single preprocessing transformation"""
    step_number: int
    name: str
    timestamp: str
    parameters: Dict[str, Any]
    statistics_before: Optional[Dict[str, Any]] = None
    statistics_after: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict:
        """Convert to dictionary"""
        return asdict(self)


@dataclass
class PreprocessingState:
    """
    Complete preprocessing state tracking

    Inspired by alphapeptstats' preprocessing_info dictionary pattern.
    Tracks ALL transformations applied to the data for reproducibility.

    SCIENTIFIC RATIONALE:
    - Reproducibility requires complete documentation of all preprocessing steps
    - Each transformation can affect downstream analysis results
    - Complete audit trail enables validation and troubleshooting
    - Allows comparison of results across different preprocessing strategies
    """

    # ========== Transformation Flags ==========
    tic_normalized: bool = False
    log2_transformed: bool = False
    scaled: bool = False
    filtered: bool = False
    imputed: bool = False

    # ========== Transformation Parameters ==========
    log_transform_pseudocount: float = 1.0
    log_transform_base: int = 2
    scaler_type: str = "RobustScaler"
    normalization_method: str = "tic"

    # ========== Filtering Parameters ==========
    min_detection_pct: float = 0.30
    min_samples: int = 5
    missing_data_method: str = "skipna"

    # ========== Data Statistics ==========
    n_glycopeptides_original: int = 0
    n_glycopeptides_filtered: int = 0
    n_samples_cancer: int = 0
    n_samples_normal: int = 0
    n_samples_total: int = 0

    # ========== Inf/NaN Tracking ==========
    inf_cleaned_count: int = 0
    nan_count_final: int = 0

    # ========== History Log ==========
    transformation_history: List[TransformationStep] = field(default_factory=list)

    # ========== Metadata ==========
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())
    version: str = "3.1.0"

    def add_transformation(self, name: str, parameters: Dict[str, Any] = None,
                          statistics_before: Dict[str, Any] = None,
                          statistics_after: Dict[str, Any] = None):
        """
        Record a transformation step

        Args:
            name: Name of the transformation (e.g., "TIC Normalization")
            parameters: Parameters used for the transformation
            statistics_before: Data statistics before transformation
            statistics_after: Data statistics after transformation
        """
        step = TransformationStep(
            step_number=len(self.transformation_history) + 1,
            name=name,
            timestamp=datetime.now().isoformat(),
            parameters=parameters or {},
            statistics_before=statistics_before,
            statistics_after=statistics_after
        )
        self.transformation_history.append(step)

    def to_dict(self) -> Dict:
        """Export as dictionary"""
        result = {
            "transformations": {
                "tic_normalized": self.tic_normalized,
                "log2_transformed": self.log2_transformed,
                "scaled": self.scaled,
                "filtered": self.filtered,
                "imputed": self.imputed,
            },
            "parameters": {
                "log_transform_pseudocount": self.log_transform_pseudocount,
                "log_transform_base": self.log_transform_base,
                "scaler_type": self.scaler_type,
                "normalization_method": self.normalization_method,
                "min_detection_pct": self.min_detection_pct,
                "min_samples": self.min_samples,
                "missing_data_method": self.missing_data_method,
            },
            "statistics": {
                "n_glycopeptides_original": self.n_glycopeptides_original,
                "n_glycopeptides_filtered": self.n_glycopeptides_filtered,
                "n_samples_cancer": self.n_samples_cancer,
                "n_samples_normal": self.n_samples_normal,
                "n_samples_total": self.n_samples_total,
                "inf_cleaned_count": self.inf_cleaned_count,
                "nan_count_final": self.nan_count_final,
            },
            "history": [step.to_dict() for step in self.transformation_history],
            "metadata": {
                "created_at": self.created_at,
                "version": self.version,
            }
        }
        return result

    def save_to_file(self, filepath: str):
        """
        Save preprocessing state to JSON file

        Args:
            filepath: Path to save the JSON file
        """
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)

        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load_from_file(cls, filepath: str) -> 'PreprocessingState':
        """
        Load preprocessing state from JSON file

        Args:
            filepath: Path to the JSON file

        Returns:
            PreprocessingState instance
        """
        with open(filepath, 'r') as f:
            data = json.load(f)

        # Create instance
        state = cls()

        # Load transformations
        trans = data.get('transformations', {})
        state.tic_normalized = trans.get('tic_normalized', False)
        state.log2_transformed = trans.get('log2_transformed', False)
        state.scaled = trans.get('scaled', False)
        state.filtered = trans.get('filtered', False)
        state.imputed = trans.get('imputed', False)

        # Load parameters
        params = data.get('parameters', {})
        state.log_transform_pseudocount = params.get('log_transform_pseudocount', 1.0)
        state.log_transform_base = params.get('log_transform_base', 2)
        state.scaler_type = params.get('scaler_type', 'RobustScaler')
        state.normalization_method = params.get('normalization_method', 'tic')
        state.min_detection_pct = params.get('min_detection_pct', 0.30)
        state.min_samples = params.get('min_samples', 5)
        state.missing_data_method = params.get('missing_data_method', 'skipna')

        # Load statistics
        stats = data.get('statistics', {})
        state.n_glycopeptides_original = stats.get('n_glycopeptides_original', 0)
        state.n_glycopeptides_filtered = stats.get('n_glycopeptides_filtered', 0)
        state.n_samples_cancer = stats.get('n_samples_cancer', 0)
        state.n_samples_normal = stats.get('n_samples_normal', 0)
        state.n_samples_total = stats.get('n_samples_total', 0)
        state.inf_cleaned_count = stats.get('inf_cleaned_count', 0)
        state.nan_count_final = stats.get('nan_count_final', 0)

        # Load history
        history = data.get('history', [])
        for step_dict in history:
            step = TransformationStep(
                step_number=step_dict['step_number'],
                name=step_dict['name'],
                timestamp=step_dict['timestamp'],
                parameters=step_dict['parameters'],
                statistics_before=step_dict.get('statistics_before'),
                statistics_after=step_dict.get('statistics_after')
            )
            state.transformation_history.append(step)

        # Load metadata
        metadata = data.get('metadata', {})
        state.created_at = metadata.get('created_at', datetime.now().isoformat())
        state.version = metadata.get('version', '3.1.0')

        return state

    def summary(self) -> str:
        """
        Generate human-readable summary

        Returns:
            Multi-line string summarizing preprocessing
        """
        lines = [
            "=" * 80,
            "PREPROCESSING SUMMARY",
            "=" * 80,
            "",
            "Data Statistics:",
            f"  Original glycopeptides: {self.n_glycopeptides_original}",
            f"  Filtered glycopeptides: {self.n_glycopeptides_filtered}",
            f"  Retention rate: {self.n_glycopeptides_filtered / max(self.n_glycopeptides_original, 1) * 100:.1f}%",
            f"  Cancer samples: {self.n_samples_cancer}",
            f"  Normal samples: {self.n_samples_normal}",
            f"  Total samples: {self.n_samples_total}",
            "",
            "Transformations Applied:",
        ]

        # Add transformation flags
        transformations = [
            ("TIC Normalization", self.tic_normalized, self.normalization_method),
            ("Log2 Transformation", self.log2_transformed, f"pseudocount={self.log_transform_pseudocount}"),
            ("Scaling", self.scaled, self.scaler_type),
            ("Detection Filtering", self.filtered, f"≥{self.min_detection_pct*100:.0f}%"),
            ("Imputation", self.imputed, ""),
        ]

        for name, applied, detail in transformations:
            status = "✓" if applied else "✗"
            detail_str = f" ({detail})" if detail else ""
            lines.append(f"  [{status}] {name}{detail_str}")

        lines.extend([
            "",
            "Data Quality:",
            f"  Inf values cleaned: {self.inf_cleaned_count}",
            f"  Final NaN count: {self.nan_count_final}",
            f"  Missing data handling: {self.missing_data_method}",
            "",
            "Transformation History:",
        ])

        for step in self.transformation_history:
            lines.append(f"  {step.step_number}. {step.name}")
            if step.parameters:
                for key, value in step.parameters.items():
                    # Format value based on type
                    if isinstance(value, float):
                        value_str = f"{value:.4f}" if value < 1 else f"{value:.2f}"
                    else:
                        value_str = str(value)
                    lines.append(f"      - {key}: {value_str}")

        lines.extend([
            "",
            "=" * 80,
            f"Created: {self.created_at}",
            f"Version: {self.version}",
            "=" * 80,
        ])

        return "\n".join(lines)

    def get_alphapeptstats_style_dict(self) -> Dict:
        """
        Export in alphapeptstats preprocessing_info dictionary style

        This format is compatible with alphapeptstats' approach to tracking
        preprocessing state, making it easy to compare results.

        Returns:
            Dictionary in alphapeptstats style
        """
        return {
            "log2_transformed": self.log2_transformed,
            "normalized": self.tic_normalized,
            "scaled": self.scaled,
            "filtered": self.filtered,
            "imputed": self.imputed,
            "normalization_method": self.normalization_method,
            "scaler": self.scaler_type,
            "log_pseudocount": self.log_transform_pseudocount,
            "n_features_original": self.n_glycopeptides_original,
            "n_features_filtered": self.n_glycopeptides_filtered,
            "n_samples": self.n_samples_total,
            "missing_data_method": self.missing_data_method,
        }


class PreprocessingTracker:
    """
    Convenience wrapper for tracking preprocessing steps

    Usage:
        tracker = PreprocessingTracker()
        tracker.record_transformation("TIC Normalization", {...})
        tracker.save("Results/preprocessing_state.json")
    """

    def __init__(self):
        self.state = PreprocessingState()

    def record_transformation(self, name: str, parameters: Dict[str, Any] = None,
                             statistics_before: Dict[str, Any] = None,
                             statistics_after: Dict[str, Any] = None):
        """Record a transformation step"""
        self.state.add_transformation(name, parameters, statistics_before, statistics_after)

    def set_sample_counts(self, n_cancer: int, n_normal: int):
        """Set sample counts"""
        self.state.n_samples_cancer = n_cancer
        self.state.n_samples_normal = n_normal
        self.state.n_samples_total = n_cancer + n_normal

    def set_glycopeptide_counts(self, n_original: int, n_filtered: int = None):
        """Set glycopeptide counts"""
        self.state.n_glycopeptides_original = n_original
        if n_filtered is not None:
            self.state.n_glycopeptides_filtered = n_filtered
        else:
            self.state.n_glycopeptides_filtered = n_original

    def mark_tic_normalized(self, method: str = "tic"):
        """Mark that TIC normalization was applied"""
        self.state.tic_normalized = True
        self.state.normalization_method = method

    def mark_log_transformed(self, pseudocount: float = 1.0, base: int = 2):
        """Mark that log transformation was applied"""
        self.state.log2_transformed = True
        self.state.log_transform_pseudocount = pseudocount
        self.state.log_transform_base = base

    def mark_scaled(self, scaler_type: str = "RobustScaler"):
        """Mark that scaling was applied"""
        self.state.scaled = True
        self.state.scaler_type = scaler_type

    def mark_filtered(self, min_detection_pct: float, min_samples: int):
        """Mark that filtering was applied"""
        self.state.filtered = True
        self.state.min_detection_pct = min_detection_pct
        self.state.min_samples = min_samples

    def increment_inf_cleaned(self, count: int):
        """Increment count of inf values cleaned"""
        self.state.inf_cleaned_count += count

    def set_final_nan_count(self, count: int):
        """Set final NaN count"""
        self.state.nan_count_final = count

    def save(self, filepath: str):
        """Save to JSON file"""
        self.state.save_to_file(filepath)

    def print_summary(self):
        """Print human-readable summary"""
        print(self.state.summary())

    def get_dict(self) -> Dict:
        """Get as dictionary"""
        return self.state.to_dict()

    def get_alphapeptstats_dict(self) -> Dict:
        """Get in alphapeptstats style"""
        return self.state.get_alphapeptstats_style_dict()
