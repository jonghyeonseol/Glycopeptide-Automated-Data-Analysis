"""
pGlyco Auto Combine - Source Package
Glycoproteomics data integration and analysis toolkit
"""

__version__ = "2.0.0"
__author__ = "pGlyco Auto Combine Team"

# Import key components for easier access
from .logger_config import setup_logging, get_logger
from .utils import (
    replace_empty_with_zero,
    to_numeric_safe,
    get_sample_columns,
    get_all_sample_columns,
    save_trace_data,
    ensure_directory,
    calculate_fold_change
)
from .config_validator import load_and_validate_config
from .data_loader import DataLoader
from .annotator import GlycanAnnotator
from .analyzer import GlycanAnalyzer
from .visualizer import GlycanVisualizer

__all__ = [
    # Core classes
    'DataLoader',
    'GlycanAnnotator',
    'GlycanAnalyzer',
    'GlycanVisualizer',

    # Configuration
    'load_and_validate_config',
    'setup_logging',
    'get_logger',

    # Utilities
    'replace_empty_with_zero',
    'to_numeric_safe',
    'get_sample_columns',
    'get_all_sample_columns',
    'save_trace_data',
    'ensure_directory',
    'calculate_fold_change',
]
