"""
Logging Configuration for pGlyco Auto Combine
Centralized logging setup to avoid conflicts
"""

import logging
import sys
from pathlib import Path
from typing import Optional

from .constants import LOG_FORMAT, LOG_LEVEL_INFO


# Global flag to ensure logging is only configured once
_logging_configured = False


def setup_logging(level: str = LOG_LEVEL_INFO,
                  log_file: Optional[Path] = None,
                  console: bool = True) -> None:
    """
    Configure logging for the entire application

    Should be called once at the start of the application.
    Subsequent calls will be ignored to prevent reconfiguration.

    Args:
        level: Logging level (INFO, DEBUG, WARNING, ERROR)
        log_file: Optional path to log file for file output
        console: Whether to output to console (default: True)
    """
    global _logging_configured

    if _logging_configured:
        return

    # Get root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    # Remove any existing handlers to avoid duplicates
    root_logger.handlers = []

    # Create formatter
    formatter = logging.Formatter(LOG_FORMAT)

    # Console handler
    if console:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(level)
        console_handler.setFormatter(formatter)
        root_logger.addHandler(console_handler)

    # File handler (optional)
    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)

    _logging_configured = True


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger instance for a specific module

    Args:
        name: Name of the logger (typically __name__)

    Returns:
        Logger instance
    """
    return logging.getLogger(name)


def reset_logging() -> None:
    """
    Reset logging configuration (mainly for testing)
    """
    global _logging_configured
    _logging_configured = False
    root_logger = logging.getLogger()
    root_logger.handlers = []
