"""
Metadata Collection Module for pGlyco Auto Combine
Implements ALCOA++ compliance for data integrity and traceability

Captures:
- Execution timestamp (ISO 8601 format)
- User identification (system username + optional researcher ID)
- System information (Python version, OS, hostname, platform)
- Git commit hash (links outputs to exact code version)
- Software environment (package versions)
- Configuration snapshot

This module ensures all outputs are attributable, contemporaneous, and traceable.
"""

import os
import sys
import platform
import socket
import subprocess
import getpass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Any
import json

from .logger_config import get_logger

logger = get_logger(__name__)


class MetadataCollector:
    """
    Singleton class for collecting execution metadata

    Ensures all outputs have complete traceability information for ALCOA++ compliance.
    """

    _instance: Optional['MetadataCollector'] = None
    _initialized: bool = False

    def __new__(cls):
        """Singleton pattern - only one instance allowed"""
        if cls._instance is None:
            cls._instance = super(MetadataCollector, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        """Initialize metadata collector (only once)"""
        if self._initialized:
            return

        logger.info("Initializing MetadataCollector for ALCOA++ compliance...")

        # Execution timing
        self.execution_start = datetime.now(timezone.utc)
        self.execution_start_str = self.execution_start.isoformat()

        # User identification
        self.system_username = getpass.getuser()
        self.researcher_id = None  # Can be set by user

        # System information
        self.hostname = socket.gethostname()
        self.platform_system = platform.system()
        self.platform_release = platform.release()
        self.platform_version = platform.version()
        self.platform_machine = platform.machine()

        # Python information
        self.python_version = sys.version
        self.python_version_short = platform.python_version()
        self.python_implementation = platform.python_implementation()

        # Git information
        self.git_commit_hash = self._get_git_commit_hash()
        self.git_branch = self._get_git_branch()
        self.git_has_changes = self._check_git_changes()

        # Package versions
        self.package_versions = self._get_package_versions()

        # Pipeline version
        self.pipeline_version = self._get_pipeline_version()

        # Execution ID (unique identifier)
        self.execution_id = self._generate_execution_id()

        self._initialized = True
        logger.info(f"Metadata collection initialized. Execution ID: {self.execution_id}")

    def _get_git_commit_hash(self) -> str:
        """Get current git commit hash"""
        try:
            result = subprocess.run(
                ['git', 'rev-parse', 'HEAD'],
                capture_output=True,
                text=True,
                check=True,
                timeout=5
            )
            return result.stdout.strip()[:8]  # Short hash
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            return "unknown"

    def _get_git_branch(self) -> str:
        """Get current git branch"""
        try:
            result = subprocess.run(
                ['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
                capture_output=True,
                text=True,
                check=True,
                timeout=5
            )
            return result.stdout.strip()
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            return "unknown"

    def _check_git_changes(self) -> bool:
        """Check if there are uncommitted changes"""
        try:
            result = subprocess.run(
                ['git', 'status', '--porcelain'],
                capture_output=True,
                text=True,
                check=True,
                timeout=5
            )
            return len(result.stdout.strip()) > 0
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            return False

    def _get_package_versions(self) -> Dict[str, str]:
        """Get versions of key packages"""
        packages = {}

        # Core data science packages
        package_names = [
            'pandas', 'numpy', 'scipy', 'scikit-learn',
            'matplotlib', 'seaborn', 'statsmodels'
        ]

        for pkg_name in package_names:
            try:
                if pkg_name == 'scikit-learn':
                    import sklearn
                    packages[pkg_name] = sklearn.__version__
                else:
                    module = __import__(pkg_name)
                    packages[pkg_name] = getattr(module, '__version__', 'unknown')
            except ImportError:
                packages[pkg_name] = 'not installed'

        return packages

    def _get_pipeline_version(self) -> str:
        """Get pipeline version from git tags or default"""
        try:
            result = subprocess.run(
                ['git', 'describe', '--tags', '--always'],
                capture_output=True,
                text=True,
                check=True,
                timeout=5
            )
            return result.stdout.strip()
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            return "v2.1.0"  # Default version

    def _generate_execution_id(self) -> str:
        """Generate unique execution identifier"""
        timestamp = self.execution_start.strftime("%Y%m%d_%H%M%S")
        return f"{timestamp}_{self.system_username}_{self.hostname}"

    def set_researcher_id(self, researcher_id: str):
        """
        Set researcher identification (optional)

        Args:
            researcher_id: Researcher identifier (e.g., email, ORCID, or name)
        """
        self.researcher_id = researcher_id
        logger.info(f"Researcher ID set: {researcher_id}")

    def get_metadata_dict(self) -> Dict[str, Any]:
        """
        Get complete metadata as dictionary

        Returns:
            Dictionary with all metadata fields
        """
        return {
            'execution': {
                'execution_id': self.execution_id,
                'start_timestamp': self.execution_start_str,
                'start_timestamp_local': self.execution_start.astimezone().isoformat(),
            },
            'user': {
                'system_username': self.system_username,
                'researcher_id': self.researcher_id or 'not specified',
            },
            'system': {
                'hostname': self.hostname,
                'platform': self.platform_system,
                'platform_release': self.platform_release,
                'platform_version': self.platform_version,
                'machine': self.platform_machine,
            },
            'python': {
                'version': self.python_version_short,
                'implementation': self.python_implementation,
                'full_version': self.python_version,
            },
            'git': {
                'commit_hash': self.git_commit_hash,
                'branch': self.git_branch,
                'has_uncommitted_changes': self.git_has_changes,
            },
            'pipeline': {
                'version': self.pipeline_version,
            },
            'packages': self.package_versions,
        }

    def get_metadata_header_lines(self, prefix: str = "#") -> str:
        """
        Get metadata formatted as comment lines for CSV/text files

        Args:
            prefix: Comment prefix (default: "#")

        Returns:
            Multi-line string with metadata
        """
        lines = []
        lines.append(f"{prefix} ========================================")
        lines.append(f"{prefix} pGlyco Auto Combine - Execution Metadata")
        lines.append(f"{prefix} ========================================")
        lines.append(f"{prefix}")
        lines.append(f"{prefix} Execution ID: {self.execution_id}")
        lines.append(f"{prefix} Generated: {self.execution_start_str}")
        lines.append(f"{prefix}")
        lines.append(f"{prefix} User: {self.system_username}")
        if self.researcher_id:
            lines.append(f"{prefix} Researcher: {self.researcher_id}")
        lines.append(f"{prefix}")
        lines.append(f"{prefix} System: {self.platform_system} {self.platform_release}")
        lines.append(f"{prefix} Hostname: {self.hostname}")
        lines.append(f"{prefix} Python: {self.python_version_short}")
        lines.append(f"{prefix}")
        lines.append(f"{prefix} Pipeline Version: {self.pipeline_version}")
        lines.append(f"{prefix} Git Commit: {self.git_commit_hash}")
        lines.append(f"{prefix} Git Branch: {self.git_branch}")
        if self.git_has_changes:
            lines.append(f"{prefix} WARNING: Uncommitted changes detected!")
        lines.append(f"{prefix}")
        lines.append(f"{prefix} Key Package Versions:")
        for pkg, version in self.package_versions.items():
            lines.append(f"{prefix}   {pkg}: {version}")
        lines.append(f"{prefix} ========================================")
        lines.append(f"{prefix}")

        return "\n".join(lines)

    def save_metadata_json(self, output_path: Path):
        """
        Save complete metadata as JSON file

        Args:
            output_path: Path to save JSON metadata
        """
        metadata = self.get_metadata_dict()

        with open(output_path, 'w') as f:
            json.dump(metadata, f, indent=2)

        logger.info(f"Metadata saved to {output_path}")

    def get_execution_duration(self) -> float:
        """
        Get execution duration in seconds (from start until now)

        Returns:
            Duration in seconds
        """
        now = datetime.now(timezone.utc)
        duration = (now - self.execution_start).total_seconds()
        return duration

    def finalize(self) -> Dict[str, Any]:
        """
        Finalize metadata collection at end of execution

        Returns:
            Complete metadata including execution duration
        """
        metadata = self.get_metadata_dict()

        # Add completion information
        completion_time = datetime.now(timezone.utc)
        duration = self.get_execution_duration()

        metadata['execution']['completion_timestamp'] = completion_time.isoformat()
        metadata['execution']['duration_seconds'] = duration
        metadata['execution']['duration_formatted'] = self._format_duration(duration)

        logger.info(f"Execution completed in {metadata['execution']['duration_formatted']}")

        return metadata

    def _format_duration(self, seconds: float) -> str:
        """Format duration in human-readable format"""
        if seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            minutes = seconds / 60
            return f"{minutes:.1f}m"
        else:
            hours = seconds / 3600
            return f"{hours:.2f}h"


# Global singleton instance accessor
def get_metadata_collector() -> MetadataCollector:
    """
    Get the global MetadataCollector instance

    Returns:
        MetadataCollector singleton instance
    """
    return MetadataCollector()
