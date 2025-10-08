"""
Data Integrity Verification Module for pGlyco Auto Combine
Implements ALCOA++ data integrity requirements

Provides:
- SHA-256 checksum calculation for files
- Data integrity verification
- Manifest generation and validation
- Tamper detection for inputs and outputs

This ensures:
- Original data is preserved and verifiable
- Outputs can be validated against checksums
- Long-term data integrity verification
"""

import hashlib
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime, timezone

from .logger_config import get_logger
from .audit_logger import get_audit_logger

logger = get_logger(__name__)


class DataIntegrityManager:
    """
    Manager for data integrity verification using SHA-256 checksums

    Implements ALCOA++ requirements for original and accurate data.
    """

    def __init__(self):
        """Initialize data integrity manager"""
        self.audit_logger = get_audit_logger()
        self.checksums: Dict[str, str] = {}

    def calculate_file_checksum(self, file_path: Path, algorithm: str = 'sha256') -> str:
        """
        Calculate checksum for a file

        Args:
            file_path: Path to file
            algorithm: Hash algorithm (default: sha256)

        Returns:
            Hexadecimal checksum string
        """
        file_path = Path(file_path)

        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")

        # Create hasher
        hasher = hashlib.sha256() if algorithm == 'sha256' else hashlib.md5()

        # Read file in chunks (memory-efficient for large files)
        chunk_size = 65536  # 64KB chunks

        with open(file_path, 'rb') as f:
            while True:
                chunk = f.read(chunk_size)
                if not chunk:
                    break
                hasher.update(chunk)

        checksum = hasher.hexdigest()

        # Store checksum
        self.checksums[str(file_path)] = checksum

        return checksum

    def calculate_directory_checksums(
        self,
        directory: Path,
        pattern: str = '*',
        recursive: bool = False
    ) -> Dict[str, str]:
        """
        Calculate checksums for all files matching pattern in directory

        Args:
            directory: Directory to scan
            pattern: File pattern (glob) to match
            recursive: Whether to scan recursively

        Returns:
            Dictionary mapping file paths to checksums
        """
        directory = Path(directory)
        checksums = {}

        # Find matching files
        if recursive:
            files = list(directory.rglob(pattern))
        else:
            files = list(directory.glob(pattern))

        logger.info(f"Calculating checksums for {len(files)} files in {directory}...")

        for file_path in files:
            if file_path.is_file():
                try:
                    checksum = self.calculate_file_checksum(file_path)
                    # Store relative path
                    relative_path = file_path.relative_to(directory.parent)
                    checksums[str(relative_path)] = checksum
                except Exception as e:
                    logger.warning(f"Failed to calculate checksum for {file_path}: {e}")

        self.checksums.update(checksums)

        logger.info(f"Calculated {len(checksums)} checksums")

        return checksums

    def verify_file_checksum(self, file_path: Path, expected_checksum: str) -> bool:
        """
        Verify file checksum against expected value

        Args:
            file_path: Path to file
            expected_checksum: Expected checksum (hex string)

        Returns:
            True if checksum matches, False otherwise
        """
        try:
            actual_checksum = self.calculate_file_checksum(file_path)

            if actual_checksum == expected_checksum:
                logger.info(f"✓ Checksum verified: {file_path.name}")
                return True
            else:
                logger.error(f"✗ Checksum mismatch: {file_path.name}")
                logger.error(f"  Expected: {expected_checksum}")
                logger.error(f"  Actual:   {actual_checksum}")
                return False
        except Exception as e:
            logger.error(f"Failed to verify checksum for {file_path}: {e}")
            return False

    def save_manifest(
        self,
        manifest_path: Path,
        checksums: Optional[Dict[str, str]] = None,
        metadata: Optional[Dict] = None
    ):
        """
        Save checksum manifest to JSON file

        Args:
            manifest_path: Path to save manifest
            checksums: Checksums to save (uses self.checksums if None)
            metadata: Additional metadata to include
        """
        manifest_path = Path(manifest_path)

        checksums_to_save = checksums if checksums is not None else self.checksums

        manifest = {
            'manifest_version': '1.0',
            'created': datetime.now(timezone.utc).isoformat(),
            'algorithm': 'SHA-256',
            'file_count': len(checksums_to_save),
            'checksums': checksums_to_save,
        }

        if metadata:
            manifest['metadata'] = metadata

        with open(manifest_path, 'w') as f:
            json.dump(manifest, f, indent=2)

        logger.info(f"Manifest saved: {manifest_path} ({len(checksums_to_save)} files)")

        # Log to audit trail
        self.audit_logger.log_file_saved(
            'checksum_manifest',
            manifest_path,
            manifest_path.stat().st_size
        )

    def load_manifest(self, manifest_path: Path) -> Dict[str, str]:
        """
        Load checksum manifest from JSON file

        Args:
            manifest_path: Path to manifest file

        Returns:
            Dictionary mapping file paths to checksums
        """
        manifest_path = Path(manifest_path)

        if not manifest_path.exists():
            raise FileNotFoundError(f"Manifest not found: {manifest_path}")

        with open(manifest_path) as f:
            manifest = json.load(f)

        checksums = manifest.get('checksums', {})

        logger.info(f"Loaded manifest: {manifest_path} ({len(checksums)} files)")

        return checksums

    def verify_manifest(
        self,
        manifest_path: Path,
        base_directory: Optional[Path] = None
    ) -> Tuple[int, int, List[str]]:
        """
        Verify all files in manifest have correct checksums

        Args:
            manifest_path: Path to manifest file
            base_directory: Base directory for resolving relative paths (optional)

        Returns:
            Tuple of (verified_count, failed_count, failed_files)
        """
        checksums = self.load_manifest(manifest_path)

        verified_count = 0
        failed_count = 0
        failed_files = []

        logger.info(f"Verifying {len(checksums)} files from manifest...")

        for file_path_str, expected_checksum in checksums.items():
            # Resolve path
            if base_directory:
                file_path = base_directory / file_path_str
            else:
                file_path = Path(file_path_str)

            # Verify checksum
            if self.verify_file_checksum(file_path, expected_checksum):
                verified_count += 1
            else:
                failed_count += 1
                failed_files.append(file_path_str)

        # Log results
        logger.info(f"Verification complete: {verified_count} verified, {failed_count} failed")

        if failed_count > 0:
            logger.error(f"Failed files: {failed_files}")

        # Log to audit trail
        self.audit_logger.log_qc_check(
            'manifest_verification',
            'PASS' if failed_count == 0 else 'FAIL',
            details={
                'verified_count': verified_count,
                'failed_count': failed_count,
                'failed_files': failed_files,
            }
        )

        return verified_count, failed_count, failed_files

    def create_input_manifest(
        self,
        dataset_dir: Path,
        output_path: Path
    ) -> Dict[str, str]:
        """
        Create manifest for input dataset files

        Args:
            dataset_dir: Directory containing input CSV files
            output_path: Path to save manifest

        Returns:
            Dictionary of checksums
        """
        logger.info("Creating input data manifest...")

        # Calculate checksums for all CSV files
        checksums = self.calculate_directory_checksums(
            dataset_dir,
            pattern='*.csv',
            recursive=False
        )

        # Save manifest
        self.save_manifest(
            output_path,
            checksums=checksums,
            metadata={
                'description': 'Input dataset checksums',
                'dataset_directory': str(dataset_dir),
            }
        )

        return checksums

    def create_output_manifest(
        self,
        results_dir: Path,
        output_path: Path
    ) -> Dict[str, str]:
        """
        Create manifest for output files

        Args:
            results_dir: Directory containing results
            output_path: Path to save manifest

        Returns:
            Dictionary of checksums
        """
        logger.info("Creating output files manifest...")

        checksums = {}

        # Calculate checksums for CSV files
        csv_checksums = self.calculate_directory_checksums(
            results_dir,
            pattern='*.csv',
            recursive=False
        )
        checksums.update(csv_checksums)

        # Calculate checksums for PNG files
        png_checksums = self.calculate_directory_checksums(
            results_dir,
            pattern='*.png',
            recursive=False
        )
        checksums.update(png_checksums)

        # Calculate checksums for TXT files
        txt_checksums = self.calculate_directory_checksums(
            results_dir,
            pattern='*.txt',
            recursive=False
        )
        checksums.update(txt_checksums)

        # Save manifest
        self.save_manifest(
            output_path,
            checksums=checksums,
            metadata={
                'description': 'Output files checksums',
                'results_directory': str(results_dir),
            }
        )

        return checksums

    def get_checksum_summary(self) -> Dict[str, any]:
        """
        Get summary of checksums

        Returns:
            Dictionary with checksum statistics
        """
        return {
            'total_files': len(self.checksums),
            'algorithm': 'SHA-256',
            'checksums': self.checksums,
        }


def calculate_checksum(file_path: Path) -> str:
    """
    Convenience function to calculate SHA-256 checksum for a file

    Args:
        file_path: Path to file

    Returns:
        Hexadecimal checksum string
    """
    manager = DataIntegrityManager()
    return manager.calculate_file_checksum(file_path)


def verify_checksum(file_path: Path, expected_checksum: str) -> bool:
    """
    Convenience function to verify file checksum

    Args:
        file_path: Path to file
        expected_checksum: Expected checksum

    Returns:
        True if checksum matches
    """
    manager = DataIntegrityManager()
    return manager.verify_file_checksum(file_path, expected_checksum)
