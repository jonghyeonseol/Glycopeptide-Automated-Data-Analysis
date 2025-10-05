"""
Custom Exceptions for pGlyco Auto Combine
Provides specific error types for better error handling and debugging
"""


class PGlycoAutoError(Exception):
    """Base exception for all pGlyco Auto Combine errors"""
    pass


# ==============================================================================
# Configuration Errors
# ==============================================================================

class ConfigurationError(PGlycoAutoError):
    """Raised when configuration is invalid or missing"""
    pass


class MissingConfigKeyError(ConfigurationError):
    """Raised when a required configuration key is missing"""
    def __init__(self, key: str):
        self.key = key
        super().__init__(f"Missing required configuration key: {key}")


class InvalidConfigValueError(ConfigurationError):
    """Raised when a configuration value is invalid"""
    def __init__(self, key: str, value, reason: str = ""):
        self.key = key
        self.value = value
        self.reason = reason
        message = f"Invalid configuration value for '{key}': {value}"
        if reason:
            message += f" - {reason}"
        super().__init__(message)


# ==============================================================================
# Data Loading Errors
# ==============================================================================

class DataLoadError(PGlycoAutoError):
    """Raised when data loading fails"""
    pass


class NoDataFilesError(DataLoadError):
    """Raised when no CSV files are found in the dataset directory"""
    def __init__(self, directory: str):
        self.directory = directory
        super().__init__(f"No CSV files found in directory: {directory}")


class MissingColumnError(DataLoadError):
    """Raised when required columns are missing from input data"""
    def __init__(self, filename: str, missing_columns: list):
        self.filename = filename
        self.missing_columns = missing_columns
        super().__init__(
            f"Missing required columns in {filename}: {', '.join(missing_columns)}"
        )


class InvalidDataFormatError(DataLoadError):
    """Raised when data format is invalid"""
    def __init__(self, filename: str, reason: str):
        self.filename = filename
        self.reason = reason
        super().__init__(f"Invalid data format in {filename}: {reason}")


class EmptyDataError(DataLoadError):
    """Raised when no valid data is loaded"""
    def __init__(self, context: str = ""):
        message = "No valid data loaded"
        if context:
            message += f": {context}"
        super().__init__(message)


# ==============================================================================
# Annotation Errors
# ==============================================================================

class AnnotationError(PGlycoAutoError):
    """Raised when annotation fails"""
    pass


class InvalidGlycanCompositionError(AnnotationError):
    """Raised when glycan composition string is invalid"""
    def __init__(self, composition: str, reason: str = ""):
        self.composition = composition
        self.reason = reason
        message = f"Invalid glycan composition: {composition}"
        if reason:
            message += f" - {reason}"
        super().__init__(message)


# ==============================================================================
# Analysis Errors
# ==============================================================================

class AnalysisError(PGlycoAutoError):
    """Raised when statistical analysis fails"""
    pass


class InsufficientDataError(AnalysisError):
    """Raised when insufficient data for analysis"""
    def __init__(self, analysis_type: str, required: int, actual: int):
        self.analysis_type = analysis_type
        self.required = required
        self.actual = actual
        super().__init__(
            f"Insufficient data for {analysis_type}: "
            f"required {required}, got {actual}"
        )


class MatrixShapeError(AnalysisError):
    """Raised when matrix dimensions are incompatible"""
    def __init__(self, expected_shape: tuple, actual_shape: tuple, context: str = ""):
        self.expected_shape = expected_shape
        self.actual_shape = actual_shape
        message = f"Matrix shape mismatch: expected {expected_shape}, got {actual_shape}"
        if context:
            message += f" - {context}"
        super().__init__(message)


class NormalizationError(AnalysisError):
    """Raised when normalization fails"""
    def __init__(self, reason: str):
        super().__init__(f"Normalization failed: {reason}")


# ==============================================================================
# Visualization Errors
# ==============================================================================

class VisualizationError(PGlycoAutoError):
    """Raised when visualization generation fails"""
    pass


class PlotGenerationError(VisualizationError):
    """Raised when a specific plot fails to generate"""
    def __init__(self, plot_type: str, reason: str):
        self.plot_type = plot_type
        self.reason = reason
        super().__init__(f"Failed to generate {plot_type}: {reason}")


class MissingVisualizationDataError(VisualizationError):
    """Raised when required data for visualization is missing"""
    def __init__(self, plot_type: str, missing_data: str):
        self.plot_type = plot_type
        self.missing_data = missing_data
        super().__init__(
            f"Missing required data for {plot_type}: {missing_data}"
        )


# ==============================================================================
# File I/O Errors
# ==============================================================================

class FileOperationError(PGlycoAutoError):
    """Raised when file operations fail"""
    pass


class OutputDirectoryError(FileOperationError):
    """Raised when output directory cannot be created or accessed"""
    def __init__(self, directory: str, reason: str):
        self.directory = directory
        self.reason = reason
        super().__init__(
            f"Cannot access/create output directory {directory}: {reason}"
        )


class TraceDataSaveError(FileOperationError):
    """Raised when trace data cannot be saved"""
    def __init__(self, filename: str, reason: str):
        self.filename = filename
        self.reason = reason
        super().__init__(f"Failed to save trace data {filename}: {reason}")


# ==============================================================================
# Validation Errors
# ==============================================================================

class ValidationError(PGlycoAutoError):
    """Raised when data validation fails"""
    pass


class SampleCountMismatchError(ValidationError):
    """Raised when sample counts don't match expected values"""
    def __init__(self, expected: int, actual: int, sample_type: str = ""):
        self.expected = expected
        self.actual = actual
        self.sample_type = sample_type
        message = f"Sample count mismatch: expected {expected}, got {actual}"
        if sample_type:
            message += f" for {sample_type}"
        super().__init__(message)


class ValueRangeError(ValidationError):
    """Raised when values are outside expected range"""
    def __init__(self, value, min_val, max_val, field_name: str = ""):
        self.value = value
        self.min_val = min_val
        self.max_val = max_val
        message = f"Value {value} outside valid range [{min_val}, {max_val}]"
        if field_name:
            message = f"{field_name}: {message}"
        super().__init__(message)
