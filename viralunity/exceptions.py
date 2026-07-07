"""Custom exceptions for ViralUnity pipeline.

Each exception carries a stable, machine-readable ``code`` so an embedding
service can surface a structured error payload (``to_dict()``) to callers
instead of only a log line.
"""


class ViralUnityError(Exception):
    """Base exception for all ViralUnity errors."""

    #: Stable, machine-readable identifier for this error class. Overridden by
    #: subclasses; may also be overridden per-instance via the ``code`` kwarg.
    code = "viralunity_error"

    def __init__(self, message: str = "", *, code: str = None):
        super().__init__(message)
        self.message = message
        if code is not None:
            self.code = code

    def to_dict(self) -> dict:
        """Return a structured, JSON-serializable representation of the error."""
        return {
            "error": type(self).__name__,
            "code": self.code,
            "message": str(self),
        }


class ValidationError(ViralUnityError):
    """Raised when validation of arguments or data fails."""

    code = "validation_error"


class ViralUnityFileNotFoundError(ViralUnityError):
    """Raised when a required file or directory is not found."""

    code = "file_not_found"


class ConfigurationError(ViralUnityError):
    """Raised when there's an error in configuration."""

    code = "configuration_error"


class SampleSheetError(ViralUnityError):
    """Raised when there's an error processing the sample sheet."""

    code = "sample_sheet_error"


class SampleConfigurationNotFoundError(ViralUnityError):
    """Raised when the sample configuration is not found."""

    code = "sample_configuration_not_found"


class Kraken2DatabaseNotFoundError(ViralUnityError):
    """Raised when the Kraken2 database is not found."""

    code = "kraken2_database_not_found"


class KronaDatabaseNotFoundError(ViralUnityError):
    """Raised when the Krona database is not found."""

    code = "krona_database_not_found"


class ReferenceNotFoundError(ViralUnityError):
    """Raised when the reference sequence file is not found."""

    code = "reference_not_found"


class PrimerSchemeNotFoundError(ViralUnityError):
    """Raised when the primer scheme file is not found."""

    code = "primer_scheme_not_found"


class AdaptersNotFoundError(ViralUnityError):
    """Raised when the Illumina adapter sequences file is not found or not provided."""

    code = "adapters_not_found"


class TaxdumpNotFoundError(ViralUnityError):
    """Raised when the NCBI taxdump directory is not found or invalid."""

    code = "taxdump_not_found"


class DiamondDatabaseNotFoundError(ViralUnityError):
    """Raised when the Diamond database or assembly summary is not found."""

    code = "diamond_database_not_found"
