"""Typed exceptions used across FRONTStr.

Keep this module free of heavyweight imports so any layer can raise these.
"""

from __future__ import annotations


class FrontstrError(Exception):
    """Base class for every FRONTStr-specific exception."""


class IngestError(FrontstrError):
    """Raised when the input file cannot be ingested (format, headers, integrity)."""


class PanelError(FrontstrError):
    """Raised for malformed panel definitions."""


class EvidenceError(FrontstrError):
    """Raised in the sequence-pileup / clustering stages."""


class InterpretationError(FrontstrError):
    """Raised when forensic classification rules cannot resolve a result."""


class ReportError(FrontstrError):
    """Raised when emitting the HTML / PDF / CSV report fails."""
