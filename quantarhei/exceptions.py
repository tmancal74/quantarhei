"""Custom exception hierarchy for quantarhei.

All public exceptions inherit from :class:`QuantarheiError` so that users can
catch any quantarhei-specific error with a single ``except QuantarheiError``
clause while still being able to handle more specific cases.
"""


class QuantarheiError(Exception):
    """Base class for all quantarhei exceptions."""


class UnitsError(QuantarheiError):
    """Raised when an unknown or incompatible unit is encountered."""


class BasisError(QuantarheiError):
    """Raised on basis/eigenbasis related errors."""


class ConfigurationError(QuantarheiError):
    """Raised when configuration cannot be created or is invalid."""


class BuildError(QuantarheiError):
    """Raised when an object must be built before use, or build fails."""


class ImplementationError(QuantarheiError):
    """Raised when functionality is not yet implemented."""
