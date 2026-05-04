"""Class representing absorption spectrum."""

from .absbase import AbsSpectrumBase


class AbsSpectrum(AbsSpectrumBase):
    """Class representing absorption spectrum.

    Inherits all functionality from ``AbsSpectrumBase``.

    Parameters
    ----------
    axis : FrequencyAxis
        Frequency axis for the spectrum.
    data : numpy.ndarray
        Spectral data array.
    """

    pass
