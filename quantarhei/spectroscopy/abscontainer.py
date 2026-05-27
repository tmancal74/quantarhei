from __future__ import annotations

from typing import Any

from ..core.saveable import Saveable
from ..exceptions import QuantarheiError


class AbsSpectrumContainer(Saveable):
    """Container for a set of absorption spectra sharing a common frequency axis.

    Parameters
    ----------
    axis : FrequencyAxis, optional
        Shared frequency axis for all stored spectra. Can be set later via
        ``set_axis``.
    """

    def __init__(self, axis: Any = None) -> None:

        self.axis = axis
        self.count = 0
        self.spectra: dict[str, Any] = {}

    def set_axis(self, axis: Any) -> None:
        self.axis = axis

    def set_spectrum(self, spect: Any, tag: Any = None) -> None:
        """Store an absorption spectrum, checking axis compatibility.

        Parameters
        ----------
        spect : AbsSpectrum
            Spectrum object to store. Its frequency axis must match the
            container's axis (or becomes the axis if none is set yet).
        tag : str or None, optional
            Key under which to store the spectrum. If ``None``, the current
            count value is used as the string key.

        Raises
        ------
        Exception
            If the spectrum's frequency axis is incompatible with the
            container's existing axis.
        """
        frq = spect.axis

        if self.axis is None:
            self.axis = frq

        if self.axis.is_equal_to(frq):
            if tag is None:
                tag1 = str(self.count)
            else:
                tag1 = str(tag)
            self.spectra[tag1] = spect
            self.count += 1
        else:
            raise QuantarheiError("Incompatible time axis (equal axis required)")

    def get_spectrum(self, tag: Any) -> Any:
        """Return the spectrum identified by tag.

        Parameters
        ----------
        tag : str or int
            Key identifying the stored spectrum. Integer tags are converted
            to strings.

        Returns
        -------
        AbsSpectrum
            The stored spectrum corresponding to ``tag``.

        Raises
        ------
        Exception
            If no spectrum is stored under ``tag``.
        """
        if not isinstance(tag, str):
            tag = str(tag)

        if tag in self.spectra.keys():
            return self.spectra[tag]
        raise QuantarheiError("Unknown spectrum")

    def get_spectra(self) -> list[Any]:
        """Return all stored spectra sorted by their string tags.

        Returns
        -------
        list
            List of stored spectra in tag-sorted order.
        """
        ven = [value for (key, value) in sorted(self.spectra.items())]
        return ven
