"""Calculation of 2D spectra from nonlinear-response containers."""

from __future__ import annotations

from typing import Any

from .. import signal_TOTL
from ..core.time import TimeAxis
from ..exceptions import QuantarheiError
from ..utils import derived_type
from .labsetup import LabSetup
from .twodcontainer import TwoDResponseContainer, TwoDSpectrumContainer


class TwoDSpectrumCalculator:
    """Calculate finite- or delta-pulse spectra from 2D responses.

    This initial implementation establishes the interface between response
    calculation and spectrum calculation. Delta pulses require no convolution,
    so :meth:`calculate` delegates to the existing impulsive conversion on the
    bootstrapped response container. Finite-pulse convolution will be added
    behind the same interface.

    Parameters
    ----------
    t1axis : TimeAxis
        Experimental coherence-time axis.
    t2axis : TimeAxis
        Experimental waiting-time axis.
    t3axis : TimeAxis
        Experimental detection-time axis.
    lab : LabSetup
        Laboratory setup containing three configured pulse shapes.
    explicit_convolution : bool, optional
        If ``True`` (default), finite pulses require an explicit convolution,
        which is not implemented yet. If ``False``, calculate the impulsive
        spectrum on the supplied axes and apply the approximate spectral-pulse
        overlay to the resulting spectra.
    """

    t1axis = derived_type("t1axis", TimeAxis)
    t2axis = derived_type("t2axis", TimeAxis)
    t3axis = derived_type("t3axis", TimeAxis)
    lab = derived_type("lab", LabSetup)

    def __init__(
        self,
        t1axis: TimeAxis,
        t2axis: TimeAxis,
        t3axis: TimeAxis,
        lab: LabSetup,
        explicit_convolution: bool = True,
    ) -> None:
        self.t1axis = t1axis
        self.t2axis = t2axis
        self.t3axis = t3axis
        self.lab = lab
        self.explicit_convolution = explicit_convolution
        self.response_container: TwoDResponseContainer | None = None

    def get_response_axes(self) -> tuple[TimeAxis, TimeAxis, TimeAxis]:
        """Return axes suitable for calculating the required responses.

        Delta pulses, and finite pulses used with the approximate overlay,
        use the experimental axes unchanged. Copies are returned so a response
        calculator cannot modify the axes owned by this calculator. Explicit
        finite-pulse convolution will extend these axes in a later
        implementation.
        """
        if not self.lab.has_delta_pulses() and self.explicit_convolution:
            raise NotImplementedError(
                "Response-axis suggestions for finite pulses are not implemented"
            )
        return (
            self.t1axis.deepcopy(),
            self.t2axis.deepcopy(),
            self.t3axis.deepcopy(),
        )

    def suggest_response_axes(self) -> tuple[TimeAxis, TimeAxis, TimeAxis]:
        """Alias for :meth:`get_response_axes`."""
        return self.get_response_axes()

    def bootstrap(self, response_container: TwoDResponseContainer) -> None:
        """Set the response container from which spectra will be calculated."""
        if not isinstance(response_container, TwoDResponseContainer):
            raise TypeError("response_container must be a TwoDResponseContainer")

        axis = response_container.axis
        if axis is None or not axis.is_equal_to(self.t2axis):
            raise ValueError(
                "Response-container waiting-time axis does not match t2axis"
            )
        self.response_container = response_container

    def calculate(self, stype: Any = signal_TOTL) -> TwoDSpectrumContainer:
        """Calculate and return a container of 2D spectra.

        Delta pulses use the established impulsive response-to-spectrum
        conversion. For finite pulses, ``explicit_convolution=False`` applies
        a post-processing spectral overlay; explicit convolution is reserved
        for a later implementation.
        """
        if self.response_container is None:
            raise QuantarheiError(
                "TwoDSpectrumCalculator must be bootstrapped before calculation"
            )
        if not self.lab.has_delta_pulses() and self.explicit_convolution:
            raise NotImplementedError("Finite-pulse convolution is not implemented")

        spectra = self.response_container.get_TwoDSpectrumContainer(stype=stype)
        if not self.lab.has_delta_pulses():
            for spectrum in spectra.spectra.values():
                spectrum.overlay_pulses(self.lab)

        return spectra
