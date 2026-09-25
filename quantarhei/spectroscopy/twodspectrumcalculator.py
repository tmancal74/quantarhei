"""Calculation of 2D spectra from nonlinear-response containers."""

from __future__ import annotations

import math
from typing import Any

import numpy

from .. import signal_NONR, signal_REPH, signal_TOTL
from ..core.managers import Manager, energy_units
from ..core.time import TimeAxis
from ..exceptions import QuantarheiError
from ..utils import derived_type
from .labsetup import LabSetup
from .twodcalculator import TwoDResponseCalculator
from .twodcontainer import TwoDResponseContainer, TwoDSpectrumContainer
from .twodresponse import TwoDResponse
from .twodspect import TwoDSpectrum


def _apply_response_window(data: numpy.ndarray) -> numpy.ndarray:
    """Apply the endpoint half-weight used before Fourier transformation."""
    ret = data.copy()
    ret[:, 0] *= 0.5
    ret[0, :] *= 0.5
    return ret


def _fourier_transform_response(data: numpy.ndarray, signal: str) -> numpy.ndarray:
    """Transform one raw response contribution using Quantarhei conventions."""
    data = _apply_response_window(data)
    if signal == signal_REPH:
        transformed = numpy.fft.fft(data, axis=1)
    elif signal == signal_NONR:
        transformed = numpy.fft.ifft(data, axis=1) * data.shape[1]
    else:
        raise QuantarheiError("Unknown 2D signal type: " + signal)
    return numpy.fft.fftshift(numpy.fft.ifft(transformed, axis=0))


def _pad_response_data(
    data: numpy.ndarray, pad: int, window: numpy.ndarray | None = None
) -> numpy.ndarray:
    """Apply optional terminal windowing and zero padding to raw response data."""
    if window is not None:
        size = int(len(window) / 2)
        data = data.copy()
        data[len(data) - size :, :] *= window[size:, None]
        data[:, len(data) - size :] *= window[None, size:]
    if pad > 0:
        data = numpy.hstack((data, numpy.zeros((data.shape[0], pad))))
        data = numpy.vstack((data, numpy.zeros((pad, data.shape[1]))))
    return data


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
    pulse_support_sigma : float, optional
        Effective half-width of a Gaussian pulse, expressed as a number of
        field-envelope standard deviations.  It is used only when suggesting
        response axes for explicit finite-pulse convolution.  The default of
        ``4.0`` retains field amplitudes down to ``exp(-8)`` of their peak.
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
        pulse_support_sigma: float = 4.0,
    ) -> None:
        if pulse_support_sigma <= 0.0:
            raise ValueError("pulse_support_sigma must be positive")
        self.t1axis = t1axis
        self.t2axis = t2axis
        self.t3axis = t3axis
        self.lab = lab
        self.explicit_convolution = explicit_convolution
        self.pulse_support_sigma = pulse_support_sigma
        self.response_container: TwoDResponseContainer | None = None
        self.response_calculator: TwoDResponseCalculator | None = None
        self.system: Any = None
        self.rwa: float | None = None
        self._response_calculator_kwargs: dict[str, Any] = {}
        self._response_bootstrap_kwargs: dict[str, Any] = {}

    def get_response_axes(self) -> tuple[TimeAxis, TimeAxis, TimeAxis]:
        """Return axes suitable for calculating the required responses.

        Delta pulses, and finite pulses used with the approximate overlay,
        use the experimental axes unchanged. Copies are returned so a response
        calculator cannot modify the axes owned by this calculator.

        For explicit finite-pulse convolution, the internal response axes are
        extended beyond the largest requested delay by the maximum difference
        between two effective pulse times.  A time-defined Gaussian is cut at
        ``pulse_support_sigma`` standard deviations of its *field envelope*;
        sampled non-Gaussian pulses use their configured temporal support.
        The extra range is rounded upward to the appropriate grid spacing.
        The waiting-time axis uses the finest supplied time step, because it
        will be integrated during the finite-pulse calculation.
        """
        if self.lab.has_delta_pulses() or not self.explicit_convolution:
            return (
                self.t1axis.deepcopy(),
                self.t2axis.deepcopy(),
                self.t3axis.deepcopy(),
            )

        margin = self._finite_pulse_delay_margin()
        t1axis = self._extend_axis(self.t1axis, margin, self.t1axis.step)
        t3axis = self._extend_axis(self.t3axis, margin, self.t3axis.step)
        quadrature_step = min(self.t1axis.step, self.t2axis.step, self.t3axis.step)
        t2axis = self._extend_axis(self.t2axis, margin, quadrature_step)
        return t1axis, t2axis, t3axis

    @staticmethod
    def _extend_axis(axis: TimeAxis, margin: float, step: float) -> TimeAxis:
        """Extend a causal response axis at its upper end on a chosen grid."""
        upper = axis.max + margin
        length = math.ceil((upper - axis.start) / step - 1.0e-12) + 1
        return TimeAxis(
            axis.start,
            length,
            step,
            atype=axis.atype,
            frequency_start=axis.frequency_start,
        )

    def _finite_pulse_delay_margin(self) -> float:
        """Return a conservative, grid-independent finite-pulse delay margin."""
        params = self.lab.saved_params
        if params is None:
            raise QuantarheiError(
                "Finite-pulse response axes require configured pulses"
            )

        half_supports: list[float] = []
        for index, pulse in enumerate(params):
            if (
                self.lab.pulse_definition_domain == "time"
                and pulse.get("ptype") == "Gaussian"
            ):
                fwhm = float(pulse["FWHM"])
                fwhm_type = pulse.get("FWHM_type", "intensity").lower()
                if fwhm_type == "intensity":
                    exponent_factor = 2.0 * numpy.log(2.0)
                elif fwhm_type in ("amplitude", "field"):
                    exponent_factor = 4.0 * numpy.log(2.0)
                else:
                    raise QuantarheiError("Unknown Gaussian FWHM_type")
                sigma = fwhm / math.sqrt(2.0 * exponent_factor)
                half_supports.append(self.pulse_support_sigma * sigma)
            else:
                half_supports.append(self._sampled_pulse_half_support(index))

        # Any relative delay involves two pulse times.  The sum of the two
        # largest supports is conservative for unequal pulse widths.
        half_supports.sort(reverse=True)
        return half_supports[0] + half_supports[1]

    def _sampled_pulse_half_support(self, index: int) -> float:
        """Return support of a non-analytic pulse around its configured center."""
        if not self.lab.has_timedomain:
            self.lab.convert_to_time()
        pulse = self.lab.pulse_t[index]
        if pulse is None:
            raise QuantarheiError("Finite pulse has no time-domain representation")
        center = self.lab.pulse_centers[index]
        return max(abs(pulse.axis.min - center), abs(pulse.axis.max - center))

    def suggest_response_axes(self) -> tuple[TimeAxis, TimeAxis, TimeAxis]:
        """Alias for :meth:`get_response_axes`."""
        return self.get_response_axes()

    @staticmethod
    def _frequency_axis(axis: TimeAxis, pad: int, rwa: float) -> Any:
        """Build an absolute-frequency axis for a raw response time axis."""
        padded = TimeAxis(axis.start, axis.length + pad, axis.step)
        padded.atype = "complete"
        frequency = padded.get_FrequencyAxis()
        frequency.data += rwa
        frequency.start += rwa
        return frequency

    @classmethod
    def convert_response(
        cls, response: TwoDResponse, stype: Any = signal_TOTL, pad: int = 0
    ) -> TwoDSpectrum:
        """Convert one raw fixed-``t2`` response slice to a 2D spectrum.

        The response retains rephasing and non-rephasing data separately until
        this point because their Fourier conventions differ.
        """
        if response.domain == "frequency":
            return response.get_TwoDSpectrum(dtype=stype)

        t1axis = response.t1axis
        t3axis = response.t3axis
        if t1axis is None or t3axis is None:
            raise QuantarheiError("Time axes of the 2D response are not set")

        requested = (signal_REPH, signal_NONR) if stype == signal_TOTL else (stype,)
        data_out: numpy.ndarray | None = None
        window = None
        if pad > 0:
            from scipy.signal import windows as sig

            window = sig.tukey(40, 1, sym=False)

        for signal in requested:
            if signal not in (signal_REPH, signal_NONR):
                raise QuantarheiError("Unknown 2D signal type: " + str(signal))
            response.set_data_flag(signal)
            raw_data = response.d__data
            if raw_data is None:
                continue
            padded_data = _pad_response_data(raw_data, pad, window)
            transformed = _fourier_transform_response(padded_data, signal)
            transformed *= padded_data.shape[0] * t1axis.step * t3axis.step
            if data_out is None:
                data_out = transformed
            else:
                data_out += transformed

        if data_out is None:
            raise QuantarheiError("Response has no data for the requested signal")

        spectrum = TwoDSpectrum()
        spectrum.set_axis_1(cls._frequency_axis(t1axis, pad, response.rwa))
        spectrum.set_axis_3(cls._frequency_axis(t3axis, pad, response.rwa))
        spectrum.rwa = response.rwa
        spectrum.set_t2(response.t2)
        spectrum.set_data(data_out, dtype=stype)
        return spectrum

    @classmethod
    def convert_response_container(
        cls, responses: TwoDResponseContainer, stype: Any = signal_TOTL
    ) -> TwoDSpectrumContainer:
        """Convert all raw response slices in a container to spectra."""
        if responses.itype not in ("ValueAxis", "TimeAxis", "FrequencyAxis"):
            raise QuantarheiError("Response container must be indexed by an axis")
        assert responses.axis is not None
        spectra = TwoDSpectrumContainer(responses.axis.deepcopy())
        pad = getattr(responses, "pad", 0)
        for value in responses.axis.data:
            response = responses.get_response(value)
            spectrum = cls.convert_response(response, stype=stype, pad=pad)
            spectra.set_spectrum(spectrum, tag=value)
        return spectra

    def bootstrap(
        self,
        sample: TwoDResponseContainer | Any,
        *,
        rwa: float | None = None,
        response_calculator_kwargs: dict[str, Any] | None = None,
        response_bootstrap_kwargs: dict[str, Any] | None = None,
    ) -> None:
        """Mount a sample or install pre-calculated responses.

        ``sample`` is normally an Aggregate or OpenSystem.  In that workflow,
        this calculator owns the experimental axes and laboratory setup, and
        creates a :class:`TwoDResponseCalculator` when :meth:`calculate` is
        called.  This keeps a usual calculation close to the experimental
        sequence: configure the laboratory, choose measurement times, mount
        the sample, and measure.

        A :class:`TwoDResponseContainer` remains accepted for advanced use:
        callers can calculate, inspect, or modify responses separately and
        then use this calculator for the spectrum step.

        Parameters
        ----------
        response_calculator_kwargs
            Keyword arguments forwarded to ``TwoDResponseCalculator`` when a
            system is supplied, e.g. relaxation settings.
        response_bootstrap_kwargs
            Keyword arguments forwarded to its ``bootstrap`` method.  ``lab``
            is supplied by this calculator and may not be overridden.
        rwa
            Optional rotating-wave reference frequency in the active energy
            units. If omitted, the system's RWA suggestion is used.
        """
        if isinstance(sample, TwoDResponseContainer):
            axis = sample.axis
            if axis is None or not axis.is_equal_to(self.t2axis):
                raise ValueError(
                    "Response-container waiting-time axis does not match t2axis"
                )
            self.response_container = sample
            self.response_calculator = None
            self.system = None
            self.rwa = None
            return

        bootstrap_kwargs = dict(response_bootstrap_kwargs or {})
        if "lab" in bootstrap_kwargs:
            raise ValueError(
                "The laboratory setup is owned by TwoDSpectrumCalculator; "
                "do not pass 'lab' in response_bootstrap_kwargs"
            )
        if "rwa" in bootstrap_kwargs:
            raise ValueError(
                "Pass 'rwa' directly to TwoDSpectrumCalculator.bootstrap(), "
                "so its active energy units are recorded at bootstrap time"
            )

        self.system = sample
        self.response_container = None
        self.response_calculator = None
        self.rwa = None if rwa is None else Manager().convert_energy_2_internal_u(rwa)
        self._response_calculator_kwargs = dict(response_calculator_kwargs or {})
        self._response_bootstrap_kwargs = bootstrap_kwargs

    def _calculate_responses(self) -> TwoDResponseContainer:
        """Calculate responses for the system installed during bootstrap."""
        if self.system is None:
            raise QuantarheiError(
                "TwoDSpectrumCalculator must be bootstrapped before calculation"
            )

        t1axis, t2axis, t3axis = self.get_response_axes()
        calculator = TwoDResponseCalculator(
            t1axis,
            t2axis,
            t3axis,
            system=self.system,
            **self._response_calculator_kwargs,
        )
        with energy_units("int"):
            rwa = self.rwa
            if rwa is None:
                rwa = self.system.get_RWA_suggestion()
            calculator.bootstrap(
                rwa=rwa, lab=self.lab, **self._response_bootstrap_kwargs
            )
        responses = calculator.calculate()
        self.response_calculator = calculator
        self.response_container = responses
        return responses

    def calculate(self, stype: Any = signal_TOTL) -> TwoDSpectrumContainer:
        """Calculate and return a container of 2D spectra.

        Delta pulses use the established impulsive response-to-spectrum
        conversion. For finite pulses, ``explicit_convolution=False`` applies
        a post-processing spectral overlay; explicit convolution is reserved
        for a later implementation.
        """
        if not self.lab.has_delta_pulses() and self.explicit_convolution:
            raise NotImplementedError("Finite-pulse convolution is not implemented")

        responses = self.response_container
        if responses is None:
            responses = self._calculate_responses()

        spectra = self.convert_response_container(responses, stype=stype)
        if not self.lab.has_delta_pulses():
            for spectrum in spectra.spectra.values():
                spectrum.overlay_pulses(self.lab)

        return spectra
