"""Reference finite-pulse convolution of third-order responses.

This module deliberately starts with direct quadrature.  It is the numerical
reference against which later interpolated or accelerated implementations must
be tested.
"""

from __future__ import annotations

import numpy

from ..core.managers import Manager
from ..core.time import TimeAxis
from ..exceptions import QuantarheiError
from .labsetup import LabSetup


class FinitePulseConvolver:
    """Direct RWA convolution of a rephasing third-order response.

    The input response is sampled as ``response[t3, t2, t1]``.  The
    implementation follows the two rephasing pulse-ordering terms of Eq.
    (4.40) in the supplied reference page.  It is intentionally a
    correctness-first reference implementation; it neither interpolates the
    response nor applies a Fourier transform.
    """

    def __init__(
        self,
        t1axis: TimeAxis,
        t2axis: TimeAxis,
        t3axis: TimeAxis,
        lab: LabSetup,
        rwa_frequency: float = 0.0,
    ) -> None:
        if lab.number_of_pulses != 3:
            raise QuantarheiError("FinitePulseConvolver requires exactly three pulses")
        if lab.has_delta_pulses():
            raise QuantarheiError("Delta pulses use the impulsive response workflow")
        if min(t1axis.start, t2axis.start, t3axis.start) < 0.0:
            raise QuantarheiError("Response time axes must start at zero or later")
        self.t1axis = t1axis
        self.t2axis = t2axis
        self.t3axis = t3axis
        self.lab = lab
        # NonLinearResponse stores coherence phases relative to this frequency.
        # The pulse factors must use the same rotating frame so their phases
        # remain slow rather than restoring an optical-frequency oscillation.
        self.rwa = Manager().convert_energy_2_internal_u(rwa_frequency)

    def convolve_rephasing(
        self, response: numpy.ndarray, t: float, T: float, tau: float
    ) -> complex:
        """Return the rephasing RWA polarization at ``(t, T, tau)``.

        ``t``, ``T`` and ``tau`` are the detection, population and coherence
        delays of the RWA expression.  The non-rephasing and double-coherence
        kernels remain separate future work because they have different field
        signatures.
        """
        expected_shape = (
            self.t3axis.length,
            self.t2axis.length,
            self.t1axis.length,
        )
        if response.shape != expected_shape:
            raise ValueError(
                "response must have shape (t3axis.length, t2axis.length, t1axis.length)"
            )

        t3 = self.t3axis.data[:, None, None]
        t2 = self.t2axis.data[None, :, None]
        t1 = self.t1axis.data[None, None, :]
        omega1, omega2, omega3 = self.lab.omega - self.rwa
        field1, field2, field3 = self.lab.get_labfields()

        e1m = numpy.conj(field1.envelope_at(t + T + tau - t3 - t2 - t1))
        first_ordering = (
            field2.envelope_at(t + T - t3 - t2)
            * field3.envelope_at(t - t3)
            * numpy.exp(
                -1j * (omega1 - omega2 - omega3) * t3
                - 1j * (omega1 - omega2) * t2
                - 1j * omega1 * t1
            )
        )
        second_ordering = (
            field3.envelope_at(t - t3 - t2)
            * field2.envelope_at(t + T - t3 - t1)
            * numpy.exp(
                -1j * (omega1 - omega2 - omega3) * t3
                - 1j * (omega1 - omega3) * t2
                - 1j * omega3 * t1
            )
        )
        prefactor = numpy.exp(
            1j * (omega1 - omega2 - omega3) * t
            + 1j * (omega1 - omega2) * T
            + 1j * omega1 * tau
        )
        volume = self.t1axis.step * self.t2axis.step * self.t3axis.step
        return complex(
            prefactor
            * numpy.sum(response * e1m * (first_ordering + second_ordering))
            * volume
        )
