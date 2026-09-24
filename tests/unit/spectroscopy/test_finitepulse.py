import numpy
import numpy.testing as npt

import quantarhei as qr
from quantarhei.spectroscopy.finitepulse import FinitePulseConvolver


def test_rephasing_reference_kernel_uses_corrected_pulse_two_argument():
    response_axes = tuple(qr.TimeAxis(0.0, 1, 1.0) for _ in range(3))
    pulse_axis = qr.TimeAxis(-5.0, 11, 1.0, atype="complete")
    lab = qr.LabSetup(nopulses=3)
    one = qr.DFunction(pulse_axis, numpy.ones(pulse_axis.length))
    ramp = qr.DFunction(pulse_axis, pulse_axis.data.astype(complex))
    lab.set_pulse_shapes(
        pulse_axis,
        (
            {"ptype": "numeric", "function": one},
            {"ptype": "numeric", "function": ramp},
            {"ptype": "numeric", "function": one},
        ),
    )
    convolver = FinitePulseConvolver(*response_axes, lab)

    # With all response times zero and zero carrier frequencies, both
    # rephasing orderings equal E2(t + T).  The handwritten correction removes
    # tau from that argument: t=0, T=1 and tau=2 therefore gives 2, not 6.
    value = convolver.convolve_rephasing(numpy.ones((1, 1, 1)), 0.0, 1.0, 2.0)
    npt.assert_allclose(value, 2.0)


def test_rephasing_reference_kernel_uses_three_dimensional_quadrature_volume():
    axes = tuple(qr.TimeAxis(0.0, 2, 1.0) for _ in range(3))
    pulse_axis = qr.TimeAxis(-10.0, 21, 1.0, atype="complete")
    lab = qr.LabSetup(nopulses=3)
    one = qr.DFunction(pulse_axis, numpy.ones(pulse_axis.length))
    lab.set_pulse_shapes(
        pulse_axis, tuple({"ptype": "numeric", "function": one} for _ in range(3))
    )
    convolver = FinitePulseConvolver(*axes, lab)

    # Eight response-grid points, two orderings, and a unit integration volume.
    value = convolver.convolve_rephasing(numpy.ones((2, 2, 2)), 1.0, 1.0, 1.0)
    npt.assert_allclose(value, 16.0)


def test_rephasing_kernel_uses_carrier_detunings_in_the_response_rwa_frame():
    axes = tuple(qr.TimeAxis(0.0, 2, 1.0) for _ in range(3))
    pulse_axis = qr.TimeAxis(-10.0, 21, 1.0, atype="complete")
    lab = qr.LabSetup(nopulses=3)
    lab.set_pulse_frequencies([10.0, 10.0, 10.0])
    one = qr.DFunction(pulse_axis, numpy.ones(pulse_axis.length))
    lab.set_pulse_shapes(
        pulse_axis,
        tuple({"ptype": "numeric", "function": one} for _ in range(3)),
    )

    # A response calculated in the omega_RWA=10 frame is slowly varying.  If
    # the pulse carrier is also 10, all residual phase factors must be unity.
    convolver = FinitePulseConvolver(*axes, lab, rwa_frequency=10.0)
    value = convolver.convolve_rephasing(numpy.ones((2, 2, 2)), 1.0, 1.0, 1.0)
    npt.assert_allclose(value, 16.0)
