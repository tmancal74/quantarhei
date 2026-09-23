from unittest.mock import Mock

import numpy.testing as npt
import pytest

import quantarhei as qr


def _delta_lab() -> qr.LabSetup:
    lab = qr.LabSetup(nopulses=3)
    pulse_axis = qr.TimeAxis(-10.0, 21, 1.0, atype="complete")
    pulse = {"ptype": "delta", "area": 1.0}
    lab.set_pulse_shapes(pulse_axis, (pulse, pulse, pulse))
    return lab


def _axes():
    return (
        qr.TimeAxis(0.0, 8, 1.0),
        qr.TimeAxis(0.0, 3, 10.0),
        qr.TimeAxis(0.0, 8, 1.0),
    )


def test_delta_pulse_representation_has_unit_area():
    lab = _delta_lab()

    assert lab.has_delta_pulses()
    for pulse in lab.pulse_t:
        npt.assert_allclose(pulse.data.sum() * pulse.axis.step, 1.0)


def test_delta_pulses_suggest_unchanged_response_axes():
    t1, t2, t3 = _axes()
    calculator = qr.TwoDSpectrumCalculator(t1, t2, t3, _delta_lab())

    response_axes = calculator.get_response_axes()

    for suggested, experimental in zip(response_axes, (t1, t2, t3)):
        assert suggested.is_equal_to(experimental)
        assert suggested is not experimental


def test_delta_calculation_wraps_existing_impulsive_conversion():
    t1, t2, t3 = _axes()
    calculator = qr.TwoDSpectrumCalculator(t1, t2, t3, _delta_lab())
    responses = qr.TwoDResponseContainer(t2axis=t2)
    expected = qr.TwoDSpectrumContainer(t2axis=t2)
    responses.get_TwoDSpectrumContainer = Mock(return_value=expected)

    calculator.bootstrap(responses)
    result = calculator.calculate(stype=qr.signal_REPH)

    assert result is expected
    responses.get_TwoDSpectrumContainer.assert_called_once_with(stype=qr.signal_REPH)


def test_calculate_requires_bootstrap():
    t1, t2, t3 = _axes()
    calculator = qr.TwoDSpectrumCalculator(t1, t2, t3, _delta_lab())

    with pytest.raises(qr.QuantarheiError, match="bootstrapped"):
        calculator.calculate()


def test_finite_pulses_are_explicitly_deferred():
    t1, t2, t3 = _axes()
    lab = qr.LabSetup(nopulses=3)
    pulse_axis = qr.TimeAxis(-50.0, 101, 1.0, atype="complete")
    pulse = {"ptype": "Gaussian", "FWHM": 15.0, "amplitude": 1.0}
    lab.set_pulse_shapes(pulse_axis, (pulse, pulse, pulse))
    calculator = qr.TwoDSpectrumCalculator(t1, t2, t3, lab)

    with pytest.raises(NotImplementedError, match="finite pulses"):
        calculator.get_response_axes()


def test_bootstrap_checks_waiting_time_axis():
    t1, t2, t3 = _axes()
    calculator = qr.TwoDSpectrumCalculator(t1, t2, t3, _delta_lab())
    other_t2 = qr.TimeAxis(0.0, 4, 10.0)
    responses = qr.TwoDResponseContainer(t2axis=other_t2)

    with pytest.raises(ValueError, match="waiting-time axis"):
        calculator.bootstrap(responses)
