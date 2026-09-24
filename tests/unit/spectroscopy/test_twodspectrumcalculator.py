from unittest.mock import Mock

import numpy.testing as npt
import pytest

import quantarhei as qr
import quantarhei.spectroscopy.twodspectrumcalculator as spectrum_module


def _delta_lab() -> qr.LabSetup:
    lab = qr.LabSetup(nopulses=3)
    pulse_axis = qr.TimeAxis(-10.0, 21, 1.0, atype="complete")
    pulse = {"ptype": "delta", "area": 1.0}
    lab.set_pulse_shapes(pulse_axis, (pulse, pulse, pulse))
    return lab


def _finite_lab() -> qr.LabSetup:
    lab = qr.LabSetup(nopulses=3)
    pulse_axis = qr.TimeAxis(-50.0, 101, 1.0, atype="complete")
    pulse = {"ptype": "Gaussian", "FWHM": 15.0, "amplitude": 1.0}
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
    calculator = qr.TwoDSpectrumCalculator(t1, t2, t3, _finite_lab())

    with pytest.raises(NotImplementedError, match="finite pulses"):
        calculator.get_response_axes()


def test_finite_pulse_overlay_suggests_unchanged_response_axes():
    t1, t2, t3 = _axes()
    calculator = qr.TwoDSpectrumCalculator(
        t1, t2, t3, _finite_lab(), explicit_convolution=False
    )

    response_axes = calculator.get_response_axes()

    for suggested, experimental in zip(response_axes, (t1, t2, t3)):
        assert suggested.is_equal_to(experimental)
        assert suggested is not experimental


def test_finite_pulse_overlay_is_applied_after_impulsive_conversion():
    t1, t2, t3 = _axes()
    lab = _finite_lab()
    calculator = qr.TwoDSpectrumCalculator(t1, t2, t3, lab, explicit_convolution=False)
    responses = qr.TwoDResponseContainer(t2axis=t2)
    expected = qr.TwoDSpectrumContainer(t2axis=t2)
    spectrum = Mock()
    expected.spectra = {t2.data[0]: spectrum}
    responses.get_TwoDSpectrumContainer = Mock(return_value=expected)

    calculator.bootstrap(responses)
    result = calculator.calculate(stype=qr.signal_REPH)

    assert result is expected
    responses.get_TwoDSpectrumContainer.assert_called_once_with(stype=qr.signal_REPH)
    spectrum.overlay_pulses.assert_called_once_with(lab)


def test_bootstrap_checks_waiting_time_axis():
    t1, t2, t3 = _axes()
    calculator = qr.TwoDSpectrumCalculator(t1, t2, t3, _delta_lab())
    other_t2 = qr.TimeAxis(0.0, 4, 10.0)
    responses = qr.TwoDResponseContainer(t2axis=other_t2)

    with pytest.raises(ValueError, match="waiting-time axis"):
        calculator.bootstrap(responses)


def test_system_bootstrap_calculates_responses_with_owned_lab_and_axes(monkeypatch):
    t1, t2, t3 = _axes()
    lab = _delta_lab()
    calculator = qr.TwoDSpectrumCalculator(t1, t2, t3, lab)

    class System:
        def get_RWA_suggestion(self):
            return 3.5

    system = System()
    responses = qr.TwoDResponseContainer(t2axis=t2)
    expected = qr.TwoDSpectrumContainer(t2axis=t2)
    responses.get_TwoDSpectrumContainer = Mock(return_value=expected)

    created = []

    class ResponseCalculator:
        def __init__(self, t1axis, t2axis, t3axis, *, system, dynamics):
            self.axes = (t1axis, t2axis, t3axis)
            self.system = system
            self.dynamics = dynamics
            self.bootstrap = Mock()
            self.calculate = Mock(return_value=responses)
            created.append(self)

    monkeypatch.setattr(spectrum_module, "TwoDResponseCalculator", ResponseCalculator)
    calculator.bootstrap(
        system,
        response_calculator_kwargs={"dynamics": "full"},
        response_bootstrap_kwargs={"pad": 4},
    )

    result = calculator.calculate()

    assert result is expected
    assert len(created) == 1
    response_calculator = created[0]
    assert response_calculator.system is system
    assert response_calculator.dynamics == "full"
    for supplied, owned in zip(response_calculator.axes, (t1, t2, t3)):
        assert supplied.is_equal_to(owned)
        assert supplied is not owned
    response_calculator.bootstrap.assert_called_once_with(rwa=3.5, lab=lab, pad=4)
    response_calculator.calculate.assert_called_once_with()
    assert calculator.response_calculator is response_calculator
    assert calculator.response_container is responses


def test_system_bootstrap_reserves_lab_for_spectrum_calculator():
    t1, t2, t3 = _axes()
    calculator = qr.TwoDSpectrumCalculator(t1, t2, t3, _delta_lab())

    with pytest.raises(ValueError, match="laboratory setup is owned"):
        calculator.bootstrap(object(), response_bootstrap_kwargs={"lab": _delta_lab()})


def test_system_bootstrap_records_explicit_rwa_in_active_units(monkeypatch):
    t1, t2, t3 = _axes()
    lab = _delta_lab()
    calculator = qr.TwoDSpectrumCalculator(t1, t2, t3, lab)
    responses = qr.TwoDResponseContainer(t2axis=t2)
    responses.get_TwoDSpectrumContainer = Mock(
        return_value=qr.TwoDSpectrumContainer(t2axis=t2)
    )

    class System:
        def get_RWA_suggestion(self):
            raise AssertionError("An explicit RWA should take precedence")

    response_calculator = Mock()
    response_calculator.calculate.return_value = responses
    monkeypatch.setattr(
        spectrum_module,
        "TwoDResponseCalculator",
        Mock(return_value=response_calculator),
    )

    with qr.energy_units("1/cm"):
        calculator.bootstrap(System(), rwa=12000.0)
    calculator.calculate()

    with qr.energy_units("1/cm"):
        expected_rwa = qr.convert(12000.0, "1/cm", "int")
    response_calculator.bootstrap.assert_called_once_with(rwa=expected_rwa, lab=lab)
