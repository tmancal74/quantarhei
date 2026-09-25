from unittest.mock import Mock

import numpy
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


def test_response_container_retains_raw_time_slices_until_conversion():
    t1, _, t3 = _axes()
    t2 = qr.TimeAxis(0.0, 1, 10.0)
    raw_rephasing = numpy.arange(t1.length * t3.length, dtype=complex).reshape(
        t3.length, t1.length
    )
    raw_nonrephasing = 1j * raw_rephasing
    response = qr.TwoDResponse()
    response.t1axis = t1
    response.t3axis = t3
    response.set_t2(t2.data[0])
    response._add_data(raw_rephasing, resolution="signals", dtype=qr.signal_REPH)
    response._add_data(raw_nonrephasing, resolution="signals", dtype=qr.signal_NONR)
    responses = qr.TwoDResponseContainer(t2axis=t2)
    responses.set_spectrum(response)

    stored = responses.get_response(t2.data[0])
    stored.set_data_flag(qr.signal_REPH)
    npt.assert_allclose(stored.d__data, raw_rephasing)
    assert isinstance(stored.t1axis, qr.TimeAxis)
    assert isinstance(stored.t3axis, qr.TimeAxis)

    spectra = qr.TwoDSpectrumCalculator.convert_response_container(responses)
    spectrum = spectra.get_spectrum(t2.data[0])
    expected_rephasing = spectrum_module._fourier_transform_response(
        raw_rephasing, qr.signal_REPH
    )
    expected_nonrephasing = spectrum_module._fourier_transform_response(
        raw_nonrephasing, qr.signal_NONR
    )
    expected = (expected_rephasing + expected_nonrephasing) * (
        t3.length * t1.step * t3.step
    )
    npt.assert_allclose(spectrum.data, expected)
    assert isinstance(spectrum.xaxis, qr.FrequencyAxis)
    assert isinstance(spectrum.yaxis, qr.FrequencyAxis)


def test_frequency_domain_mock_response_is_not_transformed_twice():
    t1, _, t3 = _axes()
    t2 = qr.TimeAxis(0.0, 1, 10.0)
    response = qr.TwoDResponse()
    response.domain = "frequency"
    response.set_axis_1(qr.FrequencyAxis(0.0, t1.length, 1.0))
    response.set_axis_3(qr.FrequencyAxis(0.0, t3.length, 1.0))
    response.set_t2(t2.data[0])
    data = numpy.ones((t3.length, t1.length), dtype=complex)
    response._add_data(data, resolution="signals", dtype=qr.signal_REPH)
    responses = qr.TwoDResponseContainer(t2axis=t2)
    responses.set_spectrum(response)

    spectrum = qr.TwoDSpectrumCalculator.convert_response_container(
        responses, stype=qr.signal_REPH
    ).get_spectrum(t2.data[0])

    npt.assert_allclose(spectrum.data, data)


def test_delta_calculation_wraps_existing_impulsive_conversion():
    t1, t2, t3 = _axes()
    calculator = qr.TwoDSpectrumCalculator(t1, t2, t3, _delta_lab())
    responses = qr.TwoDResponseContainer(t2axis=t2)
    expected = qr.TwoDSpectrumContainer(t2axis=t2)
    calculator.convert_response_container = Mock(return_value=expected)

    calculator.bootstrap(responses)
    result = calculator.calculate(stype=qr.signal_REPH)

    assert result is expected
    calculator.convert_response_container.assert_called_once_with(
        responses, stype=qr.signal_REPH
    )


def test_calculate_requires_bootstrap():
    t1, t2, t3 = _axes()
    calculator = qr.TwoDSpectrumCalculator(t1, t2, t3, _delta_lab())

    with pytest.raises(qr.QuantarheiError, match="bootstrapped"):
        calculator.calculate()


def test_finite_pulses_extend_response_axes_for_explicit_convolution():
    t1 = qr.TimeAxis(0.0, 100, 5.0)
    t2 = qr.TimeAxis(0.0, 2, 100.0)
    t3 = qr.TimeAxis(0.0, 100, 5.0)
    lab = qr.LabSetup(nopulses=3)
    pulse_axis = qr.TimeAxis(-200.0, 81, 5.0, atype="complete")
    pulse = {"ptype": "Gaussian", "FWHM": 30.0, "amplitude": 1.0}
    lab.set_pulse_shapes(pulse_axis, (pulse, pulse, pulse))
    calculator = qr.TwoDSpectrumCalculator(t1, t2, t3, lab)

    response_t1, response_t2, response_t3 = calculator.get_response_axes()

    # An intensity FWHM of 30 fs gives 4 sigma = 72.1 fs for the
    # field envelope.  The pair-delay margin is 144.1 fs, rounded to 145 fs.
    assert response_t1.start == 0.0
    assert response_t1.step == 5.0
    assert response_t1.max == 640.0
    assert response_t3.is_equal_to(response_t1)
    assert response_t2.start == 0.0
    assert response_t2.step == 5.0
    assert response_t2.max == 245.0


def test_finite_pulse_axis_margin_respects_gaussian_fwhm_convention():
    t1, t2, t3 = _axes()
    lab = _finite_lab()
    for parameters in lab.saved_params:
        parameters["FWHM"] = 30.0
        parameters["FWHM_type"] = "amplitude"
    calculator = qr.TwoDSpectrumCalculator(t1, t2, t3, lab)

    response_t1, _, _ = calculator.get_response_axes()

    # The amplitude-FWHM convention has a 101.9 fs pair margin, rounded to
    # the 1 fs t1 grid.
    assert response_t1.max == 109.0


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
    calculator.convert_response_container = Mock(return_value=expected)

    calculator.bootstrap(responses)
    result = calculator.calculate(stype=qr.signal_REPH)

    assert result is expected
    calculator.convert_response_container.assert_called_once_with(
        responses, stype=qr.signal_REPH
    )
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
    calculator.convert_response_container = Mock(return_value=expected)

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
    calculator.convert_response_container = Mock(
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
