"""Placeholder absorption for systems without molecular response machinery."""

import numpy as np
import pytest

import quantarhei as qr
from quantarhei.spectroscopy.abscalculator import LinSpectrumCalculator


class MinimalOpenSystem(qr.OpenSystem):
    def get_Hamiltonian(self):
        return qr.Hamiltonian(data=np.diag([0.0, 2.0]))


@pytest.mark.parametrize(
    "calculator_type", [LinSpectrumCalculator, qr.AbsSpectrumCalculator]
)
@pytest.mark.parametrize("raw", [False, True])
def test_empty_open_system_absorption(calculator_type, raw):
    time = qr.TimeAxis(0.0, 32, 2.0)
    calc = calculator_type(time, system=MinimalOpenSystem([0.0, 2.0]))
    with qr.energy_units("int"):
        calc.bootstrap(rwa=2.0)
        spectrum = calc.calculate(raw=raw)
        assert isinstance(spectrum, qr.AbsSpectrum)
        assert spectrum.data.shape == (32,)
        np.testing.assert_array_equal(spectrum.data, np.zeros(32))
        np.testing.assert_allclose(spectrum.axis.data, calc.frequencyAxis.data[16:48])
