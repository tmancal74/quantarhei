"""Two-level molecule 2D regression spectra, represented by a one-site aggregate."""

from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest

import quantarhei as qr
from quantarhei.utils.vectors import X

TEST_DIR = Path(__file__).parent


def calculate_tlm_spectra(underdamped=False):
    """Calculate the same bath/time-grid cases as the Chlorophyll regression."""
    t1 = qr.TimeAxis(0.0, 50, 5.0)
    t2 = qr.TimeAxis(0.0, 2, 100.0)
    t3 = qr.TimeAxis(0.0, 50, 5.0)
    with qr.energy_units("1/cm"):
        molecule = qr.Molecule([0.0, 16807.0])
        bath = qr.CorrelationFunction(
            t1,
            dict(
                ftype="OverdampedBrownian",
                reorg=140.0,
                cortime=60.0,
                T=300.0,
                matsubara=20,
            ),
        )
        if underdamped:
            # S = lambda / omega = 0.1; amplitude damping time = 3 ps.
            bath += qr.CorrelationFunction(
                t1,
                dict(
                    ftype="UnderdampedBrownian",
                    reorg=150.0,
                    freq=1500.0,
                    gamma=qr.convert(2.0 / 3000.0, "int", "1/cm"),
                    T=300.0,
                ),
            )
    molecule.set_dipole(0, 1, [1.0, 0.0, 0.0])
    molecule.set_transition_environment((0, 1), bath)
    # Direct Molecule bootstrap is not supported by the current calculator.
    system = qr.Aggregate(molecules=[molecule])
    system.build(mult=2)
    lab = qr.LabSetup()
    lab.set_pulse_polarizations(pulse_polarizations=(X, X, X), detection_polarization=X)
    calc = qr.TwoDResponseCalculator(t1, t2, t3, system=system)
    with qr.energy_units("1/cm"):
        calc.bootstrap(rwa=16807.0, pad=0, lab=lab)
    return calc.calculate().get_TwoDSpectrumContainer()


@pytest.fixture(scope="module", params=[False, True], ids=["overdamped", "underdamped"])
def tlm_spectra(request):
    return request.param, calculate_tlm_spectra(underdamped=request.param)


@pytest.mark.parametrize("t2", [0, 100])
def test_tlm_reference_spectrum(tlm_spectra, t2):
    """Compare finite, nonzero TLM spectra with stored complex reference data."""
    underdamped, container = tlm_spectra
    spectrum = container.get_spectrum(float(t2))
    assert spectrum.data.shape == (50, 50)
    assert spectrum.get_t2() == t2
    assert np.all(np.isfinite(spectrum.data))
    assert np.max(np.abs(spectrum.data)) > 0.0
    bath = "underdamped" if underdamped else "overdamped"
    reference = np.loadtxt(TEST_DIR / f"twod_tlm_{bath}_data_{t2}.dat", dtype=complex)
    # References include the current FFT integral normalization.
    npt.assert_allclose(spectrum.data, reference, rtol=1.0e-7, atol=0.0)
