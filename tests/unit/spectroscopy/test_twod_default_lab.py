"""TwoDResponseCalculator.bootstrap without an explicit LabSetup (issue #335)."""

import numpy.testing as npt

import quantarhei as qr
from quantarhei.utils.vectors import X


def _calculator():
    t1 = qr.TimeAxis(0.0, 20, 10.0)
    t2 = qr.TimeAxis(0.0, 2, 50.0)
    t3 = qr.TimeAxis(0.0, 20, 10.0)
    bath_axis = qr.TwoDResponseCalculator(t1, t2, t3).get_joint_time_axis()
    with qr.energy_units("1/cm"):
        m1 = qr.Molecule([0.0, 12000.0])
        m2 = qr.Molecule([0.0, 12300.0])
        params = dict(ftype="OverdampedBrownian", reorg=40.0, cortime=100.0, T=100)
        for m in (m1, m2):
            m.set_transition_environment(
                (0, 1), qr.CorrelationFunction(bath_axis, params)
            )
    m1.set_dipole(0, 1, [1.0, 0.8, 0.8])
    m2.set_dipole(0, 1, [0.8, 0.8, 0.0])
    agg = qr.Aggregate(molecules=[m1, m2])
    with qr.energy_units("1/cm"):
        agg.set_resonance_coupling(0, 1, 100.0)
    agg.build(mult=2)
    return qr.TwoDResponseCalculator(t1, t2, t3, system=agg)


def _spectra(lab):
    calc = _calculator()
    with qr.energy_units("1/cm"):
        calc.bootstrap(rwa=12100.0, lab=lab)
    return calc.calculate().get_TwoDSpectrumContainer()


def test_default_lab_matches_explicit_xxxx_lab():
    """Default lab is all-parallel X polarization, identical to explicit XXXX."""
    lab = qr.LabSetup()
    lab.set_pulse_polarizations(pulse_polarizations=(X, X, X), detection_polarization=X)
    default = _spectra(None)
    explicit = _spectra(lab)
    for t2 in (0.0, 50.0):
        data = default.get_spectrum(t2).data
        assert abs(data).max() > 0.0
        # same code path and inputs; only floating-point noise is allowed
        npt.assert_allclose(data, explicit.get_spectrum(t2).data, rtol=1e-12, atol=0)
