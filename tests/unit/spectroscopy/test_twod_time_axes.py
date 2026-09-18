"""Joint bath grid and delayed system construction checks."""

import numpy as np
import pytest

import quantarhei as qr


def test_joint_axis():
    calc = qr.TwoDResponseCalculator(
        qr.TimeAxis(0, 50, 5), qr.TimeAxis(0, 2, 100), qr.TimeAxis(0, 50, 5)
    )
    axis = calc.get_joint_time_axis()
    assert (axis.start, axis.step, axis.length) == (0, 5, 119)
    assert axis.data[-1] == 590


@pytest.mark.parametrize(
    "axes",
    [
        ((5, 3, 5), (0, 2, 10), (0, 3, 5)),
        ((0, 3, 5), (1, 2, 10), (0, 3, 5)),
        ((0, 3, 5), (0, 2, 7), (0, 3, 5)),
        ((0, 3, 5), (-5, 2, 10), (0, 3, 5)),
        ((0, 3, 5), (0, 2, 10), (5, 3, 5)),
    ],
)
def test_invalid_axes_rejected_in_constructor(axes):
    with pytest.raises(ValueError):
        qr.TwoDResponseCalculator(*(qr.TimeAxis(*args) for args in axes))


@pytest.mark.parametrize(
    "step,length,valid",
    [(5, 7, True), (5, 9, True), (2.5, 13, True), (5, 6, False), (10, 4, False)],
)
def test_delayed_bath_construction(step, length, valid):
    molecule = qr.Molecule([0, 1])
    molecule.set_dipole(0, 1, [1, 0, 0])
    system = qr.Aggregate(molecules=[molecule])
    calc = qr.TwoDResponseCalculator(
        qr.TimeAxis(0, 3, 5), qr.TimeAxis(0, 2, 10), qr.TimeAxis(0, 3, 5), system=system
    )
    with qr.energy_units("int"):
        bath = qr.CorrelationFunction(
            qr.TimeAxis(0, length, step),
            dict(ftype="OverdampedBrownian", reorg=0.01, cortime=60, T=300),
        )
    molecule.set_transition_environment((0, 1), bath)
    system.build(mult=2)
    lab = qr.LabSetup()
    x = [1.0, 0.0, 0.0]
    lab.set_pulse_polarizations(pulse_polarizations=(x, x, x), detection_polarization=x)
    if valid:
        calc.bootstrap(rwa=1, lab=lab)
        spectrum = calc.calculate_one(0)
        assert spectrum.data.shape == (3, 3)
        assert np.isfinite(spectrum.data).all()
        # A post-bootstrap change must be detected before response evaluation.
        system.get_SystemBathInteraction().TimeAxis = qr.TimeAxis(0, 2, 5)
        with pytest.raises(ValueError, match="cover"):
            calc.calculate_one(0)
    else:
        with pytest.raises(ValueError, match="Bath"):
            calc.bootstrap(rwa=1, lab=lab)
