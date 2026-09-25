from types import SimpleNamespace

import numpy
import numpy.testing as npt

import quantarhei as qr
from quantarhei.spectroscopy.lineshapes import voigt2D


def test_voigt2d_uses_first_axis_for_omega3_and_second_for_omega1():
    """An asymmetric grid exposes the storage order used by mock 2D spectra."""
    omega1 = numpy.array([1.0, 2.0])
    omega3 = numpy.array([10.0, 11.0, 12.0])

    data = voigt2D(omega1, 1.5, 1.0, 0.0, omega3, 11.0, 1.0, 0.0)

    assert data.shape == (len(omega3), len(omega1))
    npt.assert_allclose(data, numpy.outer(data[:, 0] / data[0, 0], data[0, :]))


def test_mock_calculator_keeps_omega3_omega1_order_on_asymmetric_axes():
    """Mock spectra use the same axis order as Fourier-transformed responses."""
    calculator = qr.MockTwoDResponseCalculator(
        qr.TimeAxis(0.0, 2, 1.0),
        qr.TimeAxis(0.0, 1, 1.0),
        qr.TimeAxis(0.0, 3, 1.0),
    )
    calculator.bootstrap(rwa=0.0)
    pathway = SimpleNamespace(
        order=0,
        relax_order=0,
        frequency=(0.0, 0.0),
        pref=1.0,
        widths=(-1.0, -1.0, -1.0, -1.0),
        dephs=(-1.0, -1.0, -1.0, -1.0),
        pathway_type="NR",
    )

    data = calculator.calculate_pathway(pathway)

    assert data.shape == (calculator.oa3.length, calculator.oa1.length)
