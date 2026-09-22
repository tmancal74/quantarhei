import unittest

"""
*******************************************************************************


    Tests of the quantarhei.qm.SystemBathInteraction class


*******************************************************************************
"""

from quantarhei.qm.liouvillespace.systembathinteraction_test import (
    TestSystemBathInteraction,
)


class SBITest(unittest.TestCase):
    """Tests for the SystemBathInteraction class"""

    def setUp(self, verbose=False):

        self.verbose = verbose

    def test_system_bath_interaction_creation(self):
        """(SystemBathInteraction) Testing creation"""
        sbi = TestSystemBathInteraction(name="dimer-2-env")

        self.assertEqual(sbi.sbitype, "Linear_Coupling")

    def test_get_temperature_and_has_temperature(self):
        """(SystemBathInteraction) Testing has_temperature() and get_temperature()"""
        sbi = TestSystemBathInteraction(name="dimer-2-env")

        self.assertTrue(sbi.has_temperature())

        T = sbi.get_temperature()
        self.assertAlmostEqual(300.0, T)

        sbi = TestSystemBathInteraction(name="dimer-2-lind")
        self.assertFalse(sbi.has_temperature())


def test_lineshape_storage_extrapolates_linearly_for_response_time_sums():
    """Response storage interpolates inside the bath grid and has a linear tail."""
    import numpy as np
    import numpy.testing as npt

    import quantarhei as qr
    from quantarhei.core.dfunction import DFunction
    from quantarhei.qm.corfunctions.correlationfunctions import c2g

    axis = qr.TimeAxis(0.0, 8, 5.0)
    with qr.energy_units("1/cm"):
        cf = qr.CorrelationFunction(
            axis,
            dict(ftype="OverdampedBrownian", reorg=30.0, cortime=60.0, T=300.0),
        )
    molecule = qr.Molecule([0.0, 1.0])
    molecule.set_transition_environment((0, 1), cf)
    system = qr.Aggregate(molecules=[molecule])
    system.build()
    sbi = system.get_SystemBathInteraction()
    storage = sbi.get_goft_storage()
    expected = DFunction(axis, c2g(axis, cf.data)).as_spline_function()
    function = storage.funcs[0]
    npt.assert_allclose(function(axis.data), expected(axis.data))
    npt.assert_allclose(function(12.5), expected(12.5))
    end = axis.data[-1]
    slope = (expected(end) - expected(end - axis.step)) / axis.step
    npt.assert_allclose(function(1000.0), expected(end) + (1000.0 - end) * slope)
    npt.assert_allclose(function(end + 1.0e-8), expected(end), atol=1.0e-8)
    arguments = np.array([[0.0, 12.5], [35.0, 1000.0]])
    npt.assert_allclose(
        function(arguments),
        expected(np.minimum(arguments, end)) + np.maximum(arguments - end, 0) * slope,
    )
    storage.create_data(reset={"t2": 100.0})
    sums = axis.data[:, None] + 100.0 + axis.data[None, :]
    npt.assert_allclose(
        storage[:, "t1+t2+t3"], (expected(end) + (sums - end) * slope)[None, :, :]
    )
