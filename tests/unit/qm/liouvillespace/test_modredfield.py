import unittest

import numpy

"""
*******************************************************************************


    Tests of the quantarhei.qm.RedfieldRelaxationTensor and
                 quantarhei.qm.RefieldRateMatrix classes


*******************************************************************************
"""


from quantarhei import (
    Aggregate,
    CorrelationFunction,
    Molecule,
    TimeAxis,
    eigenbasis_of,
    energy_units,
)
from quantarhei.qm import (
    ModifiedRedfieldRateMatrix,
    ModRedfieldRelaxationTensor,
    TDModRedfieldRelaxationTensor,
)
from quantarhei.qm.corfunctions import SpectralDensity


class TestModRedfield(unittest.TestCase):
    """Tests for the RedfieldRelaxationTensor class"""

    def setUp(self, verbose=False):

        self.verbose = verbose

        time = TimeAxis(0.0, 1000, 1.0)
        self.time = time
        with energy_units("1/cm"):
            params = {
                "ftype": "OverdampedBrownian",
                "reorg": 30.0,
                "T": 300.0,
                "cortime": 100.0,
            }

            cf1 = CorrelationFunction(time, params)
            cf2 = CorrelationFunction(time, params)
            sd = SpectralDensity(time, params)

            m1 = Molecule([0.0, 10000.0])
            m1.set_transition_environment((0, 1), cf1)
            m2 = Molecule([0.0, 10000.0])
            m2.set_transition_environment((0, 1), cf2)

            agg = Aggregate(molecules=[m1, m2])

            agg.set_resonance_coupling(0, 1, 100.0)

        agg.build()

        self.H1 = agg.get_Hamiltonian()
        self.sbi1 = agg.get_SystemBathInteraction()

        # sd.convert_2_spectral_density()
        with eigenbasis_of(self.H1):
            de = self.H1.data[1, 1] - self.H1.data[0, 0]
        self.c_omega_p = sd.at(de, approx="spline")  # interp_data(de)
        self.c_omega_m = sd.at(-de, approx="spline")  # interp_data(-de)

    def test_comparison_of_rates(self):
        """Testing that modified Redfield tensor and rate matrix agree"""
        dim = self.H1.dim
        KT = numpy.zeros((dim, dim), dtype=numpy.float64)

        RT = ModRedfieldRelaxationTensor(self.H1, self.sbi1)
        RR = ModifiedRedfieldRateMatrix(self.H1, self.sbi1)

        with eigenbasis_of(self.H1):
            for n in range(dim):
                for m in range(dim):
                    KT[n, m] = numpy.real(RT.data[n, n, m, m])

        numpy.testing.assert_allclose(KT, RR.data, rtol=1.0e-10, atol=1.0e-12)

    def test_calculation_is_forbidden_in_basis_context(self):
        """Testing that modified Redfield tensor is not created in basis context"""
        with self.assertRaises(Exception) as context:
            with eigenbasis_of(self.H1):
                ModRedfieldRelaxationTensor(self.H1, self.sbi1)

        self.assertTrue("MUST NOT be called" in str(context.exception))

    def test_td_calculation_is_forbidden_in_basis_context(self):
        """Testing that TD modified Redfield tensor is not created in basis context"""
        with self.assertRaises(Exception) as context:
            with eigenbasis_of(self.H1):
                TDModRedfieldRelaxationTensor(self.H1, self.sbi1, initialize=False)

        self.assertTrue("MUST NOT be called" in str(context.exception))

    # def test_propagation_in_different_basis(self):
    #     """(REDFIELD) Testing comparison of propagations in different bases

    #     """

    #     #LT1 = ModRedfieldRelaxationTensor(self.H1, self.sbi1, as_operators=True)
    #     LT2 = ModRedfieldRelaxationTensor(self.H1, self.sbi1, as_operators=False)

    #     time = TimeAxis(0.0, 1000, 1.0)

    #     #prop1 = ReducedDensityMatrixPropagator(time, self.H1, LT1)
    #     prop2 = ReducedDensityMatrixPropagator(time, self.H1, LT2)

    #     rho0 = ReducedDensityMatrix(dim=self.H1.dim)
    #     rho0.data[1,1] = 1.0

    #     #with eigenbasis_of(self.H1):
    #     ##if True:
    #     #    rhot1_e = prop1.propagate(rho0)

    #     with eigenbasis_of(self.H1):
    #         rhot2_e = prop2.propagate(rho0)

    #     #rhot1_l = prop1.propagate(rho0)
    #     rhot2_l = prop2.propagate(rho0)

    #     #numpy.testing.assert_allclose(rhot1_l.data, rhot1_e.data,
    #     #                              rtol=1.0e-5, atol=1.0e-12)
    #     #numpy.testing.assert_allclose(rhot2_l.data, rhot1_e.data,
    #     #                              rtol=1.0e-5, atol=1.0e-12)
    #     #numpy.testing.assert_allclose(rhot1_e.data, rhot2_e.data,
    #     #                              rtol=1.0e-5, atol=1.0e-12)
    #     numpy.testing.assert_allclose(rhot2_l.data, rhot2_e.data,
    #                                   rtol=1.0e-5, atol=1.0e-12)
