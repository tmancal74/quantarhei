import unittest

import numpy

"""
*******************************************************************************


    Tests of the quantarhei.qm.propagators.poppropagator module


*******************************************************************************
"""


from quantarhei import PopulationPropagator, TimeAxis


class TestPopulationPropagator(unittest.TestCase):
    """Tests population propagator module"""

    def test_of_population_evolution_1(self):
        """Testing population evolution matrix 2x2 starting from t = 0"""
        KK = numpy.array([[-1.0 / 100.0, 1.0 / 100.0], [1.0 / 100.0, -1.0 / 100.0]])

        U0 = numpy.eye(2)
        Ntd = 10

        t = TimeAxis(0.0, 1000, 1.0)
        prop = PopulationPropagator(t, rate_matrix=KK)

        td = TimeAxis(0.0, Ntd, 10.0)
        U = prop.get_PropagationMatrix(td)

        # analytical result

        Ucheck = numpy.zeros((2, 2, Ntd))
        Ucheck[0, 0, :] = 0.5 * (1.0 + numpy.exp(2.0 * KK[0, 0] * td.data))
        Ucheck[1, 1, :] = Ucheck[0, 0, :]
        Ucheck[1, 0, :] = 0.5 * (1.0 - numpy.exp(2.0 * KK[0, 0] * td.data))
        Ucheck[0, 1, :] = Ucheck[1, 0, :]

        for n in range(Ntd):
            numpy.testing.assert_allclose(U[:, :, n], Ucheck[:, :, n])

    def test_of_population_evolution_2(self):
        """Testing population evolution matrix 2x2 starting from t > 0"""
        KK = numpy.array([[-1.0 / 100.0, 1.0 / 100.0], [1.0 / 100.0, -1.0 / 100.0]])

        U0 = numpy.eye(2)
        Ntd = 10

        t = TimeAxis(0.0, 1000, 1.0)
        prop = PopulationPropagator(t, rate_matrix=KK)

        td = TimeAxis(2.0, Ntd, 10.0)
        U = prop.get_PropagationMatrix(td)

        # analytical result

        Ucheck = numpy.zeros((2, 2, Ntd))
        Ucheck[0, 0, :] = 0.5 * (1.0 + numpy.exp(2.0 * KK[0, 0] * td.data))
        Ucheck[1, 1, :] = Ucheck[0, 0, :]
        Ucheck[1, 0, :] = 0.5 * (1.0 - numpy.exp(2.0 * KK[0, 0] * td.data))
        Ucheck[0, 1, :] = Ucheck[1, 0, :]

        for n in range(Ntd):
            numpy.testing.assert_allclose(U[:, :, n], Ucheck[:, :, n])

    def test_jump_expansion_zero_order_matches_old_correction(self):
        """Testing that zero jump expansion is compatible with old correction"""
        KK = numpy.array([[-1.0 / 100.0, 1.0 / 100.0], [1.0 / 100.0, -1.0 / 100.0]])

        t = TimeAxis(0.0, 1000, 1.0)
        td = TimeAxis(2.0, 10, 10.0)
        prop = PopulationPropagator(t, rate_matrix=KK)

        _, cor = prop.get_PropagationMatrix(td, corrections=0)
        jumps = prop.get_JumpExpansion(td, max_order=0)

        self.assertEqual(len(jumps), 1)
        numpy.testing.assert_allclose(jumps[0], cor[0])

    def test_jump_expansion_reconstructs_constant_rate_propagator(self):
        """Testing that jump expansion reconstructs a constant-rate propagator"""
        KK = numpy.array([[-2.0 / 1000.0, 1.0 / 1000.0], [2.0 / 1000.0, -1.0 / 1000.0]])

        t = TimeAxis(0.0, 1001, 1.0)
        td = TimeAxis(0.0, 11, 20.0)
        prop = PopulationPropagator(t, rate_matrix=KK)

        U = prop.get_PropagationMatrix(td)
        jumps = prop.get_JumpExpansion(td, max_order=10)
        Uj = numpy.sum(numpy.array(jumps), axis=0)

        numpy.testing.assert_allclose(Uj, U, atol=1.0e-7)

    def test_jump_expansion_accepts_time_dependent_rates(self):
        """Testing jump expansion with a time-dependent rate matrix"""
        t = TimeAxis(0.0, 101, 1.0)
        rates = numpy.zeros((t.length, 2, 2), dtype=numpy.float64)

        k01 = 0.001 * (1.0 + t.data / t.data[-1])
        k10 = 0.0005 * numpy.ones(t.length, dtype=numpy.float64)
        rates[:, 0, 1] = k01
        rates[:, 1, 0] = k10
        rates[:, 0, 0] = -k10
        rates[:, 1, 1] = -k01

        td = TimeAxis(0.0, 6, 20.0)
        prop = PopulationPropagator(t, rate_matrix=rates)
        jumps = prop.get_JumpExpansion(td, max_order=2)

        self.assertEqual(len(jumps), 3)
        for jump in jumps:
            self.assertEqual(jump.shape, (2, 2, td.length))

        numpy.testing.assert_allclose(jumps[0][:, :, 0], numpy.eye(2))

    def test_time_dependent_propagation_matrix(self):
        """Testing propagation matrix with time-dependent rate matrix"""
        t = TimeAxis(0.0, 101, 1.0)
        rates = numpy.zeros((t.length, 2, 2), dtype=numpy.float64)
        rates[:, 0, 0] = -0.001 * (1.0 + t.data / t.data[-1])
        rates[:, 1, 1] = -0.002

        td = TimeAxis(0.0, 6, 20.0)
        prop = PopulationPropagator(t, rate_matrix=rates)
        U, jumps = prop.get_PropagationMatrix(td, corrections=0)

        numpy.testing.assert_allclose(U, jumps[0], atol=1.0e-8)
