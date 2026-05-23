import unittest

import numpy

from quantarhei import CorrelationFunction, Hamiltonian, TimeAxis, energy_units
from quantarhei.qm import (
    Operator,
    SystemBathInteraction,
    TDFoersterRateMatrix,
    TDFoersterRelaxationTensor,
)
from quantarhei.qm.corfunctions import CorrelationFunctionMatrix


class TestTDFoersterRateMatrix(unittest.TestCase):
    """Tests for the time-dependent Foerster rate matrix."""

    def setUp(self):
        time = TimeAxis(0.0, 100, 1.0)
        with energy_units("1/cm"):
            params = {
                "ftype": "OverdampedBrownian",
                "reorg": 30.0,
                "T": 300.0,
                "cortime": 100.0,
            }
            cf = CorrelationFunction(time, params)

            self.ham = Hamiltonian(data=[[0.0, 100.0], [100.0, 0.0]])

        cfm = CorrelationFunctionMatrix(time, 2, 1)
        cfm.set_correlation_function(cf, [(0, 0), (1, 1)])

        k1 = Operator(data=numpy.array([[1.0, 0.0], [0.0, 0.0]]))
        k2 = Operator(data=numpy.array([[1.0, 0.0], [0.0, 0.0]]))
        self.sbi = SystemBathInteraction([k1, k2], cfm)

    def test_rate_matrix_matches_relaxation_tensor_population_block(self):
        """Testing time-dependent Foerster rates against tensor population block."""
        rate_matrix = TDFoersterRateMatrix(self.ham, self.sbi)
        tensor = TDFoersterRelaxationTensor(self.ham, self.sbi)

        rates_from_tensor = numpy.zeros_like(rate_matrix.data)
        for aa in range(self.ham.dim):
            for bb in range(self.ham.dim):
                rates_from_tensor[:, aa, bb] = numpy.real(
                    tensor.data[:, aa, aa, bb, bb]
                )

        numpy.testing.assert_allclose(rate_matrix.data, rates_from_tensor)


if __name__ == "__main__":
    unittest.main()
