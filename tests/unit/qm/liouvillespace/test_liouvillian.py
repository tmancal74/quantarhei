import unittest

import numpy

from quantarhei import Hamiltonian
from quantarhei.qm.hilbertspace.operators import DensityMatrix
from quantarhei.qm.liouvillespace.liouvillian import Liouvillian


class TestLiouvillian(unittest.TestCase):
    """Tests for the Liouvillian superoperator"""

    def setUp(self):
        self.H2 = Hamiltonian(data=[[0.0, 1.0], [1.0, 0.2]])
        self.H3 = Hamiltonian(data=[[0.0, 1.0, 0.5], [1.0, 0.2, 0.3], [0.5, 0.3, 0.4]])

    def test_liouvillian_matches_explicit_commutator_2x2(self):
        """Liouvillian applied to rho equals [H, rho] for a 2x2 system"""
        H = self.H2
        rho = DensityMatrix(data=[[0.9, 0.05], [0.05, 0.1]])

        Li = Liouvillian(H)
        result = Li.apply(rho)

        expected = numpy.dot(H.data, rho.data) - numpy.dot(rho.data, H.data)
        numpy.testing.assert_allclose(result.data, expected, atol=1e-12)

    def test_liouvillian_matches_explicit_commutator_3x3(self):
        """Liouvillian applied to rho equals [H, rho] for a 3x3 system"""
        H = self.H3
        rho = DensityMatrix(
            data=[[0.5, 0.1, 0.05], [0.1, 0.3, 0.02], [0.05, 0.02, 0.2]]
        )

        Li = Liouvillian(H)
        result = Li.apply(rho)

        expected = numpy.dot(H.data, rho.data) - numpy.dot(rho.data, H.data)
        numpy.testing.assert_allclose(result.data, expected, atol=1e-12)

    def test_liouvillian_tensor_shape(self):
        """Liouvillian data must be (dim, dim, dim, dim)"""
        for dim in (2, 3, 5):
            data = numpy.diag(numpy.arange(dim, dtype=float))
            H = Hamiltonian(data=data)
            Li = Liouvillian(H)
            self.assertEqual(Li.data.shape, (dim, dim, dim, dim))

    def test_liouvillian_diagonal_hamiltonian_gives_zero_populations(self):
        """[H, rho] leaves diagonal of rho unchanged when H is diagonal"""
        H = Hamiltonian(data=numpy.diag([0.0, 1.0, 2.0]))
        rho_data = numpy.diag([0.5, 0.3, 0.2]).astype(complex)
        rho = DensityMatrix(data=rho_data)

        Li = Liouvillian(H)
        result = Li.apply(rho)

        numpy.testing.assert_allclose(numpy.diag(result.data), 0.0, atol=1e-12)
