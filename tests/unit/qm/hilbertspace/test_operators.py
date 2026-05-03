import unittest

import numpy

"""
*******************************************************************************


    Tests of the quantarhei.Molecule class


*******************************************************************************
"""

from quantarhei import REAL, Molecule
from quantarhei.qm import Operator, ReducedDensityMatrix


class TestReducedDensityMatrix(unittest.TestCase):
    """Tests for the Manager class"""

    def setUp(self):
        self.en = [0.0, 1.0, 2.0]
        self.m = Molecule(name="Molecule", elenergies=self.en)
        self.m.set_dipole(0, 1, [1.0, 0.0, 0.0])
        self.m.set_dipole(0, 2, [0.5, 0.0, 0.0])

        self.rho_eq = self.m.get_thermal_ReducedDensityMatrix()

    def test_excitation_by_delta(self):
        """Testing reduced density matrix excitation by delta-pulse"""
        rho_exp = numpy.zeros((3, 3), dtype=REAL)
        rho_exp[0, 0] = 1.0

        self.assertTrue(numpy.allclose(self.rho_eq._data, rho_exp))

        dd = self.m.get_TransitionDipoleMoment()
        epol = [1.0, 0.0, 0.0]

        rho_ex = self.rho_eq.excite_delta(dmoment=dd, epolarization=epol)

        rho_exp = numpy.array(
            [
                [0.00 + 0.0j, 0.00 + 0.0j, 0.00 + 0.0j],
                [0.00 + 0.0j, 1.00 + 0.0j, 0.50 + 0.0j],
                [0.00 + 0.0j, 0.50 + 0.0j, 0.25 + 0.0j],
            ]
        )

        self.assertTrue(numpy.allclose(rho_exp, rho_ex._data))

        rho_ex2 = self.rho_eq.excite_delta(
            dmoment=dd,
            epolarization=[1.0 / numpy.sqrt(2.0), 1.0 / numpy.sqrt(2.0), 0.0],
        )

        rho_exp2 = numpy.array(
            [
                [0.000 + 0.0j, 0.000 + 0.0j, 0.000 + 0.0j],
                [0.000 + 0.0j, 0.500 + 0.0j, 0.250 + 0.0j],
                [0.000 + 0.0j, 0.250 + 0.0j, 0.125 + 0.0j],
            ]
        )

        self.assertTrue(numpy.allclose(rho_ex2._data, rho_exp2))

    def test_conversion_2_populations(self):
        """Testing conversion of reduced density matrix to population vector"""
        rdm = ReducedDensityMatrix(
            data=[[0.5, 0.0, 0.1], [0.0, 0.3, 0.0], [0.1, 0.0, 0.2]]
        )

        pop = rdm.get_populations()

        self.assertTrue(numpy.allclose(pop, [0.5, 0.3, 0.2]))

    def test_saveable(self):
        """Testing reduced density matrix as Saveable"""
        rdm = ReducedDensityMatrix(
            data=[[0.5, 0.0, 0.1], [0.0, 0.3, 0.0], [0.1, 0.0, 0.2]]
        )

        # import h5py

        # with h5py.File("test_file_operators",driver="core",
        #                   backing_store=False) as fid:
        import tempfile

        with tempfile.TemporaryFile() as fid:
            rdm.save(fid)  # , test=True)
            fid.seek(0)

            rdm2 = ReducedDensityMatrix()

            rdm2 = rdm2.load(fid)  # , test=True)

        self.assertTrue(numpy.allclose(rdm2.data, rdm.data))


class TestOperatorAdd(unittest.TestCase):
    """Tests that Operator.__add__ does not mutate the left operand."""

    def test_add_returns_new_object(self):
        """op1 + op2 must not modify op1."""
        op1 = Operator(data=numpy.array([[1.0, 0.0], [0.0, 2.0]]))
        op2 = Operator(data=numpy.array([[3.0, 0.0], [0.0, 4.0]]))
        original = op1.data.copy()

        result = op1 + op2

        # op1 must be unchanged
        self.assertTrue(numpy.allclose(op1.data, original))
        # result must equal the sum
        self.assertTrue(numpy.allclose(result.data, original + op2.data))
        # result must be a distinct object
        self.assertIsNot(result, op1)

    def test_add_chained(self):
        """op1 + op2 + op3 must leave op1 and op2 unchanged."""
        op1 = Operator(data=numpy.eye(2))
        op2 = Operator(data=2.0 * numpy.eye(2))
        op3 = Operator(data=3.0 * numpy.eye(2))
        d1 = op1.data.copy()
        d2 = op2.data.copy()

        result = op1 + op2 + op3

        self.assertTrue(numpy.allclose(op1.data, d1))
        self.assertTrue(numpy.allclose(op2.data, d2))
        self.assertTrue(numpy.allclose(result.data, 6.0 * numpy.eye(2)))
