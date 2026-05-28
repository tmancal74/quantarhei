import unittest
import warnings

import numpy

import quantarhei as qr


def _aggregate():
    m1 = qr.Molecule(name="Molecule 1", elenergies=[0.0, 1.0])
    m2 = qr.Molecule(name="Molecule 2", elenergies=[0.0, 1.0])

    time = qr.TimeAxis(0.0, 100, 1.0)
    params = dict(ftype="OverdampedBrownian", reorg=20, cortime=100, T=300)
    with qr.energy_units("1/cm"):
        cf = qr.CorrelationFunction(time, params)

    m1.set_transition_environment((0, 1), cf)
    m1.position = [0.0, 0.0, 0.0]
    m1.set_dipole(0, 1, [10.0, 0.0, 0.0])

    m2.set_transition_environment((0, 1), cf)
    m2.position = [10.0, 0.0, 0.0]
    m2.set_dipole(0, 1, [10.0, 0.0, 0.0])

    agg = qr.Aggregate(name="TestAgg", molecules=[m1, m2])
    agg.set_coupling_by_dipole_dipole()
    agg.build()

    return agg


class TestOpenSystemRateMatrix(unittest.TestCase):
    """Tests for unified OpenSystem rate matrix construction."""

    def setUp(self):
        self.agg = _aggregate()

    def test_redfield_rate_matrix_aliases(self):
        """Testing unified Redfield rate matrix construction."""
        rate_matrix = self.agg.get_RateMatrix(relaxation_theory="Redfield")
        with warnings.catch_warnings(record=True) as warning_list:
            warnings.simplefilter("always")
            old_rate_matrix = self.agg.get_RedfieldRateMatrix()

        self.assertIsInstance(rate_matrix, qr.qm.RedfieldRateMatrix)
        self.assertTrue(issubclass(warning_list[0].category, DeprecationWarning))
        self.assertIn(
            "use get_RateMatrix(relaxation_theory='Redfield') instead",
            str(warning_list[0].message),
        )
        numpy.testing.assert_allclose(rate_matrix.data, old_rate_matrix.data)

    def test_foerster_rate_matrix_aliases(self):
        """Testing unified Foerster rate matrix construction."""
        rate_matrix = self.agg.get_RateMatrix(relaxation_theory="Foerster")
        with warnings.catch_warnings(record=True) as warning_list:
            warnings.simplefilter("always")
            old_rate_matrix = self.agg.get_FoersterRateMatrix()

        self.assertIsInstance(rate_matrix, qr.qm.FoersterRateMatrix)
        self.assertTrue(issubclass(warning_list[0].category, DeprecationWarning))
        self.assertIn(
            "use get_RateMatrix(relaxation_theory='Foerster') instead",
            str(warning_list[0].message),
        )
        numpy.testing.assert_allclose(rate_matrix.data, old_rate_matrix.data)

    def test_time_dependent_rate_matrices(self):
        """Testing unified time-dependent rate matrix construction."""
        redfield = self.agg.get_RateMatrix(
            relaxation_theory="Redfield", time_dependent=True
        )
        foerster = self.agg.get_RateMatrix(
            relaxation_theory="Foerster", time_dependent=True
        )

        self.assertIsInstance(redfield, qr.qm.TDRedfieldRateMatrix)
        self.assertIsInstance(foerster, qr.qm.TDFoersterRateMatrix)
        self.assertEqual(redfield.data.shape[0], self.agg.sbi.TimeAxis.length)
        self.assertEqual(foerster.data.shape[0], self.agg.sbi.TimeAxis.length)


if __name__ == "__main__":
    unittest.main()
