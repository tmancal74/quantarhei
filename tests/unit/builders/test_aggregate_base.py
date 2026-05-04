import unittest

import numpy

from quantarhei import (
    Aggregate,
    Mode,
    Molecule,
    energy_units,
)
from quantarhei.builders.aggregate_test import TestAggregate
from quantarhei.qm import ReducedDensityMatrix


def _dimer(with_modes: bool = False, hr: tuple = (0.1, 0.3)) -> Aggregate:
    m1 = Molecule(elenergies=[0.0, 1.0])
    m2 = Molecule(elenergies=[0.0, 1.1])
    m1.position = [0.0, 0.0, 0.0]
    m2.position = [10.0, 0.0, 0.0]
    m1.set_dipole(0, 1, [10.0, 0.0, 0.0])
    m2.set_dipole(0, 1, [0.0, 10.0, 0.0])
    if with_modes:
        for m, h in zip([m1, m2], hr):
            mod = Mode(0.01)
            m.add_Mode(mod)
            mod.set_nmax(0, 4)
            mod.set_nmax(1, 4)
            mod.set_HR(1, h)
    agg = Aggregate(molecules=[m1, m2])
    with energy_units("1/cm"):
        agg.set_resonance_coupling(0, 1, 100.0)
    agg.build()
    return agg


class TestElsignatures(unittest.TestCase):
    def setUp(self) -> None:
        self.agg = _dimer()

    def test_ground_state_only_when_mult_zero(self) -> None:
        sigs = list(self.agg.elsignatures(mult=0))
        self.assertEqual(len(sigs), 1)
        self.assertEqual(sigs[0], (0, 0))

    def test_single_exciton_band_has_nmono_states(self) -> None:
        sigs = list(self.agg.elsignatures(mult=1, mode="EQ"))
        self.assertEqual(len(sigs), self.agg.nmono)

    def test_lq_includes_lower_bands(self) -> None:
        sigs_lq = list(self.agg.elsignatures(mult=1, mode="LQ"))
        sigs_eq = list(self.agg.elsignatures(mult=1, mode="EQ"))
        self.assertGreater(len(sigs_lq), len(sigs_eq))

    def test_invalid_mode_raises(self) -> None:
        with self.assertRaises(Exception):
            list(self.agg.elsignatures(mult=1, mode="INVALID"))

    def test_negative_mult_raises(self) -> None:
        with self.assertRaises(Exception):
            list(self.agg.elsignatures(mult=-1))


class TestTotalNumberOfStates(unittest.TestCase):
    def setUp(self) -> None:
        self.agg = _dimer()
        self.vagg = _dimer(with_modes=True)

    def test_electronic_states_ground_plus_band(self) -> None:
        n = self.agg.total_number_of_electronic_states(mult=1)
        self.assertEqual(n, 3)

    def test_number_of_states_in_band_one(self) -> None:
        n = self.agg.number_of_states_in_band(band=1)
        self.assertEqual(n, self.agg.nmono)

    def test_total_states_with_modes_larger_than_electronic(self) -> None:
        n_el = self.vagg.total_number_of_electronic_states(mult=1)
        n_tot = self.vagg.total_number_of_states(mult=1)
        self.assertGreater(n_tot, n_el)

    def test_number_of_electronic_states_in_band(self) -> None:
        n = self.agg.number_of_electronic_states_in_band(band=1)
        self.assertEqual(n, self.agg.nmono)


class TestFcFactor(unittest.TestCase):
    def setUp(self) -> None:
        self.vagg = _dimer(with_modes=True)

    def test_fc_factor_same_state_is_nonzero(self) -> None:
        states = list(self.vagg.allstates(mult=1))
        _, state = states[self.vagg.Nb[0]]
        fc = self.vagg.fc_factor(state, state)
        self.assertGreater(abs(fc), 0.0)

    def test_fc_factor_returns_real(self) -> None:
        states = list(self.vagg.allstates(mult=1))
        _, state = states[self.vagg.Nb[0]]
        fc = self.vagg.fc_factor(state, state)
        self.assertIsInstance(fc, float)


class TestGetNearestMolecule(unittest.TestCase):
    def setUp(self) -> None:
        self.agg = _dimer()

    def test_nearest_returns_other_molecule(self) -> None:
        m0 = self.agg.monomers[0]
        nearest, dist = self.agg.get_nearest_Molecule(m0)
        self.assertIs(nearest, self.agg.monomers[1])

    def test_distance_is_positive(self) -> None:
        m0 = self.agg.monomers[0]
        _, dist = self.agg.get_nearest_Molecule(m0)
        self.assertGreater(dist, 0.0)


class TestGetMoleculeByName(unittest.TestCase):
    def setUp(self) -> None:
        self.agg = TestAggregate(name="dimer-2-env")

    def test_get_by_name_returns_correct_molecule(self) -> None:
        name = self.agg.monomers[0].name
        mol = self.agg.get_Molecule_by_name(name)
        self.assertEqual(mol.name, name)

    def test_get_index_by_name_is_consistent(self) -> None:
        name = self.agg.monomers[0].name
        idx = self.agg.get_Molecule_index(name)
        self.assertEqual(self.agg.monomers[idx].name, name)


class TestRemoveMolecule(unittest.TestCase):
    def test_remove_reduces_count(self) -> None:
        agg = TestAggregate(name="dimer-2-env")
        n_before = agg.nmono
        mol = agg.monomers[0]
        agg.remove_Molecule(mol)
        self.assertEqual(agg.nmono, n_before - 1)


class TestGetHamiltonian(unittest.TestCase):
    def setUp(self) -> None:
        self.agg = _dimer()

    def test_hamiltonian_dimension_matches_ntot(self) -> None:
        H = self.agg.get_Hamiltonian()
        self.assertEqual(H.dim, self.agg.Ntot)

    def test_hamiltonian_is_hermitian(self) -> None:
        H = self.agg.get_Hamiltonian()
        numpy.testing.assert_array_almost_equal(H.data, H.data.conj().T)


class TestDiagonalize(unittest.TestCase):
    def setUp(self) -> None:
        self.agg = _dimer()

    def test_diagonalize_sets_hd(self) -> None:
        self.agg.diagonalize()
        self.assertIsNotNone(self.agg.HD)

    def test_diagonalize_eigenvalues_are_real(self) -> None:
        self.agg.diagonalize()
        numpy.testing.assert_array_almost_equal(self.agg.HD, numpy.real(self.agg.HD))

    def test_diagonalize_idempotent(self) -> None:
        self.agg.diagonalize()
        hd_first = self.agg.HD.copy()
        self.agg.diagonalize()
        numpy.testing.assert_array_equal(self.agg.HD, hd_first)

    def test_ss_is_unitary(self) -> None:
        self.agg.diagonalize()
        SS = self.agg.SS
        numpy.testing.assert_array_almost_equal(
            numpy.dot(SS, SS.conj().T), numpy.eye(SS.shape[0])
        )


class TestConvertToGroundVibBasis(unittest.TestCase):
    def setUp(self) -> None:
        self.vagg = _dimer(with_modes=True)

    def test_convert_density_matrix_preserves_trace(self) -> None:
        H = self.vagg.get_Hamiltonian()
        rho = ReducedDensityMatrix(dim=H.dim)
        rho._data[self.vagg.Nb[0], self.vagg.Nb[0]] = 1.0
        converted = self.vagg.convert_to_ground_vibbasis(rho)
        numpy.testing.assert_almost_equal(numpy.trace(converted._data), 1.0, decimal=5)

    def test_convert_returns_same_dimension(self) -> None:
        H = self.vagg.get_Hamiltonian()
        rho = ReducedDensityMatrix(dim=H.dim)
        rho._data[self.vagg.Nb[0], self.vagg.Nb[0]] = 1.0
        converted = self.vagg.convert_to_ground_vibbasis(rho)
        self.assertEqual(converted.dim, self.vagg.Ntot)

    def test_convert_zero_operator_stays_zero(self) -> None:
        H = self.vagg.get_Hamiltonian()
        rho = ReducedDensityMatrix(dim=H.dim)
        converted = self.vagg.convert_to_ground_vibbasis(rho)
        numpy.testing.assert_array_almost_equal(
            converted._data, numpy.zeros_like(converted._data)
        )


class TestTraceOverVibrationsHamiltonian(unittest.TestCase):
    def setUp(self) -> None:
        self.vagg = _dimer(with_modes=True)

    def test_trace_hamiltonian_has_electronic_dimension(self) -> None:
        H = self.vagg.get_Hamiltonian()
        Hel = self.vagg.trace_over_vibrations(H)
        self.assertEqual(Hel.dim, self.vagg.Nel)

    def test_trace_hamiltonian_is_hermitian(self) -> None:
        H = self.vagg.get_Hamiltonian()
        Hel = self.vagg.trace_over_vibrations(H)
        numpy.testing.assert_array_almost_equal(Hel.data, Hel.data.conj().T)


class TestGetRWASuggestion(unittest.TestCase):
    def setUp(self) -> None:
        self.agg = _dimer()

    def test_rwa_suggestion_is_positive(self) -> None:
        rwa = self.agg.get_RWA_suggestion()
        self.assertGreater(rwa, 0.0)

    def test_rwa_suggestion_within_energy_range(self) -> None:
        H = self.agg.get_Hamiltonian()
        rwa = self.agg.get_RWA_suggestion()
        energies = numpy.diag(H.data)
        self.assertGreaterEqual(rwa, numpy.min(energies))
        self.assertLessEqual(rwa, numpy.max(energies))


class TestGetStateVector(unittest.TestCase):
    def setUp(self) -> None:
        self.agg = _dimer()

    def test_impulsive_state_vector_is_nonzero(self) -> None:
        sv = self.agg.get_StateVector(condition_type="impulsive_excitation")
        self.assertGreater(numpy.sum(numpy.abs(sv.data) ** 2), 0.0)

    def test_state_vector_has_correct_dimension(self) -> None:
        sv = self.agg.get_StateVector(condition_type="impulsive_excitation")
        self.assertEqual(sv.dim, self.agg.Ntot)


class TestCoupling(unittest.TestCase):
    def setUp(self) -> None:
        self.agg = _dimer()
        states = list(self.agg.allstates(mult=1))
        Nb0 = self.agg.Nb[0]
        _, self.gs = states[0]
        _, self.es1 = states[Nb0]
        _, self.es2 = states[Nb0 + 1]

    def test_coupling_between_ground_and_excited_is_zero(self) -> None:
        c = self.agg.coupling(self.gs, self.es1)
        self.assertAlmostEqual(c, 0.0)

    def test_coupling_within_exciton_band_is_nonzero(self) -> None:
        c = self.agg.coupling(self.es1, self.es2)
        self.assertGreater(abs(c), 0.0)

    def test_coupling_is_symmetric(self) -> None:
        c12 = self.agg.coupling(self.es1, self.es2)
        c21 = self.agg.coupling(self.es2, self.es1)
        self.assertAlmostEqual(c12, c21)
