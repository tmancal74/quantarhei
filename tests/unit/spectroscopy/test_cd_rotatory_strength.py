"""Rotatory strengths and CD of excitonic aggregates (issue #268).

The expected values here do not come from the site-basis matrices of the
aggregate. They come from

* the analytic homodimer result: with equal site energies E, coupling
  J > 0, dipoles d_1, d_2 at R_1, R_2 the excitons are
  |+-> = (|1> +- |2>) / sqrt(2) with E_+- = E +- J and
  c_1 c_2 = +-1/2, so the length-form rotatory strength is

      R_+- = E_+- c_1 c_2 (R_1 - R_2) . (d_1 x d_2)
           = +-(E +- J) / 2 (R_1 - R_2) . (d_1 x d_2)

  (standard point-dipole exciton result, e.g. Somsen, van Grondelle and
  van Amerongen, Biophys. J. 71 (1996) 1934, without the 1/2c prefactor that
  the calculator omits),
* CircDichSpectrumCalculator._excitonic_rot_dip, which evaluates
  sum_{k<l} c_k c_l (R_l - R_k) . (d_l x d_k) directly from monomer data,
* translation invariance, which the physical rotatory strength must have.
"""

import unittest

import numpy

from quantarhei import (
    Aggregate,
    CorrelationFunction,
    Molecule,
    TimeAxis,
    energy_units,
)
from quantarhei.exceptions import QuantarheiError
from quantarhei.spectroscopy.abscalculator import LinSpectrumCalculator
from quantarhei.spectroscopy.circular_dichroism import CircDichSpectrumCalculator

DIPS = (numpy.array([1.0, 0.0, 0.0]), numpy.array([0.3, 1.0, 0.2]))
POSS = (numpy.array([0.0, 0.0, -2.5]), numpy.array([1.0, 0.0, 2.5]))
TA = TimeAxis(0.0, 1000, 1.0)


def _dimer(
    shift=(0.0, 0.0, 0.0),
    energies=(12000.0, 12300.0),
    J=100.0,
    velocity=(),
    positions=True,
    environment=False,
):
    """Coupled dimer (energies and J in 1/cm), optionally translated."""
    mols = []
    with energy_units("1/cm"):
        for k in range(2):
            m = Molecule(elenergies=[0.0, energies[k]])
            m.set_dipole(0, 1, DIPS[k])
            if positions:
                m.position = POSS[k] + numpy.asarray(shift)
            if k in velocity:
                m.set_velocity_dipole_from_dipole()
            if environment:
                params = dict(ftype="OverdampedBrownian", reorg=20, cortime=100, T=300)
                m.set_transition_environment((0, 1), CorrelationFunction(TA, params))
            mols.append(m)
        agg = Aggregate(molecules=mols)
        agg.set_resonance_coupling(0, 1, J)
    agg.build()
    return agg


def _rotatory_strengths(agg):
    HH = agg.get_Hamiltonian()
    SS = HH.diagonalize()
    energy = numpy.diag(HH.data).copy()
    calc = LinSpectrumCalculator(TA, system=agg)
    rot = [
        calc._excitonic_rotatory_strength_fullv(SS, agg, energy, n)
        for n in range(1, HH.dim)
    ]
    return numpy.array(rot), SS, energy


# Translations are in Angstrom; scale of R_n is ~3 (int. energy * A * D^2),
# rounding errors are ~1e-15 relative, so rtol=1e-10 is a strict check.
SHIFTS = ([10.0, 0.0, 0.0], [0.0, 7.0, -3.0], [-25.0, 40.0, 13.0])


class TestLengthForm(unittest.TestCase):
    def test_analytic_homodimer(self):
        E, J = 12000.0, 100.0
        agg = _dimer(energies=(E, E), J=J)
        rot, _, _ = _rotatory_strengths(agg)
        geom = numpy.dot(POSS[0] - POSS[1], numpy.cross(DIPS[0], DIPS[1]))
        with energy_units("1/cm"):
            Eint = agg.convert_energy_2_internal_u(numpy.array([E - J, E + J]))
        expected = numpy.array([-Eint[0] / 2, Eint[1] / 2]) * geom
        numpy.testing.assert_allclose(rot, expected, rtol=1e-10)

    def test_translation_invariance(self):
        ref, _, _ = _rotatory_strengths(_dimer())
        self.assertTrue(numpy.all(numpy.abs(ref) > 1.0))
        for shift in SHIFTS:
            rot, _, _ = _rotatory_strengths(_dimer(shift=shift))
            numpy.testing.assert_allclose(rot, ref, rtol=1e-10)

    def test_agrees_with_circular_dichroism_calculator(self):
        # CircDichSpectrumCalculator uses the same point-dipole model without
        # the transition-energy factor: R_n = E_n * rot_dip_n.
        agg = _dimer(shift=SHIFTS[1])
        rot, SS, energy = _rotatory_strengths(agg)
        cdc = CircDichSpectrumCalculator(TA, system=agg)
        ref = [
            (energy[n + 1] - energy[0]) * cdc._excitonic_rot_dip(SS, agg, n)
            for n in range(2)
        ]
        numpy.testing.assert_allclose(rot, ref, rtol=1e-10)

    def test_cd_spectrum_translation_invariant_and_nonzero(self):
        spectra = []
        for shift in ([0.0, 0.0, 0.0], SHIFTS[2]):
            agg = _dimer(shift=shift, environment=True)
            calc = LinSpectrumCalculator(TA, system=agg)
            with energy_units("1/cm"):
                calc.bootstrap(rwa=12150.0)
            spectra.append(calc.calculate()["CD"].data)
        scale = numpy.max(numpy.abs(spectra[0]))
        self.assertGreater(scale, 0.0)
        numpy.testing.assert_allclose(spectra[1], spectra[0], atol=1e-10 * scale)

    def test_missing_position_gives_no_extrinsic_cd(self):
        rot, _, _ = _rotatory_strengths(_dimer(positions=False))
        numpy.testing.assert_array_equal(rot, 0.0)


class TestVelocityForm(unittest.TestCase):
    def test_analytic_homodimer(self):
        # v_a = -i E d_a: R_+- = +-E^2 / (2 E_+-) (R_1 - R_2) . (d_1 x d_2)
        E, J = 12000.0, 100.0
        agg = _dimer(energies=(E, E), J=J, velocity=(0, 1))
        rot, _, _ = _rotatory_strengths(agg)
        geom = numpy.dot(POSS[0] - POSS[1], numpy.cross(DIPS[0], DIPS[1]))
        with energy_units("1/cm"):
            Es = agg.convert_energy_2_internal_u(numpy.array([E, E - J, E + J]))
        expected = numpy.array([-1.0, 1.0]) * Es[0] ** 2 / (2 * Es[1:]) * geom
        numpy.testing.assert_allclose(rot, expected, rtol=1e-10)

    def test_translation_invariance(self):
        ref, _, _ = _rotatory_strengths(_dimer(velocity=(0, 1)))
        self.assertTrue(numpy.all(numpy.abs(ref) > 1.0))
        for shift in SHIFTS:
            rot, _, _ = _rotatory_strengths(_dimer(shift=shift, velocity=(0, 1)))
            numpy.testing.assert_allclose(rot, ref, rtol=1e-10)

    def test_extrinsic_moment_equals_intrinsic_magnetic_dipole(self):
        # Moving R_a x v_a from the position into an intrinsic magnetic
        # dipole m_a (molecules without positions) must give the same
        # rotatory strength; this pins the relative sign of RRv and RRm.
        ref, _, _ = _rotatory_strengths(_dimer(velocity=(0, 1)))
        agg = _dimer(velocity=(0, 1), positions=False)
        for k, mol in enumerate(agg.monomers):
            v = mol.dvmoments[0, 1, :]
            mol.set_magnetic_dipole(0, 1, numpy.cross(POSS[k], v))
        agg.rebuild()
        rot, _, _ = _rotatory_strengths(agg)
        numpy.testing.assert_allclose(rot, ref, rtol=1e-10)

    def test_mixed_velocity_dipoles_raise(self):
        agg = _dimer(velocity=(1,))
        with self.assertRaises(QuantarheiError):
            _rotatory_strengths(agg)


if __name__ == "__main__":
    unittest.main()
