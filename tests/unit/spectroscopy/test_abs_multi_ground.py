import unittest

import numpy

from quantarhei import (
    AbsSpectrumCalculator,
    CorrelationFunction,
    Mode,
    Molecule,
    TimeAxis,
    energy_units,
)


class TestAbsMultiGroundState(unittest.TestCase):
    """Regression test for #272: _spect_from_dyn with multiple ground states.

    A molecule with a vibrational mode having nmax(0) > 1 produces multiple
    ground-state levels. The _spect_from_dyn function must iterate over ALL
    ground states, not just the first one.
    """

    def test_multi_ground_state_absorption_from_dynamics(self):
        time = TimeAxis(0.0, 500, 1.0)

        with energy_units("1/cm"):
            mol = Molecule(elenergies=[0.0, 12000.0])
            mol.set_dipole(0, 1, [0.0, 1.0, 0.0])
            params = dict(ftype="OverdampedBrownian", reorg=30, cortime=100, T=300)
            cf = CorrelationFunction(time, params)
            mol.set_transition_environment((0, 1), cf)

            mod = Mode(frequency=800.0)
            mol.add_Mode(mod)
            mod.set_nmax(0, 3)
            mod.set_nmax(1, 3)
            mod.set_HR(1, 0.1)

        mol.set_electronic_rwa([0, 1])

        prop = mol.get_ReducedDensityMatrixPropagator(
            time, relaxation_theory="stR", time_dependent=True
        )

        abs_calc = AbsSpectrumCalculator(time, system=mol)
        abs_calc.bootstrap(prop=prop)

        abs_alt = abs_calc.calculate(from_dynamics=True, alt=True)
        abs_single = abs_calc.calculate(from_dynamics=True, alt=False)

        with energy_units("1/cm"):
            y_alt = abs_alt.data
            y_single = abs_single.data

        max_val = numpy.max(numpy.abs(y_single))
        self.assertGreater(max_val, 0.0)

        diff = numpy.max(numpy.abs(y_alt - y_single))
        rdiff = diff / max_val
        self.assertLess(rdiff, 0.01)

        ham = mol.get_Hamiltonian()
        n_ground = ham.rwa_indices[1]
        self.assertGreater(n_ground, 1)


if __name__ == "__main__":
    unittest.main()
