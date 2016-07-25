# -*- coding: utf-8 -*-

from tests.build.test_molecules import TestMolecule
from tests.build.test_molecules import TestMoleculeVibrations

t1 = TestMolecule()
t1.setUp()

t1.test_get_Hamiltonian()

t2 = TestMoleculeVibrations()
t2.setUp()

t2.test_molecule_with_vibrations_1()
t2.test_thermal_density_matrix_0_temp()
#t2.test_thermal_density_matrix_finite_temp()
t2.test_thermal_density_matrix_finite_temp_nondiag()
