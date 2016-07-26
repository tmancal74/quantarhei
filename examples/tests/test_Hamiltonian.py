# -*- coding: utf-8 -*-
"""
    Quantarhei units tests to be run in console

"""
from tests.qm.hilbertspace.test_hamiltonian import TestHamiltonian

t = TestHamiltonian()
t.setUp(verbose=True)

t.test_units_management()

