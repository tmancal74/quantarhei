# -*- coding: utf-8 -*-

from tests.qm.hilbertspace.test_operators import TestReducedDensityMatrix


t = TestReducedDensityMatrix()
t.setUp()
t.test_excitation_by_delta()
t.test_conversion_2_populations()


